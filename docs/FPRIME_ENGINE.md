# F_p (prime field) engine — design

Extend WDSat to natively solve systems of polynomial equations over a prime
field F_p, the natural target of Weil-descent ECDLP reductions. The existing
GF(2) engine (cnf/xorset/xorgauss) is left untouched; this is a parallel
engine selected by input format, reusing the search driver pattern and the
copy-on-write journal/undo design.

## 1. Scope

- Input: polynomial system over F_p, monomial degree <= 4 (same cap as ANF).
- Prime p < 2^24 (Weil descent only targets small fields). Arithmetic via
  64-bit products + `% p`; fast paths later if needed.
- Deliverable semantics: a solution (assignment of the original variables to
  F_p values) or UNSAT. The printed model uses variable IDs in file order.
- Out of scope for now: binary fields F_{2^k} (that is the existing engine),
  p >= 2^24 (correct but slow), ECDLP instance generation (see §9.4).

## 2. Input format (`.fp` files)

```
p fp <n_vars> <n_eqs> <p>
e <c1> <t1> <c2> <t2> ... 0            (polynomial equation, see term syntax)
c <c1> <t1> ...                        (constant part folded into the lhs)
```

- Term syntax: `v` (variable, degree 1) or `.d a b` / `.d a b c` / `.d a b c d`
  (degree 2-4 monomial). Coefficients are non-negative ints reduced mod p.
- Equation: sum c_k * term_k == 0 (mod p).
- The reader pre-enumerates the **term lattice**: every sub-monomial of every
  monomial gets its own term ID (1..n_terms), exactly like the GF(2) engine's
  monomial variables. Original variables keep IDs 1..n_vars; term IDs above
  are internal. `__MAX_ID__` bounds n_terms (767 today; raise for larger
  systems, see §8).

## 3. Representation

- `fp_p` (uint32_t) global prime; `fp_neg`, `fp_inv` helpers
  (inverse by Fermat for general p; table for p <= 2^16).
- Row = coefficient vector: `uint32_t row[n_terms]` + `uint32_t row_const`.
  Row width = n_terms+1 words, computed once at load (`fp_row_words`), so
  journal copies only the meaningful width (same trick as `xg_row_words`).
- State per term i:
  - `fp_row[i]` — current expression of term i as a polynomial in terms
    (after substitutions), or a constant row when i is assigned.
  - `fp_alive[i]` — row i is an active equation (GF(2): `xorgauss_equivalent`).
  - `fp_val[i]` / `fp_undef[i]` — assigned value in F_p or undef
    (GF(2): `xorgauss_assignment` ±).
  - `fp_degree[i]` — number of currently unassigned factors (monomials only).
- `fp_terms_of_var[v][...]` — for each original variable v: the monomial IDs
  containing it, each with the ID of its remaining-factor term (precomputed
  from the lattice; GF(2): `monomials_to_column`).

Memory (n_terms = N, row = (N+1) words):

| item | size |
|---|---|
| row matrix (N rows) | N*(N+1)*4B |
| values/flags/degrees | ~9N B |
| journal | dynamic arena, §7 |

N=767: matrix ~2.4 MB. N=4096: ~67 MB (config change required).

## 4. Arithmetic

- `fp_mul(a,b) = (uint64_t)a*b % p` (p < 2^24 => product < 2^48, safe).
- `fp_div(a,b) = fp_mul(a, fp_inv(b))`; b==0 only via checked code paths.
- All reductions happen at the row level (after whole-row ops), not per
  cell, where possible: row add/sub/mul-scalar produce uint64 intermediates
  reduced once per cell.

## 5. Propagation

### 5.1 Assignment of term v to value a (`fp_set_true` analog: `fp_assign(v, a)`)

Worklist-driven, mirroring `xorgauss_set_true`/`xorgauss_infer`:

1. If v is already assigned: check consistency, else conflict (false).
2. If row v is a reduction of another term t (already "equivalent"): check a
   against the evaluation of t, conflict on mismatch, else propagate v's
   assignment to t and continue.
3. Otherwise (fresh row): set row v to the constant a; for every alive row i
   with nonzero coefficient c at term v:
   - **a != 0**: row i <- row i - (c/a) * row v. row v is now constant, so
     this removes the v term and folds the rest in.
   - **a == 0**: row i's v coefficient is just zeroed.
   - After the update, classify row i:
     - constant row -> assign the term of row i to that constant;
     - single nonzero term t with coefficient c -> assign t to const/c;
     - else leave alive.
4. Push new assignments onto the worklist (cascade).

### 5.2 Monomial collapse (degree bookkeeping)

When original variable v is assigned a:
- **a == 0**: every monomial m containing v is assigned 0 (cascade);
  `fp_degree[m] = 0`.
- **a != 0**: for each monomial m containing v with remaining-factor term r
  (r may be a variable, a pre-enumerated lower monomial, or 1 = "none"):
  - m is replaced by `a * r` everywhere: this is the F_p generalization of
    `xorgauss_replace` — substitute term m by `a*r` in all alive rows (same
    "for every alive row containing m" scan), then mark m reduced.
  - `fp_degree[m]` drops by 1; when it reaches 0 the monomial equals the
    product of all absorbed values (already constant by then).

The term lattice being static means `a * r` is always an existing term (or a
pure scalar times an existing term) — no dynamic term creation, same property
the GF(2) engine relies on.

### 5.3 Linear core

Degree-1 equations are propagated by the same row-reduction machinery (a
linear row is just a row with one nonzero per... no: a linear row has several
term-1 entries). Gaussian-elimination-style: when a row reduces to
const + c*t, the term t is immediately assigned. Full linear systems are
therefore solved with no search at all (Phase 1 deliverable).

## 6. Search

DFS over original variables in a static order (identity first; heuristics
later). At each node:
- pick the next unassigned variable v in order;
- if the linear core has already determined v (row v reduced to a constant)
  skip (propagation consumed it);
- otherwise branch over values. Value order: 0 first (zero collapses the most
  terms), then 1, 2, ... — cheap to compute, effective pruning.
- Stop condition: all original variables assigned and all rows satisfied.

Domain size is p: p-ary branching is only tractable because (a) the linear
core eliminates most terms, (b) the remaining free terms in Weil-descent
ECDLP systems are few (the log's digit count), (c) zero-branching collapses
monomials. This matches the GF(2) engine's profile (small DFS + strong
propagation). If a system has many free terms, the solver reports it is
searching p-ary — a documented limitation, not a bug.

## 7. Undo / journal

Generalize the XG journal 1:1:
- Segment per search node (breakpoint macro pattern).
- Copy-on-write for rows: first modification of a row in a segment saves
  (row, alive flag, degree, val/undef state) into the segment.
- Per-row segment mark **stack** (the nested-undo bug in the GF(2) journal is
  already understood; implement the stack from the start).
- Assignment history: (term, value, was_undef) triples instead of signed
  literals.
- Journal storage: dynamic arena (the static worst case 52 x N x (N+1) words
  is ~123 MB at N=767 with 4-byte coeffs — too big for BSS). Arena of
  (term_id, flags, row-copy) entries with per-segment offset list; 256 MB cap
  with a clear error.

## 8. Config

- `__MAX_ID__`: keep 767 for Phase 1-2 (ECC Weil-descent systems fit: k=8
  descent => 8 log digits + O(10^2) monomial terms). Raise via config when a
  real instance demands it; matrix cost is quadratic.
- `__MAX_DEGREE__ = 4` unchanged.
- New: `__FP_MAX_P__ = 1<<24`, journal arena cap.

## 9. Verification

Independent reference: Python checker (`wdsat_ref_fp.py`): parse `.fp`,
evaluate polynomials mod p, (a) verify any model, (b) brute-force for small
systems (n_vars <= 8, p <= 31).

1. **Linear**: random linear systems over F_p (n up to 200, p random prime
   <= 10^6) — verdict must match numpy-free rational-free Gaussian
   elimination; solutions verified by substitution.
2. **Nonlinear, brute-forceable**: random sparse degree <= 3 systems, n <= 7,
   p <= 13, planted solutions — verdict + model must match exhaustive search.
3. **Property**: model of every SAT result re-evaluated by the checker.
4. **ECDLP smoke test (Phase 3+)**: small curve over F_p with a point of
   small order, known log; descend via the formal-group log equations
   (generated in Sage as a throwaway script, instance committed as a
   benchmark). The engine must recover the log.

## 10. Roadmap

| phase | content | exit criterion |
|---|---|---|
| 0 | this design | review |
| 1 | `fp.c` arithmetic, `.fp` reader + term lattice, linear solve (Gaussian elimination, RREF, free vars := 0), Python checker (`wdsat_ref_fp.py`) | §9.1 green — done: 40/40 randomized linear trials + hand tests (p=2..10^5) |
| 2 | monomial collapse + substitution + worklist cascade | §9.2 green |
| 3 | DFS + journal/undo (arena) + value ordering | §9.2 green, models verified |
| 4 | ECDLP instance generation + first real instance | §9.4 green |
| 5 | perf (row packing, Barrett reduction, row-width scaling) | timing report |

GF(2) benchmarks must stay byte-identical throughout (separate code paths,
shared only config).
