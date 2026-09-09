# Binary Semaev / SAT / factor-base review

This is an implementation review and small-field diagnostic, not a
Crypto Autoresearcher measurement-admitted result or a cryptanalytic advance.
The code was developed in an isolated WDSat worktree from `83130eb`. The
Crypto Autoresearcher source was read, not modified. No hypothesis status,
experiment ledger, or claimed attack complexity was changed.

## Recommended direction

Start with a subgroup-compatible affine factor space `c + V`, not automatically
the polynomial span containing 1. Add the necessary trace constraints to SAT,
and choose V using actual subgroup-point counts, decomposition coverage,
encoding size, verified relation yield and eventual independent matrix rank.
For Koblitz curves, test Frobenius-invariant affine spaces when their dimension
exists. For general binary curves, field squaring need not preserve the curve;
no Frobenius orbit compression is inferred merely from an invariant x-set.

## Mathematical review

For the ordinary binary model

`y^2 + xy = x^3 + A*x^2 + B`, `B != 0`,

the third summation polynomial is

`S3(x,y,z) = (xy+xz+yz)^2 + xyz + B`.

The generator uses this formula with either two summands or a finite-auxiliary
chain for three summands. An algebraic root is not a complete relation
certificate: lift each coordinate in the actual field, recover signs, check
subgroup membership, and verify the point sum. This implementation rejects and
blocks unsuccessful factor-coordinate projections and continues SAT search.
It retains repeated summands and uses a separate exhaustive group-law oracle.
`exhausted_polynomial_models` describes the chosen encoding only; the finite
chain excludes infinity intermediates and no general completeness theorem for
all representation/symmetry choices is asserted.

Every point in an odd-order subgroup is a double. The doubling formula gives
`x(2P) = lambda^2 + lambda + A`, hence `Tr(x) = Tr(A)`. On `c+V` this becomes
one linear Boolean constraint per factor coordinate. It is a necessary filter,
not sufficient in general for rational lifting or subgroup membership.
For odd n, `Tr(1)=1`; translate trace-zero direction spaces by 1 when compatible
with the target curve. The tool also keeps unfiltered and random controls.

An affine space is squaring-invariant exactly when its direction space is
invariant and `c^2+c` belongs to that direction space. Under a normal basis,
Frobenius is a cyclic shift. For odd n, `T^n-1` is square-free, and invariant
linear-space dimensions are sums of degrees of its irreducible factors.
The code constructs the actual factor kernels, then checks rank and closure.
At n=7 the available dimensions are 0,1,3,4,6,7. At n=17 they are
0,1,8,9,16,17: there is no invariant direction space of dimension 6.
For n=163, the order of 2 modulo 163 is 162, giving dimensions 0,1,162,163.
An affine translation does not repair a missing invariant direction dimension.

Frobenius compression is established prior art. On the target prime subgroup,
`phi(P) = [lambda]P`; compressed relation coefficients must therefore include
signed powers of lambda. Applying Frobenius to one relation does not produce
new independent rows. Orbit-size reductions must account for stabilizers.
The laboratory counts exact signed orbit sizes on closed Koblitz factor bases;
it does not implement the compressed relation matrix or discrete-log extraction.
[Galbraith, Granger, Merz and Petit, On Index Calculus Algorithms for Subfield Curves](https://pure.royalholloway.ac.uk/files/39395430/Accepted_Manuscript.pdf)
(retrieved during this review).

The SAT literature motivates retaining parity constraints, introducing shared
product variables, and controlling branching on factor-coordinate variables.
Its propagation result and symmetry routine depend on the authors' symmetric
S4 representation. They cannot simply be transferred to a different S3 chain.
Published WDSat/Magma timing ratios are not reproduced here.
[Trimoska, Ionica and Dequen, A SAT-Based Approach for Index Calculus on Binary Elliptic Curves](https://pmc.ncbi.nlm.nih.gov/articles/PMC7334981/)
and [Parity (XOR) Reasoning for the Index Calculus Attack](https://arxiv.org/abs/2001.11229)
(retrieved during this review).

Semaev's proposed asymptotic improvement relies on solving assumptions; the
presence of cubic Boolean equations alone establishes no subexponential SAT
complexity. This review makes no such transfer.
[Semaev, New algorithm for the discrete logarithm problem on elliptic curves](https://arxiv.org/abs/1504.01175)
(abstract retrieved during this review).

## Code findings and changes

1. **Incorrect result after initial conflict, fixed.** `wdsat_solve` discarded
   the return from `wdsat_infer_unitary`. Input `p cnf 1 2` followed by unit
   clauses `1 0` and `-1 0` returned assignment `0`. The solver now reports
   UNSAT immediately.
2. **Partial models from a partial `-g` list, fixed.** With three variables,
   `x 2 3 0`, and `-g 1`, the old result was `022` (2 means unassigned).
   A branching prefix now appends all omitted unary variables. The symmetry
   terminal case likewise finishes unresolved auxiliary variables by search.
3. **Broken `-m` option, fixed.** Its getopt specification lacked a colon,
   despite consuming `optarg`. The option now requires an argument, and
   `getopt` is stored in int rather than char. Basic dimension checks were added.
4. **Unsafe ordering combinations, guarded.** Duplicate/out-of-range branching
   IDs are errors. Symmetry requires canonical first blocks and excludes `-t`.
   These checks cannot prove input symmetry. The laboratory does not enable `-b`.
5. **Generator optimization, implemented.** Sparse GF(2^n)-coefficient Boolean
   monomials use bit-mask unions for multiplication. Squaring squares only field
   coefficients and immediately applies Boolean idempotence. CNF shares AND
   gates and uses linear-size XOR chains. ANF output keeps native parity rows.
6. **Harness integration findings, not modified.**
   `/Volumes/SSD990/crypto-autoresearcher/src/semaev_tree.py` restricts the
   factor space to a polynomial basis, puts u variables before x variables,
   and defaults B to the field generator rather than 1 (so its default is not
   a Koblitz curve). Its planted-target function samples without replacement
   and can return the identity despite the downstream need for an affine x.
   These choices must be explicit before sharing inputs with the new encoder.

Existing limitations remain: WDSat uses compile-time capacity constants and a
legacy input parser; this patch is not general parser hardening. Initialization
still uses the Gaussian module even without `-x`, so the flags do not define
a completely Gaussian-free ablation. UNSAT retains process status 0. The ANF
export handles constant contradictions using opposing unit clauses to avoid
the legacy empty-XOR assumption. Generic solver failures and timeouts are
distinct from UNSAT in the laboratory.

## Retained diagnostic results

All runs below use GF(2^7), dimension 3, m=3, seed 20260909 and eight sampled
nonidentity targets. Exact coverage enumerates all three-fold group sums with
repetition, independently of SAT. Timing fields are diagnostic only: small
process-startup-dominated runs and concurrent activity preclude a speedup claim.

On the Koblitz curve A=B=1, the group order is 142 and subgroup order is 71:

| Factor space | Subgroup points | Exact coverage out of 70 | Verified sampled targets |
|---|---:|---:|---:|
| Polynomial linear | 6 | 28 | 2/8 |
| Random linear | 6 | 20 | 1/8 |
| Random affine | 6 | 38 | 7/8 |
| Trace-one affine control | 10 | 58 | 6/8 |
| Invariant linear, either kernel | 0 | 0 | 0/8 |
| Invariant affine, either kernel | 14 | 70 | 8/8 |

The two winning direction bases, in polynomial-bit representation for the
recorded field modulus, are `[98,76,48]` and `[2,4,16]`, both with offset 1.
Each corresponding factor base is one signed Frobenius orbit of size 14.
This is a toy candidate selection; one compressed column in this tiny subgroup
does not establish an asymptotic advantage.

The trace filter preserved every sampled relation outcome and reduced rejected
factor-coordinate projections from 99 to 0 across the Koblitz comparison.
The eight sampled targets are insufficient to rank general success probabilities;
the exact coverage column is the stronger quantity for this toy field.

On the general binary control A=1, B encoded as 3, the group/subgroup orders are
also 142/71. The trace-one affine control has 10 subgroup points and covers
64/70 targets (7/8 sampled verified), versus 38/70 and 40/70 for the two
field-invariant affine spaces. Squaring is not a curve endomorphism here;
the report deliberately omits compressed-orbit columns.

Artifacts are under `diagnostics/binary-sat-review-20260909/`: paired Koblitz
reports, a general binary report, and a native solver cross-check. Each raw
artifact directory retains its generator source, variable maps, original ANF,
every attempted CNF, solver stdout/stderr, command, exit status and input hash.
The archive contains the same directories for durable reproduction.

## Validation and continuation

- Six C-solver regression tests include 100 random ANF systems in four flag
  configurations, each compared to exhaustive truth tables.
- Seven Sage tests check descent against direct field arithmetic, affine and
  invariant spaces, trace predicates, both export polarities against SAT
  solvers, and actual point-sum/subgroup certificates.
- 864 WDSat invocations on 216 generated initial ANF systems agree with retained
  CaDiCaL SAT/UNSAT statuses. Every returned WDSat model is complete and satisfies
  the ANF under a separate evaluator. No `-b` assumption is used here.

The next useful scaling step is a frozen, independently reviewed comparison at
n=17: dimension-6 trace-compatible affine controls versus feasible dimension-8
invariant affine directions, charging their different SAT variable counts and
factor-base sizes. A separate arm can evaluate normal-basis Frobenius-related
factor spaces or orbit unions when the desired invariant dimension is absent.
Do not attribute permutation symmetry or full-orbit compression to those
encodings without checking the actual constraints and coefficient map.
Measure encoding cost, SAT/UNSAT/timeout mix, memory, valid relations per second,
matrix rank growth and final extraction cost against an automorphism-aware rho
baseline. Larger execution and scientific promotion remain separate work.
