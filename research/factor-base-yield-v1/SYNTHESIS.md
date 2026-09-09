# Factor-base yield result for binary Koblitz Semaev index calculus

## Outcome and evidence class

This campaign found an independently reviewed, reproducible **bounded
implementation breakthrough** in factor-base relation yield on the tested
Koblitz curve over `GF(2^19)`. It is a measured small-field result, not a
cryptanalytic break, asymptotic improvement, novelty finding, standardized-curve
result, or comparison against automorphism-aware Pollard rho.

With curve

`E: y^2 + xy = x^3 + x^2 + 1`,

field modulus `z^19 + z^5 + z^2 + z + 1`, and prime subgroup order
`r = 262543`, fix three summands, dimension `k=6`, Frobenius schedule
`[0,1,2]`, exact factor-membership clauses, the Gaussian WDSat binary, and the
same stopping/rank rules. Replacing the fixed trace basis

`[2,4,8,16,32,64]`

with the membership-selected sparse basis

`[8,16,128,1024,4096,32768]`

at affine offset 1 produced the following independently verified differences:

| Quantity | Fixed trace-k6 | Selected sparse-k6 |
|---|---:|---:|
| Coordinate assignments | 64 | 64 |
| Valid subgroup x assignments | 23 | 45 |
| Valid-assignment density | 35.9375% | 70.3125% |
| Actual factor points including signs | 46 | 90 |
| Sign-folded columns | 23 | 45 |
| Exact membership clauses, three blocks | 123 | 57 |
| Exact reachable nonzero targets | 80,038 / 262,542 | 245,738 / 262,542 |
| Exact target coverage | 30.4858% | 93.5995% |

The selected directions are polynomial exponents `[3,4,7,10,12,15]`.
Membership-only enumeration of all 18,564 six-direction candidates, caching
31,180 distinct x inputs, reproduces this winner without a subgroup logarithm
table. The runtime computes the Frobenius eigenvalue as the matching root of
`T^2 - T + 2` modulo `r`; it obtains `lambda = 84450` using only public point
equality.

## Full relation extraction

Three streams used different public `Q` values and secrets. Every arm continued
until a relation dependency recovered the secret and `[d]P=Q` passed.

| Stream | Arm | Attempts | SAT relations | UNSAT | Solver s | Recorded per-target pipeline s |
|---|---|---:|---:|---:|---:|---:|
| 202609107 | Fixed | 90 | 24 | 66 | 103.177878 | 124.157804 |
| 202609107 | Sparse | 44 | 43 | 1 | 44.943814 | 59.726330 |
| 202609108 | Fixed | 93 | 23 | 70 | 92.995857 | 110.961086 |
| 202609108 | Sparse | 45 | 42 | 3 | 37.118965 | 47.775270 |
| 202609109 | Fixed | 68 | 22 | 46 | 65.058583 | 78.314258 |
| 202609109 | Sparse | 39 | 38 | 1 | 30.213126 | 39.309332 |

Median attempts were 90 versus 44. The ratio of solver-only medians is
`2.505346x`; the ratio of recorded per-target pipeline medians is `2.322563x`.
The selected basis used fewer attempts and less recorded time in every paired
stream despite having more matrix columns.

Solver time measures each WDSat subprocess. Per-target pipeline time also
includes formula/membership construction, artifact writes, model validation,
reachability lookup, witness/rewrite checks and rank insertion. The separately
recorded 13.793538-second validation setup is common to the arms and does not
represent complete program setup. Whole-program wall time was not measured.

## Independent validation

The independent addendum verdict is **CONCUR** for the narrow claim above. Its
checker imports no producer modules and verified:

- 379 independently regenerated mixed ANF/CNF formulas;
- all 192 SAT assignments, lifts, subgroup memberships and signed point sums;
- all 187 UNSAT targets against complete point-addition reachability sets;
- every Frobenius coefficient rewrite and all six relation matrices at every
  prefix;
- secrets 27183, 16180 and 14142 at precisely the recorded first useful
  dependency;
- all 1,937 files against Git snapshot `3c5e675`;
- a byte-identical rebuild of the native solver, SHA-256
  `64ef817176c5ebf0dc4623d20790b88efb8ced65920471cb5ad153ffada1adbc`.

The independent report is `independent-review-n19-ablation.md`, with compact
checker and build artifacts beside it. The curated pull request omits the
producer result and thousands of raw generated files. See `EVIDENCE.md` for the
immutable evidence branch and exact commit containing the complete retained
tree, including `n19-same-dimension-ablation-20260909-a/result.json`.

## What changed in the implementation

The work also repaired WDSat correctness and made the comparison executable:

1. Initial propagation contradictions now return UNSAT rather than a false
   assignment.
2. Partial branching prefixes now complete every model; unsafe symmetry/order
   combinations are rejected.
3. The enhanced parser accepts mixed ANF and general CNF clauses.
4. A right-sized `s3-chain-19` build profile supports these descended systems.
5. The generator supports distinct affine factor spaces per summand, exact
   membership clauses, retained variable maps and point certificates.
6. Frobenius-image relations rewrite onto one column set with verified signed
   powers of `lambda`.
7. Full relation matrices recover the public discrete logarithm rather than
   stopping at polynomial roots or coverage estimates.

Build and core validation commands are:

```sh
make -C src
python3 -m unittest discover -s tests -p 'test_binary_solver.py' -v
DOT_SAGE=/private/tmp/codex-sage-factor-base sage -python tests/check_binary_factor_base.py
DOT_SAGE=/private/tmp/codex-sage-factor-base sage -python tests/check_yield_screen.py
DOT_SAGE=/private/tmp/codex-sage-factor-base sage -python tests/check_yield_controls.py
make -C src s3-chain-19
```

## Prior art and remaining boundary

SAT solving for binary summation-polynomial systems, affine/symmetry-breaking
factor spaces, dimension tuning, and Frobenius-aware relation rewriting are all
prior art. The experiment establishes the behavior of this selected sparse
basis under the frozen implementation. It does not establish novelty of the
combination. Relevant primary sources are:

- Galbraith and Gebregiyorgis, *Summation polynomial algorithms for elliptic
  curves in characteristic two*, ePrint 2014/806.
- Trimoska, Ionica and Dequen, *A SAT-Based Approach for Index Calculus on
  Binary Elliptic Curves*, ePrint 2019/313.
- Galbraith, Granger, Merz and Petit, *On Index Calculus Algorithms for
  Subfield Curves*, ePrint 2020/1315 / SAC 2020.

The intervention is the basis choice. Density, sumset geometry, exact clause
shape and solver behavior are consequences that move together; the experiment
does not prove density alone causes the full timing gain. Stronger claims need
held-out curves and field representations, equal-density controls, larger
degrees, complete setup and program-wall accounting, repeated machine-level
timing with uncertainty rules, and comparison with an automorphism-aware rho
baseline.
