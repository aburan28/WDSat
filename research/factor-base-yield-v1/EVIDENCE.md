# Evidence layout

This pull request keeps the implementation and compact reviewed records in the
main diff. The full immutable research branch, including every generated mixed
ANF/CNF input, factor filter, variable map, solver stdout/stderr and intermediate
run, is published at:

`aburan28/WDSat`, branch `evidence/binary-sat-review-20260909`, commit
`55d55b2620d768d9f7c78dcd8990a0689533c1d0`.

The producer result is
`n19-same-dimension-ablation-20260909-a/result.json` on that evidence commit.
The independent compact checker result and byte-identical build receipt are
included in this pull request beside the review. To rerun the checker against
all raw files, check out the evidence commit and run:

```sh
DOT_SAGE=/private/tmp/reviewer-n19-ablation-sage \
  sage -python research/factor-base-yield-v1/independent-review-n19-ablation-check.py
```

The evidence branch is intentionally separate so the code-review diff is not
expanded by thousands of generated files and large custody archives.
