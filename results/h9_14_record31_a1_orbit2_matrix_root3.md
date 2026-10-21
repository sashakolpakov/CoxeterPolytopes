# Record 31, `a=1`, source orbit 2: matrix root class 3 closed

## Outcome

Canonical matrix root class 3 is completely excluded from the compact Coxeter
problem.  Its 19 depth-two classes were exhausted in 19 disjoint singleton
fragments.  Exactly five literal ridge-valid targets occur, all in depth-two
class 14, and they remain five distinct `C2 x S9` target orbits.  All five
pass the connected mod-two homology-sphere audit, but production-order and
hostile-order unrestricted compact Coxeter-label joins both leave zero
terminal states.

This closes root class 3 only.  Together with the earlier root-1 and root-2
results, three of the six canonical root classes of source orbit 2 are now
closed.  Root classes 4 through 6 remain open, so neither source orbit 2 nor
the full Record-31 `a=1` schedule is yet excluded.  Nothing here concerns
noncompact finite-volume chambers.

## Exact partition

Every fragment terminated normally under a 0.60 GiB watchdog.  Each archive
pins the generated C++ source, compiled binary, watchdog, stderr digest,
watchdog disposition, elapsed time, and peak RSS.

| depth-two class | nodes | pure 120-facet targets | literal ridge targets | peak GiB | seconds |
|---:|---:|---:|---:|---:|---:|
| 1 | 1,766,599 | 9 | 0 | 0.082 | 37.714 |
| 2 | 3,757,715 | 0 | 0 | 0.165 | 85.084 |
| 3 | 866,390 | 19 | 0 | 0.040 | 15.394 |
| 4 | 891,233 | 16 | 0 | 0.042 | 15.911 |
| 5 | 4,706,131 | 0 | 0 | 0.247 | 91.247 |
| 6 | 1,322,808 | 1 | 0 | 0.069 | 26.328 |
| 7 | 4,676,314 | 12 | 0 | 0.236 | 105.555 |
| 8 | 3,037,300 | 1 | 0 | 0.148 | 59.719 |
| 9 | 2,121,105 | 1 | 0 | 0.087 | 39.141 |
| 10 | 660,987 | 1 | 0 | 0.027 | 11.512 |
| 11 | 4,943,472 | 35 | 0 | 0.239 | 99.260 |
| 12 | 4,010,976 | 1 | 0 | 0.179 | 83.027 |
| 13 | 2,107,396 | 5 | 0 | 0.094 | 40.671 |
| 14 | 2,882,783 | 877 | 5 | 0.099 | 47.259 |
| 15 | 2,245,567 | 60 | 0 | 0.079 | 35.996 |
| 16 | 2,906,860 | 216 | 0 | 0.121 | 46.829 |
| 17 | 1,559,497 | 9 | 0 | 0.069 | 33.061 |
| 18 | 802,477 | 7 | 0 | 0.038 | 13.611 |
| 19 | 708,044 | 19 | 0 | 0.036 | 12.211 |

The ranges are disjoint and cover exactly `1,...,19`.  Fragment node counts,
pure-target counts, and literal multiplicities are nonadditive under restarted
deeper quotient tables and are not combined into a theorem-facing total.

The five ridge targets canonicalize to five families.  Direct blocker/ridge
replay gives `{5:120}` and `{2:540}` for each.  All five have the scheduled
face vector

```text
(14,87,319,770,1281,1477,1142,540,120)
```

and pass every connected mod-two link audit.  The two unrestricted label
joins independently leave terminal state count zero; the hostile replay uses
a different support order and reconstructs the 108 order-four and 420
order-five classical Lanner placements.

## Artifacts and limits

- fragments:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_d2_*.json`;
- exact canonical union:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_target_union.json`;
- production downstream:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_downstream.json`;
- hostile downstream:
  `data/enumerations/h9_14_record31_a1_orbit2_matrix_root3_independent_downstream.json`;
- assembler:
  `scripts/assemble_h9_14_record31_a1_orbit2_matrix_partition.py`;
- verifier:
  `scripts/verify_h9_14_record31_a1_orbit2_matrix_root3.py`.

The generated engine SHA-256 is
`0bf66b0bd22827f90e4034ba80539eb1cc57b712cdd59e40070b4150e074ddfc`.
The canonical-union archive has SHA-256
`1b9100bea75c2017535487f8140dffaa416cb0d2d779b7bc0a4df0da9ebafc00`;
the production and hostile downstream archives have SHA-256 values
`191376eda54fbd448c0423fcb08d34b6f5882b8e170c9b837929d76f58fd63fe`
and
`7bc697e1e9ce300de967ad85633e1b8637c9d161189f1ee9c53054c7282e76cd`.

This is a terminal exact root-class theorem, not a timeout promotion, a Gram
realization, a full source-orbit theorem, or a noncompact statement.
