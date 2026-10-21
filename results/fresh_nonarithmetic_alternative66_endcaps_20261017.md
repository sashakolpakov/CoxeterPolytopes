# Five wall-66 alternatives: endpoint comparison

This is a final bounded numerical construction diagnostic, not a new exact chamber or a nonexistence theorem. Generation and fresh saved replay both pass.

The [source](../scripts/fresh_nonarithmetic_alternative66_endcaps_20261017.py) and [data](fresh_nonarithmetic_alternative66_endcaps_20261017.json) compare the five native-admissible wall-66 recipes with labels `(2,2,4)`, `(3,2,3)`, `(3,2,4)`, `(4,2,3)`, `(4,2,4)` on the fixed frame `(58,7,33)`. They preserve the exact predecessor [63-wall control](fresh_nonarithmetic_release49_62_replacement_20261017.md) and its [complete intrinsic boundary](fresh_nonarithmetic_release49_62_boundary_20261017.md). Only the first alternative has that independent exact boundary certificate; the other four are compared numerically here.

For each alternative, the two endpoint frames are `(33,30,66)` and `(58,60,66)`. A cap must yield finite corners with the first and second anchors, and with the first and third anchors. All seven allowed A2 contact profiles are retained. The complete finite rank-at-most-five diagram classification reduces finite labels to 2–5 except for isolated I2(m) components. A symbolic label 6 records each entire m≥6 tail; it is not a cutoff. The inherited H4 support admits no nonorthogonal finite extension. The old wall60 has its nonzero A2 profile, which is included in the corner tests and all native product offsets.

| 66 labels | Left low / tail profiles | Right low / tail profiles | Accepted left / right low caps |
|---|---:|---:|---:|
| 2,2,4 | 38 / 2 | 16 / 2 | 28 / 0 |
| 3,2,3 | 50 / 4 | 9 / 0 | 42 / 0 |
| 3,2,4 | 38 / 2 | 9 / 0 | 24 / 0 |
| 4,2,3 | 50 / 4 | 5 / 0 | 41 / 0 |
| 4,2,4 | 38 / 2 | 5 / 0 | 23 / 0 |

There are 516 real low-label branches. All 158 accepted left caps pass the numerical native relations and common-center test. Each has a separate positive-definite support census and a full projective triple-intersection diagnostic; each retains 10 intrinsic open ends. Thus none is a strict boundary improvement. All future right-end branches have wall59 as a blocker. The [release59 follow-up](fresh_nonarithmetic_alternative66_release59_20261017.md) records its deletion costs rather than silently dropping that constraint.

Every one of the 16 high-label families has 33 continuous cosine samples, including the limiting value 2. These samples are diagnostics, not discrete Coxeter normals or an exhaustive interval proof. No high-tail nonexistence claim is made. No low paired candidate survives because the right pool is empty.

One left cap changes the exterior graph substantially despite its unchanged open-end count. The [changed-chain follow-up](fresh_nonarithmetic_alternative66_chain_patch_20261017.md) uses that adjacency change. No new geometry in this note is promoted as exact, and the finite-volume construction goal remains open.

Verification: generation session93829 and fresh saved replay49758 both exited0. Frozen predecessor files are unchanged.
