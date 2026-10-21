# 1,399 local arithmetic quadrilateral cores for dimension13

The completed B1 dotted attachment queue contains1,399 singular stopping
states with at most15 roots. Every one now passes a separate exact local
geometry certificate: its common-white face is a compact quadrilateral
with four prescribed finite vertices, and its reflections preserve an
integral module for an admissible Q(sqrt5) form. The ambient chamber may
have infinite volume.

The [source](../scripts/fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.py)
and [certificate](fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.json)
reconstruct every field Gram from its finite and dotted labels. They check
the exact normalized products, totally positive root norms, nonpositive
physical off-diagonal entries, and a nonzero radical with alternating
black signs. Every adjacent corner is identified with separately checked
finite Coxeter components. A positive corner of rank n−2, a strictly
dotted indefinite pair, and the radical force inertia(n−2,1,1).

Deleting black root3 gives a full basis. Its conjugate Gram is checked
positive definite by exact Q(sqrt5) elimination, proving conjugate
inertia(n−1,0,1). All15,610 simple-cycle tests are algebraically integral.
Individual Cartan entries are not assumed integral: the finite simple-path
module construction supplies the invariant integral module from cyclic
integrality.

The [local realization lemma](fresh_singular_quad_local_realization_20261017.md)
proves existence of a strict timelike center, the four actual finite
corners, and completeness of their face cone. Its product-coordinate
description excludes additional ideal or exterior rays. A separate
read-only internal mathematical review of this lemma passed.

| Number of core roots | Certified cores |
|---|---:|
| 4 | 2 |
| 5 | 17 |
| 6 | 31 |
| 7 | 60 |
| 8 | 88 |
| 9 | 112 |
| 10 | 158 |
| 11 | 143 |
| 12 | 164 |
| 13 | 187 |
| 14 | 214 |
| 15 | 223 |

The [spectator extension](fresh_sqrt5_rank13_dotted_b1_spectators_20261017.md)
adds every finite orthogonal white forest of rank15−n, placing these
quadrilaterals in ambient dimension13. Completeness here is restricted
to the selected singular states of the B1 dotted queue. This does not
establish a complete dimension13 face domain or a new dimension bound.
The parent queue's full attachment replay remains a separate check.

Final generation and fresh saved replay passed, including every exact
matrix, corner, radical and cycle check. Syntax, source manifests and
local report links are checked separately. The package is frozen.

```sh
python3 scripts/fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017.py --verify
```
