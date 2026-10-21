# Complete cone of the exact 36-wall H9 chamber

The [exact clipping program](../scripts/fresh_nonarithmetic_escape35_20261016_cone.py) and [saved census](fresh_nonarithmetic_escape35_20261016_cone.json) give the complete cone of the [36-wall chamber](fresh_nonarithmetic_escape35_20261016.md): **100 finite rays, no ideal rays, 1,068 hyperideal rays and no opposite-sheet timelike or null rays**. A separately certified strict spacelike escape proves infinite volume.

The input is the frozen [complete 32-wall cone](fresh_h9_joint_pair_20261015_two_extension.json), with its [independent full-cone review](fresh_h9_joint_pair_20261015_two_review.md). The program reconstructs every input ray from its literal active covectors, orients and normalizes it exactly, checks every initial wall sign and full active set, and verifies its norm classification. All 972 input rays are recovered: 90 finite and 882 hyperideal.

Four exhaustive halfspace clips then add walls 32–35 in order:

| Number of walls | Cut-away old rays | Crossing two-faces | Finite rays | Hyperideal rays | All rays |
|---:|---:|---:|---:|---:|---:|
| 32 | — | — | 90 | 882 | 972 |
| 33 | 32 | 112 | 92 | 960 | 1,052 |
| 34 | 5 | 25 | 95 | 977 | 1,072 |
| 35 | 145 | 167 | 98 | 996 | 1,094 |
| 36 | 49 | 123 | 100 | 1,068 | 1,168 |

There are no zero-pairing old rays at these four clipping stages. The increasing number of exterior rays is an important limitation: although the latest wall closes a compact quadrilateral and two old finite-edge ends, these additions have not closed the chamber globally.

## Completeness of each clip

The inherited cone is pointed and full-dimensional. The retained strict timelike point H preserves full dimension after every clip. For each negative/positive ray pair, intersect their active-wall sets. The pair spans a two-dimensional face precisely when the complete common-zero face contains exactly those two rays. The checker tests this using the full incidence bitsets of the current complete ray list, rather than assuming that a common eight-wall support is sufficient. The eight-common-wall preliminary filter is necessary because a two-face of a full-dimensional ten-dimensional pointed cone has active covector rank eight.

If the negative ray is v, the positive ray is w, and the new wall covector is l, the new boundary ray is

\[
l(w)v-l(v)w.
\]

Both coefficients are positive. Therefore every old inequality is retained, and its old active set is exactly the intersection of the two old active sets. The new wall pairing is checked to vanish exactly. This also proves that no unrecorded old active inequality appears through cancellation. Every new ray of a halfspace clip either was an old ray or arises from such a crossing two-face; thus the procedure preserves completeness inductively. The full active set uniquely identifies each extreme ray of this pointed cone, and duplicate masks are checked not to arise.

Each newly created ray is classified by exact norm and future orientation. The final finite active sets agree with all 100 vertices of the separate local finite-corner update. The saved data contain every final active set, every crossing provenance record, intermediate counts and a digest of the exact coordinate reconstruction. They do not require trusting an unstored numerical hull.

## Consequences for the next construction step

The full incidence list gives the following ray counts on the recent facets:

| Facet | Finite rays | Hyperideal rays |
|---:|---:|---:|
| 32 | 2 | 70 |
| 33 | 3 | 22 |
| 34 | 3 | 164 |
| 35 | 2 | 121 |

Of the 1,068 hyperideal rays, 761 lie on none of facets 33, 34 or 35. The other 307 lie on exactly one of those facets, split as 22, 164 and 121; no extreme ray lies on two of them. The latest three caps therefore do not yet connect their new boundary regions. This is a complete-cone diagnosis, not a proof that these regions cannot be connected by further compatible walls. A next cap or cap set should be assessed against multiple remaining regions and the full cone, rather than only one saved escape.

## Arithmetic and scope

Coordinates use eight formal coefficients over the frozen octic field, for the three positive radicals introduced by walls 33–35. Multiplication uses squarefree bit masks. Sign tests use outward rational interval bounds, with an exact recursive conjugate-norm fallback. Only base-field divisions occur while normalizing input rays; no new-extension division, irreducibility assertion or approximate rank test is used.

```sh
python3 scripts/fresh_nonarithmetic_escape35_20261016_cone.py
python3 scripts/fresh_nonarithmetic_escape35_20261016_cone.py --verify
```

Generation and fresh saved replay pass. The arithmetic implementation is separate from the local cap's nested-pair operations but shares its frozen base field and literal new normals. This is not an independent complete-cone review. The input 32-wall cone was independently reviewed; the four exhaustive exact clips and the final census are newly certified here. This source, data and report are final and frozen.
