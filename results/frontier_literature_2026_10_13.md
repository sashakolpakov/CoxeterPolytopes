# Literature and new nonarithmetic face search, 2026-10-13

No new-dimensional hyperbolic Coxeter polytope has been established by this
session. The repository is explicitly a frontier-search project, and its
primary compact target remains dimension nine. Existing high-dimensional
positive entries carefully distinguish a new algebraic realization from a
new abstract polytope.

## Literature checked

- Compact examples remain known through dimension eight; the general upper
  bound is 29. Finite-volume examples are known through dimension 19 and in
  dimension 21; the general upper bound is 995. These are distinct frontiers.
  See the [Felikson–Tumarkin maintained page](https://www.maths.dur.ac.uk/users/anna.felikson/Polytopes/polytopes.html).
- The noncompact arithmetic case is sharper: dimension 20 and dimensions
  above 21 are excluded by Esselmann. Consequently a first noncompact
  finite-volume example in dimension 20 would have to be nonarithmetic, and a search among
  integral Lorentzian lattices cannot find one. See
  [Belolipetsky, introduction and §6](https://arxiv.org/pdf/1506.03111).
- An important update to the repository's August literature checkpoint is
  [McLeod, arXiv:2608.02894v2](https://arxiv.org/html/2608.02894v2), dated
  2026-10-05. It explicitly acknowledges that facets need not be Coxeter
  polytopes and searches all 304 types surviving the Lannér obstruction,
  with DRAT checks. Thus the repository's objection to version one is not
  automatically an objection to version two. This session has not replayed
  those external certificates. The paper claims uniqueness of the known
  compact H6/10 chamber, not a new dimension.
- [Guo–Ma–Zang–Zheng, arXiv:2607.14715](https://arxiv.org/html/2607.14715v1)
  constructs infinite sets of pairwise incommensurable noncompact examples
  in dimensions 4, 5, 6, 7, and 9. This does not extend the dimension frontier.
  Allcock's older infinite-population construction through dimension 19
  should not be confused with an infinite population of *commensurability
  classes* or with new nonarithmetic dimensions.

## A different constructive experiment

[Vinberg's 2015 paper, §§3–8](https://sfb701.math.uni-bielefeld.de/preprints/sfb14047.pdf)
glues the 37-wall chamber of `I_(18,1)` to the 24-wall chamber of
`II_(17,1) + <2>` along congruent orthogonal facets. The common facet is the
19-wall chamber of `II_(17,1)`. The resulting finite-volume noncompact H18
chamber is nonarithmetic. Certain Coxeter faces give examples in dimensions
14 and at most 12; the sufficient face criterion excludes components of
types `A_n` and `D5`. Therefore testing *all* faces, including these omitted
types, is a legitimate additional experiment. It might populate new
dimensions within the nonarithmetic subclass, but would not by itself give
a first finite-volume example in those dimensions.

No implementation of this mixture was found by searching the repository.
The new script
[`explore_vinberg_h18_mixture_faces.py`](../scripts/explore_vinberg_h18_mixture_faces.py)
uses the repository's exact `I_(18,1)` root enumerator and the archived
rank-19 K3 wall Gram matrix, avoiding transcription of the paper's diagram.
It matches the common facet by its exact Gram graph, chooses an 18-vector
basis, and constructs every cross-parent pairing in `Q(sqrt(2))`.

The resulting Gram matrix has 40 walls and exact inertia `(18,1,21)`.
All original Coxeter angle tests pass. One two-cycle is `12+8sqrt(2)`;
together with noncompactness this independently excludes arithmeticity.
The two source chambers are established known examples, so the mixture is
an exact reconstruction of Vinberg's construction, not a new example.

Every one of the 40 facets was screened by an exact Schur complement:

| H17 facet outcome | Count |
|---|---:|
| Exact Gram match to known 22-wall `I_(17,1)` chamber | 3 |
| Exact Gram match to known 19-wall `II_(17,1)` chamber | 13 |
| Forbidden Coxeter angle | 24 |

All 19 facets crossing the gluing seam fail the angle gate. The 16 positive
facets lie in one arithmetic parent. Thus the apparently promising H17
outputs are independently identified rediscoveries.

The script also exhausts seam-crossing codimension-two, -three, and -five
faces. These are natural H16, H15, and H13 candidates. This is not a complete
search of all nonarithmetic faces: a Coxeter face contained in an arithmetic
parent is only guaranteed to be quasi-arithmetic, and can in principle be
properly quasi-arithmetic. Such faces require separate arithmeticity tests.
The completed exact search gives:

| Face dimension | Seam subsets checked | Coxeter survivors |
|---|---:|---:|
| 17 | 19 | 0 |
| 16 | 171 | 0 |
| 15 | 969 | 0 |
| 13 | 11,628 | 0 |

Every subset in these rows is elliptic, and every rejection has an explicit
forbidden-angle witness in the archive. This closes 12,787 seam-crossing
face candidates of this one known H18 chamber.

For each selected elliptic subset `S`, its facets are exactly the elliptic
one-node extensions; their intrinsic Gram matrix is
`G_ij - G_iS G_SS^-1 G_Sj`.

The only finite-angle squares possible in `Q(sqrt(2))` are
`0, 1/4, 1/2, 3/4, (2+sqrt(2))/4`, corresponding to labels
`2,3,4,6,8`. This follows from the degree and field of
`cos(2*pi/m)`. Numerical arithmetic only suggests a candidate forbidden
pair. Every rejection is checked in the exact field, and a complete exact
pair scan is used if no numerical witness is proposed.

The independent checker
[`verify_vinberg_h18_mixture_faces.py`](../scripts/verify_vinberg_h18_mixture_faces.py)
uses its own `Fraction` Gaussian elimination and quadratic sign routine,
without importing the search or its arithmetic implementation. It compares
the witness domain to every eligible seam subset, so missing records do not
silently count as exclusions. Its complete replay returned
`PASS_INDEPENDENT_RATIONAL_REPLAY` for all 12,787 witnesses. Both scripts
also pass Python compilation.

Run from the repository root:

```sh
PYTHONPATH=src:. python3 scripts/explore_vinberg_h18_mixture_faces.py --seam-codimensions 2 3 5
python3 scripts/verify_vinberg_h18_mixture_faces.py
```

The generated archive is
[`vinberg_h18_mixture_facets_2026_10_13.json`](../data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json).
This experiment only concerns faces of this specific H18 mixture. It
cannot exclude other polytopes, arithmetic fields, lattices, gluings,
dimensions, or deformation constructions.
