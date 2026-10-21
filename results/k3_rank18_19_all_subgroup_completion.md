# Full subgroup closure for the rank-18 and rank-19 K3 wall stocks

## Result

Every subgroup class of the exact weighted wall-diagram groups has now been
enumerated for the two remaining high-rank finite-automorphism K3 chambers.

| Picard rank | complete walls | exact wall group | element classes | subgroup classes |
|---:|---:|---|---:|---:|
| 18 | 19 | $C_2$ | 2 | 2 |
| 19 | 24 | $S_3$ | 3 | 4 |

The proper subgroup classes and their fixed dimensions are:

| source | subgroup | conjugacy-class size | fixed rank | hyperbolic dimension | rank gate |
|---|---|---:|---:|---:|---|
| rank 18 | $C_2$ | 1 | 10 | 9 | below target |
| rank 19 | $C_2$ | 3 | 11 | 10 | retained |
| rank 19 | $C_3$ | 1 | 7 | 6 | below target |
| rank 19 | $S_3$ | 1 | 5 | 4 | below target |

Thus the complete non-involution part of the denominator contributes no
fixed section above $H^9$.  The sole retained class is the transposition
class already frozen as the rank-19 quotient-nef construction.  There is no
new hyperbolic Coxeter polytope and no new algebraic realization.

Together with the separate rank-15--17 completion, this closes **all subgroup
classes of all five complete rank-15--19 K3 wall stocks**.

## Direct theorem-grade replay of the sole survivor

Kondō realizes the full $S_3$ action on the rank-19 K3 surface.  Its three
transpositions are conjugate by realized surface automorphisms and give three
labelled copies of the same fixed nef cone.  A representative has fourteen
wall orbits.  The orbit $(23,24)$ has square $8$, so its hyperplane misses
hyperbolic space; the other thirteen orbit sums are the facets.

The exact checks give:

| invariant | value |
|---|---:|
| hyperbolic dimension | 10 |
| facets | 13 |
| finite vertices | 34 |
| ideal vertices / cusps | 4 |
| hyperideal extreme rays | 0 |
| elliptic rank-9 edge diagrams | 208 |

All thirteen inequalities have explicit strict facet witnesses: their
equality rays span rank 10 and a separately recorded exact ray pairs strictly
positively with the wall.  The nonright Coxeter relations are

$$
9\times m=3,\qquad 2\times m=4,\qquad 2\times m=\infty.
$$

There are no ultraparallel facet pairs.  Exact double description gives the
$34+4+0$ ray count above.  The literal Vinberg--CoxIter audit checks all 208
elliptic edges and obtains

$$
(2,0)^{137}(1,1)^{66}(0,2)^5.
$$

Every edge therefore has exactly two finite-or-ideal endpoints.  The four
full parabolic neighborhoods are exactly the four ideal rays, there are no
bad edges, and the chamber is finite-volume and noncompact.

## Integral divisor lattice

Every wall-diagram symmetry lifts from the published curve basis to an
integral unimodular Lorentz isometry of Néron--Severi.  For a transposition,
the saturated fixed lattice has

$$
\operatorname{rank}=11,\qquad
|\det|=512,\qquad
\operatorname{SNF}=1^2 2^9,
$$

with signature $(1,10)$.  All thirteen spacelike orbit normals are integral
reflecting vectors in this saturated lattice, and their coordinate matrix
has Smith diagonal $1^{11}$.  Hence they generate the entire fixed lattice,
not a hidden finite-index sublattice.

Because the full $S_3$ action is geometric, pullback identifies this fixed
nef cone with the divisor cone of the normal projective quotient by a
transposition.  This algebraic realization is already the one recorded in
the frozen quotient-nef theorem; the present result closes its surrounding
subgroup denominator rather than producing a second realization.

## Exact deduplication and bounded lineage

The new certificate compares the primitive facet Gram entry by entry with
all three archived transposition folds.  Explicit facet permutations identify
all three, and the representative agrees exactly with the previously frozen
rank-19 quotient-nef Gram.

Every other positive in the local lineage ledger is excluded by the exact
dimension or facet-count invariant.  The only additional frozen 13-facet
$H^{10}$ comparator is the classical $f_2^{10}$ chamber.  Exact edge-labelled
Coxeter-graph isomorphism rejects it; independently, $f_2^{10}$ has 29 finite
vertices and 3 cusps, versus 34 and 4 here.

The pinned CoxIter snapshot contains four ordinary $H^{10}$/13-facet records,
and its index-two-doubling archive contains four more.  Exact relation-matrix
isomorphism finds no match among all eight.  This is a bounded corpus result,
not a worldwide priority theorem.  In particular it does not overturn the
status of the chamber as an already frozen/published K3 quotient-nef
construction.

## Reproducibility

- exact exhaustion and geometry certificate:
  [`k3_rank18_19_all_subgroup_completion.json`](../data/enumerations/k3_rank18_19_all_subgroup_completion.json)
- certifier:
  [`certify_k3_rank18_19_all_subgroup_completion.py`](../scripts/certify_k3_rank18_19_all_subgroup_completion.py)
- regression tests:
  [`test_k3_rank18_19_all_subgroup_completion.py`](../tests/test_k3_rank18_19_all_subgroup_completion.py)
- earlier direct involution archive:
  [`k3_high_rank_all_involution_folds.json`](../data/candidates/k3_high_rank_all_involution_folds.json)
- published quotient-nef input:
  [`k3_quotient_nef_slices.json`](../data/known_examples/k3_quotient_nef_slices.json)

Replay with

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_k3_rank18_19_all_subgroup_completion.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_k3_rank18_19_all_subgroup_completion.py
```

The closure quantifier is the full exact weighted wall group of each of these
two complete published chambers.  It is not a classification of all rank-18
or rank-19 K3 surfaces or all hyperbolic Coxeter polytopes.
