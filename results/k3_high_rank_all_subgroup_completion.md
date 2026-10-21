# Full subgroup closure for the rank-15--17 K3 wall stocks

## Result

The full weighted-diagram subgroup search is now closed for the three
published complete K3 nef chambers of Picard ranks 15, 16, and 17.

| Picard rank | complete walls | exact weighted group | element classes | subgroup classes | proper classes of fixed rank at least 11 |
|---:|---:|---|---:|---:|---:|
| 15 | 21 | $S_3\times C_2$, order 12 | 6 | 10 | 3 |
| 16 | 19 | $C_2\times C_2$, order 4 | 4 | 5 | 1 |
| 17 | 19 | $C_2$, order 2 | 2 | 2 | 0 |

Thus all 12 element conjugacy classes and all 17 subgroup conjugacy classes
have been exhausted.  Exactly four proper subgroup classes survive the fixed
rank gate:

| source | subgroup class | fixed space | complete fixed cone | disposition |
|---|---|---:|---:|---|
| rank 15 | $C_2$, class size 3 | $H^{13}$ | 19 facets | already certified involution fold |
| rank 15 | $C_3$, class size 1 | $H^{12}$ | 17 facets | already certified triality fold |
| rank 15 | $S_3$, class size 1 | $H^{12}$ | 17 facets | new subgroup class, but exactly the same cone as the $C_3$ fold |
| rank 16 | $C_2$, class size 1 | $H^{14}$ | 17 facets | already certified involution fold |

There is therefore **no new hyperbolic Coxeter polytope** in the newly opened
part of this denominator.  The only new high subgroup is the rank-15 $S_3$,
and it does not define a new fixed cone.

## Exact $C_3=S_3$ fixed-cone theorem

Let $C_3$ be the triality subgroup and let $S_3$ be the high fixed-rank
subgroup of the exact rank-15 weighted wall group.  The computation proves

$$
C_3\subset S_3,
$$

and both groups have precisely the same 17 wall orbits:

$$
1,2,\ldots,13,
(14\ 15\ 16),(17\ 18\ 19),20,21.
$$

Their orbit-sum Gram matrices agree entry by entry and have rank 13.  Since
the 21 published walls are the complete nef-wall stock and span the divisor
space, equality of these orbit partitions gives equality of the invariant
wall span, invariant half-space system, and hence the complete fixed cone.
The extra transpositions in $S_3$ impose no additional wall orbit and no
additional fixed direction.

Consequently the existing strict-facet and volume certificate transfers
without rerunning a weaker shortcut:

- all 17 inequalities are genuine facets with strict relative-interior
  witnesses;
- the orbit Gram has Coxeter relation counts
  $116\times m=2$, $12\times m=3$, $2\times m=6$, three parallel pairs,
  and three ultraparallel pairs;
- exact cone enumeration has 93 finite, 3 ideal, and 0 hyperideal extreme
  rays;
- the literal Vinberg/CoxIter audit checks all 623 elliptic rank-11 edge
  diagrams, obtains endpoint histogram
  $(2,0)^{496}(1,1)^{124}(0,2)^3$, and finds no bad edge.

The chamber is finite-volume and noncompact.  Its abstract lineage is already
closed: it is Scharlau--Kirschmer `List12/Gram12` entry 6, with the displayed
17 walls forming the full reflection chamber (index one).

## What happened to every other subgroup

In rank 15, the remaining subgroup fixed dimensions are at most seven:
the central $C_2$, the other $C_2$ class, $C_2^2$, $C_6$, the second $S_3$,
and the full $S_3\times C_2$ all fail the rank-above-nine target before any
expensive geometry is needed.  In particular, the cyclic order-six class has
fixed hyperbolic dimension six; its nine wall orbits include an isotropic and
a timelike orbit sum, in agreement with the earlier independent audit.

In rank 16, two of the three $C_2$ classes have fixed hyperbolic dimensions
seven and eight, and the full Klein four group has dimension seven.  In rank
17, the sole nonidentity subgroup has fixed hyperbolic dimension eight.

The identity classes reproduce the ambient published chambers and are not
counted as proper casts.

## Algebraic-realization boundary

Every weighted-diagram permutation was lifted in the leading published curve
basis and checked to be an integral unimodular Lorentz isometry of the
Néron--Severi lattice.  That is enough to define the divisor-wall fixed cone;
it is not, by itself, a Torelli certificate for an automorphism of a selected
K3 surface.

For the general rank-15 surface, Kondō realizes only the central $C_2$ graph
action on Néron--Severi.  Neither the $C_3$ triality nor the high $S_3$ is
asserted here to be a surface-automorphism action.  Special-period realization
remains unassessed and is unnecessary for the abstract Coxeter-polytope
statement.  The two involution realization statuses are merely referenced
from their frozen archive and are not reclassified here.

## Reproducibility

- exact certificate:
  [`k3_high_rank_all_subgroup_completion.json`](../data/enumerations/k3_high_rank_all_subgroup_completion.json)
- certifier:
  [`certify_k3_high_rank_all_subgroup_completion.py`](../scripts/certify_k3_high_rank_all_subgroup_completion.py)
- tests:
  [`test_k3_high_rank_all_subgroup_completion.py`](../tests/test_k3_high_rank_all_subgroup_completion.py)
- referenced involution certificate:
  [`k3_high_rank_involution_folds_rho15_17.json`](../data/enumerations/k3_high_rank_involution_folds_rho15_17.json)
- referenced triality certificate:
  [`k3_rank15_order3_h12_fold.json`](../data/enumerations/k3_rank15_order3_h12_fold.json)

Replay with

```bash
PYTHONPATH=src:. python3 \
  scripts/certify_k3_high_rank_all_subgroup_completion.py --check-archive
PYTHONPATH=src:. pytest -q \
  tests/test_k3_high_rank_all_subgroup_completion.py
```

The quantifier is deliberately bounded: this closes every subgroup of the
full exact weighted wall groups of these three complete published chambers.
It is not a classification of all K3 divisor lattices, all surface
automorphisms, or all hyperbolic Coxeter polytopes.
