# Characteristic-3 Fermat quartic: exact Frobenius-fold no-go in H^11

## Result

The arithmetic Frobenius-fixed slice of the complete Kondo--Shimada
induced chamber is **not a Coxeter polyhedron**.  This is a local and exact
obstruction: two walls individually proved to be genuine facets have
normalized squared product `3/8`, which is not `cos(pi/m)^2` for any integer
`m >= 2`.  There are 864 such forbidden certified-facet pairs in total.
Consequently no Vinberg finite-volume test is appropriate for this slice.

This is a negative result for one particularly canonical high-rank cast,
not a negative result for the supersingular K3 stock or for all subgroups of
its chamber stabilizer.

## Exact input and fixed space

The surface is the Fermat quartic in characteristic 3, a supersingular K3
surface of Artin invariant 1.  Its Neron--Severi lattice in the published
basis has rank 22, determinant `-9`, and signature `(1,21)`.  The chamber
has the complete `112 + 648 + 5184 = 5944` wall list from the assignments
`FQWW112`, `FQWW648`, and `FQWW5184`.

The published Frobenius matrix `FQAF` satisfies, exactly,

- `FQAF^2 = I`;
- `FQAF * FQGramS * FQAF^T = FQGramS`;
- `FQh0 * FQAF = FQh0`; and
- `rank ker(FQAF-I) = 12`.

Thus its fixed positive cone is hyperbolic 11-space.

## Restricted walls

The contragredient action of Frobenius was derived directly on all 5944
published wall covectors.  Its orbit census is:

| wall stock | fixed walls | two-cycles | wall orbits |
|---|---:|---:|---:|
| 112 | 8 | 52 | 60 |
| 648 | 0 | 324 | 324 |
| 5184 | 48 | 2568 | 2616 |
| total | 56 | 2944 | 3000 |

Orbit sums restrict the inequalities to the fixed space.  After identifying
positive proportional sums, the 3000 orbits give 2316 distinct wall rays;
684 coincidences are recorded rather than silently counted as facets.  The
rays span rank 12, and their exact Gram inertia is `(1,11,2304)`.

## Genuine-facet certificate

For a restricted normal `n`, let

`p_n = h0 - (<h0,n>/<n,n>) n`.

Whenever `p_n^2 > 0`, `<p_n,n> = 0`, and `<p_n,n'> > 0` for every other
distinct restricted wall ray `n'`, the point `p_n` lies in the relative
interior of that wall in the hyperbolic chamber.  Therefore the wall is a
genuine facet.  Exact integer arithmetic certifies this condition for 384
distinct rays.  This number is an exact certified lower bound, not a claim
that the remaining rays have all been globally reduced.

One explicit forbidden pair is formed by restricted rays 1 and 128.  Their
primitive covectors are

```
(-2,2,1,1,2,-2,1,1,1,0,1,0,0,0,1,0,0,0,2,1,0,1)
(1,2,0,1,2,1,1,0,1,1,0,0,2,2,2,0,1,1,0,2,0,2).
```

Using `FQGramS^-1`, their squares are `-4` and `-8/3`, and their product is
`2`.  Hence

`2^2 / ((-4)(-8/3)) = 3/8`.

Their strict relative-interior witnesses have squares `5` and `10`,
respectively.  Across all 384 certified facets, the forbidden ratios are:

- `3/8`: 576 pairs;
- `1/16`: 288 pairs.

Both are rational values strictly between 0 and 1 but outside the Coxeter
list `0, 1/4, 1/2, 3/4, 1` for rational squared cosines of submultiples of
pi.  This proves the no-go without assuming that the 384 facets are the
complete facet set.

## Scope and replay

The computation covers exactly the Frobenius involution `FQAF`.  It does
not claim a census of conjugacy classes in the full chamber stabilizer, and
surface-automorphism realization is not used as a gate.  No finite-volume
claim is made because the Coxeter-angle gate already fails.

- Script: `scripts/search_ss_fermat_quartic_char3_folds.py`
- Exact JSON: `data/enumerations/ss_fermat_quartic_char3_frobenius_fold.json`
- Tests: `tests/test_ss_fermat_quartic_char3_frobenius_fold.py`
- Source provenance: `data/external/kondo_shimada_fermat_quartic_char3/README.md`
