# Independent hostile audit: global silver Lambda22 nonreflectivity

Terminal verdict: **PASS_GLOBAL_NONREFLECTIVITY**.

For

`M=<-(sqrt(2)-1)> + Lambda22` over `Z[sqrt(2)]`,

the full integral reflection subgroup has infinite index and its complete
reflection chamber has infinite hyperbolic volume.  This is a global theorem
for the displayed host, not a bounded-prefix conclusion.

## Independent arithmetic

The replay imports no production silver module.  It reparses the rank-22
Gram, checks the explicit map on all 23 basis vectors, reconstructs its
integral inverse, and verifies the full form identity.  Its Lorentz-plane
block has determinant one and trace `10+8sqrt(2)>2`, so the map has infinite
order and preserves the chosen time component.

The eight complete mirror-norm classes are replayed separately at time
coordinate zero.  From `y=u+sqrt(2)v` one has
`u^2+2v^2=Tr(q)/2`; the ordinary/S2/S3/S6 minima exclude all eight cases.
Thus the basepoint is regular.

## Separator bound

If an oriented reflective root separates `v0` from `g*v0`, endpoint signs
give

`(d,y) >= ((4-sqrt(2))/2)x`.

Cauchy then gives

`x^2/q <= 136+96sqrt(2)=(8+6sqrt(2))^2`,

and hence `H<=56+40sqrt(2)<150`.  The same endpoint inequality forces

`(z,r)>=(10-7sqrt(2))x/6>0`, for `z=(1,d/3)`.

The H150 frontier is not merely trusted by hash.  This script independently
enumerates all 58 arithmetic controls in a proved Minkowski coefficient box,
all even binary Gram candidates, their exact binary minimum-four test and
pairing gates, and every gated PSD extension by `d`.  It reproduces 18
retained levels, 72 necessary types, and exactly seven post-H14 nonnegative
types.  Clean-room modular kernels and PARI shell enumeration prove minima
`min(S2,S3,S6)=(12,24,96)`, the critical S2 norm-12 pairing spectrum
`{0,2,4}`, and the S2 mod-8 residue classes.  These data eliminate all
seven threats; in particular, the complete nonnegative extension sets in the
two shell-pairing cases are singletons `(18,12)` and `(22,14)`.

Thus every actual root through H150 pairs strictly negatively with `z`, and
no mirror separates the regular endpoints.  Since the integral isometry
permutes all integral-reflection mirrors, it stabilizes their common open
chamber.

The restriction is a hyperbolic translation of positive length.  A small
ball in the chamber has infinitely many pairwise disjoint translates inside
that chamber, proving infinite volume.  Since the ambient admissible
arithmetic orthogonal group has finite covolume, the reflection subgroup has
infinite index.  The conjugate form is definite, so the host is anisotropic
and its ambient arithmetic quotient is cocompact.

The production artifact compared after the independent derivation is
`data/enumerations/silver_lambda22_global_nonreflectivity.json`, SHA-256 `db24c44646bb4e46b31a0d9c168635b752d4d6a133f845e6c95dffca3259accf`.

This does **not** prove field-wide nonexistence over `Q(sqrt(2))` or exclude
other `H22` hosts.
