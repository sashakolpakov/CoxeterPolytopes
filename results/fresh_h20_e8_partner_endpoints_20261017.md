# Endpoint profiles of the E8 triangular branch

The three new open edges of the proposed162-wall E8 branch have **159 fixed finite endpoint profiles**, plus one explicitly unbounded dihedral profile per edge. The finite-diagram census is exact. The ambient branch and the parameter-interval computations remain numerical; this package does not certify a new exact wall or finite-volume chamber.

All four source/data pairs below pass generation and fresh saved replay. They are final for this package:

| Calculation | Source | Data |
|---|---|---|
| Complete finite endpoint profiles and conic intervals | [source](../scripts/fresh_h20_e8_partner_endpoints_20261017.py) | [data](fresh_h20_e8_partner_endpoints_20261017.json) |
| Coincidence with an old facet | [source](../scripts/fresh_h20_e8_partner_endpoints_20261017_bridges.py) | [data](fresh_h20_e8_partner_endpoints_20261017_bridges.json) |
| Coordinated partner rotation | [source](../scripts/fresh_h20_e8_partner_endpoints_20261017_coordinated.py) | [data](fresh_h20_e8_partner_endpoints_20261017_coordinated.json) |
| Coordinated E8 boost and rotation | [source](../scripts/fresh_h20_e8_partner_endpoints_20261017_boost.py) | [data](fresh_h20_e8_partner_endpoints_20261017_boost.json) |

The input is the frozen [E8 triangular-cluster recipe](fresh_h20_e8_triangle_cluster_20261016.md). Its target support is

`T=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87,151]`, of type `A10+A4+A2+A1`.

The numerical target H3 has24 rays,10 finite and14 exterior, with maximum squared chart radius2.22840934760. Its eight open finite edges are the old five `(95,142),(95,157),(142,152),(145,152),(145,157)` and the three new partner edges `(159,160),(159,161),(160,161)`. No full ambient vertex census is claimed.

## Complete finite endpoint domains

For an actual finite endpoint, `T` plus the two edge walls plus the new wall must be spherical of rank20. A connected component containing the new wall contains whole old irreducible components: deleting the new vertex cannot cut an old connected component. Thus it suffices to delete each vertex from every finite irreducible diagram of rank at most20, match the remainder to a union of whole old components, and retain all labelled isomorphisms. The script includes A, B, D, E, F, H and the disconnected new A1 case. All roots use square2; labels therefore determine the normalized Gram, including the B/F bonds without a root-length ambiguity.

| Edge | Old rank19 support | Fixed profiles |
|---|---|---:|
|159,160|A12+A4+A2+A1|57|
|159,161|B12+A4+A2+A1|45|
|160,161|A12+A4+A2+A1|57|

The sum is159. Each edge also permits a new root with its only nonright contact at the old isolated A1 root151, of arbitrary label `m>=3`. Its projection cost is `2 cos²(pi/m)` and residual square `2 sin²(pi/m)`. This single family incorporates all I2 labels, including small labels omitted from the fixed list. No other unbounded label family is possible in a finite connected diagram of rank at least3. Every fixed profile saves its exact support Gram, inverse, contact cost, positive residual norm, and deletion/isomorphism witnesses.

## Analytic conics, numerical inequalities

After removing the prescribed positive T projection, each fixed endpoint problem has two linear constraints in a Lorentz4 space. Its square equation is a hyperbola. On each branch the residual normal has the form

`r(k)=p+k a+k^(-1)b`, with `k>0`.

Every old full product includes its T-projection offset and is consequently `A_j k+B_j/k+C_j`. Requiring every unprescribed product to be dotted is a quadratic inequality in k. The script intersects all such intervals and a sufficient common-center inequality. Exactly one branch of each of the159 fixed profiles has a nonempty interval. Three representative points per interval are saved; those samples do not exhaust all geometry scores. Every sampled cap leaves the same maximum exterior radius, so none removes the dominant open region.

Allowing one designated actual old H3 facet to have a nondotted product, while every other unprescribed product stays dotted, gives no parameter interval for any fixed profile. This is only a numerical interval result in that restricted contact regime. If all remaining gates are relaxed to acuteness, **54 intervals survive** with an intrinsically positive-definite new/old pair. These open relaxations explicitly prevent a general obstruction claim. The three unbounded151 families and simultaneous finite-angle coincidences remain outside this computation.

The simple polar of the three partner planes also cannot close the region while retaining the E8 wall158. Their residual Gram is Lorentzian of rank3, and its orthogonal line is precisely the positive line spanned by the158 projection. A positive multiple is obtuse to158. A negative multiple imposes the opposite halfspace on T, destroying the interior of that face. This is a statement about this polar construction, not about other endpoint profiles.

## Coordinated alternatives

The rotation experiment keeps the E8 wall fixed, rotates its B3 partner triple, and asks a new cap to meet both an old open edge and one partner edge. The two old finite labels range over2 through5, subject to a positive-definite corner. There are24 old-edge patterns and ten endpoint profiles right to T. Eliminating the rotation gives480 quartic equations. All960 numerical positive real solutions are retained and tested; none is compatible. This is an alternative branch, not an extension retaining the original partner positions.

The second experiment permits the E8 boost to vary while retaining its original ideal direction. It uses all551 old-compatible caps in the frozen attachment scout, including nonright T profiles. The partner-angle equations give a quadratic in the cap's product with the moving timelike center, followed by a quadratic for the boost. There are4,324 matching profile problems and8,648 positive solutions, with no accepted cluster. The largest solved boost is10.8801077228, whereas the numerical old-acuteness interval for this E8 direction begins at21.3726111269; requiring all old nonsupport contacts dotted raises the threshold to63.5515144250. Thus even acuteness fails before the remaining gates for every solved branch.

These equations avoid a sampled rotation grid, but their coefficients and root isolation are numerical. The cap bank, finite label range and fixed ideal direction remain restrictions. They prove no all-parameter exclusion. The useful outputs are the complete finite endpoint diagrams, explicit continuous families, and a well-defined failure of these coordinated alternatives.
