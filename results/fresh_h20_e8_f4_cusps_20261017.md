# Six explicit walls complete both F4 cusp links

The [exact162-wall chamber](fresh_h20_e8_exact_cluster_20261017.md) has two deficient links `affineF4+E8+A6`. **Six compatible norm2 walls at common height6 complete both to `affineF4+affineE8+affineA7`.** The resulting chamber has168 actual facets, at least seven proper cusps and at least188 finite vertices. It cuts the previously saved defective-cusp escape.

It still has infinite volume. There are no new finite vertices, the target H3 face is unchanged, and a fresh strict spacelike escape is certified. Moreover, an all-real height argument shows that merely lowering these heights, making them unequal, or reversing the missing Euclidean direction cannot make these particular completion profiles meet that H3 face.

The [six-wall source](../scripts/fresh_h20_e8_f4_cusps_20261017.py) and [data](fresh_h20_e8_f4_cusps_20261017.json) pass generation and fresh saved replay. The [height-gate source](../scripts/fresh_h20_e8_f4_cusps_20261017_height_gate.py) and [data](fresh_h20_e8_f4_cusps_20261017_height_gate.json) also pass generation and fresh saved replay. All are final. These are primary certificates, with the exact162 reconstruction and the [conditional completion formulas](fresh_h20_e8_cusp_completion_20261017.md) explicitly inherited.

## Literal six-wall construction

Use the two null vectors nu0 and nu1, their E8 highest root theta, A6 endpoint weights omega_f,omega_l, and transverse vectors z from the conditional certificate. Each z has square8/7 and is orthogonal to the19 incident roots. For each cusp add

`-theta+6nu`, `-omega_f+z+6nu`, `-omega_l-z+6nu`.

The fixed normalization of nu is the positive F4 kernel in the inherited data. Their squares are2. Within each triple, the E8 root is right to the other two, and the A7 endpoints have mutual product−1. All nine cross-cusp pairs are strictly dotted.

The script derives a sufficient common integer bound8 from all old-wall, center and cross-cusp inequalities, then checks every integer height from1 upward. Height6 is the first certified height making every old nonincident and cross-cusp contact strictly dotted while preserving the inherited strict center. Heights1 through5 have saved failing products. This is a statement about common integer heights in that strict-dotted regime, not a minimum over arbitrary real or unequal heights.

All987 new/prior gates are checked:867 dotted,112 right and eight label3. The same strict future center remains inside all168 halfspaces; its orthogonal facet projections certify every wall as actual. Both completed null vectors have exactly22 incident roots. Their positive relations are the original F4 relation giving nu, the E8 relation giving6nu, and the A7 relation giving12nu. The Euclidean component ranks are4,8,7, totaling19, and the normal span is `nu-perp` of dimension20. Any recession direction in the horospherical link pairs nonpositively with every root; the positive relations force every such product to vanish, and the rank condition forces the direction to be zero. Thus each link is bounded. All five previously proper cusp rays have strictly negative products with the six new walls.

No new finite vertex is possible. A finite vertex containing a new root must use only the22 roots of that root's completed cusp diagram, since every other wall is dotted to it. A positive-definite subset must omit at least one root from each of its three affine components, leaving at most19 roots. A finite H20 vertex requires20. Roots from different new triples cannot share a finite vertex because their products are dotted. Old finite vertices persist by the common-center argument already proved for the162 chamber, so the inherited lower bound remains188.

## What was closed and what remains

The second new wall at cusp0 has product approximately `+0.000278925056725` with the old strict escape. Thus the six-wall batch really removes the open region represented by that vector; it is not merely a formal link construction.

At cusp0 the new walls are dotted to T roots3 and10; at cusp1 they are dotted to T roots2 and9. None can meet the hyperbolic T-face, so that entire H3 face and its compact triangular boundary face remain unchanged.

A new exact escape is obtained from rational coefficients in the residual basis formed by the T projections of old roots `[95,142,145,31]`:

`(-2638991277/1250000000, 5464275411/2500000000, -16560862307/10000000000, 14082107/9765625)`.

Subtract `1/64` times the sum of the17 T fundamental weights. The resulting vector has square approximately `1.21031093653`, negative product with the common center, and strictly negative product with every wall. Its rational seed was suggested by a numerical exterior point, but all certificate signs are outward exact interval checks. This proves infinite volume without relying on a numerical cone census.

## Lower heights cannot make these profiles close the H3 face

The separate gate uses only the literal21-root core. For a wall `r=b+hnu`, write J_T for the positive T inverse Gram and d for its T products. Its projection to `T-perp` has square

`q(h)=2-d(h)^T J_T d(h)`.

A wall that meets the hyperbolic T-face must have positive projected square. The gate proves the opposite throughout each admissible real height domain.

For the E8 affine wall `r=-theta+hnu`, coherence forces `h>0`: pairing `theta+r=hnu` with a strict future center has a negative left side and a negative nu product. Its product with old root10 at cusp0, or9 at cusp1, is exactly−h. Those old roots and r have square2, so a nonzero acute Coxeter product has magnitude at least1. Hence `h>=1`. At this lower endpoint, q and its derivative are strictly negative, as is its quadratic coefficient. Therefore q remains negative for every larger h.

For either A6 endpoint wall, every possible quotient direction is `+z` or `-z`; adding a multiple of nu is absorbed into the arbitrary real height. The gate checks both orientations and both endpoints at both cusps: eight exact cases. In each case a specified old core wall forces `h>=h0` by acuteness. The shifted polynomial `q(h0+t)` has all three coefficients strictly negative, so q is negative for every `t>=0`. The data contain the exact algebraic coefficients, the gating old root and outward sign certificates.

This covers every real height, including unequal heights, for the specified norm2 affineE8 plus endpoint-affineA7 completion profiles. It does not classify other affine-link diagrams, changed old-wall sets or the other finite endpoint profiles. The completed cusps provide real local progress, while further H3 closure must use a different incidence pattern rather than another height choice in this family.
