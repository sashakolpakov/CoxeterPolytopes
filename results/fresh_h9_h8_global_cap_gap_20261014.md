# Scalar H8 attachment conditions remain feasible in a new cap direction

An exact all-dotted cap direction satisfies the available projected-norm and diameter necessities for attaching the fixed compact 11-facet H8 source to the mixed-D4 chamber. Its only low projected norm is old wall **8**, with

`kappa_8² = a = (sqrt(5)-1)/4`.

Every other normalized old wall has `kappa² > 4a`. All twelve projected centers have pairwise distance smaller than the source diameter, while the old balls' union has diameter larger than the source. Thus these scalar necessary conditions do not themselves give a universal exclusion over cap directions. **No isometric source placement, full ball coverage, or finite-volume completion is constructed or ruled out for this cap.**

The [exact producer](../scripts/fresh_h9_h8_global_cap_gap_20261014.py) and [certificate](fresh_h9_h8_global_cap_gap_20261014.json) use frozen Fraction field arithmetic and outward-rounded 128-bit intervals. Generation and saved replay pass. The construction is outside the earlier frozen checkpoint.

## Literal cap recipe

Let H be the strict center of the frozen [mixed-D4 chamber](fresh_h9_mixed_d4_coupling_20261014.json), and set

`V=(3.346178,4.2047,3.038572,1.952451,.909731,.150188,-1.149919,-.592165,-.912382,-.750588)`.

All decimals here mean exact rationals. Define `z=H+t(V-H)`, where t is the unique root in

`1.00000001124723266977347773466025 < t < 1.00000001124723266977347773466026`

of `(z,n_8)²=(2+a)z²`. The certificate records the exact quadratic coefficients over the existing quartic field, opposite endpoint signs, and a derivative with constant nonzero sign on the interval. This specifies an exact algebraic cap; the interval is not a numerical acceptance tolerance.

The cap norm is approximately `0.0005128229644`. Its pairing with every old normal and the old strict center is negative, and all twelve old-wall relations are strictly dotted. Its projection of H gives a strict timelike point on the new facet. Thus the cap itself is a valid coherent added wall.

The other projected norm squares are all greater than `sqrt(5)-1`; the smallest is approximately `1.67`. The center diameter is approximately `2.34`, and the ball pair indexed `(7,8)` has union diameter approximately `4.90`. These decimals are explanatory only: all accepted comparisons use exact algebraic cosh bounds. The source diameter is reconstructed from all 41 vertices and 820 pairs, with unique extremal pair `(32,33)` and exact cosh `17+8sqrt(5)`.

## A general necessary restriction on covering pairs

For normalized old normals of square two, let `kappa_i²` be the negative squared norm of their projections to the cap hyperplane. Write `c_ij=-(n_i,n_j)`. The ball radius is `asinh(sqrt(2)/kappa_i)`, and the center-distance cosh is

`(sqrt((kappa_i²+2)(kappa_j²+2))+c_ij)/(kappa_i*kappa_j)`.

Both radii and the center distance decrease as either projected norm increases. For any finite-angle old pair, `c_ij<=phi`. If both projected norms are at least `sqrt(5)-1`, the resulting ball-union diameter is at most approximately `4.166966`, strictly less than the source diameter `4.245100`. The exact checker verifies this inequality without inverse hyperbolic functions.

Consequently, if the source's diameter endpoints are covered by two such old balls, their pair must either involve a low projected center or correspond to one of the three old dotted pairs:

`(4,11), (5,6), (8,11)`.

The [separate two-anchor gate](fresh_h9_h8_corner_anchors_20261014.md) sharpens the low-center count to at most one. The cap above satisfies that stronger count as well. Because its sole low norm is exactly a, any compatible placement must put its normalized wall-8 projected center on source vertex 0, 32, or 33 with minimum contact scale one. This gives three concrete corner-anchored placement problems, rather than another unrestricted cap-direction search.

Remaining work is to orient the source around one of those three anchors, impose every old/source sidewall Coxeter relation, and check the entire required coverage. Passing the present scalar relaxations is not evidence that those further conditions will pass.
