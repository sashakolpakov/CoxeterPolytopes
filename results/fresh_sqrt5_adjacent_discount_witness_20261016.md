Class 4202 is a real local arithmetic quadrilateral, and the preserved candidate does not pay it. For retained facet 4 and common spectator forest (I_2(5)^2+A_{11}), its payment is exactly

\[
1-\delta_{\mathrm{adj}}+w_{14612}+w_{24603}.
\]

Consequently this face requires (w_{14612}+w_{24603}\geq\delta_{\mathrm{adj}}). Both bonuses are zero in the preserved numerical candidate. Interpreting its saved decimal coefficients as exact rationals gives payment (1820357345943253/5000000000000000\), about 0.36407, below the required one. This is a reproducible candidate failure, not an assertion that the saved numerical vector was an exact feasible primal.

The core has eight normals. Its nonright bonds are (05=13=5), (23=25=37=45=3), and (67=4). Scale root 6 to square four and every other root to square two. The Gram entries then lie in \(\mathbb Q(\sqrt5)\): the 5-bonds have product \(-\phi\), the 3-bonds product \(-1\), and the 4-bond product \(-2\). The full core Gram has inertia \((6,1,1)\), while its conjugate has inertia \((7,0,1)\). Quotienting its one-dimensional relation yields an admissible Lorentz form of dimension seven, hence a local chamber in \(H^6\).

The certificate supplies literal quotient coordinates and a strict future center of square \((225-113\sqrt5)/2<0\). White roots \(4,5,6,7\) are positive definite. Each of the four adjacent black pairs gives a future finite intersection point satisfying all remaining inequalities strictly. An exhaustive six-pair cone check finds precisely those four feasible rays and no others. Thus the entire white-supported two-face is a compact quadrilateral.

Adding orthogonal spherical blocks \(I_2(5)^2+A_{11}\) supplies 15 more white roots and embeds this face into a local chamber in \(H^{21}\), with retained facet 4. A literal strict center remains timelike; the four face vertices have exactly 21 active normals. The 23-normal Gram has physical inertia \((21,1,1)\) and conjugate inertia \((22,0,1)\). All 529 Cartan coefficients are in \(\mathbb Z[\phi]\), so the reflections preserve their admissible integral root module.

The unique baseline-one corner is the adjacent black pair in B6. The other corners are \(H_3+A_1+I_2(4)\), \(H_3+B_3\), and \(A_3+A_1+I_2(4)\), before adding common spectators. The current H3-center bonus misses the H3–B3 pair: its H3 black root is the 5-edge leaf, not the center. A possible additional channel is that leaf paired with the B3 3-edge endpoint; its budget needs a new complete check. No repair or improved bound is asserted here.

The [script](../scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py) and [result](fresh_sqrt5_adjacent_discount_witness_20261016.json) passed generation and a fresh saved replay. All geometry uses a standalone Fraction implementation of \(\mathbb Q(\sqrt5)\), separate from the catalogue arithmetic. Payment compilation uses the frozen proposed-rule graph helper and the preserved coefficient snapshot. This establishes an actual compact face of a local arithmetic chamber; it does **not** establish a compact or finite-volume ambient polytope, nor prove that the face must occur in one.
