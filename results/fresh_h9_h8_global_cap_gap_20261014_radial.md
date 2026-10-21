# Uniform radial reduction for H8 attachments

For every all-dotted cap direction, the fixed compact H8 source can be covered by the twelve old-wall balls only in the following branches:

* all twelve normalized projected norms satisfy `kappa²>=sqrt(5)-1`; or
* the unique lower norm belongs to old polar wall **11**, its center is source vertex **32 or 33**, and its exceptional-ray contact scale is **1 or sqrt(2)**.

At scale sqrt(2), ball 4 must cover the opposite E8 vertex and ball 8 must cover source vertex 0. At scale 1, both assignments of these two balls remain possible. No existence or exclusion is asserted for these remaining branches.

The [producer](../scripts/fresh_h9_h8_global_cap_gap_20261014_radial.py) and [exact certificate](fresh_h9_h8_global_cap_gap_20261014_radial.json) verify the old pair graph, source-corner distances and every comparison with outward-rounded 128-bit algebraic intervals. Generation and saved replay pass. The exact field and interval implementations are reused from frozen certificates. This is a primary proof certificate, not a separate independent arithmetic audit.

The covering condition is necessary for the raw attachment formed from the retained twelve old walls, the cap, and the fixed source's perpendicular sidewalls. Further nonperpendicular walls could add covering balls; they are outside this radial-cover conclusion. The result is not a nonexistence theorem for arbitrary completions retaining an H8 facet.

## Exact constants and monotonicity

Normalize each old normal to square two, and put `a=(sqrt(5)-1)/4`, `B=4a`. The [source gap lemma](fresh_h9_h8_facet_cover_20261014.md) and [two-anchor gate](fresh_h9_h8_corner_anchors_20261014.md) imply that every projected norm is at least a and at most one is below B. A lower-norm center is one of source vertices 0, 32, 33, scaled by a contact s with `kappa²=a*s²`. Since `1<=s<2`, the first allowed values are `1,sqrt(2),phi`, followed by larger finite-angle magnitudes.

The certificate reconstructs all three source vectors, normalizing their minimum nonzero contact to one. Their squared norms are -a. Source distances obey

`cosh D_A=9+4sqrt(5)` for pairs `(0,32),(0,33)`, and
`cosh D_B=17+8sqrt(5)` for pair `(32,33)`.

Thus `D_A≈3.57964082` and `D_B≈4.24510025`.

Let x and y be the projected squared norms of two centers and C the negative product of their normalized old roots. Then

`cosh d = sqrt(1+2/x)*sqrt(1+2/y) + C/sqrt(x*y)`,

and the second ball's radius is `asinh(sqrt(2/y))`. All terms decrease as x or y increases because C is nonnegative. Consequently the greatest possible radial reach from a low center to a high ball is obtained at the minimum permitted norms.

For `x>=a`, `y>=B` and `C<=1`, that reach is at most approximately **3.54606370**, strictly below D_A. The low center's own ball has radius at most approximately **1.66346**, also below D_A. A high ball has diameter at most approximately **2.12255012**, strictly below D_A.

## Two remote vertices eliminate low walls 0–10

The exact old-root graph of pairs with `C>1` has only the four edges

`(0,1), (4,11), (5,6), (8,11)`.

The two exceptional source vertices other than the low center are each at distance at least D_A from it, and their mutual distance is also at least D_A. The own ball and every `C<=1` ball miss each remote vertex by the preceding radial bound. If the low wall is any index 0–10, it has at most one neighbor with `C>1`. That single high ball cannot cover both remote vertices because its diameter is smaller than D_A. This excludes every such low anchor, for every allowed scale and every cap direction.

## Low wall 11

Only balls 4 and 8 can cover the two remote exceptional vertices. Their old roots are right, and their projected norms are at least B. Their total union diameter is strictly below D_B. Therefore wall 11 cannot be anchored at source vertex 0, whose two remote vertices are 32 and 33, separated by D_B.

Suppose instead it is anchored at vertex 32 or 33. If `s>=phi`, even the larger old product `C_(4,11)`, with the other projected norm at B, gives radial reach at most approximately **4.22289563**, below D_B. Thus every scale at least phi is excluded. The contact alphabet then leaves only s=1 or s=sqrt(2).

At s=sqrt(2), ball 8's maximum radial reach is approximately **3.92485**, below D_B. Hence ball 4 must cover the opposite E8 vertex; since no single high ball can cover both remote vertices, ball 8 must cover vertex 0. At s=1, the certificate leaves both assignments open.

This reduces the low-center case to four corner/scale branches and six possible two-ball assignments. The no-low-center branch also remains open. The bound is uniform over cap directions and has no label cutoff or numerical-search dependence.
