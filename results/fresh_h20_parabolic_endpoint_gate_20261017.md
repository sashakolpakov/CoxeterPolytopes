# Continuous heights for the parabolic endpoint profiles

The five parabolic profiles whose walls can meet the target H3 space cannot supply neighboring facets to one another in the stated endpoint family. This holds at every allowed height: all seven cross-edge pairs are strictly dotted, and the three pairs on the same edge are obtuse. This package does not yet exclude mixed finite/parabolic endpoint pairs.

[Source](../scripts/fresh_h20_parabolic_endpoint_gate_20261017.py) and [certificate](fresh_h20_parabolic_endpoint_gate_20261017.json) passed generation and a fresh saved replay. They use the complete [20-profile census](fresh_h20_parabolic_endpoint_profiles_20261017.md), the [exact 162-wall cluster](fresh_h20_e8_exact_cluster_20261017.md), and the [168-wall completion](fresh_h20_e8_f4_cusps_20261017.md). All statements concern coherent acute extensions retaining those 168 walls, with every unprescribed old contact other than 158 strictly dotted. The pairing with 158 need only be nonpositive for the bounds below; thus the height intervals are continuous relaxations of the Coxeter domain.

For an edge support W, let J be its spherical inverse Gram, p=-WJd, and u the projection of wall 158 perpendicular to W. Set sigma=u^2>0, e=u/sqrt(sigma), and let V be the future unit finite vertex perpendicular to W and 158. The two future null rays in W-perp can be normalized as (V-e)/sqrt(sigma) and (V+e)/sqrt(sigma). The first, denoted nu, pairs -1 with 158; the second pairs +1 and is excluded by the old halfspace. Exact coefficient identities check these assertions. Outward rational intervals certify that nu pairs strictly negatively with every old168 wall outside W, and negatively with the inherited strict future center. It is therefore the actual unbounded end of the old edge, with exactly W incident among the retained walls.

A parabolic profile has p^2=2. Every coherent nondegenerate extension with this profile is r=p+h nu. Its positive affine relation is r+WJd=h nu. Pairing that relation with a common strict future center gives h>0; h=0 would be a positive linear dependence among normals. Write x=-r.158, off=p.158, so h=x+off and r(x)=p+off nu+x nu. The new root has square2 for every x. The unwanted null ray cannot occur because a positive affine relation would force that ray to be an old chamber ideal point, whereas it violates wall158.

Only the five profiles with positive projected norm can meet the hyperbolic T-space. For them the certificate proves the following necessary lower bounds:

| Edge/profile | Continuous lower bound on x | Old dotted gate supplying it |
| --- | --- | --- |
| 159,160 / affine A12 | 29/20 | 145 |
| 159,161 / affine C13 | 0 | acuteness with158 |
| 159,161 / affine C12 | 31/100 | 160 |
| 159,161 / affine B12 | 31/100 | 160 |
| 160,161 / affine A12 | 81/100 | 159 |

For each nonzero bound L, the selected old product has negative slope in x and is strictly greater than the dotted threshold at x=L. Hence every admissible cap satisfies x>L. Using the closed relaxed interval x>=L is safe. The coefficient vectors and exact interval witnesses are saved, and no assumed fixed center is used to reject a height.

For two different edges, put x=L+t and y=M+u. Each product r(x).s(y)+2 is bilinear in nonnegative t,u. All four coefficients are strictly negative, certified by outward rational intervals, for all seven pairs. The three profiles on the middle edge have constant mutual products sqrt(2), 1, sqrt(2), independent of their heights, so no two are acute. Two different heights of the same profile have product2 and also cannot coexist as distinct acute walls. The fifteen zero-projected-square profiles do not meet T and cannot add target-face sides.

Dottedness persists after perpendicular projection to T: the inverse spherical Gram is entrywise nonnegative, so subtracting the support correction makes the mutual product at least as negative and cannot increase the diagonal squares. Thus these parabolic caps do not meet each other inside the target hyperbolic space. This is a restriction on this specific endpoint family, not a general obstruction to ideal vertices or cusp completion in H20. Additional finite old contacts, different retained walls, and mixed finite/parabolic neighbors are separate questions.
