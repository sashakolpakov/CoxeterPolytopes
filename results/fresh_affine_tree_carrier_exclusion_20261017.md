# Large affine-tree carriers cannot occur in an unpaid finite quadrilateral star

**Final internally reviewed all-rank lemma.** Let an actual quadrilateral two-face have four finite vertices and zero payment at all four corners under Vinberg's original finite weights. Suppose its ambient acute Coxeter diagram has no connected affine subdiagram of order at most17. In the private/shared decomposition, neither H−U nor H−V contains a connected affine B/C/D diagram. Together with the [cycle-carrier exclusion](fresh_affine_cycle_carrier_exclusion_20261017.md), this excludes every affine subdiagram of these two complements. This is a geometric transfer ingredient, not by itself a dimension bound.

The [checker](../scripts/fresh_affine_tree_carrier_exclusion_20261017.py) and [data](fresh_affine_tree_carrier_exclusion_20261017.json) verify the complete finite resistance-role reduction and the two continuation obstructions. The resistance range is derived, not a rank cutoff. All ranks and arbitrary remaining supporting whites are handled in the proof below.

## Actual roles and the first contact with an affine carrier

Use the [finite-quadrilateral audit's](fresh_simple_finite_quad_vinberg_transfer_20261017.md) notation: W is the supporting white diagram, the opposite black pairs are {u1,u2} and {v1,v2}, and each W+ui+vj is positive definite. All four blacks lie in the same connected component H. The private diagram U contains the u-roots and private u-white components, with V defined symmetrically. Zero payment forces ui.vj=0. These actual-face and private/shared facts do not require simple ideal vertices elsewhere.

Suppose A is an affine tree contained in H−U. It contains both v-roots: a nonspherical diagram omitting either opposite black pair would lie in an actual finite corner. All its other roots are white. Its order m is at least18, so it is a classical B/C/D affine tree. The argument for H−V is symmetric.

For every white root z outside A, and also for z=u1,u2, both (A−v1)+z and (A−v2)+z are spherical. We first classify the relevant contacts with A, before using any shared-path continuation. All roots have norm2.

There are at most two contacts. For any two A-neighbors x,y of z, the unique x-to-y path in A together with z is a graph cycle. Both marked v-roots must lie on that path, or the cycle would survive a spherical deletion. Three contacts are impossible: their three pairwise tree paths cannot all contain the same two distinct nodes. Equivalently, removing the edge segment between the two marks separates opposite sides; three contact nodes cannot be pairwise on opposite sides. Thus a two-contact profile has both marks on the x-to-y path.

We will apply the following label reduction to an opposite black root z=ui, or to a white z reached by a white path from a u-root. Any nonzero contact to A then has label3 or4. To prove this, let x be the contacted A-node. At least one of the two marked deletions leaves x connected within A to the retained v-root: on a tree, choose the deletion on the other side of the attachment's projection to the v1-to-v2 path. The connected finite corner component containing the high bond then also contains the retained v and the u-root, either directly if z=ui or along the stated white path. A label5 forces that component to be H3/H4 and gives positive Vinberg payment. A label at least6 cannot support the continuation beyond an isolated finite rank-two component; in the direct marked-pair case it would itself give positive payment. The [finite-weight interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md) records these published finite classification/payment inputs.

A sole label3 or4 contact is impossible by the all-rank small-affine leaf witnesses in the [private-tree lemma](fresh_private_affine_tree_exclusion_20261017.md). Those witnesses concern A+z alone and do not assume A is private. Consequently every relevant nonzero profile has exactly two contacts, each of label3 or4.

## Complete two-contact classification from grounded resistance

Normalize the positive kernel of A by nu=1 on main nodes, nu=1/sqrt2 at a doubled terminal and nu=1/2 at a fork leaf. The scaled Gram is the tree Laplacian with resistance1 on ordinary/doubled edges and resistance2 on fork-leaf edges. Let x,y be the two contact nodes, let L be their total resistance distance, and write t_x=−z.x, t_y=−z.y, so t_x^2,t_y^2 belong to {1,2}.

For a marked root j on the path, put s_j=R(x,j). Grounding at j separates the x- and y-sides, so the cross Green entry is zero. At an endpoint the corresponding contact is deleted. The exact Schur cost is therefore

    q_j = a s_j + b(L−s_j),
    a=t_x^2 nu_x^2,  b=t_y^2 nu_y^2,
    q_j<2 for both distinct marked roots.

Every coefficient is at least1/4. Hence L<8. Opposite end decorations have resistance m−1≥17 and cannot occur. Up to exchanging endpoints, only the following roles remain:

| Contact roles | Path positions available to the marks | Complete surviving cases |
|---|---|---|
| main/main | 0,1,...,L | L=1, both contacts3; affine A2 |
| fork leaf/main | 0,2,3,...,L | Both contacts3, L=3 or4, marks L−1,L; affine A3 or A4 |
| sibling fork leaves | 0,2,4 | Both contacts3 give affine A3; exactly one contact4 forces its leaf and the fork center to be the marks |
| doubled terminal/main | 0,1,...,L | L=1 gives the terminal triangle, both contact nodes marked; L=2 forces both contacts3 and marks1,2 |

Here each displayed small affine cycle is the induced old contact path plus z. If both sibling contacts are4 there is no strict solution. In the fork/main row the main contact must be3: a doubled main contact leaves at most one admissible marked position. With a simple main contact, fork contact3 gives only L=3,4; fork contact4 gives no pair of positions. In the doubled-terminal row, a doubled main contact again leaves at most one position. Terminal contact4 gives only L=1; terminal contact3 gives L=1,2.

These assertions are the elementary strict inequalities for the affine function q(s). The checker enumerates all908 role/contact/marked-position possibilities with L<8 and records the11 surviving oriented rows. This is a complete derived domain, independent of m, not a sampled range of large trees. Every row not already a forbidden small affine diagram has at least one contact at a marked v-root.

It follows immediately that both opposite black roots ui are orthogonal to A: each is right to both v-roots, whereas every surviving nonzero contact profile would touch at least one of them.

## A first shared connector cannot have a white A-contact

Choose a shortest path from the set {u1,u2} to A in H. There is one because H is connected. Since neither u-root contacts A, its last outside node z is white and has a preceding node w. The path has no internal black nodes: any other u could be used as its start, and either v already belongs to A. By shortestness every predecessor, in particular w, is right to every root of A. The node w can be an initial u-root. The path supplies the label reduction used above for z.

The complete role list has just two non-small-affine profiles with a white A-contact. Neither admits this continuation w.z<0.

**Fork switch.** Let p be a fork center with leaves a,c; the marked roots are a,p. The extra root has z.a=−sqrt2 and z.c=−1. In the finite corner deleting p, the four nodes z,a,c,w form the doubled three-arm star with bonds4,3,and the nonright w.z bond. At its smallest possible label3 this is affine B3, with positive kernel (1,1/sqrt2,1/2,1/2) in that order. Any stronger bond makes the kernel value negative. Thus the corner cannot be positive definite.

**Two-step doubled terminal.** Let x−4−p−3−y be the initial path from a doubled terminal, with marks p,y. The extra root has simple contacts to x and y. Delete p. The old diagram separates into the isolated x and a long finite B or D tail starting at y, according as the remote decoration is doubled or forked. The new z joins x,y,w, so at the weakest continuation w−3−z it is the new fork center with leaves x,w. Together with the unchanged remote decoration this is exactly affine B or D, of order m+1. A stronger continuation bond again makes its positive-kernel value negative. This contradiction uses positive definiteness of an actual corner, so the affine obstruction is allowed to be large; it is not an illicit application of the order17 prohibition.

The only remaining first-connector profile is consequently a terminal triangle: A is affine B or C, its marked roots are a doubled terminal a and adjacent main node b, z.b=−1, and z.a=−1 or−sqrt2. In particular this first connector has no white A-contact.

## One terminal connector forces the entire carrier to be private

The diagram A−a is a large finite B or D with long endpoint b. Any other outside root r that contacts A is white, since both u-roots were already proved orthogonal to A. The diagram (A−a)+r is therefore a principal subdiagram of an actual finite corner, as is (A−a)+r+z. If r contacts a surviving A-root, the complete finite one-node extension classification at these ranks forces the sole surviving contact to be b by label3. Finite exceptional types are too small; a finite B/D extension can only prolong the long end.

But z already prolongs b. Two such roots r,z produce affine B or D when r.z=0, and a negative perturbation of its positive kernel when r.z<0. Thus no other r can contact any node of A−a. The only remaining possible A-contact is at a.

A sole label3 or4 contact at that doubled terminal creates the small affine F4 or C2 witness from the private-tree lemma. A label at least6 fails already with the additional nonzero a.z contact in the finite corner retaining a and deleting b. For label5, that same corner contains r,a,z and the original white path from z to its initial u-root. Here r is white, since both outside u-roots were already proved orthogonal to A. The connected finite component containing the5 bond would therefore contain both marked roots a and u, of type H3/H4, and have positive payment; otherwise it is already nonspherical. This excludes every further A-neighbor. The path is valid in that corner because its internal nodes are white and its initial u is retained.

Therefore z is the unique outside neighbor of A, and it touches only its two v-roots. All white nodes of A are private v-whites: their only exits from A pass through a v-root, not through a white path toward U. Conversely, any further private v-white component would require another outside A-neighbor at a v-root. The unique such neighbor z is shared, since it has the chosen white path to a u-root. Hence there is no additional private v-component and V=A exactly.

This contradicts the frozen [private affine-tree exclusion](fresh_private_affine_tree_exclusion_20261017.md). It proves that H−U has no affine B/C/D carrier; symmetry proves the same for H−V.

## Verification and exact scope

Generation and a fresh saved replay pass the complete 908 role cases, 11 surviving rows, the grounded cost identities, the exact fork-switch affine kernel, and six literal remote-decoration gluing checks at orders18,19,23. Those last finite checks corroborate the transparent all-rank B/D graph identity. Syntax, six direct source hashes and all local links pass. Separate construction-agent and root mathematical/source reads pass the support and path-specific label reductions, every resistance role, both continuation obstructions and final all-neighbor uniqueness. Neither review claims a duplicate arithmetic replay. Source, data and report are final and frozen.

The proof does not assert the false general statement that every external root with two spherical deletions contacts only black nodes. The fork switch is an explicit counterexample to that statement. What is used here is the actual shared path, its first-contact predecessor, zero finite-corner payment, and the rigidity supplied by one terminal connector. No arithmetic ground field, label grid, bound on extra supporting whites, or hypothetical ambient realization of a local Gram matrix is assumed.

Together with the cycle-carrier lemma and the absence of affine diagrams of order at most17, the result removes all affine subdiagrams from the two complements. Applying the rest of the published finite-quadrilateral payment argument and combining it with triangle/ideal-corner payments and curvature requires a separate composition audit. This report alone makes no ambient dimension claim.
