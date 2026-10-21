# Private affine trees are excluded from an unpaid finite quadrilateral star

**Theorem.** Consider the actual wall star of a quadrilateral two-face with four finite vertices in a finite-volume acute Coxeter polytope. Assume all four corners have zero payment under Vinberg's original finite-vertex weights, and the ambient diagram contains no connected affine subdiagram of order at most17. In the private/shared decomposition below, neither private diagram can be affine. In particular the large affine B/C/D alternatives are excluded. This does not yet establish the absence of affine subdiagrams in the larger indefinite complements H−U and H−V, and gives no new dimension bound.

The [checker](../scripts/fresh_private_affine_tree_exclusion_20261017.py) and [data](fresh_private_affine_tree_exclusion_20261017.json) corroborate the local identities and diagrams. The argument is all-rank; its completeness does not come from the three orders sampled by the checker.

## Actual-star hypotheses and fixed roles

Write the opposite black pairs as {u1,u2} and {v1,v2}, and the supporting white diagram as W. Each W+ui+vj is positive definite. In the connected indefinite component H of the star, the Gram matrix has negative index one and corank one, and all four black roots belong to H. These are actual-face facts proved in the [finite-quadrilateral transfer audit](fresh_simple_finite_quad_vinberg_transfer_20261017.md); their proofs use the compactness of this two-face, not simplicity of other ideal vertices. We use that geometry here, not that report's all-simple affine exclusion.

Zero corner payment implies ui.vj=0 for every i,j. A private u-white component has no white path to a v-root. Let U consist of u1,u2 and their private white components; define V symmetrically. The other white components are shared. Thus U and V are orthogonal, and any shared root adjacent to U contacts only u1,u2. In particular no shared root contacts a private white node of U. A shared root has a white path ending at a v-root; the first root after it on such a path is outside U. This is the private/shared decomposition of Vinberg §6.5, not an arbitrary partition of an abstract Gram matrix.

Any connected affine component of U must contain both u-roots: otherwise it lies in one of the positive-definite corners. Every remaining private component is attached to a u-root, so an affine U is connected. Under the order restriction it has at least18 nodes. The finite/affine classification leaves only classical A/B/C/D types. The cycle case is already excluded by the [private affine-cycle lemma](fresh_private_affine_cycle_exclusion_20261017.md), which uses only the weaker exclusion of affine diagrams of order at most9. We treat a tree U below.

The finite classification and payment facts are the published inputs: a large connected finite diagram containing a doubled bond is B; a large simply-laced finite diagram is A or D; a finite diagram with a bond5 is H2/H3/H4, and every marked pair in its connected component receives positive Vinberg payment; a bond of finite label at least6 can only occur in an isolated rank-two finite component. See Vinberg §6.1–3 and the [frozen finite-weight interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md), with the [original paper](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt). No arithmetic-field hypothesis is used.

## A complete classification of a shared U-neighbor

Choose a shared white root z with a contact to U. Such a root exists because H is connected and all native cross-black contacts vanish. Its contacts with U are confined to u1,u2. Every diagram (U−ui)+z is spherical.

First, a nonzero z.ui contact has label3 or4. A label5 would place ui and the v-end of the shared path in one H3/H4 finite corner component, giving positive payment. A label at least6 admits no connected finite continuation of that path beyond its rank-two component. This argument uses the actual shared path and zero payment; the two spherical deletions alone would not justify the label restriction.

Second, z cannot contact only one U-root. A simple or doubled leaf on a classical affine tree of order at least18 always gives an induced connected affine diagram of order at most9. The following local witnesses cover all positions; their unused continuations are omitted when taking the induced subdiagram.

| Attachment | Small induced affine witness |
|---|---|
| Simple leaf at an ordinary main-path node with simple directions of lengths at least2 and5 | E8 with arms1,2,5 |
| Simple leaf within the two-node neighborhood of a doubled end | The doubled three-arm star B3, or its one-edge extension B4 |
| Simple leaf at the doubled terminal itself | F4 path with successive labels3,4,3,3 |
| Simple leaf at a fork center | D4 four-leaf star |
| Simple leaf at a fork leaf | E8: at the old fork, the arms have lengths1,2,5 |
| Doubled leaf at a main node with two old simple neighbors | B3 star with labels4,3,3 |
| Doubled leaf at a node incident to an old doubled edge | C2 path with labels4,4 |
| Doubled leaf at a fork leaf | B4: keep its center, sibling leaf, first main-path node, the old leaf and the new doubled leaf |

For the ordinary simple-leaf row, the available all-simple path through a main node has at least m−3 edges; if both directions have length at least2, one has length at least5. The listed near-end cases cover the remaining positions. At a fork leaf, the far main direction is long enough to take five edges. These are uniform arguments for m≥18. In the potentially easy-to-miss doubled fork-leaf case, with fork p, leaves a,c, first main node t and new bond a−4−z, the induced five-node diagram on (p,a,c,t,z) has positive kernel (sqrt2,sqrt2,sqrt2/2,sqrt2/2,1). Thus it really is a forbidden small affine diagram, even though a different deletion can leave a long B path.

Consequently z contacts both marked roots. Normalize the positive cusp kernel nu of U by assigning value1 to main/fork-center nodes, 1/sqrt2 to doubled terminals, and1/2 to fork leaves. Then diag(nu) G_U diag(nu) is a weighted tree Laplacian: edge resistance is1, except that fork-leaf edges have resistance2. If R_ij is the resistance along the unique path from i to j, the grounded Green identity gives

    (G_(U−j)^−1)_ii = nu_i^2 R_ij.

Indeed the grounded Laplacian Green column at i is (R_kj+R_ij−R_ki)/2. Multiplication by the Laplacian gives e_i−e_j, proving the identity at every rank. Write t_i=−z.ui, so t_i^2 is1 or2. The two spherical deletions require

    t_i^2 nu_i^2 R_ij < 2,       t_j^2 nu_j^2 R_ij < 2.

These inequalities exhaust the following possibilities:

* Two main nodes must be adjacent and both contacts simple. Their triangle with z is affine A2, forbidden.
* Two fork leaves must be siblings, both contacts simple. Leaves at opposite forks have resistance m−1≥17 and fail the inequalities; siblings have resistance4. Together with their center and z the surviving sibling case is affine A3, forbidden.
* A main node and a fork leaf have resistance at least2, already violating the strict main-node inequality.
* Terminals at opposite ends have resistance m−1≥17, violating a terminal inequality.
* The only remaining case is a doubled terminal a and its adjacent main node b. The contact z.b is simple, while z.a has label3 or4.

Thus affine D is impossible, and an affine B or C private diagram has this precise terminal triangle. No distance-two payment shortcut is being used: long classical finite diagrams can have unpaid distance-two marked pairs.

## Exactly one connector and an exact Schur reduction

Delete the marked doubled terminal a. The connected diagram U−a is a large finite B or D, with b its long endpoint. Each shared U-neighbor extends b by a simple bond, by the preceding classification. There cannot be two distinct such roots z,z'. If z.z'=0, the induced U−a+z+z' is respectively affine B or D. If z.z'<0, its product against that same positive affine kernel is strictly negative. Either case contradicts the positive definiteness of the corner containing b,z,z'. Hence z is the unique root outside U contacting U.

The block B0=U+z is nonsingular Lorentzian. Grounding at a first eliminates the positive-definite block U−a and leaves

    [[0,−h], [−h,c]],       h>0.

Here h>0 follows from the strictly negative pairing of z with the positive null combination for U. In the literal terminal normalization h=t_a+sqrt2 and c=1, but these values are unnecessary. The determinant is −h^2 and the inverse entry indexed by z is zero.

Let X=H\(U+z). Every X-root is orthogonal to U; its only possible contact with B0 is at z. Therefore the Schur correction on X is identically zero: its Gram matrix after eliminating B0 is exactly G_X. The negative index one and corank one of H imply that G_X is positive semidefinite with corank one.

There is exactly one connected affine component A of X; all other components are positive definite. Both v-roots belong to A. Otherwise the affine component would survive one v-deletion and be contained in an actual finite corner. The no-small-affine hypothesis gives |A|≥18, so A is again classical. Since H is connected, A has a contact to z.

This step is the reason one cannot simply assume that an affine U forces an affine V. It proves the required second affine component using the full actual-star corank, without an affine-order cutoff on H−U or H−V.

## The second component becomes private and contradicts a finite corner

The finite diagram U−a+z is a large B or D with long endpoint z. In either actual corner obtained by deleting one v-root, any X-root attached to z must extend this long endpoint by a single simple bond. Moreover at most one such neighbor can survive: two would again create an affine B/D subdiagram, or a negative perturbation of one. The whole finite extension cannot be exceptional at these ranks.

Suppose z contacts a white root of A. That root survives both v-deletions. No other A-contact can coexist with it: any second white contact survives both deletions, and a contact to either v-root survives the other deletion. Hence z would be a single simple leaf on A. If A is a tree, the preceding leaf witnesses give a forbidden small affine subdiagram; if A is an A-cycle, a simple leaf gives induced E8 with arms1,2,5. Thus z has no white contact to A.

It must contact both v-roots, by simple bonds. A sole v-contact would be the same forbidden simple leaf. Any other X-component would have to attach to z to keep H connected, but such a contact survives both v-deletions alongside one of the two v-contacts. This is impossible by the preceding finite-extension restriction. Therefore X=A. Its whites are private v-whites: within H their only possible paths to U would pass through a v-root, because they have no edge to z. Thus X is precisely the private V.

If V were an A-cycle, the already proved private-cycle lemma would exclude it. If it is a tree, apply the complete two-contact classification above. V has a marked doubled terminal a' and adjacent marked main node b', and z.b'=−1. In the actual finite corner deleting a,a' and retaining b,b', the induced diagram

    (U−a) + z + (V−a')

is two large finite B/D arms joined by simple bonds at their long ends. Their remote decorations give exactly affine C for B/B, affine B for B/D or D/B, and affine D for D/D. This contradicts the corner's positive definiteness. The possible label4 on z.a or z.a' is irrelevant in this corner, since those terminals were deleted.

This proves the private-tree exclusion, and together with the inherited private-cycle lemma excludes every affine private U or V under the stated hypotheses.

## Verification and limits

Generation and a fresh saved replay pass: 1,731 exact grounded Green-column checks, 363 explicit small-affine leaf witnesses, 6,924 two-contact tests, 12 terminal-gluing graph identities, and the symbolic zero-Schur-correction identity. The affine witness catalogue has explicit positive kernels, verified as weighted Laplacians. The probes use orders18,19,23; the report supplies the all-rank geometry and classification. Syntax, five direct source hashes and all local links pass. Separate construction-agent and root mathematical/source reads pass the all-rank leaf and resistance classification, connector uniqueness, literal Schur reduction, both-deletion and connectivity arguments, and final terminal gluing. Neither review claims an independent replay of the arithmetic probes. Source, data and report are final and frozen.

The result is deliberately limited to private diagrams in the actual unpaid finite-quadrilateral decomposition. An arbitrary external root with two spherical marked deletions can contact white roots of a large affine diagram; that broader extension statement is false without further shared-path restrictions. In particular this note does not silently infer that every affine subdiagram of H−U or H−V is private. Those complement-diagram issues remain separate. Frozen predecessors, the unrestricted dimension163 proof and the all-simple dimension29 proof are unchanged.
