# Complete continuous endpoint domains for the retained H3 face

The actual face on (S=\{0,1,2,3,5,6\}=H_4+A_2) has13 finite vertices and7 exterior rays. Its complete cone has30 combinatorial one-faces. Exactly28 have positive-definite eight-wall support and thus meet hyperbolic space:15 have two finite endpoints, nine have one open ideal end, and four have two open ideal ends. A compact completion must terminate all17 open ends.

The [primary checker](../scripts/fresh_nonarithmetic_h3_endpoint_domains_20261017.py) and [saved certificate](fresh_nonarithmetic_h3_endpoint_domains_20261017.json) passed generation and fresh saved replay. They reduce every possible terminal wall of these edges to **296 finite contact profiles plus14 all-label tails**, with explicit continuous normal formulas and all55 native old-wall pairing formulas. This is a complete domain reduction, not a solved compatibility or closure problem. No new cap is accepted, and no finite-volume or noncompletion claim is made.

The domain certificate reuses the exact physical normal recipes and outward rational interval arithmetic from the [all-profile E55 gate](fresh_nonarithmetic_a2_escape_partition_20261017.md). It has not received an independent arithmetic review. The [independently reviewed complete55 cone](fresh_nonarithmetic_nonright_s6_20261016_cone_review.md) supplies the complete ray inventory; every relevant edge and Gram test is reconstructed here.

## Why these are unavoidable families

The two roots0,1 of the retained (H_4) have angle \(\pi/5\). If the (S)-face had an ideal point in a finite-volume ambient completion, those reflections would belong to its Euclidean cusp link. Their product has irrational linear trace, incompatible with the full translation lattice of a finite-volume cusp. Thus this particular three-face must be compact in any finite-volume completion.

An old geodesic edge persists under acute additions until it is terminated. At each formerly open end of a compact completed edge, at least one added normal together with the eight edge-support normals has positive-definite Gram. Such a wall is in one of the domains below. This condition is necessary; paying all17 old endpoint requirements alone would not prove that no new open region is created.

The30 raw one-faces are obtained directly from the complete common-zero ray sets. Each has exactly two extreme rays and eight common old walls. The (H_4) block is orthogonal to the other four roots in every case. Exact leading-minor tests distinguish the28 positive-definite supports from the two other one-faces. The13 supports with an open end are listed in the table.

| Two extra edge walls | Finite profiles | Isolated-A1 tails |
|---|---:|---:|
|7,30|1|0|
|7,33|7|0|
|7,49|7|0|
|8,44|36|2|
|8,52|36|2|
|30,33|36|2|
|30,44|36|2|
|30,47|36|2|
|30,49|9|0|
|30,52|13|0|
|33,49|7|0|
|34,47|36|2|
|34,52|36|2|
|Total|296|14|

## Complete contact alphabet, including every high label

Normalize all normals to square two. A terminal wall must have positive Schur residual after projection onto its positive-definite eight-root support. The inverse (H_4) block is positive and every diagonal exceeds two; therefore every terminal wall is right to (H_4).

Let (G) be the remaining four-root Gram on nodes ((5,6,a,b)), and let (J=G^{-1}\ge0). The exact certificate verifies the following for every support:

* Every isolated (A_1) slot has (J_{ii}=1/2).
* Every other slot has (J_{ii}\ge2/3).
* Every inverse diagonal is at least (1/2).

For contact magnitudes (c_i=2\cos(\pi/m_i)\ge0), positive residual requires (c^TJc<2). A contact of label at least six on a nonisolated slot already costs at least ((2/3)\cdot3=2). A label at least six on an isolated (A_1) costs at least (3/2); any second nonzero contact costs at least another (1/2), with no negative cross terms available. Therefore every label at least six must occur alone on an isolated (A_1).

All remaining profiles are the vectors of labels2–5 with (c^TJc<2). Slots with inverse diagonal at least two are first forced right. The finite exact enumeration records every accepted vector, each zero-residual vector, the rejected count, and all inverse entries. Zero residual would not give a finite endpoint and is excluded from this necessary compact-completion domain.

Each isolated-A1 tail is the full infinite family

\[
m\ge6,\qquad c=2\cos(\pi/m),\qquad \rho=2-c^2/2>0.
\]

The certificate also supplies the continuous relaxation \(\sqrt3\le c<2\). No high label is replaced by a sample or silently discarded.

For example, the edge ((7,49)) has support (H_4+A_4), with (A_4) path5–6–7–49. Its seven profiles are: entirely right; one label3 at any path node; or one label4 at either endpoint. The edge ((7,30)) has support (H_4^2), so its sole profile is entirely right.

## Exact continuous normals and all retained constraints

For an ordered eight-root support (F), let (G_F) be its Gram. The checker constructs four dual vectors (d_i\in\operatorname{span}(F)), corresponding to its final four roots, so that their products with the support are the appropriate coordinate unit vectors. With negative contact vector (-c), put

\[
p=-\sum_i c_i d_i,\qquad \rho=2-c^TJc>0.
\]

The orthogonal complement of (F) is a Lorentz plane. Project the retained future center into it and normalize to obtain (k_0^2=-1). Project an explicitly recorded coordinate vector, remove its component along (k_0), and normalize to obtain (k_1^2=1), (k_0\cdot k_1=0). Nonzero inverse pivots, negative center-projection norm and positive spatial-projection norm are all certified by rational intervals.

Every square-two normal with those support contacts is exactly

\[
r=p+s\sqrt\rho\bigl(\cosh u\,k_1+\sinh u\,k_0\bigr),
\qquad s\in\{-1,1\},\quad u\in\mathbb R.
\]

For each of the55 retained normals, the certificate stores

\[
R_{ji}=\langle d_i,n_j\rangle,\qquad
b_j=\langle k_1,n_j\rangle,\qquad
d_j'=\langle k_0,n_j\rangle.
\]

Thus its full native contact is

\[
\langle r,n_j\rangle=-R_jc+s\sqrt\rho
\bigl(b_j\cosh u+d_j'\sinh u\bigr).
\]

Every retained constraint is included, including the two nonright-S6 ambient cutters. The remaining Coxeter requirement is that each value belong to

\[
\{0\}\cup\{-2\cos(\pi/m):m\ge3\}\cup(-\infty,-2].
\]

The final interval includes both parallel and dotted pairs; further cusp obstructions may reject some of them. No finite contact outside an assumed coefficient field is rejected here.

Common future interior is also expressed exactly. The saved actual old finite vertex has support (S\cup\{7,8,11\}). Under all55 acuteness inequalities, the condition (r\cdot v<0) is necessary by finite-vertex persistence. It is sufficient because a sufficiently small perturbation of (v) toward the retained strict center then satisfies every old inequality and the new one strictly. Its pairing formula is stored alongside the old-wall formulas. This permits an arbitrary new common center rather than requiring the original selected center to survive.

With (t=e^u>0), every acuteness inequality becomes a quadratic inequality in (t). Mutual finite contacts between two endpoint families become equations in their two parameters. Solving those complete continuous domains, retaining unspecified finite labels as unresolved when necessary, is the next task. The present certificate removes the finite-label cutoff from that task; it does not report a successful coordinated closure.
