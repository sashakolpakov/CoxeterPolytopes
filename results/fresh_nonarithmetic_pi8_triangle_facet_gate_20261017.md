# The new r8 facet needs more than one closing side

**Exact generation and fresh saved replay PASS.** No single additional Coxeter wall can close the actual two-face supported by
\[
W=[0,1,2,3,5,6,55]=S_6+[r_8]
\]
as a compact triangle with its retained sides7 and30. The statement includes every finite label, arbitrary coefficient fields, and nonright contacts to the A2 factor or to \(r_8\). It does not exclude a polygon with two or more additional sides.

The [source](../scripts/fresh_nonarithmetic_pi8_triangle_facet_gate_20261017.py) and [saved certificate](fresh_nonarithmetic_pi8_triangle_facet_gate_20261017.json) bind the [58-wall tetrahedral-face construction](fresh_nonarithmetic_pi8_tetrahedron_completion_20261017.md) and the previously reviewed [complete endpoint-profile method](fresh_nonarithmetic_h3_endpoint_domains_20261017.md). This is a new primary exact certificate using the frozen physical arithmetic and finite-profile helper, not an independent implementation of those dependencies.

## Complete contact domain

The existing vertex \(W+[7,30]\) has positive rank9 diagram \(H_4+H_4+A_1\). A proposed third side \(w\) of a compact triangle must give finite endpoints on both retained sides. Thus
\(W+[7,w]\) and \(W+[30,w]\) must both be positive definite.

In a finite diagram no root can attach to the H4 factor. After omitting that fixed factor, the two endpoint supports are \(A_3+A_1\) and \(A_2+A_1+A_1\). Their complete contact domains are reconstructed from their exact inverse Grams. Nonnegative inverse entries make the projection cost monotone in all contact magnitudes. A nonisolated inverse diagonal is at least2/3, so a label at least6 already costs at least2 and is impossible. At an isolated A1 the diagonal is1/2; a label at least6 costs at least3/2, and any second nonzero contact adds at least1/2. Hence every high-label possibility is a single contact to an isolated A1.

Matching the two endpoint domains gives exactly:

- 46 finite profiles, in the order \([5,6,55,7,30]\), with all labels in2–5;
- one tail with its only nonright contact to55;
- four tails with label2,3,4,or5 to7, arbitrary label \(m\ge6\) to30, and all other contacts right.

These are necessary domains for an actual compact triangle, not an arbitrary finite-label cutoff.

## Exact rejection

The nine fixed roots \(W+[7,30]\) span a positive subspace with a negative one-dimensional orthogonal complement. For a prescribed contact row, write
\[
w=p\pm\sqrt{p^2-2}\,k,\qquad k^2=-1,
\]
where \(p\) is the uniquely determined positive-space projection. Four of the46 finite profiles have \(p^2<2\), so no norm2 normal exists. Every branch of the other42 profiles is strictly obtuse to a retained old wall. The certificate saves a strict rational interval blocker for each of the84 roots.

For the sole55-tail, \(p=-c\,r_8/2\), with \(c=2\cos(\pi/m)<2\). Thus \(p^2=c^2/2<2\), again giving no real norm2 root.

For each30-tail, write \(p=a d_u+c d_v\), where \(a\in\{0,1,\sqrt2,\varphi\}\) and \(\sqrt3\le c\le2\). The squared residual \(\rho=p^2-2\) is strictly positive throughout this entire closed interval. On the opposite kernel branch, the product with old4 is strictly positive: its linear projection term is positive and the future kernel has negative product with4.

On the future branch, the product with old49 is
\[
L(c)+k_{49}\sqrt{\rho(c)},\qquad L(c)>0,\quad k_{49}<0.
\]
The quadratic \(L(c)^2-k_{49}^2\rho(c)\) has three strictly positive degree2 Bernstein coefficients on \([\sqrt3,2]\), all stored as outward rational intervals. Consequently the product with49 is strictly positive throughout. This excludes every finite label in each tail, without adjoining its cosine field or truncating its label.

The certified compact tetrahedron from the preceding construction lies on the different support \(H_4+I_2(8)\) and remains valid. The present obstruction concerns the original S6 face's new r8 facet. The next closure problem requires at least two additional inducing sides, or a changed retained-wall construction.

Final generation76350 and fresh saved replay18096 both exited0. Syntax, direct source hashes, and report links passed before freezing. No claim of ambient finite-volume completion or general noncompletion is made.
