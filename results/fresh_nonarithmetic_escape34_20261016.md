# An escape-cutting cap closes a tetrahedral three-face in H9

A new exact cap on the retained 34-wall chamber cuts its saved strict escape and closes **three** old open finite edges. Together with old vertex 1, its three new vertices form a complete compact tetrahedral three-face. The resulting chamber has 35 actual facets, **98 finite vertices and 263 finite edges**.

A different exact strict spacelike escape remains. The chamber still has **infinite volume**, and its number of open finite-corner germs increases from 341 to 356. Compact closure of this particular three-face is the geometric gain; there is no finite-volume example or new dimension claim.

The [producer](../scripts/fresh_nonarithmetic_escape34_20261016.py) and [explicit certificate](fresh_nonarithmetic_escape34_20261016.json) are separate from the frozen [34-wall predecessor](fresh_nonarithmetic_wall3_cap_20261016.md) and its [independent review](fresh_nonarithmetic_wall3_cap_20261016_review.md).

## Choice of support

The previous escape was constructed from the exterior ray active on

\[
[0,1,2,3,5,6,7,12,30].
\]

A numerical inspection of its eight-root sub-supports suggested two complete geodesic faces. The successful cap uses the exactly positive-definite support

\[
S=[0,1,2,3,5,6,7,12]
\]

as its right support, with a \(\pi/3\) contact to old wall 8. These are new endpoint conditions, distinct from the earlier right-to-\(T+31\) family. A preliminary finite-label numerical screen was used only to select this exact recipe; it supplies no nonexistence statement, completeness claim, or acceptance gate.

## Exact normal and all old contacts

Let Q be the inherited physical form, with old normals \(n_i\), common strict timelike center H, and previous escape E. The nine covectors indexed by \(S+[8]\) have rank nine. Use the kernel k and the unique p satisfying

\[
(p,n_i)=0\ (i\in S),\qquad (p,n_8)=-1,\qquad (p,k)=0
\]

in the literal normalization saved in the certificate. With

\[
D=24+\frac{72\sqrt5}{5},\qquad r=p-\sqrt D\,k,
\]

the exact identities are

\[
p^2=14+6\sqrt5,\qquad
k^2=-\frac{5+\sqrt5}{16},\qquad (p,k)=0,
\qquad p^2+Dk^2=2.
\]

Thus \(r^2=2\). Both quadratic branches are checked. The chosen branch has \((r,H)<0\) and \((r,E)>0\); the other branch is obtuse to retained walls.

All 34 new pairings are certified: eight right contacts, the \(\pi/3\) contact with wall 8, and 25 strictly dotted contacts, including the preceding cap 33. Every unspecified contact is strictly dotted. No assumed finite-angle alphabet or numerical fit is used to accept a wall. Direct physical projection of H to each wall provides 35 future timelike relative-interior facet witnesses, strict at every other wall.

The arithmetic uses a formal square-root extension over the predecessor's formal extension of the frozen octic field. The new encoding \([a,b]\) means \(a+b\sqrt D\). Norms, products and signs are evaluated exactly; opposing-sign terms are compared by \(a^2-Db^2\), recursively in the previous field. All divisions are by nonzero old-field scalars. No irreducibility or minimal-degree claim is required.

## Complete finite-corner update and tetrahedral closure

The complete finite-neighbor list of new wall 34 is \(S\cup\{8\}\). All nine possible eight-element subsets are tested exactly. Three yield positive-definite nine-root vertex supports:

| Vertex | Incident walls |
|---:|---|
| 95 | 0,1,2,3,5,6,7,12,34 |
| 96 | 0,1,2,3,5,6,8,12,34 |
| 97 | 0,1,2,3,5,7,8,12,34 |

Their physical points are future timelike, with exactly the listed zero pairings and all other wall pairings strictly negative. All 95 old finite vertices are strictly preserved. The independently reviewed complete predecessor census and the complete new-facet census therefore give exactly 98 finite vertices.

The six new finite edges are all pairs among vertices \(1,95,96,97\). In particular, the three edges from old vertex 1 close three old open germs. The finite-edge total is 263, and

\[
9\cdot98-2\cdot263=356
\]

open finite-corner germs remain. This does not count all ideal or exterior rays.

The four vertices have common spherical support

\[
A=[0,1,2,3,5,12].
\]

In \(A^\perp\), their boundary walls are \([6,7,8,34]\). Their projected Gram matrix, with \(\phi=(1+\sqrt5)/2\), is

\[
B=\begin{pmatrix}
3/2&-1&0&0\\
-1&2&-\phi&0\\
0&-\phi&2&-1\\
0&0&-1&2
\end{pmatrix}.
\]

Every principal three-by-three submatrix is positive definite, B is nonsingular, and every entry of \(B^{-1}\) is strictly negative. The checker constructs the four dual rays \(v_j\) with

\[
(v_j,\bar n_i)=-\delta_{ij},\qquad v_j^2=(B^{-1})_{jj}<0.
\]

All four rays are future oriented and satisfy every retained inequality. Since the four projected boundary covectors form a basis of the Lorentz four-space, their halfspace cone is exactly the positive cone of these four rays. Its Klein section is a compact tetrahedron. The other retained walls contain all four vertices, so they do not cut it. This proves that the entire A-face is this compact tetrahedron, not merely that four finite points exist.

This is a compact three-face of the ambient Coxeter chamber. No claim that its intrinsic projected dihedral angles themselves form a Coxeter tetrahedron is made.

## New escape and a stopping rule for this support

The old 32-wall exterior ray at saved index 8 has active support

\[
[0,1,2,4,5,6,7,9,17].
\]

Reconstruct it as V and normalize \((V,H)=-1\). It has positive norm and nonpositive pairings with all 35 walls. The exact point

\[
E_{\mathrm{new}}=V+H/8192
\]

has positive norm, negative pairing with H, and strictly negative pairings with all 35 walls. It proves infinite volume by an open portion of the future ideal boundary inside all halfspaces. No complete exterior-ray census is needed for this conclusion.

The checker also projects \(E_{\mathrm{new}}\) orthogonally to \(A^\perp\). The projection is future timelike, strict at the four boundary inequalities, and is explicitly reconstructed as

\[
E_A=\sum_{j=0}^{3}c_jv_j,\qquad
c_j=-(E_A,\bar n_j)>0.
\]

Consequently **no further coherent acute wall perpendicular to A can cut this escape**, regardless of its field or contact labels. Here coherent means that the proposed extension has a common strict future timelike interior point.

For completeness, the finite-vertex preservation argument allows that common interior point h to differ from H. At any old finite vertex v, its positive-definite active Gram matrix G has nonnegative inverse. Write

\[
h=a v+\sum b_i n_i,\qquad a>0,\qquad
b=G^{-1}((h,n_i))\le0.
\]

An added acute normal u has \((u,n_i)\le0\) and \((u,h)<0\), hence \((u,v)<0\). It is therefore strictly negative on all four \(v_j\), and on their positive combination \(E_A\). If \(u\perp A\), then \((u,E_{\mathrm{new}})=(u,E_A)<0\). Further escape-cutting construction must change this support or use nonright contacts to it.

This stopping rule concerns one fixed six-root support. It does not exclude arbitrary completions or further nonperpendicular caps.

```sh
python3 scripts/fresh_nonarithmetic_escape34_20261016.py
python3 scripts/fresh_nonarithmetic_escape34_20261016.py --verify
```

Generation and fresh saved replay pass. The [independent review](fresh_nonarithmetic_escape34_20261016_review.md) also passes generation and saved replay, separately reconstructing the new normal, all34 contacts, all35 facet witnesses, the complete finite-corner update, the compact tetrahedron and the renewed escape with its positive vertex decomposition. Its formal extension implementation is separate from this producer and shares the explicitly bound frozen octic arithmetic. No frozen predecessor has been modified. This report is final and frozen.
