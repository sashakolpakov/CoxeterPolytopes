# Root review: homogeneous marked-D quadrilateral obstruction

21 October2026. **Full written mathematical PASS** for
`results/finite57_homogeneous_large_d_quadrilateral_schur_obstruction_20261021.md`,
SHA `3ef4130fc30f557425fdb06144342054e019c6eedda3e54a197c532b158b10ac`.
Root read all174 lines and rederived the actual role assignment, Schur
matrix, determinant and rank contradiction. This was manual proof/source
inspection only; no enumeration, test, symbolic worker, compilation or
Lean run. Root did not author the theorem. The separate full other-author
review is7de3b872; inherited geometric and payment inputs remain explicit.

The hypothesis supplies the COMPLETE rank n-2 white diagram
A_(r-4)+2A1+Z, and requires Z orthogonal to every black wall. It is
stronger than seeing one marked D at one corner. The one-contact long
end and three-contact center distinguish the black roles independently
of a temporary corner labeling. With r>=7 they are nonadjacent, and
opposite sides retain identical contacts with the fixed white path and
leaves. The resulting adjacent black products are zero. For opposite
actual walls the acute/native gap gives each of g,h either zero or at
least one; no arbitrary formal face angle is assumed Coxeter.

Writing m=r-4>=3, the path inverse has endpoint diagonal m/(m+1)
and cross1/(m+1); BOTH isolated leaves contribute1/2 to the center
pair. Eliminating the white block therefore gives long diagonals
(m+2)/(m+1), center diagonals1/(m+1), adjacent entries-1/(m+1),
and opposite entries-g-m/(m+1), -h-(2m+1)/(m+1). These match the
source exactly. The opposite differences have eigenvalues2+g,2+h.
The remaining sum block has diagonals2/(m+1)-g and
-2m/(m+1)-h and cross-2/(m+1).

Its determinant numerator is

    (m+1)*g*h+2*m*g-2*h-4.

At g=0 this is strictly negative. At g>=1 it is
h*((m+1)*g-2)+2*m*g-4>=2*m-4>0. The native gap leaves no
intermediate g. Thus the Schur complement is nonsingular. The core
has rank r+2 and the unchanged orthogonal positive Z adds n-r,
forcing rank n+2 in an actual Lorentz space of dimension n+1. This is
impossible. The proof does not need a claim that any hypothetical radical
would be a feasible cusp; it proves that the required radical is absent.

The budget limits are stated correctly. At ambient57, the current
2D17+D23 cost remains372. Its fixed D17+D23 inert cost258.5 leaves
only105.5 at the active component if a364 ceiling is sought. A competing
A17 still costs113.25, so nonhomogeneity does not supply a cheap recipient
or a once-only allocation. The older guarded marked-node savings have a
different scope. This new theorem alone deletes no coupon, proves no
universal363.75 allowance and does not lower the bound57. A global
reallocation across all nonhomogeneous corners remains open.
