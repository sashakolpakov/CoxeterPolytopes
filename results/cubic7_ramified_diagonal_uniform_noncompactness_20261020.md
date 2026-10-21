# A uniform obstruction for the ramified-seven diagonal family

**Candidate written theorem, pending independent review.** For
\(t=2\cos(2\pi/7)\), \(\beta=2-t\), \(\alpha=t\beta\), the full
integral reflection chamber of
\(T_n=\langle-\alpha\rangle\oplus I_n\) is not compact for any \(n\ge3\).
No new computation, published-classification assumption or Pell census is
used. The proof closes an entire literal family, not all cubic forms.

The imported own proofs are
[the norm/facet reduction](cubic7_ramified_diagonal_norm_and_facet_reduction_20261020.md)
(`34cef661bd0333b633d4cf0912b096a273635d3bff54db8fa358b7222385e87a`),
[the first edge](cubic7_ramified_diagonal_complete_first_edge_20261020.md)
(`6e8057bee97987b2bd231b07dcf113f769e66c86345d3f99743eca5e7b5f94ed`),
and [the second edge and finite symmetry](cubic7_ramified_diagonal_second_edge_and_finite_symmetry_20261020.md)
(`e4724ccce3557bf719586db0f7a3ae50bee5adb84da36d6a0408ded07c369c0d`).
In particular the complete primitive norm pool is {1,2}, and the full
coordinate norm-one facet is T_(n−1). Set
\(u=t^2+t-1\), \(a=t-1\), \(X=a^{-1}=t^2+2t\).

## 1. The first two proofs extend to dimension three

Those two sources conservatively state n≥5. Here is the missing low-rank
inventory, with no change to any scalar bound or exclusion in them.
At the first edge in dimension three, the retained system is B2, with
root lengths two and one and label4. Its only finite connected rank-three
extension is B3, adding a norm-two root at the long end. H3 has label5,
whose product is not in K; A3 does not contain label4. Adding a norm-one
root at the short end would require an inner product −1/2, impossible in
the integral lattice. Thus the exact same first-edge inventory applies.

At the second edge the white roots are the orthogonal pair
\(r^2=2\), \((-e_3)^2=1\). A new norm-one root is orthogonal to the old
norm-one root, and its product with r is zero or −1. A new norm-two root
has product zero or −1 with the old norm-one root, and product zero,
−1 or −u with r. If both products are nonzero, only the label3/label4
path B3 is finite; the label7/label4 path is not. This is exactly the
same seven-row table in the second-edge source. The finite-product
arguments and every algebraic interval exclusion in that source are
unchanged. Therefore its actual firstness, full new B3 germ and chamber
involution G hold for T3 itself.

In coordinates (time,1,2,3), the fixed subspace of G is
\[
\operatorname{Fix}(G)=\operatorname{span}(M,e_3),\qquad
M=(X;t^2+t,t+1,0),\quad M^2=-2t.
\]
G fixes coordinate3 and is timelike point inversion on the first three
coordinates. It stabilizes the FULL chamber C containing the B3 germ at
\(p=(1;0,0,0)\), with walls
\(\rho_1=-e_1+e_2\), \(\rho_2=-e_2+e_3\), \(\rho_3=-e_3\).

## 2. A full A3 finite germ cannot occur

This is a local statement about T3, not a census. Every A3 simple root has
canonical norm two: for a label3 pair, the squared inner product is
\(q_iq_j/4\); canonical norms (1,1) would require the nonintegral product
1/2, and (1,2) would require \(1/\sqrt2\notin K\). Thus only (2,2)
is possible along the connected A3 chain. Work at the inert prime2.
The lattice is unimodular there, since alpha is a unit at2. Its even-norm
hyperplane modulo2 is
\[
E=\{(z_0,z_1,z_2,z_3):tz_0+z_1+z_2+z_3=0\}.
\]
Indeed \(-\alpha\equiv t^2\pmod2\), and Frobenius is injective in
\(O/(2)=\mathbf F_8\). A lift of its bilinear annihilator is
\(c=(t^{-1};1,1,1)\). Its norm is
\[
c^2=3-\alpha/t^2=4-2/t,
\]
which has 2-adic valuation exactly one.

Suppose the root span S of an A3 germ is primitive at2. Its mod2 image
has dimension three, consists of even-norm vectors, and hence equals E.
Its primitive orthogonal vector v is therefore a unit multiple of c
modulo2. Writing \(v=\lambda c+2w\), with lambda a local unit, gives
\(v^2\equiv\lambda^2c^2\pmod4\), of valuation one. On the other hand,
primitive codimension-one embedding in a unimodular lattice gives
\(v_2(v^2)=v_2(\det S)=v_2(4)=2\), a contradiction. For this determinant
identity, choose w0 with \((v,w_0)=1\); the Schur complement in the basis
S,w0 is \(1/v^2\), so \(\det S/v^2\) is a unit.

If S is not primitive at2, its integral saturation has index2: the square
of the index divides detS=4. The discriminant module of A3 is cyclic
O/(4), so this is the unique index2 overlattice. For standard A3 chain
roots \(s_1,s_2,s_3\), it contains
\(v_0=(s_1+s_3)/2\). The three vectors
\(v_0,v_0+s_2,v_0-s_3\) are orthogonal norm-one vectors. Thus the
saturated positive lattice is I3, and its full norm-one/norm-two zero
system is B3 by the imported trace argument. Consequently an actual
full finite germ containing an A2 chain is B3 or A2+A1, never full A3.
H3 is excluded by its label5.

## 3. What compactness would force at the remaining edge

Assume C compact. The edge from p retaining the actual A2 walls
rho1,rho2 must end at a finite vertex. An isolated new norm-one or
norm-two root would have the form \((x;y,y,y)\). Modulo beta its norm
would give \(3y^2=1\) or2 in F7. Neither is possible: 1/3=5 and2/3=3
are nonsquares. Thus A2+A1 is excluded. Section2 excludes full A3, so
the endpoint has full B3 and a norm-one new simple wall b.

After canonical unit normalization the only two possibilities are
\[
b=(x;y+1,y,y),\qquad\text{or}\qquad
b=(x;y+1,y+1,y).
\]
They have \(x>0\), and their edge-intersection parameter is
\(\alpha x/Z\), where respectively \(Z=3y+1\) or \(3y+2\).
Their norm equation is \(Z^2-3\alpha x^2=1\). Acuteness with the old
actual wall \(\rho_3=-e_3\) gives y≥0. In both cases Z>1; in the first
case this also gives y>0, whereas the second argument needs only y≥0.

The second possibility cannot be first. Put
\(\eta=(x;y+1,y+1,y+1)\). Its norm is \(2(y+1)\), and its reflection
has integral matrix
\[
s_\eta=\begin{pmatrix}
Z&-x\mathbf1^T\\
\alpha x\mathbf1&I-(y+1)J
\end{pmatrix},\qquad Z=3y+2,
\]
where J is the three-by-three all-ones matrix. The identity
\(\alpha x^2=(y+1)(3y+1)\) verifies this formula directly.
Thus eta defines a mirror of the FULL integral arrangement. But it cuts
the same edge at \(\alpha x/(Z+1)<\alpha x/Z\), contrary to firstness.
This uses an actual integral reflection, not an assumed quotient map.

We are left with \(b=(x;y+1,y,y)\), \(Z=3y+1\).
Let R reverse coordinates1 and3. Define
\[
H=\begin{pmatrix}
Z&-x\mathbf1^T\\
\alpha x\mathbf1&-yJ-R
\end{pmatrix}.
\]
It is an integral future-preserving isometry and H²=I. For a direct Gram
check, \(S=-yJ-R\) satisfies
\(S\mathbf1=-Z\mathbf1\),
\(S^2=I+\alpha x^2J\), and \(Z^2-3\alpha x^2=1\).
It sends rho1 to rho2, rho2 to rho1, and rho3 to b. Hence it identifies
the full actual germs at p and the new B3 vertex in the SAME chamber;
therefore H stabilizes C. Its fixed plane is
\[
\operatorname{Fix}(H)=\operatorname{span}(N,e_1-e_3),\qquad
N=(Z+1;\alpha x,\alpha x,\alpha x),\quad
N^2=-2\alpha(Z+1)<0.
\]

## 4. The two involutions cannot preserve a compact chamber

No possible positive x,y above gives a common hyperbolic fixed point of
G and H. To see this, a common fixed vector has the form
\(w=\lambda M+\mu e_3\). Membership in Fix(H) requires
\(w_1+w_3=2w_2\). If lambda is zero, the fixed-plane formula forces w=0.
Otherwise it gives
\(\mu=\lambda\beta(t+1)\), so w is proportional to
\[
W=(X;t^2+t,t+1,\beta(t+1)).
\]
There may be a further time-coordinate restriction, but it cannot help:
\[
W^2=-2t+\beta^2(t+1)^2=(t-1)(2t-1)>0.
\]
Every nonzero common fixed vector is spacelike. Thus G,H have no common
fixed point in hyperbolic space.

A finite group preserving the future sheet has a common hyperbolic fixed
point: sum the orbit of a future timelike vector and normalize it.
Therefore the group generated by G,H is infinite. But the isometry group
of a compact full-dimensional hyperbolic polytope is finite, acting
faithfully on its finite spanning vertex set. Compactness and local
finiteness of the integral arrangement ensure C has finitely many walls
and vertices. This is the contradiction.

Hence C is not compact. The proved full coordinate-facet identity then
excludes compact full chambers of T_n for every n≥3.

To conclude nonreflectivity we explicitly use the standard arithmetic
lattice and compactness theorem: the integral orthogonal group of an
admissible form is an arithmetic lattice, and is cocompact when the form
is anisotropic over K. This is an EXTERNAL theorem input, as in
Borel–Harish-Chandra, *Arithmetic subgroups of algebraic groups*,
Annals of Mathematics75 (1962), §§7,11–12,
[primary paper](https://doi.org/10.2307/1970210).
Here anisotropy follows immediately from positive definiteness at either
other real embedding. If its full integral reflection subgroup were
cofinite, covolume multiplicativity would give it finite index in that
cocompact arithmetic lattice, hence a compact chamber, already excluded.
Thus the same argument rules out a cofinite full integral reflection
group. The arithmetic theorem is not claimed to have been reproved or
checked by an exact worker here.

The actual compact (4,2,4,2) two-face from the reviewed edge proofs remains
consistent with the result. No new higher-dimensional example or general
dimension ceiling is claimed; this is a uniform obstruction for this
specific ramified cubic diagonal family.
