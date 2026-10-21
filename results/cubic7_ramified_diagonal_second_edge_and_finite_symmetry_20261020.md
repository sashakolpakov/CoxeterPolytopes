# A second complete edge and a finite chamber symmetry in T_n

This paper argument uses the complete first edge
[6e8057be](cubic7_ramified_diagonal_complete_first_edge_20261020.md)
(`6e8057bee97987b2bd231b07dcf113f769e66c86345d3f99743eca5e7b5f94ed`)
and the separately reviewed norm/facet reduction
(`34cef661bd0333b633d4cf0912b096a273635d3bff54db8fa358b7222385e87a`).
There is no new run.
Use the same \(t,\beta,\alpha,u\), and put
\(a=t-1\), \(X=a^{-1}=t^2+2t\).
Assume \(n\ge5\). The actual vertex
\(Q=(t+1;\alpha u,0,\ldots)\) has walls
\(r=(u;t+1,0,\ldots),\alpha_2,\ldots,\alpha_n\).
Release \(\alpha_2\), retaining \(r+B_{n-2}\).

## Inventory and finite upper

An eventual finite simple wall has the form
\(w=(x;y,z,\epsilon,0,\ldots)\), with norm \(q\), where
\(\epsilon=0\) or 1 is its attachment to the long end of \(B_{n-2}\).
Write \((r,w)=-b\). The complete possibilities are

| role | q | epsilon | b | R=2q−2epsilon−b² |
|---|---:|---:|---|---|
| orthogonal short | 1 | 0 | 0 | 2 |
| orthogonal long | 2 | 0 | 0 | 4 |
| B extension, orthogonal to r | 2 | 1 | 0 | 2 |
| I2(4), orthogonal to B | 1 | 0 | 1 | 1 |
| I2(3), orthogonal to B | 2 | 0 | 1 | 3 |
| I2(7), orthogonal to B | 2 | 0 | u | beta |
| whole B_n bridge | 2 | 1 | 1 | 1 |

A norm-one/norm-two finite pair has only the nonzero product one: the
three conjugate squared products are at most two, one strictly, and the
trace census in the first-edge proof excludes every other algebraic
integer. Two norm-two roots have labels 3 or 7: their product is
\(2\cos(\pi/m)\in K\), whose real cyclotomic degree divides three;
the only degree-three alternatives are \(m=7,9\), and the latter field
has discriminant81 rather than49. The orthogonal case has \(b=0\).
No finite connected extension through an I2(7) is possible. The only
connected bridge through the B component is B_n. When \(B_{n-2}=B_3\),
an apparent F4 extension would require an inner product −1/2 between
two canonical norm-one roots, impossible in the integral lattice.
These observations prove completeness of the table, not a bounded list
of guessed contacts.

Let
\[
x_0=t^2+3t+2,\qquad y_0=(t+1)^2=X+1,\qquad
s=(x_0;y_0,y_0,1,0,\ldots).
\]
Direct substitution gives \(s^2=2\), \((s,r)=-1\),
\((s,\alpha_3)=-1\), with every other retained product zero.
Thus it supplies a finite B_n upper, with a collision strictly in the
forward edge ray. The explicit future endpoint is given below.

## Exact finite interval

Put \(K=2x+ub\), \(Z=2(t+1)z\), and \(A=K/X^2\in O\).
The contact and norm equations give
\[
y=\beta x-\frac b{t+1},\qquad
Z^2-2\alpha K^2=2(t+1)^2R,\qquad
(Q,w)=-\frac{\alpha K}{t+1}.
\]
Choose \(K,Z>0\). At the upper, \(K_0=X^2\) and
\(Z_0=2(t+1)^3\). A collision no later than the upper gives
\(A^2\le R\). At the other two embeddings,
\[
A_i^2\le -\frac{a_i^3}{\beta_i}R_i.
\]
Moreover \(A\equiv b\pmod2\): one has
\(K\equiv ub\) and \(X^2\equiv u\pmod2\), and \(u\) is a unit.
The exact identities
\[
\operatorname{Tr}(a^3/\beta)=-7,\qquad
0<a^3/\beta<1/48
\]
give, for constant R, \(\operatorname{Tr}(A^2)<385R/48\).

For R=1, \(b=1\), trace at most eight and \(A\equiv1\pmod2\)
leave only \(A=1\). This follows directly from the signed trace-eight
list in the first-edge proof; its other elements are not1 modulo2.

For R=2, \(b=0\), write \(A=2B\). The trace bound gives
\(\operatorname{Tr}(B^2)\le4\), hence \(B=0,\pm1\).
Positive \(A\) and \(A^2\le2\) exclude them all.

For R=4, again \(A=2B\), and \(\operatorname{Tr}(B^2)\le8\).
The positive physical possibilities \(B\le1\) are
\(1,-v_0,v_0+1,u-1\), where \(v_0=t^2-2\).
The second-conjugate bound is
\(B_2^2\le-(v_0-1)^3/(2-v_0)<45/32\).
It excludes \(-v_0\) and \(u-1\), whose second conjugates have squared
magnitudes \(u^2>3\) and \((t+1)^2>4\).
For \(B=v_0+1\), reconstruction requires \(z^2=2X^2\), impossible
since its field norm is8. For \(B=1\), one obtains a root at exactly the
upper: it is \(r+2\sigma\), where
\(\sigma=(x_0;y_0,y_0,0,\ldots)\).

The R=3 row is excluded before any square test. Modulo \(\beta\),
\(t=2\), \(t+1=3\), and \(b=1\) force \(y=2\) in \(\mathbf F_7\).
Its norm equation then requires \(z^2=5\), a nonsquare.

For R=beta, \(b=u\), the conjugate bounds give
\(\operatorname{Tr}(A^2)\le25\), since
\(\operatorname{Tr}(a^3)=-25\) and \(\beta+a^3<1\).
The congruence \(A\equiv u\pmod2\) makes all three power-basis
coefficients odd. In the trace form, their hexagonal term is at most five,
so the last two coefficients are both \(\pm1\).
Together with \(0<A\le\sqrt\beta<1\), the only possibilities are
\(2-u\) and \(1+t-t^2\).
The former has second conjugate \(2+t>3\), violating
\(A_2^2\le-(v_0-1)^3<27/8\).
For the latter, reconstruction gives
\(x=X,y=t+1,z^2=\beta X^2\), impossible since its field norm is7.

Thus there is no earlier wall. At the upper the R=1 rows give \(s\)
and \(\sigma\), but \(\sigma=s-e_3\) is a positive combination of
\(s\) and the retained walls, because \(-e_3\) is such a combination.
The R=4 root \(r+2\sigma\) is redundant too. The unique new simple wall
is s. Local finiteness on the finite edge segment proves actual firstness.

## Full germ and the honest finite-symmetry boundary

Define \(T=(t+1)^3=x_0+y_0-1\). On the coordinates
\((\text{time},1,2)\), set
\[
G=\begin{pmatrix}
T&-(x_0+u)&-x_0\\
\alpha(x_0+u)&-x_0&-y_0\\
\alpha x_0&-y_0&-y_0
\end{pmatrix},
\]
and fix all other coordinates. To check its Gram identity directly, put
\(M=(X;t^2+t,t+1)\), with \(M^2=-2t\). On this three-space,
\(G=-I-M(M,\cdot)/t\), the timelike point inversion.
Its last two columns are
\(e_3-r-s\) and \(e_3-s\), respectively, so they are orthogonal
norm-one vectors. Its first column is future timelike of norm
\(-\alpha\) and orthogonal to them. All entries lie in O; the Gram
identity makes its determinant a unit, so this is an integral isometry.
It sends \(\alpha_1\mapsto r\), \(\alpha_2\mapsto s\), and fixes
\(\alpha_i\) for \(i\ge3\). Therefore its first column is the actual
new vertex and its FULL zero germ is B_n, by invariance of the full root
arrangement and the already established full germ at p.

The map identifies actual chamber germs in the same full chamber, hence
stabilizes that chamber. But it is only an involution:
\(Gr=\alpha_1\), \(Gs=\alpha_2\), and it fixes the remaining simple
walls, whence \(G^2=I\). Its trace is \(n-3\).
Thus M is its fixed timelike vector in the displayed three-space.

The two proved edges and their images under G close an actual compact
two-face with angle sequence \((\pi/4,\pi/2,\pi/4,\pi/2)\).
Indeed projection perpendicular to the common white roots on coordinates
\(3,\ldots,n\) gives normals \(\alpha_1,-e_2,r,\sigma\), of norms
\(2,1,2,1\). Their successive products are \(-1,0,-1,0\).
This is compatible with a compact higher-dimensional chamber and is NOT
an infinite-order symmetry or a nonreflectivity certificate. The present
route has a rigorous finite compactness gate and two exact edges, but no
compact H9 example or full-host exclusion yet.
