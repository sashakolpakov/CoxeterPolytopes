# Exact maximal K3 eigen-intervals carry infinitely many fields

**Status (2026-08-10): proved rank-four, hyperbolic-dimension-one
counterfamily.**  Finite volume, complete restricted-mirror saturation,
reflection maximality in the declared primitive full-trace host, exact
prescribed eigen-nef return, and even a primitive ample polarization of fixed
square do **not** force a finite field set if hyperbolic dimension one is
allowed.

For every integer \(a\ge2\), put

\[
 d=a-1,\qquad
 t=a+\sqrt{a^2-1},\qquad
 k=\mathbf Q(t),\qquad
 A=\mathbf Z[t].
\tag{0.1}
\]

There is a marked minimal projective K3 surface \(X_a\), an integral
self-adjoint \(A\)-action on its rank-four Picard lattice, and a prescribed
arithmetic Coxeter interval

\[
 P_a\subset\mathbf H^1
\tag{0.2}
\]

such that:

1. \(P_a=\mathbf P(\operatorname {Nef}(X_a)\cap E_t^+)\) exactly;
2. the two endpoint reflections preserve the full primitive trace lattice;
3. their infinite-dihedral group contains every reflection in the complete
   restricted K3 \((-2)\)-mirror arrangement;
4. this group is reflection-maximal among \(A\)-integral \(k\)-linear
   reflections of the declared trace host;
5. \(P_a\) has finite length and minimal ground field \(k\);
6. \(X_a\) has a marked primitive ample class \(H_a\) with \(H_a^2=14\); and
7. infinitely many of the fields \(k\) are distinct.

The maximality assertion is deliberately relative to the declared integral
full-trace order and lattice.  No maximality across arbitrary overorders or
commensurable overlattices is claimed.  The result is also genuinely
one-dimensional: it does not settle the required \(n\ge2\) field theorem.

## 1. The diagonal trace lattice

Use the rank-four lattice and companion action

\[
 G_a=2\begin{pmatrix}-1&0\\0&1\end{pmatrix}
       \oplus
       2\begin{pmatrix}1-2a&-1\\-1&-1\end{pmatrix},
\qquad
 \tau_a=
 \begin{pmatrix}0&-1\\1&2a\end{pmatrix}
 \oplus
 \begin{pmatrix}0&-1\\1&2a\end{pmatrix}.
\tag{1.1}
\]

In the original coordinates \((x_1,x_2,y_1,y_2)\), make the unimodular
change

\[
 (Y,X,U,Z)=(x_2,x_1,y_1+y_2,y_1).
\tag{1.2}
\]

Then

\[
 G_a(v,v)=-2q_d(v),\qquad
 q_d=-Y^2+X^2+U^2+2dZ^2.
\tag{1.3}
\]

Thus the K3 roots are exactly the integral vectors

\[
 \mathcal R_d=\{r\in\mathbf Z^4:q_d(r)=1\}.
\tag{1.4}
\]

The \(t\)-eigenspace of \(\tau_a\) is

\[
 E_t=\{(-tA,A,(1-t)B,B):A,B\in\mathbf R\}.
\tag{1.5}
\]

Its restricted form is

\[
 q_d|_{E_t}=-(t^2-1)A^2+2d(t+1)B^2
 =(t+1)\bigl(-(t-1)A^2+2dB^2\bigr).
\tag{1.6}
\]

At the conjugate embedding \(t'=t^{-1}\), the form inside parentheses is
positive definite.  Hence

\[
 h_a=\operatorname {diag}(-(t-1),2d)
\tag{1.7}
\]

is the admissible binary \(k\)-form on the full trace module \(A^2\).

Write \(s=B/A\).  The positive interval in this chart is

\[
 |s|<R,\qquad
 R^2={t^2-1\over2d(t+1)}
     ={t\over t-1}>1,
\tag{1.8}
\]

where

\[
 (t-1)^2=2dt.
\tag{1.9}
\]

## 2. A complete root-free interval

For \(0\le s\le1\), take the future-pointing representative

\[
 p_s=(t,-1,(t-1)s,-s)\in E_t.
\tag{2.1}
\]

Two roots cut its endpoints:

\[
 \alpha_0=(0,0,-1,0),\qquad
 \alpha_1=(1,-1,1,0),
\tag{2.2}
\]

with

\[
 q_d(p_s,\alpha_0)=-(t-1)s,\qquad
 q_d(p_s,\alpha_1)=(t-1)(s-1).
\tag{2.3}
\]

### Proposition 2.1

The open interval \(0<s<1\) contains no restricted K3 root mirror.  Its
closure meets the complete restricted arrangement exactly in the two
endpoint mirrors (2.2).

### Proof

The roots of height \(Y=0\) are
\(\pm(0,1,0,0)\) and \(\pm(0,0,1,0)\).  Choose the sector \(X<0<U\).
For every other root mirror, choose its representative

\[
 r=(h,-m,u,z),\qquad h>0,\quad m,u\ge0.
\tag{2.4}
\]

Changing the signs of \(X\) and \(U\) only weakens the corresponding
inequality on this sector.  The root equation is

\[
 -h^2+m^2+u^2+2dz^2=1.
\tag{2.5}
\]

At \(s=0\),

\[
 q_d(p_0,r)=-th+m<0,
\tag{2.6}
\]

because integrality and (2.5) give \(m\le h\).

At \(s=1\), put \(A_0=h-u\ge0\).  If \(A_0=0\), then
\(m^2+2dz^2=1\), so \(m=1,z=0\), and

\[
 q_d(p_1,r)=1-u\le0,
\tag{2.7}
\]

with equality only for \(\alpha_1\).  If \(A_0\ge1\), (2.5) gives

\[
 m^2+2dz^2=1+A_0^2+2A_0u.
\tag{2.8}
\]

Weighted Cauchy--Schwarz yields

\[
 (m-2dz)^2\le(2d+1)(m^2+2dz^2).
\tag{2.9}
\]

Since \(t>2d+1\),

\[
\begin{aligned}
 &(u+tA_0)^2
 -(2d+1)(1+A_0^2+2A_0u)\\
 &\quad =
 u^2+2A_0u(t-2d-1)
 +A_0^2(t^2-2d-1)-(2d+1)>0.
\end{aligned}
\tag{2.10}
\]

Equations (2.8)--(2.10) imply

\[
 m-u-2dz<tA_0,
\tag{2.11}
\]

and hence

\[
 q_d(p_1,r)=m-u-2dz-tA_0<0.
\tag{2.12}
\]

Every root inequality is linear in \(s\).  Equations (2.6)--(2.12) make it
strict on \(0<s<1\), with only the two declared endpoint equalities.
\(\square\)

## 3. Integral endpoint reflections

On the \(A\)-module \(A^2\) with form (1.7), define

\[
 S_0=\begin{pmatrix}1&0\\0&-1\end{pmatrix},
\qquad
 S_1=
 \begin{pmatrix}
 2t-1&-2(t-1)\\
 2t&-(2t-1)
 \end{pmatrix}.
\tag{3.1}
\]

Direct multiplication gives

\[
 S_i^2=I,\qquad S_i^th_aS_i=h_a.
\tag{3.2}
\]

The fixed points of \(S_0,S_1\) in \(\mathbf H^1\) are respectively
\(s=0\) and \(s=1\).  Their entries lie in \(A\), so both preserve the
primitive full trace lattice.  In the integral basis \((1,t)\) of each
summand, multiplication by \(t\) is the companion matrix \(C_a\), and

\[
 S_1=
 \begin{pmatrix}
 2C_a-I&-2C_a+2I\\
 2C_a&-2C_a+I
 \end{pmatrix}
\tag{3.3}
\]

is a literal integral \(4\times4\) isometry of \(G_a\), commuting with
\(\tau_a\).  The same is immediate for \(S_0=I_2\oplus(-I_2)\).

Let \(W_a=\langle S_0,S_1\rangle\).  It is infinite dihedral, with

\[
 T_a=S_1S_0=
 \begin{pmatrix}
 2t-1&2(t-1)\\
 2t&2t-1
 \end{pmatrix}.
\tag{3.4}
\]

It tessellates the whole hyperbolic line by translates of \(P_a=(0,1)\).
Because \(S_0,S_1\) are integral isometries commuting with the field action,
they preserve \(E_t\) and the complete restricted K3 root arrangement.
Proposition 2.1 therefore gives:

### Corollary 3.1 (restricted saturation)

Every complete restricted K3 root mirror is a \(W_a\)-translate of one of
the two endpoint mirrors.  Hence its reflection belongs to \(W_a\), and the
complete restricted arrangement has one \(W_a\)-orbit of chambers.

Indeed, move any restricted mirror point into \(\overline P_a\); it cannot
land in the interior by Proposition 2.1.

## 4. Reflection maximality in the declared host

### Proposition 4.1

The group \(W_a\) is maximal among reflection groups preserving the
\(A\)-lattice \((A^2,h_a)\).

### Proof

Any extra integral field reflection has a mirror in some \(W_a\)-tile.
Conjugating by \(W_a\), suppose that mirror lies strictly between \(0\) and
\(1\).  Its product with \(S_0\) is an orientation-preserving integral
isometry \(M\) translating by less than \(2\ell_a\), where \(\ell_a\) is the
length of \(P_a\).

Every determinant-one isometry of the binary diagonal form (1.7) has shape

\[
 M=\begin{pmatrix}P&Q\\R&P\end{pmatrix},
\qquad
 (t-1)Q=2dR,\qquad P^2-QR=1.
\tag{4.1}
\]

At the distinguished embedding,

\[
                         1<P<2t-1.
\tag{4.2}
\]

At the conjugate embedding the form is positive definite, so

\[
                              |P'|\le1.
\tag{4.3}
\]

Write \(P=m+nt\), \(m,n\in\mathbf Z\).  Since \(P>1\ge P'\), one has
\(n>0\).  Conditions (4.2)--(4.3) exclude \(n\ge2\): for \(n=2\), the
upper bound forces \(m\le-2\), whereas the conjugate lower bound forces
\(m\ge-1\); for \(n\ge3\), the two intervals are disjoint because
\((n-2)t^2>n\).  Thus \(n=1\), and (4.3) leaves only

\[
                         P=t\quad\hbox{or}\quad P=t-1.
\tag{4.4}
\]

From (4.1),

\[
                    Q^2={2d(P^2-1)\over t-1}.
\tag{4.5}
\]

For \(P=t\), this says \(Q^2=2d(t+1)\).  If
\(Q=q_0+q_1t\), coefficient comparison gives

\[
 q_0^2-q_1^2=2d,\qquad q_1(q_0+(d+1)q_1)=d.
\tag{4.6}
\]

With \(r=q_0/q_1\), equations (4.6) force
\((r-1)^2=2(d+2)\).  Writing its integral square root as \(2s\) gives
\(d=2s^2-2\) and then

\[
 s-1=q_1^2s\quad\hbox{or}\quad s+1=q_1^2s,
\tag{4.7}
\]

both impossible.

For \(P=t-1\), equation (4.5) becomes

\[
                         Q^2=1+(2d-1)t.
\tag{4.8}
\]

Its constant coefficient would give \(q_0^2-q_1^2=1\), hence
\(q_1=0,q_0=\pm1\), contradicting the nonzero \(t\)-coefficient in
(4.8).  Thus no interior integral reflection exists. \(\square\)

The proof is internal to the declared order \(A=\mathbf Z[t]\).  It neither
uses nor claims maximality after changing the order or lattice.

## 5. Finite volume and the exact ground field

For the form (1.6), the hyperbolic coordinate is \(s/R\).  Hence

\[
 \ell_a=\operatorname {arctanh}{1\over R}
       =\operatorname {arctanh}\sqrt{{t-1\over t}},
\qquad
 \cosh^2\ell_a=t.
\tag{5.1}
\]

The interval has finite positive length.  The cyclic Cartan product of its
two ultraparallel endpoint mirrors is

\[
 {4(n_0.n_1)^2\over n_0^2n_1^2}
 =4\cosh^2\ell_a=4t.
\tag{5.2}
\]

Therefore its minimal Vinberg field is exactly
\(\mathbf Q(4t)=k\).

## 6. Exact polarized K3 return

Let \(e_i,f_i\) be three hyperbolic-plane bases in \(\Lambda_{K3}\), and let
\(v\in E_8(-1)\) be a primitive root.  The block basis embeds by

\[
\begin{aligned}
 x_1&\longmapsto e_1-f_1,&
 x_2&\longmapsto e_2+f_2,\\
 y_1&\longmapsto e_3+(2-2a)f_3+v,&
 y_2&\longmapsto v.
\end{aligned}
\tag{6.1}
\]

Its pullback Gram is \(G_a\).  After replacing the third image by its
difference with the fourth, the rows \(e_1,e_2,e_3,v\) give a
determinant-one minor, so the embedding is primitive.  A generic period in
the signature-\((2,16)\) complement gives a projective K3 surface with exact
Picard lattice \(L_a\).  Choose the
ambient complete root chamber containing \(p_{1/2}\), and use a K3 Weyl
re-marking to carry it to the ample chamber.  Conjugate
\(\tau_a,S_0,S_1,\alpha_0,\alpha_1\) by the same marking.  Proposition 2.1
then gives the prescribed equality

\[
 \boxed{
 P_a=\mathbf P\bigl(\operatorname {Nef}(X_a)\cap E_t^+\bigr).}
\tag{6.2}
\]

No hidden projected root cuts its interior.

Every returned integral NS operator is induced by a rational
divisor-product correspondence.  If \(D_i,D_i^\vee\) are intersection-dual
bases of \(\operatorname {NS}(X_a)_{\mathbf Q}\), then

\[
       Z_M=\sum_iD_i^\vee\times M(D_i)
\tag{6.3}
\]

acts as \(M\) on \(\operatorname {NS}_{\mathbf Q}\).  In particular this
applies to the field generator \(\tau_a\) and to the two trace reflections.

There is even a uniform ample polarization.  Before re-marking, take

\[
 h=(3,-1,1,0)
\tag{6.4}
\]

in the \((Y,X,U,Z)\) coordinates.  It is primitive and

\[
                         G_a(h,h)=14.
\tag{6.5}
\]

For a sector root (2.4),

\[
 q_d(h,r)=-3h_r+m+u<0,
\tag{6.6}
\]

because

\[
 m+u\le\sqrt{2(m^2+u^2)}
       \le\sqrt{2(h_r^2+1)}<3h_r.
\tag{6.7}
\]

It is also strict on the two height-zero sector walls.  Thus \(h\) lies in
the ambient chamber interior, and its re-marking \(H_a\) is primitive ample
with square \(14\).

This bounded polarization still does not control the arithmetic packet.
In the original block basis,

\[
 h=(-1,3,0,1),
\qquad
                         h.\tau_a(h)=32a-10,
\tag{6.8}
\]

while

\[
 \det L_a=-32(a-1),\qquad
 \operatorname {disc}(x^2-2ax+1)=4(a^2-1).
\tag{6.9}
\]

Maximal exact return and bounded polarization therefore fail to bound the
correspondence height, the Picard discriminant, or the field discriminant.
A finite-field theorem still needs one of those normalized integral bounds.

## 7. Ambient contrast: these lattices are not 2-reflective

The finite-volume theorem concerns the distinguished eigenline, not the
ambient hyperbolic three-space.  In fact the complete ambient K3 root
chamber has infinite volume for every \(d\).

For every \(k\ge1\), define

\[
 r_k=(dk^2,0,dk^2-1,k).
\tag{7.1}
\]

Then \(q_d(r_k)=1\).  Let \(H=dk^2\) and

\[
 b_k=(1+H^2,0,H(H-1),Hk)
      =(1,0,0,0)+Hr_k.
\tag{7.2}
\]

Thus \(b_k\in r_k^\perp\) and \(q_d(b_k)=-1-H^2<0\).  For any sector root
\(s=(h,x,u,z)\), put \(A_0=h-u\ge0\).  A direct completion of squares gives

\[
\begin{aligned}
 2A_0\bigl(h-Hq_d(r_k,s)\bigr)
 &=(H+1)(x^2-1)
   +{A_0^2(H^2+1)\over H+1}\\
 &\quad
   +2d(H+1)
      \left(z-{HkA_0\over H+1}\right)^2.
\end{aligned}
\tag{7.3}
\]

For \(A_0=0\) the same quantity is positive directly.  If \(A_0\ge1\),
the right side is nonnegative; in the only delicate case
\(A_0=1,x=0\), the integer

\[
 (H+1)z-Hk\equiv k\pmod {H+1}
\tag{7.4}
\]

has absolute value at least \(k\), because \(H+1\ge2k\).  Equality occurs
only for \(s=r_k\).  Hence \(b_k\) supports a genuine chamber facet, and a
small tangent perturbation makes the two height-zero inequalities strict.
The chamber has infinitely many distinct facets \(r_k^\perp\).

All simple-root inner products are integral, so finite Coxeter exponents are
only \(2\).  Infinitely many facets give an infinite-dimensional
\(\mathbf F_2\)-abelianization of the ambient root group; it is not finitely
generated.  A finite-covolume hyperbolic lattice is finitely generated.
Therefore the ambient root chamber has infinite volume.

It is not reflection-maximal in \(O(L_a)\), either.  The vectors

\[
 (0,1,-1,0),\qquad(0,0,0,1)
\tag{7.5}
\]

have \(q_d\)-norms \(2\) and \(2d\), and their integral reflections swap
\(X,U\) and change the sign of \(Z\), respectively.  They are not K3
\((-2)\)-root reflections.  This ambient failure does not affect the
restricted maximality proved in Proposition 4.1.

## 8. Infinitely many fields and exact scope

The fields \(\mathbf Q(\sqrt{a^2-1})\) are infinite: a fixed squarefree field
contributes only \(O(\log X)\) Pell first coordinates \(a\le X\), while the
family contains every \(a=2,\ldots,X\).

Consequently, any field-finiteness statement whose hypotheses allow
\(\mathbf H^1\) and host-relative maximality is false even after adding:

- a minimal exact-Picard K3 carrier;
- a primitive full trace lattice;
- integral self-adjoint correspondence action;
- finite-volume exact prescribed eigen-nef return;
- complete restricted-root saturation;
- reflection maximality in the declared trace host; and
- a primitive ample polarization of fixed degree \(14\).

What remains open is materially stronger:

- hyperbolic dimension \(n\ge2\);
- absolute maximality across overorders or commensurable hosts; or
- a normalized bound on the complete correspondence/operator packet.

The exact replay is supplied by the
[certifier](../scripts/certify_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py),
[archive](../data/enumerations/rank4_k3_exact_maximal_eigen_interval_infinite_fields.json),
[tests](../tests/test_rank4_k3_exact_maximal_eigen_interval_infinite_fields.py),
and [hostile scope audit](rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md).
