# Canonical escape still does not supply simultaneous effectivity

**Status (2026-08-10): proved exact \(I_{1,9}\) obstruction; proved that the
packet is not reflection-maximal on its natural host; maximal plus canonical
escape nonrealizability remains OPEN.**

This note pushes the P-first realizability gate beyond the rank-three
canonical-escape obstruction.  It gives a complete finite-volume Coxeter
packet for which

\[
                         K\notin P
\]

and every wall has the numerical adjunction color of a rational
\((-1)\)- or \((-2)\)-curve, but simultaneous effectivity is impossible.
It then exhibits the exact extra reflecting root which prevents this packet
from being a maximal counterexample.

No surface classification or chamber census is used.  The only finite
polyhedral calculation is an exact extreme-ray replay for this one displayed
packet.

## 1. The integral Hamming packet

Use

\[
 L=I_{1,9}=\langle H,E_1,\ldots,E_9\rangle,\qquad
 x.y=x_0y_0-\sum_{j=1}^9x_jy_j,
\]

and mark

\[
 K=-3H+E_1+\cdots+E_9,\qquad F=-K=3H-E_1-\cdots-E_9.  \tag{1.1}
\]

Thus \(K\) is characteristic, \(K^2=F^2=0=10-\operatorname{rk}L\), and
\(F\) is primitive.

Define eight mutually orthogonal roots:

\[
\begin{array}{ll}
A_1=E_1-E_2,&A_2=E_3-E_4,\\
A_3=E_5-E_6,&A_4=E_7-E_8,\\
A_5=2H-E_1-\cdots-E_6,&
A_6=2H-E_1-E_2-E_3-E_4-E_7-E_8,\\
A_7=2H-E_1-E_2-E_5-E_6-E_7-E_8,&
A_8=2H-E_3-\cdots-E_8.
\end{array}                                             \tag{1.2}
\]

They satisfy

\[
 A_i^2=-2,\qquad F.A_i=0,\qquad A_i.A_j=0\quad(i\ne j).
\]

Put

\[
                         \beta_i=F-A_i.                 \tag{1.3}
\]

Then

\[
\beta_i^2=-2,\quad A_i.\beta_i=2,\quad
A_i.A_j=A_i.\beta_j=\beta_i.\beta_j=0\quad(i\ne j).     \tag{1.4}
\]

Let \(\mathcal H\subset\mathbf F_2^8\) be the explicit binary code generated
by

\[
\begin{pmatrix}
0&0&0&1&0&1&1&1\\
0&0&1&0&1&0&1&1\\
0&1&0&0&1&1&0&1\\
1&0&0&0&1&1&1&0
\end{pmatrix}.                                          \tag{1.5}
\]

It has one word of weight \(0\), fourteen of weight \(4\), one of weight
\(8\), and minimum distance \(4\).  For \(c\in\mathcal H\), define

\[
 E_c=E_9-\frac12\sum_{i=1}^8c_iA_i+\frac{|c|}{4}F.      \tag{1.6}
\]

The sixteen vectors (1.6) are integral.  Direct calculation gives

\[
\begin{aligned}
 E_c^2&=-1,&F.E_c&=1,\\
 E_c.A_i&=c_i,&E_c.\beta_i&=1-c_i.                     \tag{1.7}
\end{aligned}
\]

The 32 proposed walls are

\[
 \{A_i,\beta_i:1\le i\le8\}\ \cup\
 \{E_c:c\in\mathcal H\}.                               \tag{1.8}
\]

Their canonical pairings are exactly the physical ones:

\[
 K.A_i=K.\beta_i=0,\qquad K.E_c=-1.                    \tag{1.9}
\]

Every wall reflection is integral on \(L\): for a square \(-2\) wall it is
\(x\mapsto x+(x.C)C\), and for a square \(-1\) wall it is
\(x\mapsto x+2(x.C)C\).

The ten walls

\[
\begin{gathered}
A_1,A_2,A_3,A_4,\ E_{00000000},E_{00010111},
E_{00101011},\\
E_{01001101},E_{10001110},E_{00111100}
\end{gathered}                                         \tag{1.10}
\]

have coordinate determinant \(-1\).  Hence the wall lattice is the whole
unimodular lattice \(L\), not a proper root sublattice.

## 2. It is a complete finite-volume Coxeter chamber

All distinct-wall pairings are nonnegative.  Their Coxeter-order histogram
is

\[
\begin{array}{c|ccc}
\text{order}&2&4&\infty\\ \hline
\text{unordered pairs}&352&128&16.
\end{array}                                             \tag{2.1}
\]

The vector

\[
                 h=(5,-2,-1,-2,-1,-2,-1,-2,-1,-1)     \tag{2.2}
\]

satisfies

\[
                 h^2=4,\qquad h.C=1
\]

for every one of the 32 walls.  Thus it is a strict future interior point.

An exact rational double-description calculation gives 2306 primitive
extreme rays:

\[
\begin{array}{c|ccc}
\text{square}&0&1&2\\ \hline
\text{number}&226&1040&1040.
\end{array}                                             \tag{2.3}
\]

The incidence histogram is

\[
2080\text{ finite rays on 9 walls},\quad
224\text{ ideal rays on 12 walls},\quad
2\text{ ideal rays on 16 walls}.                       \tag{2.4}
\]

Every inequality supports a rank-nine facet and contains exactly 670
extreme rays.  All extreme rays are future timelike or future isotropic;
there is no hyperideal direction.  Hence (1.8) defines a finite-volume
hyperbolic Coxeter polyhedron \(P\), and the displayed side reflections
have \(P\) as a fundamental chamber.

The sorted primitive-ray hash and the full exact replay are in the companion
certificate.

## 3. Canonical escape holds

The vector \(F=-K\) is an ideal point of \(P\):

\[
 F.A_i=F.\beta_i=0,\qquad F.E_c=1.                     \tag{3.1}
\]

On the other hand,

\[
                              K.E_c=-1,
\]

so \(K\notin P\).  This packet therefore passes the raw canonical-escape
gate which the rank-three counterexample fails.

## 4. Simultaneous effectivity is impossible

### Theorem 4.1

There is no smooth projective rational surface \(S/\mathbf C\), with a
marking \(\operatorname{NS}(S)\cong L\) carrying \(K_S\) to \(K\), on which
all sixteen classes \(A_i,\beta_i\) are represented by distinct irreducible
curves.  A fortiori, the 32-wall packet cannot be returned as the full nef
chamber of such a surface.

### Proof

Adjunction would make every displayed square \(-2\) curve smooth rational.
For each \(i\),

\[
                           D_i=A_i+\beta_i=F=-K_S.      \tag{4.1}
\]

Thus the eight \(D_i\) are distinct effective anticanonical divisors.
For \(i\ne j\), all four intersections between their components vanish by
(1.4).  Consequently \(D_i\) and \(D_j\) have disjoint supports.

A rational surface has \(\operatorname{Pic}^0(S)=0\), so equality in
\(\operatorname{NS}(S)\) is linear equivalence.  Any two of the disjoint
divisors \(D_i,D_j\) therefore span a base-point-free pencil in \(|F|\).
It defines a morphism

\[
                              f:S\longrightarrow\mathbf P^1.
\]

The class \(F\) is primitive, so the Stein factorization has connected
fibres.  Moreover

\[
 p_a(F)=1+\frac{F^2+K_S.F}{2}=1.
\]

Generic smoothness over \(\mathbf C\) makes \(f\) a genus-one fibration.

Each \(D_i\) is a reduced fibre consisting of two smooth rational curves
with intersection number \(2\).  They meet in at most two distinct points,
so

\[
                         e_{\rm top}(D_i)\ge2.          \tag{4.2}
\]

The eight displayed fibres alone therefore contribute at least \(16\) to
the Euler number.  But Noether's formula on a rational surface gives

\[
                 e_{\rm top}(S)=c_2(S)=12-K_S^2=12.    \tag{4.3}
\]

All other singular genus-one fibres have nonnegative Euler contribution.
Equations (4.2)--(4.3) contradict one another.  \(\square\)

This is an incidence/effectivity obstruction, not a lattice or
Riemann--Roch obstruction: each wall separately has the correct numerical
adjunction data.

## 5. The exact missing reflection

The packet is not reflection-maximal.  Take

\[
\begin{aligned}
r
 &=E_{00111100}-E_{00000000}\\
 &=F-\frac12(A_3+A_4+A_5+A_6)\\
 &=H-E_5-E_7-E_9.                                     \tag{5.1}
\end{aligned}
\]

Then

\[
                         r^2=-2,\qquad F.r=0,           \tag{5.2}
\]

so \(s_r(x)=x+(x.r)r\) is an integral reflection of \(L\).  More
decisively,

\[
                              h.r=0.                   \tag{5.3}
\]

Thus the mirror \(r^\perp\) contains the strict interior point (2.2).
A reflection mirror belonging to the Coxeter group \(W(P)\) cannot meet
the interior of its fundamental chamber.  Hence \(s_r\notin W(P)\).
Because \(O^+(L)\) is discrete,

\[
                         \langle W(P),s_r\rangle
\]

is a proper discrete reflection overgroup.  This proves nonmaximality
without identifying or classifying the overgroup.

There is no legitimate host repair:

1. the wall vectors already span \(L\) unimodularly by (1.10);
2. the ordinary integral host is therefore \(L\);
3. the weak-log doubled host has the same underlying abelian group, and
   scaling the form does not change its integral isometry group; and
4. \(s_r\) preserves both hosts.

Restricting maximality to norm-\(2\) roots of the doubled form would simply
discard valid square \(-2\) curve reflections.  It is not the maximality
notion needed for exact nef return.

In horospherical coordinates at the cusp, put

\[
                 t=x.F,\qquad y_i=x.A_i.
\]

The sixteen cusp walls give the box

\[
                              0\le y_i\le t,
\]

while (5.1) cuts it by

\[
                         t-\frac12(y_3+y_4+y_5+y_6)=0. \tag{5.4}
\]

The barycentre \(y_i=t/2\) lies on (5.4).  The missing reflection is exactly
a Hamming glue reflection which subdivides the \(A_1^8\) cusp box.

## 6. What maximality really adds: cusp-root saturation

### Proposition 6.1 (host-reflection saturation)

Let \(M\) be an integral Lorentzian host and \(P\) a finite-volume Coxeter
fundamental chamber whose reflection group \(W(P)\) is maximal among
discrete reflection subgroups of \(O^+(M)\).  If a negative vector
\(r\in M\otimes\mathbf Q\) has reflection \(s_r\in O^+(M)\), then

\[
                            r^\perp\cap\operatorname{int}P=\varnothing.
\tag{6.1}
\]

At every ideal cusp \(F\), the same statement says that no host-integral
reflecting hyperplane through \(F\) may cut the relative interior of a
horospherical cusp link.

### Proof

If the intersection in (6.1) were nonempty, \(s_r\) could not belong to
\(W(P)\), because mirrors of a Coxeter group do not cross a fundamental
chamber interior.  But \(O^+(M)\) is discrete, so
\(\langle W(P),s_r\rangle\) would be a proper discrete reflection
overgroup, contrary to maximality.  The cusp statement follows by taking a
horosphere sufficiently far into the cusp.  \(\square\)

This gives a genuine extra P-first incidence gate:

> a reflection-maximal packet must be saturated not merely in its displayed
> wall span, but against every reflecting vector of the adapted host,
> including glue roots visible only on cusp links.

For the Hamming packet, (5.1) is an explicit failure certificate.
Proposition 6.1 is necessary, but it does not produce a surface or prove
simultaneous curve effectivity.

### Corollary 6.2 (anticanonical cusp gate)

If, in addition, \(F=-K\) is a primitive ideal cusp of \(P\), then the
projected parabolic wall complex in

\[
                         F^\perp/\mathbf ZF
\]

is saturated against every negative class whose affine reflection preserves
the adapted host.  Thus a box of anticanonical decompositions
\(F=A_i+\beta_i\) is incompatible with relative maximality whenever a glue
reflection cuts the box interior.  This is the precise extra incidence
condition forced by maximality in the present \(K\)-escaping lane.

## 7. Scoped conclusion

What is now proved is:

\[
\boxed{
K\notin P+\text{finite volume}+\text{physical adjunction colors}
\ \not\Longrightarrow\ \text{simultaneous effectivity}.
}
\]

What is also proved is that this particular obstruction cannot be promoted
to a maximal counterexample: the root (5.1) gives a proper integral
reflection overgroup.

The following remain open:

- a nonrealizable packet satisfying both relative reflection maximality and
  \(K\notin P\);
- a theorem that those two conditions alone force a rational carrier;
- simultaneous effectivity from raw lattice/chamber data; and
- universal or exhaustive P-first surface construction.

The exact replay is
[scripts/certify_rational_pfirst_k_escape_hamming_effectivity_no_go.py](../scripts/certify_rational_pfirst_k_escape_hamming_effectivity_no_go.py).
