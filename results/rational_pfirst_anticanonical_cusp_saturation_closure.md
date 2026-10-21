# Reflection saturation closes the Hamming box into a realizable flex-pencil chamber

**Status (2026-08-10): proved exact rank-ten anticanonical-cusp closure,
relative reflection maximality, and an explicit rational realization.
Non-anticanonical cusps and other ranks remain OPEN.**

The 32-wall Hamming packet is not a maximal nonrealizability example.
Saturating its natural host-reflection arrangement replaces the eight
incompatible anticanonical decompositions by one affine nine-component
fibre.  The resulting ten-facet simplex:

1. is maximal among reflection groups inside \(O^+(I_{1,9})\);
2. has the physical canonical colors \(1\times(-1)+9\times(-2)\);
3. uses Euler budget \(10\le12\); and
4. is the full nef chamber of an explicit nine-fold blow-up of a flex cubic
   pencil.

Everything below is proved from the displayed lattice, an elementary
reflection descent, and a direct blow-up construction.  No reflection,
elliptic-surface, or chamber classification is invoked.

## 1. Exact alcove selected inside the Hamming chamber

Keep

\[
L=I_{1,9},\qquad
K=(-3,1,\ldots,1),\qquad
F=-K=(3,-1,\ldots,-1)
\]

and the Hamming strict point

\[
h=(5,-2,-1,-2,-1,-2,-1,-2,-1,-1).
\]

The missing Hamming mirrors pass through \(h\), so choose the deterministic
generic perturbation

\[
\begin{aligned}
v&=(0,1,2,3,4,5,6,7,8,9),\\
g&=100h+v.
\end{aligned}                                           \tag{1.1}
\]

It has square \(51915\), and its minimum pairing with the 32 old walls is
\(91\).  Thus \(g\) is still strictly inside the Hamming chamber.

Signed permutations and reflection in
\(H-E_1-E_2-E_3\) reduce \(g\), in ten exact steps, to

\[
(409,-193,-108,-106,-104,-98,-96,-94,-92,-91).         \tag{1.2}
\]

Pulling the standard simple roots back through those ten reflections gives
the following alcove walls:

\[
\begin{array}{lll}
r_0=E_9,\\
r_1=-H+E_3+E_5+E_7,&
r_2=E_1-E_3,&
r_3=E_3-E_5,\\
r_4=H-E_1-E_2-E_3,&
r_5=E_2-E_4,&
r_6=E_4-E_6,\\
r_7=E_6-E_8,&
r_8=E_8-E_9,&
r_9=E_5-E_7.
\end{array}                                             \tag{1.3}
\]

Their exact colors are

\[
\begin{array}{c|cc}
&r_i^2&K.r_i\\ \hline
i=0&-1&-1\\
1\le i\le9&-2&0.
\end{array}                                             \tag{1.4}
\]

They form a determinant-\(-1\) basis of \(L\).  Every old Hamming wall is
a nonnegative integral combination of the \(r_i\); only
\(E_{00000000}=r_0\) remains a simple facet.  Hence the simplex (1.3) is
contained in the old chamber, while the other 31 old mirrors lie beyond
its walls or inside its positive-root arrangement.  The archived exact
expansion table makes this containment independently replayable.

The point \(h\) lies on the common codimension-eight face

\[
 h.r_0=h.r_1=1,\qquad h.r_i=0\quad(2\le i\le9).         \tag{1.5}
\]

It is not an alcove interior point after saturation.

## 2. The saturated simplex

The Gram matrix of (1.3) has:

- one edge of order \(4\), between \(r_0\) and \(r_8\);
- eight order-\(3\) edges among \(r_1,\ldots,r_9\); and
- 36 orthogonal pairs.

The nine square-\(-2\) nodes form a connected tree with eight edges.
The ten dual vertex rays have squares

\[
                   0,0,2,6,5,4,3,2,1,1.               \tag{2.1}
\]

Thus every vertex is finite or ideal and the Coxeter simplex has finite
volume.  Its two primitive ideal vertices are

\[
\begin{aligned}
F&=(3,-1,-1,-1,-1,-1,-1,-1,-1,-1),\\
Q&=(2,-1,0,-1,0,-1,0,-1,0,0).
\end{aligned}                                           \tag{2.2}
\]

The Hamming point is exactly

\[
                              h=F+Q.                    \tag{2.3}
\]

Only \(F\) is characteristic: all its coordinates are odd, whereas \(Q\)
is not characteristic.

The primitive null relation on the nine cusp walls is

\[
F=
2r_1+4r_2+6r_3+5r_4+4r_5+3r_6+2r_7+r_8+3r_9.          \tag{2.4}
\]

If the \(r_i\) are effective curves, the reduced support of (2.4) consists
of nine rational components meeting along a tree with eight transverse
edges.  Its Euler number is therefore

\[
                              2\cdot9-8=10.             \tag{2.5}
\]

This is the precise replacement for the Hamming box's eight two-component
fibres, whose impossible minimum contribution was \(16\).

## 3. Maximality without a reflection classification

### Theorem 3.1

The side-reflection group of the simplex (1.3) contains every integral
hyperbolic reflection of \(L=I_{1,9}\).  Consequently it has no proper
reflection overgroup inside \(O^+(L)\).

### Proof

Because \(L\) is unimodular, a primitive vector \(u\) whose reflection is
integral has divisibility one.  Integrality of

\[
                         s_u(x)=x-\frac{2(x.u)}{u^2}u
\]

therefore forces \(u^2\mid2\).  A reflecting vector has negative square, so
after primitive normalization

\[
                              u^2=-1\quad\text{or}\quad-2. \tag{3.1}
\]

It is enough to work with the standard isometric simple system

\[
E_9,\quad E_i-E_{i+1}\ (1\le i\le8),\quad
H-E_1-E_2-E_3.                                         \tag{3.2}
\]

The first nine reflections generate every signed permutation of the
\(E_i\).  Up to sign and such a permutation, write any vector in (3.1) as

\[
                     u=(a,-m_1,\ldots,-m_9),
\quad
a\ge0,\quad m_1\ge\cdots\ge m_9\ge0.                   \tag{3.3}
\]

Put \(q=-u^2\in\{1,2\}\) and \(s=m_1+m_2+m_3\).  If \(a>0\) and
\(s\le a\), then

\[
\begin{aligned}
\sum_{i=1}^9m_i^2
&\le m_1^2+m_2^2+7m_3^2\\
&\le(m_1+m_2+m_3)^2\\
&\le a^2,
\end{aligned}                                           \tag{3.4}
\]

contrary to \(\sum m_i^2=a^2+q\).  The middle inequality follows by
expanding and using \(m_1,m_2\ge m_3\).

Thus \(s>a\).  Also each \(m_i\le a\), unless \(a=0\), because
\((a+1)^2>a^2+2\).  Hence \(s\le3a\).  If \(s=3a\), equality in the norm
equation forces

\[
q=2,\quad a=m_1=m_2=m_3=1,\quad m_4=\cdots=m_9=0,
\]

which is already the last root in (3.2).  Otherwise reflection in
\(H-E_1-E_2-E_3\) changes the time coefficient to

\[
                              a'=2a-s,
\]

with \(|a'|<a\).  Change the overall sign if necessary and repeat.  The
positive integer \(a\) strictly descends.

At \(a=0\), the norm equation leaves one unit \(m_i\) when \(q=1\), or two
unit \(m_i\)'s when \(q=2\).  Signed permutations carry these respectively
to \(E_9\) and \(E_8-E_9\).  Thus every integral reflection is generated by
(3.2), and therefore by its pullback (1.3).  \(\square\)

This proves the exact relative maximality required by the P-first return
argument.  It is stronger than bounded root search or diagram saturation.

## 4. Anticanonical-cusp rigidity in rank ten

### Theorem 4.1

Let \(P'\) be a finite-volume integral Coxeter chamber in
\(\mathbb P(I_{1,9}\otimes\mathbf R)\).  Assume:

1. its side-reflection group is maximal among reflection subgroups of
   \(O^+(I_{1,9})\); and
2. the primitive characteristic vector \(F=-K\) is an ideal cusp of \(P'\).

Then \(P'\) is conjugate to the simplex (1.3) by an integral reflection
isometry fixing \(F\).  In particular, its anticanonical cusp link has the
nine-wall relation (2.4), and its reduced affine support has Euler number
\(10\).

### Proof

Adjoining any integral reflection to the side group remains discrete inside
\(O^+(L)\).  Relative maximality and Theorem 3.1 therefore identify the side
group with the full integral reflection group.  Its fundamental chambers
are translates of (1.3).

Under such a translation, the cusp \(F\) must map to one of the two ideal
rays (2.2).  Isometries preserve characteristic vectors.  Since \(F\) is
characteristic and \(Q\) is not, the image must be the primitive ray \(F\)
itself.  The translating reflection isometry fixes \(F\), proving the
claim.  \(\square\)

The physical coloring now follows automatically:

\[
K.r_0=-1,\qquad K.r_i=0\ (i\ge1).
\]

A twice-primitive physical \((-4)\) wall \(B\) would have \(K.B=2\), hence

\[
                              F.B=-2<0.
\]

It cannot occur in a chamber containing the cusp \(F\).  Thus the
rank-ten anticanonical-cusp lane allows exactly one \((-1)\) facet and nine
\((-2)\) facets after full host saturation.

This is the promised cusp-link budget theorem: in this lane, reflection
saturation forces the affine fibre incidence to consume \(10\), not \(16\),
of the rational elliptic \(c_2=12\) budget.

## 5. Direct rational realization by a flex pencil

Consider the cubic pencil

\[
        Y^2Z=X^3+XZ^2+\lambda Z^3                    \tag{5.1}
\]

on \(\mathbf P^2\).  The point \(O=[0:1:0]\) is a flex of every finite
member and \(Z=0\) is its flex line.  The member at
\(\lambda=\infty\) is the triple line \(3(Z=0)\).

In local coordinates \(x=X/Y,\ z=Z/Y\) at \(O\), the pencil ideal is
generated by

\[
                         z^3,\qquad z-x^3-xz^2.         \tag{5.2}
\]

The first three blow-ups follow the flex direction.  Direct substitution
reduces (5.2), up to local units and a regular change of the second
coordinate, to

\[
                              (x^6,w).                  \tag{5.3}
\]

Six further successive point blow-ups resolve (5.3).  Thus exactly nine
blow-ups make the pencil base-point-free.

With the chronological total-transform marking, the reduced components of
the triple-line fibre are

\[
 E_i-E_{i+1}\ (1\le i\le8),\qquad
 H-E_1-E_2-E_3,                                        \tag{5.4}
\]

and the last exceptional curve \(E_9\) is a section.  Their multiplicities
are exactly those in (2.4).

The finite discriminant of (5.1), up to a nonzero scalar, is

\[
                              4+27\lambda^2.            \tag{5.5}
\]

It has two simple roots over \(\mathbf C\).  The corresponding cubics are
irreducible nodal curves, so (5.4) is the only reducible fibre.

Let \(S\) be the resulting blow-up.  It is rational, \(K_S=-F\), and
\(F\) is the fibre class.  The lattice generated by

- \(F\);
- the section \(E_9\); and
- the eight fibre components not meeting \(E_9\)

has coordinate determinant \(-1\).  It is therefore all of
\(\operatorname{NS}(S)\).  The section quotient is trivial, so \(E_9\) is
the unique section.

Finally, let \(C\) be any irreducible curve and put \(d=F.C\).  Since
\(K_S=-F\), adjunction gives

\[
                         C^2=2p_a(C)-2+d.               \tag{5.6}
\]

If \(C\) is horizontal and negative, (5.6) forces \(d=1\) and
\(p_a(C)=0\), so \(C\) is a section and hence equals \(E_9\).  If \(C\) is
vertical and negative, it is one of the nine displayed fibre components.
Thus these ten curves are all the negative curves of \(S\).

The positive-cone/Lorentz argument now shows that their ten half-spaces give
the entire nef cone.  Consequently the standard simplex, and by the
reflection marking also the pulled-back simplex (1.3), is genuinely
realized:

\[
                        P_{\rm sat}=\mathbf P(\operatorname{Nef}S).
\tag{5.7}
\]

The Euler accounting is exact:

\[
                  10+1+1=12=c_2(S),
\]

from the affine fibre and the two nodal fibres.

## 6. Consequence for the P-first program

For the Hamming packet, host-reflection saturation has a definitive
outcome:

\[
\boxed{\text{the obstruction collapses to a maximal, realizable
anticanonical chamber}.}
\]

The eight impossible decompositions do not survive as facets.  They become
non-simple positive roots after the Hamming cusp box is subdivided.

The stronger theorem proved here is rank- and cusp-specific:

> in \(I_{1,9}\), a relatively reflection-maximal chamber with primitive
> characteristic cusp \(F=-K\) has one section wall, one affine
> nine-component fibre support, and Euler contribution \(10\le12\).

Still open are:

- maximal \(K\)-escaping chambers for which \(-K\) is not a cusp;
- analogous cusp-budget theorems in other ranks or nonunimodular adapted
  hosts;
- automatic effectivity without an explicit blow-up or incidence
  construction; and
- any universal surface-existence theorem from bare maximality and
  \(K\notin P\).

The exact replay is
[scripts/certify_rational_pfirst_anticanonical_cusp_saturation_closure.py](../scripts/certify_rational_pfirst_anticanonical_cusp_saturation_closure.py).
