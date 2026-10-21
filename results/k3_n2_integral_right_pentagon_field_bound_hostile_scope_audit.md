# Hostile audit: the general K3 right-pentagon field bound

**Verdict (2026-08-11): PASS, general right-pentagon scope.**  The theorem
removes both the power coupling and any imposed mixed-height bound.  Its
automatic arithmetic input is exactly the integrality of the five facet
reflections on a finite projective returned order-module.  The exact K3
same-chamber, completeness, saturation, and declared-host reflection
maximality hypotheses remain present.  Maximality is retained but unused:
the five Gram identities and signed-place norm estimates already close this
incidence type.

## 1. Normalization and endpoint audit

For each nonadjacent facet pair, the squared normalized correlation

\[
       R_{ij}={h(e_i,e_j)^2\over h(e_i,e_i)h(e_j,e_j)}
\]

is an element of \(k\), without choosing square roots of the normal lengths
inside \(k\).  Compact-right-pentagon ultraparallelity gives
\(\sigma_0(R_{ij})>1\).  At the definite embedding, strict
Cauchy--Schwarz gives \(0<\sigma_1(R_{ij})<1\).

Both inequalities are strict.  The pairing is a nonzero field element at
the distinguished place, so it cannot vanish after applying another field
embedding.  The two normals have nonzero wedge over \(k\), so they cannot
become proportional at the definite place.  Thus, for every shifted
coordinate \(x_i=R_i-1\),

\[
                     (x_i)_0>0,\qquad -1<(x_i)_1<0.
\tag{1.1}
\]

## 2. Five-minor and recurrence audit

The exact principal four-by-four Gram determinants are

\[
\begin{aligned}
 CE-C-D-E+1,&\quad -A+BE-B-E+1,\\
 BD-B-C-D+1,&\quad AD-A-D-E+1,\\
 AC-A-B-C+1.&
\end{aligned}
\]

They all vanish because the normals live in a rank-three field space.  With

\[
 (x_0,x_1,x_2,x_3,x_4)=(A-1,B-1,C-1,D-1,E-1),
\]

the five equations are exactly

\[
                         x_{i-1}x_{i+1}=x_i+1
                         \qquad(i\bmod5).
\tag{2.1}
\]

In particular,

\[
 x_1=AC-A-C,\qquad
 x_3={C\over AC-A-C},\qquad
 x_4={A\over AC-A-C}.
\tag{2.2}
\]

These agree with the owner's formulas after writing \(A=x_0+1\) and
\(C=x_2+1\).  Independent exact replay verifies all five minors and all
five cyclic recurrence equations.

## 3. Quarter-integrality and norm audit

For rank \(m=3\), the product of two facet reflections satisfies

\[
       \operatorname {tr}_k(s_is_j)=-1+a_{ij}a_{ji}.
\]

The left side is the module trace of an endomorphism of the finite
projective order-module.  Hence each pair product

\[
                   a_{ij}a_{ji}=4R_{ij}=4(x_i+1)
\]

belongs to \(\mathcal O\).  The automatic statement is
\(4x_i\in\mathcal O\), not \(x_i\in\mathcal O\).

By (1.1), the algebraic integer \(y_i=4x_i\) has one positive value and a
conjugate in \((-4,0)\).  Therefore its norm is a negative nonzero integer,
and

\[
 1\le|N(y_i)|=(y_i)_0|(y_i)_1|<4(y_i)_0.
\]

There is no endpoint or zero-norm exception.  Dividing by four only after
this inequality gives

\[
                                  (x_i)_0>{1\over16}.
\tag{3.1}
\]

## 4. Moving-denominator audit

Put \(z=x_{i+2}\).  Two cyclic recurrences give

\[
 x_{i+1}=x_i z-1,
 \qquad
 x_{i+3}={z+1\over x_i z-1}.
\]

The denominator is positive at the distinguished place.  Using
\(x_{i+3}>1/16\) gives

\[
 x_i<16+{17\over z}.
\]

Since \(z>1/16\), every coordinate satisfies

\[
                                  1/16<x_i<288.
\tag{4.1}
\]

The index is arbitrary, so this includes
\(x_1=AC-A-C\).  The proof neither assumes \(A\ge C\) nor loses a cyclic
case.

## 5. Pair generation and discriminant audit

The integral pair cycle \(\alpha=4A\) has

\[
                   4<\alpha_0<1156,
             \qquad 0<\alpha_1<4.
\]

Thus \(\alpha\) is nonrational: a rational element has equal values under
both embeddings.  Since \(k\) is quadratic, every pair cycle is
nonrational and generates \(k\).  No length-five-only generation gap
remains.

The integral trace and norm obey

\[
                       T\le1159,\qquad N\le4623.
\]

Consequently

\[
        \operatorname {Disc}(\mathbf Z[\alpha])
          =T^2-4N\le1159^2-4=1,343,277.
\]

The field discriminant divides this order discriminant.  This proves the
field discriminant bound \(1,343,277\), hence a genuine finite field set,
without an artificial correspondence height.

## 6. Extra-integral corollary and poison audit

If \(A,C\in\mathcal O\), then \(A-1,C-1\) have conjugates in \((-1,0)\)
and negative nonzero integral norms.  Their distinguished values exceed
one.  Combining this with the quarter-integral bound on \(D-1,E-1\) gives

\[
 A_0,C_0<34,\qquad
 |\operatorname {Disc}(k)|\le1152.
\]

The exact trace--norm enumeration contains 528 rows and 336 distinct
squarefree radicands.  These are 336 candidate squareclasses, not 336
realized fields.

The \(\mathbf Q(\sqrt5)\) poison packet is exact as a cyclic coordinate
calculation.  All five recurrences hold, every \(4x_i\) is integral, and
the directed acute five-cycle \(-408-664\phi\) is integral.  Nevertheless
its first three \(x_i\) are not algebraic integers and its first
distinguished coordinate is about \(0.309<1\).  It therefore invalidates
only the division-by-four step and the sharper \(1152/336\) deduction.  The
report correctly does not promote this formal Gram packet to a K3 lattice,
an integral Cartan factorization, or a maximal chamber.

## 7. Surviving scope

The theorem closes the compact right-angled pentagon incidence with
independent moduli.  It does not settle arbitrary \(n=2\) incidence,
noncompact chambers, or polygons with different Gram graphs.  It also does
not derive a canonical mixed K3 correspondence degree from maximality.
The remaining frontier is to force a bounded recurrence packet from a
general exact maximal eigen-nef chamber, or to control its short chordless
cycle coordinates without the pentagon identities.  No reflective-lattice
classification is used.
