# Finite rational vertices: the exact discrepancy barrier and the exceptional-tail obstruction

**Status (2026-08-10): proved, carrier-relative.** This note isolates what
finite-vertex surface singularities can and cannot contribute to the desired
dimension theorem. It uses only adjunction, exact intersection matrices, and
the elementary inverse-positivity of a negative-definite curve matrix, plus
the standard nef-big contraction criterion for the explicit tail walls. It
does not use Nikulin--Vinberg, ABSW, arithmetic volume, or a classification
of finite Coxeter diagrams. All surfaces and contractions below are complex
and projective.

The positive result is uniform. Given the normal \(\mathbb Q\)-Gorenstein
contraction of an SNC finite-vertex packet away from the three hostile
incidences H1, H2, and H3, the contracted points are klt, with every log
discrepancy at least \(1/2\). The negative result is equally important: this
local bound does **not** bound the dimension. Arbitrarily long chains of
exceptional point blowups occur as actual finite vertices of nef cones, have
strictly positive discrepancies, and have finite crystallographic local
reflection type. A global theorem must therefore control the numerical orbit
or contraction of this exceptional tail; another local singularity estimate
cannot finish the proof.

## 1. Physical packet and discrepancy system

Let \(S\) be a smooth rational surface and let
\(\mathcal C=\{C_1,\ldots,C_r\}\) be an SNC collection of smooth rational
curves with negative-definite intersection matrix \(Q=(C_i.C_j)\). In the
application, these are the physical facets through one finite vertex of a
literal full-nef Coxeter chamber. Each curve has one of the three types

\[
 B^2=-4,\qquad R^2=-2,\qquad E^2=-1.                 \tag{1.1}
\]

Adjunction gives

\[
 K_S.B=2,\qquad K_S.R=0,\qquad K_S.E=-1.             \tag{1.2}
\]

Suppose there is a projective birational morphism \(f:S\to Y\) to a normal
\(\mathbb Q\)-Gorenstein surface whose exceptional locus is exactly
\(\bigcup_iC_i\), and that \(f\) is an isomorphism off this locus. Thus \(S\)
is a log resolution of the contracted points. Write

\[
 K_S=f^*K_Y+\sum_i a_iC_i.                            \tag{1.3}
\]

The discrepancy vector is the unique exact solution

\[
 Qa=k,\qquad k_i=K_S.C_i.                             \tag{1.4}
\]

Put \(A=-Q\). Because \(Q\) is negative definite and distinct effective
curves have nonnegative intersection, \(A\) is a positive-definite
M-matrix. Hence every entry of \(A^{-1}\) is nonnegative.

For completeness, this last sign is elementary. If \(Ax=y\geq0\), write
\(x=w-z\), where \(w,z\geq0\) have disjoint supports and \(z\) is the
negative part of \(x\). If \(z\ne0\), then

\[
 0\leq z^Ty=z^TAw-z^TAz<0:
\]

the first term on the right is nonpositive because the off-diagonal entries
of \(A\) are nonpositive, while the second is strictly positive by positive
definiteness. This contradiction gives \(x\geq0\), hence
\(A^{-1}\geq0\) entrywise.

## 2. One barrier vector proves the uniform klt theorem

Assume the weak-incidence conditions on this packet, for every pair of
distinct square-minus-four curves and every displayed \(R\) and \(E\):

\[
 B_i.B_j=0,\qquad B_i.R=0,\qquad \sum_iB_i.E\leq2.    \tag{2.1}
\]

Define a vector \(u\) by

\[
 u_B=\frac12,\qquad u_R=u_E=0.                        \tag{2.2}
\]

There is no diagram classification hidden here. Row-by-row substitution
gives

\[
 (Au-k)_B=0,\qquad (Au-k)_R=0,\qquad
 (Au-k)_E=1-\frac12\sum_iB_i.E\geq0.                  \tag{2.3}
\]

Since \(Aa=-k\), equations (2.2)--(2.3) imply

\[
 A(a+u)=Au-k\geq0.
\]

Multiplication by the entrywise-nonnegative matrix \(A^{-1}\) proves

\[
                  a_B\geq-\frac12,\qquad
                  a_R\geq0,\qquad a_E\geq0.           \tag{2.4}
\]

Thus the contracted points of \(Y\) are klt. Put \(A_i^{\log}=1+a_i\).
On this resolution every log discrepancy is at least \(1/2\), and this also
controls every divisor on every higher resolution. Indeed, blowing up a point
which lies on \(s=0,1,2\) components of an SNC exceptional divisor creates a
divisor \(F\) with

\[
 a_F=1+\sum_{p\in C_i}a_i,
 \qquad
 A_F^{\log}=2+\sum_{p\in C_i}a_i.
\]

For \(s=0\) this is \(2\); for \(s=1\) it is
\(A_i^{\log}+1\); and for \(s=2\) it is
\(A_i^{\log}+A_j^{\log}\). None is smaller than \(1/2\), and point blowups
preserve the SNC condition, so induction covers every divisorial valuation
over a surface. Under the strict convention for \(\epsilon\)-klt, the class
is therefore uniformly \(\epsilon\)-klt for every fixed
\(0<\epsilon<\frac12\), but not uniformly \(1/2\)-klt. Equality occurs
already when the negative section of \(\mathbb F_4\) is contracted, so the
constant cannot be increased in this normalization.

This is the local AG content behind the global weak-incidence orbifold
theorem. It is stronger than checking a finite catalogue: it applies to
every negative-definite physical packet of arbitrary rank satisfying
(2.1).

## 3. H1--H3 are the sharp failures of this barrier

Within the SNC and \(\mathbb Q\)-Gorenstein contraction scope of Section 1,
the three minimal hostile packets give exact boundary values:

| failure | physical matrix | discrepancy vector |
|---|---|---|
| H1: \(B_1.B_2=2\) | \([-4,2;2,-4]\) | \((-1,-1)\) |
| H2: \(B.R=2\) | \([-4,2;2,-2]\) | \((-1,-1)\) |
| H3: \(E.B_1=E.B_2=E.B_3=1\) | central \((-1)\), three \((-4)\) leaves | \((-2,-1,-1,-1)\) |

H1 and H2 are log canonical but not klt; here intersection number two must be
realized by two distinct transverse intersection points in order to remain
SNC. H3 is worse than log canonical. These are the minimal numerical points
where the three rows in (2.3) cease to be nonnegative. The statement is about
the physical intersection failures, not merely the abstract support graph.

The SNC qualification is essential. If the intersection number two in H1 or
H2 is instead an ordinary tangency, the same matrix still gives discrepancies
\((-1,-1)\) on \(S\), but \(S\) is not a log resolution. Blowing up the
tangency produces a divisor of discrepancy \(1-1-1=-1\) and a triple point;
blowing up that triple point produces discrepancy
\(1-1-1-1=-2\). That non-SNC realization is not log canonical. Thus the
matrix alone cannot justify the lc label outside the stated hypotheses.

## 4. Why the local theorem does not bound dimension

Take \(X_0=\mathbb P^2\), blow up a point \(p_1\), choose \(p_2\) anywhere
on the first exceptional curve, and, for \(j\geq3\), choose \(p_j\) on the
newest exceptional curve but away from its intersection with the preceding
exceptional curve. Write
\(\pi:X_m\to\mathbb P^2\) for the composite. The exceptional locus is the
SNC chain

\[
 R_1-\cdots-R_{m-1}-E_m,\qquad R_i^2=-2,\quad E_m^2=-1. \tag{4.1}
\]

Its physical intersection matrix is negative definite and unimodular. Its
leading principal determinants have the required alternating signs; the
full determinant is \((-1)^m\). Solving (1.4) gives

\[
                         (a_1,\ldots,a_m)=(1,2,\ldots,m). \tag{4.2}
\]

These packets are actual finite vertices of nef cones, not just formal
matrices. If \(L=\pi^*\mathcal O_{\mathbb P^2}(1)\), then \(L\) is nef,
\(L^2=1\), and \(L.C=0\) for an irreducible curve \(C\) exactly when \(C\)
is one of the exceptional components. The matrix in (4.1) is nonsingular and
\(\rho(X_m)=m+1\), so its \(m\) orthogonality equations cut out the positive
ray \(\mathbb R_{>0}L\). Each component is a genuine wall: if \(C^2<0\)
and \(H\) is ample, then

\[
 H_C=H+\frac{H.C}{-C^2}C
\]

is nef, has \(H_C.C=0\), and has positive intersection with every other
irreducible curve. Moreover
\(H_C^2=H^2+(H.C)^2/(-C^2)>0\). After clearing denominators,
\(H_C|_C\) is the trivial degree-zero line bundle on \(C\simeq\mathbb P^1\),
so the standard nef-big surface contraction criterion makes \(H_C\)
semiample and contracts exactly \(C\). Hence this is a genuine facet, and
\([L]\) is a finite vertex of the projectivized nef cone.

In the doubled physical lattice \(q=2(\,\cdot\,)\), reflection in every
component is integral: the reflection coefficient is \(-x.R\) for a
square-minus-two component and \(-2x.E\) for the square-minus-one end.
Adjacent \(R\)-walls have Coxeter order three and the final \(R\)-\(E\) pair
has order four. Hence for \(m\geq2\) the local root diagram is the finite
crystallographic chain with one double edge (and for \(m=1\) it is \(A_1\)).
This says nothing about whether the *entire* nef cone is a Coxeter chamber.

The discrepancies in (4.2) are increasingly positive. Moreover
\(a^TQa=a^Tk=-m\), so the canonical square drops by \(m\) under the blowups
and the discrepancy correction \(-a^TQa\) restores exactly \(m\).
Contraction returns the same canonical square. Nothing in the local
singularity, Coxeter-angle, or bounded-discrepancy packet detects the length
\(m\).

Therefore the implication

\[
 \text{uniformly klt finite vertices}\Longrightarrow
 \text{uniform dimension bound}
\]

is false. This is not a counterexample to a global dimension theorem: the
construction does not assert that the completed full-nef chamber is
finite-volume or reflection-maximal. It proves that those global hypotheses
must do indispensable work.

## 5. Exact remaining bridge in the H-free lane

For the project-native route, the missing assertion is now narrow:

> **Exceptional-tail bridge.** In a literal finite-volume full-nef carrier,
> reflection maximality must force a uniformly bounded exceptional tail, or
> force a finite numerical operator orbit of a geometric contraction of that
> tail.

The second alternative is designed to feed the already proved equivariant
minimalization theorem: a finite numerical orbit makes the full trace block
orthogonal to the exceptional lattice and descends it to a minimal carrier.
Bare local klt singularities do not establish this orbit condition. Bare
reflection maximality also cannot be declared to establish it without a new
argument; redundant blowups of a fixed lower-dimensional carrier show why
the distinction between a full-dimensional literal carrier and an external
cohomological subpacket matters.

This note proves no log-MMP core bound in the hostile lanes. It proves only
that the SNC H1 and H2 prototypes are lc but not klt, that the minimal H3 star
is non-lc, and that in the H-free lane the local discrepancy theorem cannot
control exceptional-tail length. H1--H3 therefore still require separate
global arguments. For the H-free lane, the bridge above, rather than another
census of vertex diagrams, is the next theorem required for a global
ABSW-type conclusion.

## 6. Reproducibility

The exact row reductions, the three sharp hostile packets, and the first
thirty-two tail lengths are frozen in
scripts/certify_rational_finite_vertex_discrepancy_tail_gate.py and
data/enumerations/rational_finite_vertex_discrepancy_tail_gate.json.
The arbitrary-length formula (4.2) is proved symbolically above; the finite
replay is a regression certificate, not the logical source of the theorem.
