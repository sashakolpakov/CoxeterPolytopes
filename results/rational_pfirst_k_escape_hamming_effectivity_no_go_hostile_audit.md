# Hostile audit: the \(I_{1,9}\) Hamming effectivity obstruction

**Audit status (2026-08-10): PASS with the scope stated in the source.
The packet passes canonical escape and fails simultaneous effectivity, but
it is explicitly not reflection-maximal.**

This audit checks
results/rational_pfirst_k_escape_hamming_effectivity_no_go.md independently
of any surface or chamber classification.

## 1. Lattice and wall checks

The characteristic vector

\[
K=(-3,1,\ldots,1)
\]

has square zero in \(I_{1,9}\).  Its negative
\(F=(3,-1,\ldots,-1)\) is primitive isotropic.  The eight displayed
\(A_i\) have Gram \(-2I_8\) and are orthogonal to \(F\).  Therefore
\(\beta_i=F-A_i\) satisfy

\[
\beta_i^2=-2,\quad A_i.\beta_i=2,
\]

with every cross-pairing for distinct indices equal to zero.

The four binary generators give exactly 16 codewords, with weight
enumerator

\[
                         1+14z^4+z^8.
\]

Substitution in

\[
E_c=E_9-\frac12\sum c_iA_i+\frac{|c|}{4}F
\]

returns integral vectors.  Exact dot products verify

\[
E_c^2=-1,\quad F.E_c=1,\quad
E_c.A_i=c_i,\quad E_c.\beta_i=1-c_i.
\]

Thus all adjunction colors and every reflection formula are integral.
A displayed ten-wall minor has determinant \(-1\), so the walls span the
whole host \(L\).

## 2. Coxeter and finite-volume checks

The strict vector

\[
h=(5,-2,-1,-2,-1,-2,-1,-2,-1,-1)
\]

has square \(4\) and pairing \(1\) with every wall.  All off-diagonal wall
pairings give precisely Coxeter orders \(2,4,\infty\).

Fresh exact double description returns 2306 extreme rays.  After primitive
integral normalization their square histogram is

\[
0^{226},\quad1^{1040},\quad2^{1040}.
\]

Every ray has positive time coordinate, so all lie in the same closed
future Lorentz cone.  Every one of the 32 inequalities is active on 670
rays spanning a rank-nine hyperplane.  Hence no listed wall is redundant,
there is no hyperideal extreme direction, and the finite-sided chamber has
finite hyperbolic volume.  Sorting makes the ray hash independent of the
double-description output order.

These checks establish that the side-reflection group is discrete with the
displayed polyhedron as a fundamental chamber.  They do not establish
maximality.

## 3. Canonical escape

Exact pairings give

\[
F.A_i=F.\beta_i=0,\qquad F.E_c=1.
\]

Thus \(F=-K\) is an ideal cusp.  Since \(K.E_c=-1\), \(K\) violates sixteen
wall inequalities and is not in the chamber.  There is no orientation
ambiguity because the strict future vector has all positive wall slacks.

## 4. Effectivity contradiction

Assume the sixteen square \(-2\) walls are irreducible curves on a marked
smooth rational surface over \(\mathbf C\).  Adjunction gives arithmetic
genus zero, so they are smooth rational curves.

Each

\[
D_i=A_i+\beta_i
\]

is an effective divisor of class \(F=-K_S\).  Different \(D_i\) have no
common component and total intersection zero.  Since local intersection
multiplicities of distinct irreducible curves are nonnegative, their
supports are actually disjoint.

On a rational surface \(\operatorname{Pic}^0=0\), so the \(D_i\) are
linearly equivalent.  Two disjoint members give two sections with no common
zero and hence a base-point-free pencil.  Primitivity of \(F\) rules out a
nontrivial disconnected Stein factor, and adjunction gives \(p_a(F)=1\).
Over \(\mathbf C\), generic smoothness therefore gives a genus-one
fibration.

Every displayed \(D_i\) is a reduced fibre made from two rational curves
whose intersection number is two.  Its topological Euler number is at least
two: the components meet at either one or two underlying points.  Singular
fibres of a genus-one fibration have nonnegative Euler number.  Hence eight
such fibres force

\[
e(S)\ge16.
\]

Noether's formula gives \(e(S)=12-K_S^2=12\).  This is a genuine
contradiction.  It uses simultaneous effectivity and incidence; it does
not arise from individual Riemann--Roch calculations.

## 5. Nonmaximality and host scope

For \(c=00111100\),

\[
r=E_c-E_0=H-E_5-E_7-E_9
\]

has square \(-2\), is orthogonal to \(F\), and has integral reflection on
\(L\).  Directly,

\[
                              h.r=0.
\]

The missing mirror therefore crosses the strict chamber interior.  It
cannot be a mirror of the side-reflection group, since mirrors of a Coxeter
group avoid its fundamental-chamber interior.  Adjoining it inside the
discrete group \(O^+(L)\) gives a proper discrete reflection overgroup.

The determinant-one wall subset closes the obvious host loophole.  The
minimal wall-generated host is already the unimodular lattice \(L\), which
has no proper integral overlattice for the original form.  The usual
weak-log doubled host has the same underlying abelian group and the same
isometries.  The missing reflection survives both.

This does not rule out an artificially chosen nonstandard host with a
different maximality problem.  Such a choice would not certify the
all-integral-reflection maximality used by the P-first exact-return theorem.

## 6. Saturation lemma

The source's host-reflection-saturation proposition is formal and sound.
If an integral reflecting mirror crossed the interior of a relatively
reflection-maximal chamber, adjoining that reflection would give a proper
discrete reflection overgroup.  At a cusp, an affine reflecting hyperplane
which cuts the relative interior of a sufficiently deep horospherical link
also cuts the hyperbolic chamber interior.

For this packet, the cusp coordinates obey

\[
0\le y_i\le t,
\qquad
r^\perp:\quad t-\frac12(y_3+y_4+y_5+y_6)=0.
\]

The second hyperplane contains the box barycentre.  Hence the \(A_1^8\)
Hamming cusp is visibly unsaturated.

For a primitive anticanonical cusp \(F=-K\), the quotient
\(F^\perp/\mathbf ZF\) is exactly the finite Euclidean lattice on which this
local saturation must be checked.  The source therefore obtains a real
maximality consequence, but only in the anticanonical-cusp lane; bare
\(K\notin P\) does not by itself assert that \(-K\) is a cusp.

## 7. Final audit boundary

The following implications are proved:

\[
\begin{aligned}
&K\notin P+\text{physical integral walls}
   \not\Rightarrow \text{simultaneous effectivity},\\
&\text{relative reflection maximality}
   \Rightarrow \text{host-reflection-saturated cusp links}.
\end{aligned}
\]

The audit does not claim:

- a maximal plus \(K\)-escape counterexample;
- that cusp saturation implies effectivity;
- that maximality constructs a rational surface; or
- exhaustiveness of any P-first compiler.

Those questions remain OPEN.
