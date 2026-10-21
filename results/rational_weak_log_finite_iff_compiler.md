# `RLogComp`: a finite iff/`EMPTY` compiler for maximal weak-log rational carriers

**Status (2026-08-11): proved over \(\mathbf C\), with an executable exact
route certificate.**  This theorem closes the *decision* problem for one
precise rational Property-4 species.  Starting with a finite-sided
finite-volume arithmetic reflection chamber, it constructs a finite list of
explicit constructible cluster/Hilbert loci.  Their union is nonempty if and
only if the chamber has a complex rational weak-log carrier with simultaneous
prime facet curves and literal full-nef return.  Each locus has a terminating
Gr\"obner/quantifier-elimination nonemptiness test, so the output is either a
finite carrier-stack presentation or `EMPTY`.

The result does **not** assert that the output is nonempty for every raw
maximal chamber.  In particular, it turns Property 4 into a finite algebraic
geometry calculation; it does not prejudge that calculation.  The adapted
`P10` packet shows why `EMPTY` is a necessary possible answer.

The two new ingredients are:

1. in the big branch, the single universal system \(|192D|\) is enough to
   test the otherwise circular nef premise; and
2. in the nonzero isotropic branch, no semiampleness or period test is
   needed: the affine cusp facets produce an effective multiple of \(D\),
   and this proves nefness directly.

Here

\[
 B=\sum B_i,\qquad \Delta=B/2,\qquad
 D=-(K_S+\Delta).
\]

The exact arithmetic replay is
[`certify_rational_weak_log_finite_iff_compiler.py`](../scripts/certify_rational_weak_log_finite_iff_compiler.py),
with its frozen
[`archive`](../data/enumerations/rational_weak_log_finite_iff_compiler.json)
and
[`tests`](../tests/test_rational_weak_log_finite_iff_compiler.py).

## 1. Exact target category

Let \(P\subset\mathbb H^n\), \(n\ge2\), be a finite-sided finite-volume
Coxeter chamber.  Assume that its side group is maximal among arithmetic
reflection overgroups in the real hyperbolic isometry group.  Absolute
maximality, or maximality among all discrete reflection overgroups, is more
than sufficient.

A **weak-log rational return** of \(P\) is the following package over
\(\mathbf C\).

1. A smooth projective rational surface \(S\), an oriented marking of its
   full N\'eron--Severi lattice \(L\), and distinct geometrically integral
   facet divisors \(C_1,\ldots,C_f\).
2. Every \(C_i\) is a smooth rational curve of square \(-1,-2\), or \(-4\).
   The square-minus-four curves are precisely a disjoint boundary
   \(B=B_1+\cdots+B_b\).
3. With

   \[
       M=\langle L,B_1/2,\ldots,B_b/2\rangle,
       \qquad q_M=2q_L,
   \]

   every side reflection preserves \(M\), and the oriented halfspaces of
   the \(C_i\) reproduce the labelled chamber \(P\).
4. The class \(D=-(K_S+B/2)\) is nef.
5. The return is literal:

   \[
                         P=\mathbf P(\operatorname {Nef}S).
   \]

The adjective *full* is essential: \(L=\operatorname {NS}(S)\), not a
proper trace slice or a selected divisor sublattice.  Geometric integrality
in item 1, together with the numerical adjunction rows, makes every facet
curve smooth and rational: an integral curve of arithmetic genus zero on a
smooth surface has no normalization genus or delta defect.

## 2. The raw host list is finite

Normalize the Coxeter normal lines by \(u_i^2=-1\).  In a weak-log return the
physical square of the curve on the \(i\)-th line is one of

\[
                         C_i^2\in\{-1,-2,-4\}.             \tag{2.1}
\]

There are therefore only \(3^f\) physical Gram candidates:

\[
             C_i.C_j=\sqrt{(-C_i^2)(-C_j^2)}\,(u_i.u_j).   \tag{2.2}
\]

Reject a candidate unless all entries are integral and the displayed classes
span the Lorentz space.  Finite volume guarantees full real span for every
actual return.  Let \(R=\sum_i\mathbf ZC_i\).  Any full Picard lattice of a
return satisfies

\[
                         R\subseteq L\subseteq R^*.         \tag{2.3}
\]

The discriminant group \(R^*/R\) is finite, so its subgroups enumerate a
finite list of integral overlattices.  Retain only unimodular lattices of
signature \((1,n)\).  This is exhaustive because the N\'eron--Severi lattice
of a smooth complex rational surface is unimodular.

This is a rational-form test, not merely a real-inertia test.  Every proposed
Picard host must be \(\mathbf Q\)-isometric, up to the one common positive
similarity, to the full physical wall space.  Equivalently, constructing
\(L\) in (2.3) takes place inside that same rational quadratic space.  In
particular rational isotropy, determinant square class, and local Hasse data
cannot be changed by declaring a new integral host.

The facet adjunction equations

\[
                 K.C_i=-2-C_i^2                            \tag{2.4}

\]

determine \(K\in L_{\mathbf Q}\) uniquely, since the facets span.  Retain
only \(K\in L\) which is characteristic and satisfies

\[
                         K^2=10-\operatorname {rk}L.        \tag{2.5}

\]

The \((-4)\)-rows determine \(B\); test their orthogonality, form \(M\),
and test every reflection on \(M\) by exact lattice membership.  Thus no
canonical class, physical scale, boundary label, or Picard host is smuggled
into the geometric stage.

If the minimal Vinberg field of \(P\) is not \(\mathbf Q\), this list is
automatically empty: a literal rational full-nef return has rational curve
normals and rational cyclic Cartan products.  This is a rejection theorem,
not a construction of a nonrational carrier.

## 3. The birational list is finite and exhaustive

Fix one surviving numerical host.  If a weak-log return exists, every
\((-1)\)-curve spans a negative extremal Mori ray and hence its orthogonal
hyperplane is a facet of the complete nef cone.  It is therefore one of the
displayed \(C_i\).

Run every sequential contraction supported on the finite facet list.  If
\(E_1,\ldots,E_s\) have already been contracted and \(A=(E_i.E_j)\), the
pullback of the current image of a candidate \(C\) is

\[
        e_C=C-EA^{-1}E^{\mathsf t}GC.                      \tag{3.1}

\]

It is an allowed next pivot exactly when

\[
                         e_C^2=K.e_C=-1.                   \tag{3.2}

\]

Boundary rows are allowed to become pivots later in a program; only their
original curves, not their images on intermediate models, are required to
have square \(-4\).  This small enlargement avoids assuming that a
minimalization can be chosen boundary-preserving.

The recursion has depth at most \(\rho-1\) and finite branching.  It stops at
\(\mathbf P^2\) or at a Hirzebruch surface \(\mathbb F_e\).  Conversely,
any smooth rational return admits such a sequence.  A curve contracted later
has negative strict transform on the original surface, hence was already a
displayed facet.  Thus the finite recursion loses no return.

There is one necessary endpoint finiteness check: the marked pair
\((\operatorname {NS}(\mathbb F_e),K_{\mathbb F_e})\) alone remembers only
the parity of \(e\), not its magnitude.  If \(e\ge2\), however, the negative
section \(\Sigma\) has an irreducible strict transform \(\widetilde\Sigma\)
on \(S\), with

\[
                 \widetilde\Sigma^2=-e-t<0,               \tag{3.3}
\]

where \(t\ge0\) is the number of blowup centres lying on its successive
strict transforms.  Literal full-nef return makes
\(\widetilde\Sigma^\perp\) a displayed facet, so the physical alphabet
forces \(e+t\in\{1,2,4\}\).  Since \(\mathbb F_1\) is not minimal, the only
ruled endpoints are

\[
                         e\in\{0,2,3,4\}.                 \tag{3.4}
\]

Thus the endpoint surface list itself is finite; this conclusion would not
follow from the numerical pair \((\operatorname {NS},K)\) alone.

For each endpoint program \(\tau\), reverse the contractions.  The resulting
proper/free/satellite proximity pattern defines a finite-type ordered cluster
scheme \(T_\tau\), covered by finitely many explicit affine blowup charts.
This is the abstract arbitrary-input version of PB-Atlas Stages A--B2a.

## 4. Simultaneous physical curves form a finite constructible locus

Let \(\mathcal S_\tau\to T_\tau\) be the universal iterated blowup.  Every
facet vector gives a relative line bundle \(\mathcal L_i\).  Take the
relative Hilbert scheme of divisors in that class and its geometrically
integral locus.  Their fibre product over \(T_\tau\) is denoted

\[
                          \operatorname {Prime}_\tau(P).    \tag{4.1}
\]

It is constructible and of finite type.  Equivalently, on the finite PB
charts it is cut out by the actual jet matrices, exact-multiplicity opens,
and factorization-avoidance conditions.  This is precisely where the
bi-anticanonical section algebra enters: a formal class or expected jet rank
is not substituted for a divisor-prime section.

All curves occur on the **same** cluster point.  Their intersections are
then their prescribed lattice intersections.  Distinct prime divisors of
intersection zero are disjoint, so the \((-4)\) rows give a disjoint smooth
boundary and \((\mathcal S_\tau,B/2)\) is klt.

## 5. The nef gate is finite in all three numerical cases

Numerically require \(D\) to lie in the proposed closed chamber and require
\(D^2\ge0\).  The remaining geometric nefness test splits as follows.

### 5.1 The strong-Coble case: \(D=0\)

There is no nefness debt.  Simultaneous prime facets and
\(B=-2K\) already give the strong-Coble package.  Reflection maximality then
returns the complete nef cone.

### 5.2 The nonzero cusp case: \(D^2=0\)

This case does not need a semiample multiple.  The following lemma is the
key new finite certificate.

### Lemma 5.1 (positive cusp divisor)

Let \(P\) have finite volume and let \(0\ne D\) be a rational null vector in
its homogeneous wall cone.  Then a subset \(J\) of the facets incident to
\([D]\) and positive rational numbers \(a_j\) satisfy

\[
                         \lambda D=\sum_{j\in J}a_jC_j,
                 \qquad \lambda>0.                         \tag{5.1}
\]

The relation is found by finite exact linear algebra.

#### Proof

The finite-volume wall-cone lemma makes every null ray in the wall cone an
extreme ray, hence an ideal vertex.  Its active facet normals span
\(D^\perp\).  Choose a minimal subset spanning that hyperplane.  Its Gram
matrix is negative semidefinite; after changing sign it is a symmetric
singular M-matrix because Coxeter off-diagonal intersections are
nonnegative.  Exactly one irreducible block is singular.  The kernel of an
irreducible singular M-matrix is generated by a strictly positive vector.
The corresponding linear combination of its facet normals spans the radical
of \(D^\perp\), namely \(\mathbf RD\).  Orient it toward \(D\).  All data
are rational, giving (5.1); clearing denominators makes the coefficients
integral. \(\square\)

On \(\operatorname {Prime}_\tau(P)\), clear denominators in (5.1) and put

\[
                         F=mD=\sum_{j\in J}b_jC_j.          \tag{5.2}

\]

This is an effective divisor.  Each component in its support has
\(F.C_j=mD.C_j=0\), while every other irreducible curve meets each prime
component nonnegatively.  Hence \(F\), and therefore \(D\), is nef.

This explains the exact role of the isotropic period obstruction.  The
normal character can decide whether \(D\) is semiample, but semiampleness is
not required for exact nef return.  The effective affine cusp divisor proves
the required nefness even when its normal character is non-torsion.

### 5.3 The big case: \(D^2>0\)

Put \(\Delta=B/2\) and \(L_D=2D=-2K-B\), which is Cartier.  On any actual
return, \(D\) is nef and big and \((S,\Delta)\) is klt.  Fujino's effective
base-point-free theorem applies with

\[
       \dim S=2,\qquad a=0,\qquad
       aL_D-(K+\Delta)=D.                                 \tag{5.3}
\]

Its explicit multiplier is

\[
 m(2,0)=2^{2+1}(2+1)!(\lceil0\rceil+2)=96.                \tag{5.4}

Consequently every actual return satisfies

\[
                         |96L_D|=|192D|
                    \quad\text{base-point-free}.           \tag{5.5}

Conversely, (5.5) immediately proves that \(D\) is nef.  Thus (5.5) is an
iff-safe *fixed* finite gate: necessity comes from the theorem and
sufficiency from global generation.  The source is Osamu Fujino,
[Effective base point free theorem for log canonical pairs, Theorem 1.1 and
Remark 1.2](https://www.math.nagoya-u.ac.jp/~fujino/effective4.pdf).

On a finite affine cover of \(\mathcal S_\tau\), the complete system in
(5.5) is computed by finite module/jet matrices.  After the standard
cohomology-and-base-change flattening stratification, the support of the
cokernel of the evaluation map

\[
 \pi^*\pi_*\mathcal O(192D)\longrightarrow\mathcal O(192D) \tag{5.6}
\]

is closed over each stratum.  Proper projection to the stratum and taking
the complement produce a finite constructible base-point-free locus.

## 6. Exact-return theorem and the compiler

For each \(\tau\), intersect \(\operatorname {Prime}_\tau(P)\) with the
appropriate nef locus of Section 5 and call the result
\(\mathcal R_\tau(P)\).  There are finitely many \(\tau\).

### Theorem 6.1 (`RLogComp` iff/`EMPTY` theorem)

Under the hypotheses of Section 1,

\[
 P\text{ has a complex weak-log rational return}
 \quad\Longleftrightarrow\quad
 \coprod_\tau\mathcal R_\tau(P)(\mathbf C)\ne\varnothing. \tag{6.1}
\]

For every point on the right,

\[
                         P=\mathbf P(\operatorname {Nef}S). \tag{6.2}

\]

#### Proof

Necessity of the finite host and MMP lists was proved in Sections 2--3.
Every actual prime facet tuple gives a point of (4.1).  In the big case it
passes (5.5) by Fujino's fixed bound; in the null case Lemma 5.1 supplies
(5.2); and the zero case is tautological.  Thus every return appears on the
right.

Conversely, a point of \(\mathcal R_\tau(P)\) constructs the marked smooth
rational surface, all simultaneous prime facet curves, the disjoint
\((-4)\) boundary, the adapted integral host, and actual nefness of \(D\).
The maximal weak-log exact-return theorem then classifies every hidden
negative curve into the adapted integral reflection alphabet.  Arithmetic
reflection maximality saturates all those mirrors, and the ample orientation
gives (6.2), including ideal rays. \(\square\)

All schemes and constructible sets above are defined by finitely many
polynomials over a computable number field.  Geometric integrality can be
implemented by the finite factorization strata (or relative Hilbert/Chow
maps), and base-point-freeness by (5.6).  Gr\"obner elimination plus the
Rabinowitsch trick decides whether each constructible stratum has an
algebraically closed point.  Therefore the procedure terminates with either
a finite constructible carrier-stack presentation or `EMPTY`.  This is a
mathematical compiler theorem; the current PB-Atlas implementation realizes
the numerical and plane-chart front end but not yet the enormous generic
degree-\(192D\) elimination.

## 7. `P10`, `P11`, and the limits of maximality

The theorem separates three meanings that the bare mirror arrangement does
not separate.

1. **Decorated adapted `P10`.**  The requested physical class is
   \(B=2(-K)\), so \(D=0\).  This packet is finite-volume and reflection
   maximal, yet its simultaneous prime-facet locus is empty over
   \(\mathbf C\).  On the flex carrier the same class is the nonreduced
   divisor \(2F\); on a strong-Coble blowup a prime \(-2K\) can exist, but
   not simultaneously with the requested complete `P10` facet stock.  This
   proves that maximality cannot force the `PRIME` gate.
2. **Uncoloured `P10` and `P11`.**  Their flex carriers interpret the special
   mirror as the primitive irreducible genus-one curve \(F=-K\), of square
   \(-1\) and \(-2\), respectively.  These are exact positive full-nef
   returns and give the sharp anticanonical-facet bound \(n\le11\), but they
   are not objects of the smooth-rational-facet weak-log target in Section 1.
3. **Ordinary smooth-rational colouring.**  With only
   \((-1,-1)\) and \((-2,0)\) adjunction rows, the existing P-first theorem
   is stronger: finite volume and relative maximality force \(n\le9\) and
   the standard flex-cluster carrier directly.

The broader direct full-NS theorem \(n\le11\) remains a theorem about every
already constructed carrier with the declared negative-wall alphabet.  In
hostile incidence lanes where no nef weak-log class exists, maximality alone
does not make the physical curves effective and does not give a finite
hidden-curve test.  `RLogComp` does not silently claim to solve that larger
problem.

### 7.1 Raw-maximal chambers with every rational/log return empty

There is a particularly sharp rational-field example.  The compact
arithmetic \((2,4,6)\) triangle is absolutely reflection-maximal and has
one rational Gram

\[
 G=
 \begin{pmatrix}
 -12&6&0\\
 6&-6&3\\
 0&3&-2
 \end{pmatrix}.                                           \tag{7.1}
\]

Its rational form is anisotropic.  Indeed, an isotropic vector would give a
nonzero rational, hence primitive integral, solution of

\[
 6x^2-6xy+3y^2-3yz+z^2=0.                                \tag{7.2}
\]

Reduction modulo \(3\) first forces \(3\mid z\).  Write \(z=3z_1\),
divide (7.2) by \(3\), and reduce again.  The equation becomes

\[
                         2x^2+xy+y^2=0\pmod3,              \tag{7.3}
\]

whose only solution is \(x=y=0\).  Thus \(3\) divides all three
coordinates, contradicting primitivity.

By contrast, the Picard lattice of a complex rational surface of rank three
is \(I_{1,2}\), which is rationally isotropic.  Hence the `HOST` list for
this raw triangle is empty.  There is also a completely geometric version:
every rank-three rational surface has a nef isotropic pullback of a ruling
after minimalization, so its projectivized full nef chamber has an ideal
point; it cannot equal a compact triangle.  Thus even

```text
k(P)=Q + compact arithmetic + absolute reflection maximality
```

does not force a rational/log carrier.  See the companion
[`(2,4,6) raw-maximal rational no-go`](compact_246_raw_maximal_property4_counterexample.md).

There is a second undecorated, completely P-native negative example.  The
compact arithmetic \((2,5,6)\) triangle is absolutely reflection-maximal and
has minimal Vinberg field

\[
                              \mathbf Q(\sqrt5).            \tag{7.4}
\]

Every literal full-nef rational-surface chamber has minimal Vinberg field
\(\mathbf Q\), because its facet normals are rational curve classes and all
cyclic Cartan products are rational.  Therefore the \((2,5,6)\) triangle has
**no** rational or rational-log full-nef carrier, regardless of how its
mirrors are physically coloured.  This answers the branchwise no-go question
without relying on a proposed canonical marking.

The companion
[`(2,5,6) K3 PHYS-empty theorem`](k3_compact_maximal_256_phys_empty.md)
proves more: even the primitive full-trace/full-Picard K3 compiler is empty,
by an inert-prime square-class obstruction.  Neither statement excludes a
proper projected eigenslice in a larger Picard lattice or another carrier
species.  What they exclude is the inference

```text
raw reflection maximality => nonempty physical carrier compiler.
```

## 8. Consequences and exact scope

Whenever `RLogComp` is nonempty, the existing AG transfer gives

\[
                         n\le25,\qquad k(P)=\mathbf Q,      \tag{8.1}
\]

and only finitely many relatively reflection-maximal weak-log chambers admit
such carriers.  If \(B=0\), the sharper weak-del-Pezzo conclusion is
\(n\le9\).

What is now closed is the logical arrow

```text
fixed raw maximal rational P
    -> finite HOST/MMP/PRIME/NEF constructible calculation
    -> exact carrier stack or EMPTY
    -> literal same-P return and AG bounds.
```

What is not closed is universal nonemptiness of that stack, nor the global
disjunction with K3/eigen-nef, Enriques, genus-one-facet, and hostile direct
rational species.  The theorem therefore advances Property 4 from an
unbounded search to a terminating decision problem on the full weak-log
species, without converting `iff/EMPTY` into a false essential-surjectivity
claim.

## Replay

```bash
python3 scripts/certify_rational_weak_log_finite_iff_compiler.py \
  --check-archive

python3 -m pytest -q \
  tests/test_rational_weak_log_finite_iff_compiler.py
```
