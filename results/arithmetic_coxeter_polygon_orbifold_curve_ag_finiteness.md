# Arithmetic Coxeter polygons from algebraic orbifold curves

## Exact inverse carrier and a genuine curve-theoretic finiteness theorem

**Status.**  This is an unconditional two-dimensional theorem.  It is not a
proof of finiteness in all hyperbolic dimensions.  Unlike the weak blow-up or
quadric-arrangement encodings, however, the algebraic object is used
essentially: the orientation double of a Coxeter polygon is a genus-zero
arithmetic orbifold curve, and the finiteness theorem for genus-zero
arithmetic Fuchsian/Shimura curves bounds the possible ambient arithmetic
groups.  No ABSW, Nikulin reflective-lattice finiteness theorem, or Vinberg
dimension theorem is used.

The argument is deliberately stated for marked finite-area polygons; compact
polygons are the proper subcase and ideal vertices become cusps.

## 1. Categories

Let \(\mathsf {Cox}^{\max}_{2,\mathrm {arith}}\) have as objects marked
finite-area arithmetic Coxeter polygons \(P\subset\mathbb H^2\) whose full
side-reflection group

\[
                         W_P<\operatorname {Isom}(\mathbb H^2)
\]

is maximal among discrete reflection groups.  Morphisms are marked hyperbolic
isometries.

Let \(\mathsf {OrbCurv}^{\mathbf R}_{0,\mathrm {arith}}\) have as objects

\[
        (\overline {\mathcal X}_{\mathbf R},D_{\rm cusp},c,\mathcal H),
\]

where:

1. \(\overline {\mathcal X}_{\mathbf R}\) is a proper real algebraic
   Deligne--Mumford curve and \(D_{\rm cusp}\) is its marked reduced cusp
   divisor; the complement has complex analytification equal to an arithmetic
   hyperbolic orbifold;
2. the coarse compactified curve has geometric genus zero;
3. \(c\) denotes the real descent involution on the complex analytification,
   induced by the orientation-reversing coset; and
4. \(\mathcal H\) is the choice of the mirrored half/component together with
   its cyclic facet marking.

The last item is finite real framing.  It prevents the complex curve alone
from forgetting which real half is the original chamber.
Morphisms are holomorphic or antiholomorphic (equivalently, real-semilinear
or dianalytic) isomorphisms of the orbifold log curves preserving
\(D_{\rm cusp}\), the descent involution, the chosen half, and its cyclic
marking.  Thus orientation-reversing marked isometries are also carried
functorially.

## 2. The functor and exact return

Write

\[
                    \Gamma_P=W_P\cap\operatorname {Isom}^+(\mathbb H^2).
\]

Define

\[
 F_2(P)=
 \left(\overline {[\mathbb H^2/\Gamma_P]}_{\mathbf R},
 D_{{\rm cusp},P},c_P,\mathcal H_P\right).
                                                        \tag{2.1}
\]

The analytic orbifold, its compactification, and its real involution are
algebraic by orbifold Riemann existence, GAGA, and real descent.  When the
group is congruence this is the usual Shimura/modular orbifold curve.  The
initial subgroup \(\Gamma_P\) need not be congruence, so no canonical Shimura
model is asserted for it.  The marked cusp divisor is retained because the
proper coarse curve alone does not distinguish an ideal vertex.

> **Theorem 2.1 (exact orbifold-curve carrier).**  The assignment \(F_2\) is
> functorial on marked isometries.  There is a return functor \(R_2\) on its
> image satisfying
>
> \[
>                         R_2F_2(P)\cong P.              \tag{2.2}
> \]
>
> It recovers the cyclic side list, finite labels, cusps, complete chamber,
> and the hyperbolic metric, not only the abstract Coxeter presentation.

**Proof.**  Topologically \(P\) is a disk.  Doubling it across its mirrored
boundary gives a sphere.  A vertex of angle \(\pi/m\) becomes a cone point of
order \(m\), and an ideal vertex becomes a puncture.  Hence the compactified
coarse curve of \([\mathbb H^2/\Gamma_P]\) has genus zero.

The quotient by \(c_P\), with the selected half \(\mathcal H_P\), is exactly
the mirrored orbifold \(P\).  Orbifold uniformization gives its unique
complete curvature-\((-1)\) metric.  The cone order at a doubled vertex is
\(m\), so the returned interior angle is \(\pi/m\); punctures return ideal
vertices.  The cyclic real marking returns the complete side order and
therefore the same polygon.  Conjugating \(W_P\) by a marked isometry
conjugates \(\Gamma_P\), \(c_P\), and the chosen half, proving functoriality.
\(\square\)

The real framing is essential.  The bare complex coarse curve is often just
\(\mathbf P^1\), and therefore cannot remember the polygon.

## 3. The algebraic-curve finiteness input

The theorem used here is D. Long, C. Maclachlan, and A. Reid,
[*Arithmetic Fuchsian Groups of Genus Zero*](https://web.ma.utexas.edu/users/areid/genus0_final.pdf),
Theorem 1.1:

> there are finitely many conjugacy classes of **congruence** arithmetic
> Fuchsian groups whose orbifold curves have genus zero; consequently the
> same holds for maximal arithmetic Fuchsian groups of genus zero.

Their proof is a theorem about arithmetic orbifold/Shimura curves, not a
reflective-lattice classification.  Maximal arithmetic Fuchsian groups are
congruence.  The uniform automorphic spectral gap and the Zograf inequality
for a genus-\(g\) orbifold curve bound the area in terms of \(g\); for
\(g=0\), the displayed bound in that paper is

\[
                    \operatorname {Area}\le {128\pi\over3}.  \tag{3.1}
\]

Together with finiteness of arithmetic Fuchsian groups of bounded area this
gives the stated finite list.  Equivalently, this can be viewed through the
linear gonality philosophy for congruence modular/Shimura curves: a
genus-zero curve has gonality one, while congruence complexity forces
gonality to grow.  The spectral/congruence input is indispensable; the stack
\(\mathcal M_{0,r}\) by itself is not finite.

Thus this result is genuinely mediated by an algebraic-curve invariant
(genus, or gonality), but it is not advertised as a proof from BAB/KSBA
alone.

## 4. Pushback to maximal reflection polygons

We first record the precise congruence bridge needed below.

> **Lemma 4.0 (maximal full lattice congruence bridge).**  If
> \(\Lambda<\operatorname {Isom}(\mathbb H^2)\) is a maximal arithmetic
> lattice, then \(\Lambda^+\) contains a principal congruence subgroup of a
> reference arithmetic Fuchsian group.  In particular \(\Lambda^+\) is a
> congruence arithmetic Fuchsian group.

**Proof.**  Choose a faithful representation of the defining \(k\)-group and
a reference \(\mathcal O_k\)-lattice \(L_0\).  Since \(\Lambda\) is
arithmetic,

\[
 K=\operatorname {core}_{\Lambda}
       (\Lambda\cap\operatorname {Stab}(L_0))
\]

has finite index in \(\Lambda\), is normal there, and preserves \(L_0\).
For coset representatives \(\lambda_i\) of \(\Lambda/K\), put
\(L=\sum_i\lambda_iL_0\).  If
\(\lambda\lambda_i=\lambda_jk\) with \(k\in K\), then
\(\lambda(\lambda_iL_0)=\lambda_jL_0\); hence \(L\) is stabilized by
\(\Lambda\).  Its full arithmetic stabilizer is a
discrete lattice containing \(\Lambda\); maximality gives

\[
                  \Lambda=\operatorname {Stab}(L).       \tag{4.0}
\]

Choose a nonzero ideal \(\mathfrak a\) with
\(\mathfrak aL_0\subset L\subset\mathfrak a^{-1}L_0\).  A sufficiently deep
principal congruence subgroup, for example at an ideal contained in
\(\mathfrak a^2\), acts trivially enough on both inclusions to stabilize
\(L\).  Take it from the outset inside the norm-one,
orientation-preserving reference Fuchsian group.  For
\(g\equiv1\pmod {\mathfrak a^2}\), one has
\((g-1)L\subset\mathfrak aL_0\subset L\), and the same calculation for
\(g^{-1}\) gives equality.  By (4.0) this principal subgroup lies in
\(\Lambda^+\).  Under
\(\operatorname {Spin}(q)\simeq\operatorname {SL}_1(B)\), deep principal
levels in this faithful orthogonal model and in the quaternion-order
definition are cofinal; this is congruence in the sense used by
Long--Maclachlan--Reid.  This is the standard
lattice-stabilizer proof of the maximal-arithmetic congruence lemma. \(\square\)

> **Theorem 4.1 (AG-mediated finiteness in \(\mathbb H^2\)).**  There are
> finitely many marked-isometry classes of maximal finite-area arithmetic
> Coxeter polygons.  The proof factors as
>
> \[
> P\longmapsto F_2(P)
> \longmapsto\text{genus-zero arithmetic orbifold curve}
> \longmapsto\text{finite ambient list}
> \longmapsto P.                                        \tag{4.1}
> \]

**Proof.**  Let \(W=W_P\) and choose a maximal arithmetic lattice
\(\Lambda<\operatorname {Isom}(\mathbb H^2)\) containing \(W\).  The subgroup
of \(\Lambda\) generated by all reflections of \(\Lambda\) is conjugacy
invariant, discrete, and contains \(W\).  Reflection maximality therefore
forces it to equal \(W\).  In particular

\[
                         W\triangleleft\Lambda.          \tag{4.2}
\]

Moreover, \(\Lambda/W\) is the finite symmetry group of the finite-area
fundamental polygon \(P\).  This is the elementary maximal-lattice lemma used
in reflection-group arguments; it does not invoke a finiteness theorem.

Put \(\Lambda^+=\Lambda\cap\operatorname {Isom}^+(\mathbb H^2)\) and
\(\Gamma=W^+\).  Equation (4.2) gives a finite orbifold quotient whose
compactified coarse curves form

\[
 \overline X_\Gamma^{\rm coarse}\simeq\mathbf P^1
       \longrightarrow \overline X_{\Lambda^+}^{\rm coarse}.
\]

A finite quotient of \(\mathbf P^1\) again has genus zero.  Lemma 4.0 says
that \(\Lambda^+\) is congruence.  The first clause of Theorem 1.1 of
Long--Maclachlan--Reid therefore leaves only finitely many conjugacy classes
of \(\Lambda^+\).

The full lattice \(\Lambda\) lies in the normalizer of \(\Lambda^+\), whose
quotient by \(\Lambda^+\) is finite.  Hence only finitely many full
\(\Lambda\) occur.  For fixed \(\Lambda\), there are finitely many conjugacy
classes of reflections: they are the finitely many codimension-one orbifold
singular strata in a finite-sided fundamental orbifold.  Since \(W\) is the
subgroup generated by all reflections of \(\Lambda\), it is actually
determined by \(\Lambda\).  Exact return in Theorem 2.1 now gives only
finitely many polygons. \(\square\)

## 5. What this does and does not solve

This theorem meets the requested logical pattern in dimension two:

```text
raw maximal arithmetic polygon
  -> exact real-framed algebraic orbifold curve
  -> independent genus-zero arithmetic-curve finiteness
  -> finite list of the original polygons.
```

It is not the desired all-dimensional theorem.  The reason the construction
works in dimension two is structural: \(\operatorname {SO}(1,2)^\circ\) is
of Hermitian type, so its arithmetic quotients are algebraic curves.  For
\(n>2\), \(\operatorname {SO}(1,n)/\operatorname {SO}(n)\) is not Hermitian;
the hyperbolic quotient is not itself an algebraic Shimura variety.  The
determinant-line type-IV construction embeds it only as a totally real,
infinite-index slice, so higher-dimensional reflection maximality does not
currently imply an AG boundedness theorem on the ambient Shimura stack.

The theorem therefore supplies a successful base case and isolates the
remaining global debt precisely: construct a bounded algebraic avatar for
the non-Hermitian higher-dimensional quotients, or prove that the existing
quadric-arrangement or mixed K3/Fourier--Mukai avatar lands in a uniformly
bounded, cohomologically discrete locus.
