# Higher-dimensional Coxeter orbifolds are not algebraic quotient varieties

**Status (2026-08-12): exact real-analytic orbifold quotient, but a terminal
obstruction to extending the genus-zero algebraic-curve proof beyond
\(\mathbf H^2\).**  The obstruction occurs on an actual compact,
reflection-maximal arithmetic chamber in \(\mathbf H^6\).  It is independent
of ABSW, Nikulin, and reflective-lattice finiteness.

The conclusion is deliberately narrow.  It does not obstruct the global
quadric-arrangement functor or an external K3/Shimura carrier.  It proves
that the quotient which makes the two-dimensional argument work cannot
itself be the required higher-dimensional complex algebraic DM object.

## 1. What the quotient remembers exactly

Let \(P\subset\mathbf H^n\) be a finite-volume Coxeter chamber, let \(W\)
be its side-reflection group, and put

\[
                         \Gamma=W^+.
\]

The differentiable orbifold

\[
                   \mathcal O_P=[\mathbf H^n/\Gamma]                 \tag{1.1}
\]

is an exact carrier.  Its orientation-reversing involution has a selected
half isometric to \(P\); the cyclic ridge inertia recovers every finite
Coxeter angle, and the hyperbolic developing map recovers the complete
metric chamber.

On coarse spaces, (1.1) is the double of \(P\) along its mirror facets.
After adding its finite set of cusp points,

\[
       |\overline{\mathcal O}_P|\simeq \operatorname {Dbl}(\overline P)
       \simeq S^n.                                                   \tag{1.2}
\]

For compact \(P\), no end compactification is needed.  Formula (1.2) is the
higher-dimensional analogue of the genus-zero coarse curve in dimension
two.  The decisive difference is that (1.1) is a **real** orbifold and not a
complex algebraic orbifold when \(n>2\).

## 2. A compact reflection-maximal arithmetic \(\mathbf H^6\) control

Start with the known compact arithmetic Bugaenko-double chamber
\(P_{6,10}\).  Its reflection group \(W_0\) is a cocompact arithmetic
lattice in \(\operatorname {Isom}(\mathbf H^6)\).

There is a maximal discrete reflection overgroup \(W_{\max}\supset W_0\),
without using a reflection-group finiteness theorem.  Indeed, if a discrete
group \(R\) contains the lattice \(W_0\), then the covering of locally
symmetric orbifolds shows that \([R:W_0]\) is finite; otherwise the
finite-volume quotient by \(W_0\) would have infinitely many positive-volume
sheets over the quotient by \(R\).  Thus \(R\) is itself a lattice and

\[
 [R:W_0]={\operatorname {vol}(\mathbf H^6/W_0)
                 \over\operatorname {vol}(\mathbf H^6/R)}.           \tag{2.1}
\]

The Kazhdan--Margulis lower covolume bound in fixed dimension therefore
bounds the set of positive integers

\[
 \mathcal I=\{[R:W_0]:R\supset W_0\text{ is a discrete reflection group}\}.
                                                                         \tag{2.2}
\]

This nonempty bounded subset of \(\mathbf Z_{>0}\) has a maximum.  Choose
\(W_{\max}\) realizing it.  A strict discrete reflection overgroup
\(R\supsetneq W_{\max}\) would satisfy

\[
 [R:W_0]=[R:W_{\max}][W_{\max}:W_0]>[W_{\max}:W_0],
\]

contradicting maximality in (2.2).  Since \(W_{\max}\) is a finite
extension of \(W_0\), the two groups are commensurable; hence cocompactness
and arithmeticity pass to \(W_{\max}\).  Its chamber \(P_{\max}\) is
therefore a compact absolutely reflection-maximal arithmetic source object
in \(\mathbf H^6\).

Put \(\Gamma_{\max}=W_{\max}^+\).  A torsion-free normal finite-index
subgroup \(\Gamma_0\) exists by Selberg's lemma, and

\[
                    M=\mathbf H^6/\Gamma_0                            \tag{2.3}
\]

is a closed real-hyperbolic six-manifold.

## 3. The direct algebraization obstruction

> **Theorem 3.1 (no complex-algebraic quotient avatar).**  There is no
> smooth proper complex algebraic DM stack \(\mathscr X\) whose analytic
> orbifold is
>
> \[
>             \mathscr X^{\rm an}\simeq
>             [\mathbf H^6/\Gamma_{\max}].                            \tag{3.1}
> \]
>
> In particular, the maximal arithmetic quotient in Section 2 has no
> smooth proper projective-DM (and hence no smooth proper KSBA-orbifold)
> avatar, and it is not a Baily--Borel arithmetic quotient.  This does not
> exclude a singular KSBA pair or another external algebraic carrier that
> merely records the Coxeter data.

**Proof.**  If (3.1) existed, the normal subgroup \(\Gamma_0\) would define
a finite representable etale analytic cover of \(\mathscr X^{\rm an}\).
The Riemann existence theorem for finite etale covers of proper DM stacks
algebraizes it to \(X_0\to\mathscr X\).  Because \(\Gamma_0\) is
torsion-free, \(X_0\) has no orbifold stabilizers.  It is a smooth proper
complex algebraic space with

\[
                X_0^{\rm an}\simeq M,
                \qquad \pi_1(X_0)=\Gamma_0.                           \tag{3.2}
\]

Chow's lemma for algebraic spaces, followed by resolution, gives a smooth
projective variety birational to \(X_0\).  Topological fundamental group is
a birational invariant for smooth complete complex varieties (equivalently,
blowing up a smooth center does not change it).  Thus \(\Gamma_0\) would be
a Kahler group.

Carlson--Toledo prove that the fundamental group of a closed
real-hyperbolic \(n\)-manifold, for \(n>3\), is not a Kahler group.  Applied
to (2.3), this contradicts (3.2). \(\square\)

For a polarized avatar there is an even more elementary coarse obstruction.
The coarse space in (1.2) is \(S^6\), whose \(H^2\) vanishes, whereas every
positive-dimensional projective variety has a nonzero ample class.  The
stack argument above is stronger because it survives passage to the
torsion-free cover and sees the orbifold fundamental group rather than only
the coarse sphere.

The same harmonic-map theorem also blocks a proposed projective branched
cover.  If \(Y\) is compact Kahler and \(N\) is a closed real-hyperbolic
\(n\)-manifold with \(n>2\), Carlson--Toledo's rank theorem makes every map
\(Y\to N\) trivial on homology in degrees greater than two.  In equal real
dimension it has degree zero.  Hence no finite projective branched cover of
nonzero degree can turn (2.3) into an algebraic quotient.

This is not merely the parity observation that an odd-dimensional real
orbifold cannot itself be complex analytic.  It excludes the first even
higher-dimensional control where complex dimension three would be
topologically possible.

## 4. Why the apparent escape routes do not give AG finiteness

### 4.1 Real algebraic models are too flexible

Theorem 3.1 concerns complex algebraization of the quotient orbifold.  A
closed hyperbolic manifold can certainly be represented as the real locus
of some real algebraic variety: the Nash--Tognoli theorem does this for
every closed smooth manifold.  Thus the predicate

```text
the quotient has some real algebraic model
```

places no restriction on dimension, volume, field, or facet data.  Such a
model is noncanonical, and its complex points do not equal the hyperbolic
quotient.  Adding the developing map, Gram local system, and selected real
half can restore exact return, but then all Coxeter information resides in
unbounded decorations rather than in the algebraic model.

This explains the difference between (1.1) and the existing exact
quadric-arrangement or orthogonal-Shimura carriers: those are honest complex
AG objects with a marked **totally real fixed component**.  They do not
claim that the fixed component is an algebraic subvariety.

### 4.2 Coarse sphere topology and low cohomology are constant

The coarse compactification (1.2) has the homology of a sphere for every
Coxeter chamber, regardless of its number of facets, angles, field, or
volume.  The singular cell stratification contains the missing data, but
its complexity is exactly what needs to be bounded.

There is a sharp higher-dimensional control.  The integral \(\mathbf H^{10}\)
redoubling family \(P_N\) has

\[
 f(P_N)=N+12,\qquad
 [W(P_1):W(P_N)]=N,\qquad
 \operatorname {vol}(P_N)=N\operatorname {vol}(P_1).                 \tag{4.1}
\]

The seed is a finite-covolume subgroup of an integral orthogonal group, so
it and every finite-index descendant are arithmetic.  Yet every oriented
coarse cusp compactification is still \(S^{10}\).  These descendants are
nonmaximal and are **not** a counterexample to boundedness of the maximal
class.  Their exact role is narrower: coarse topology, arithmeticity, and
low cohomology alone do not produce a bound.  A new theorem would have to
turn **reflection maximality** into a bound on the singular stratification
or covolume.

Stack cohomology does not repair this.  Since \(\mathbf H^n\) is
contractible, the orbifold homotopy type of (1.1) is the classifying-space
type of \(\Gamma\), not that of its coarse sphere.  Calling those varying
group-cohomology classes algebraic does not put them in a bounded family.

### 4.3 Congruence gives no Baily--Borel compactification

The H2 proof works because

\[
        \operatorname {SO}(1,2)^\circ/\operatorname {SO}(2)
\]

is a Hermitian symmetric domain.  For \(n>2\),

\[
        \operatorname {SO}(1,n)^\circ/\operatorname {SO}(n)          \tag{4.2}
\]

is not Hermitian.  Even if a maximal arithmetic lattice is congruence,
Baily--Borel does not apply to (4.2).  Borel--Serre instead supplies a real
manifold-with-corners compactification.  It has no canonical ample bundle,
Hilbert polynomial, KSBA pair, or algebraic degree to which projective AG
boundedness can be applied.  For the compact \(\mathbf H^6\) control there
is not even a cusp boundary on which a different compactification could
change the conclusion.

Embedding (4.2) as a totally real fixed locus in an orthogonal or unitary
Shimura variety does produce an exact **external** carrier.  The Coxeter
group is then infinite index in the ambient arithmetic group, and real
reflection maximality does not bound the ambient polarization or special
divisor.  This separate obstruction is proved in the
[orthogonal-Shimura audit](orthogonal_shimura_maximal_reflection_pushback_no_go.md).

### 4.4 Mirror doubling is not an algebraic branched cover

The map

\[
        \mathbf H^n/W^+\longrightarrow \mathbf H^n/W=P              \tag{4.3}
\]

branches along the mirror facets, which have real codimension one.  A
holomorphic branch locus is a complex-analytic subset and therefore has
even real codimension at each smooth generic point; it cannot be a real
hypersurface.  (For a singular normal target a finite map may be
quasi-etale and branch only in higher complex codimension, so purity in
complex codimension one is not asserted here.)  Thus (4.3) is a real
reflection double, not a complex algebraic branched cover.

In the oriented quotient, a generic ridge stabilizer is cyclic and is
locally compatible with root-stack inertia.  This local resemblance does
not provide a global complex structure; Theorem 3.1 rules it out on the
compact \(\mathbf H^6\) control.  A double cover or root stack of a separate
complexified hyperplane arrangement is an external arrangement carrier,
not an algebraization of (4.3).

### 4.5 Rigidity isolates a fixed group, not the global class

Mostow rigidity implies that, once a lattice \(\Gamma<SO(1,n)\) with
\(n\ge3\) is fixed, its hyperbolic realization is isolated.  This is not a
finite-type theorem for all Coxeter inputs.  The presentation, dimension,
field, covolume, number of mirrors, and integral structure all vary.
Character-variety or local-system moduli become finite type only after
fixing precisely those data.

The natural known replacement for genus-area in higher dimension is
analytic: conformal volume plus a congruence spectral gap bounds
hyperbolic covolume.  That is the ABSW route.  It is intentionally not used
here and cannot be relabelled as an AG theorem.

## 5. Exact verdict and the missing theorem

The quotient-orbifold route has the following terminal ledger.

| proposed step | result |
|---|---|
| \(P\mapsto[\mathbf H^n/W^+]\) as a real-analytic orbifold | exact |
| coarse oriented compactification | \(S^n\), hence genus zero only when \(n=2\) |
| complex algebraic DM uniformization for all maximal inputs | impossible by Theorem 3.1 |
| arbitrary real algebraic model | exists too generally to give a bound |
| Baily--Borel after congruence | unavailable because (4.2) is non-Hermitian |
| projective branched-cover/root-stack replacement | not the quotient; direct nonzero-degree cover is obstructed |
| coarse topology or low cohomology bounds volume/facets | false without a new maximality theorem, by (4.1) |
| new bound on \(n\), \([k:\mathbf Q]\), facets, or chambers | not obtained |

Because direct algebraization is impossible, a successful AG proof must use
an external exact carrier.  The minimal missing bridge is therefore:

> **Open theorem QVOL-AG (algebraic volume comparator).**  Construct for
> every absolutely reflection-maximal arithmetic chamber an exact
> real-framed complex algebraic carrier \((\mathscr A_P,Z_P,C_P)\) and a
> polarized algebraic invariant \(D(P)\) such that
>
> 1. the carriers and decorations with \(D(P)\le C\) form one finite-type
>    stack;
> 2. the coarse-sphere/mirror condition together with reflection maximality
>    gives an absolute bound \(D(P)\le C\); and
> 3. \(D(P)\) bounds hyperbolic covolume and a Gram-generating arithmetic
>    packet, hence dimension, field degree, and facets.

No current Baily--Borel, Borel--Serre, KSBA, K-moduli, root-stack,
character-variety, or real-algebraization theorem supplies clause 2 or 3.
The conformal-volume/spectral theorem supplies an analytic substitute, but
that is exactly the excluded classical reflection-finiteness route.

Thus the H2 orbifold-curve theorem is a genuine exceptional isomorphism,
not the first case of a direct quotient-AG proof in every dimension.

## 6. Primary sources and exact controls

1. J. A. Carlson and D. Toledo,
   [*Harmonic mappings of Kahler manifolds to locally symmetric spaces*](https://www.numdam.org/item/PMIHES_1989__69__173_0.pdf),
   Publ. Math. IHES 69 (1989), 173--201: the rank/homology obstruction and
   the non-Kahler-group consequence for compact real-hyperbolic lattices.
2. W. L. Baily, Jr. and A. Borel,
   [*Compactification of arithmetic quotients of bounded symmetric domains*](https://annals.math.princeton.edu/1966/84-3/p11),
   Ann. of Math. 84 (1966), 442--528: the Hermitian/bounded-domain
   algebraic compactification theorem.
3. A. Borel and J.-P. Serre,
   [*Corners and arithmetic groups*](https://eudml.org/doc/139559),
   Comment. Math. Helv. 48 (1973), 436--483: the real
   manifold-with-corners compactification.
4. A. Tognoli,
   [*Su una congettura di Nash*](https://www.numdam.org/item/ASNSP_1973_3_27_1_167_0/),
   Ann. Scuola Norm. Sup. Pisa 27 (1973), 167--185: real algebraization of
   compact smooth manifolds.
5. The compact arithmetic control and its exact index-two lineage are in
   [`p6_10_bugaenko_double_lineage_and_classification_audit.md`](p6_10_bugaenko_double_lineage_and_classification_audit.md).
6. The unbounded arithmetic coarse-sphere control is
   [`k3_s3_h10_redoubling_family.md`](k3_s3_h10_redoubling_family.md).
7. I. Agol, M. Belolipetsky, P. Storm, and K. Whyte,
   [*Finiteness of arithmetic hyperbolic reflection groups*](https://arxiv.org/abs/math/0612132),
   Groups Geom. Dyn. 2 (2008), 481--498: cited only to identify the excluded
   conformal-volume/spectral substitute, not used in any proof above.
