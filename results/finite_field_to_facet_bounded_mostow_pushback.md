# From finitely many Vinberg fields to chamber finiteness: the exact Mostow reduction

**Status (2026-08-12): unconditional downstream theorem and sharp remaining
gap.  No reflective-lattice finiteness theorem of Nikulin or
Agol--Belolipetsky--Storm--Whyte is used.**

## 1. Outcome

Suppose a class of finite-volume arithmetic Coxeter chambers in dimensions
at least three has only finitely many Vinberg ground fields and has bounded
dimension.  The remaining global problem is smaller than a simultaneous
facet/height/bad-reduction/lattice problem.

The following single additional bound is enough:

\[
                    \#\operatorname {Facets}(P)\le F.                \tag{1.1}
\]

Indeed, a finite field set gives a finite alphabet of finite Coxeter edge
labels.  A facet bound therefore gives finitely many Coxeter matrices.  A
fixed Coxeter matrix gives a fixed marked abstract Coxeter group, and
Mostow--Prasad rigidity says that this group has at most one faithful,
discrete, finite-covolume reflection realization in a fixed
\(\mathbf H^n\), up to conjugacy and the finite choice of sides of its
mirrors.  Thus (1.1) gives finitely many exact normalized Grams and chambers.

In particular, once (1.1) is known, there is **no independent dotted-height
or bad-prime theorem left to prove**.  Those bounds follow by taking a
maximum over the resulting finite chamber set.  For a fixed arithmetic
reflection representation, Jordan--Zassenhaus also leaves only finitely many
native invariant lattices up to similarity.

What is not proved is (1.1) for absolutely maximal arithmetic reflection
groups.  Fixed fields, bounded dimension, stable arrangement moduli,
integrality, and quotient-stack finite type do not imply it.  Exact
redoubling families keep the field, dimension, arithmetic lattice and root
alphabet fixed while their facets and dotted heights tend to infinity.
They are proper finite-index reflection subgroups, so absolute reflection
maximality is precisely the hypothesis that must be converted into (1.1).

## 2. Finite fields give a finite Coxeter alphabet

For an intersecting pair of facets with angle \(\pi/m\), the length-two
Vinberg cycle is

\[
 c_m=4\cos ^2(\pi/m)=2+2\cos(2\pi/m).                                \tag{2.1}
\]

It lies in the minimal Vinberg ground field \(k(P)\).  For \(m\ge3\),

\[
 \mathbf Q(c_m)=\mathbf Q(\zeta_m+\zeta_m^{-1}),\qquad
 [\mathbf Q(c_m):\mathbf Q]={\varphi(m)\over2}.                      \tag{2.2}
\]

The case \(m=2\) is rational and is treated separately.  Consequently, if
\([k(P):\mathbf Q]\le D\), then every finite edge label satisfies

\[
                             \varphi(m)\le2D.                         \tag{2.3}
\]

Since \(\varphi(m)\to\infty\), only finitely many \(m\)'s occur.  Thus a
finite field set is stronger than needed here: a uniform field-degree bound
already gives a finite alphabet

\[
              \mathcal A_D=\{m\ge2:\varphi(m)\le2D\}\cup\{\infty\}. \tag{2.4}
\]

The one symbol \(\infty\) deliberately covers both parallel and
ultraparallel pairs.  Their metric distinction is not part of the abstract
Coxeter presentation; it will be recovered by rigidity.

## 3. Fixed presentation gives at most one complete realization

Let \(M=(m_{ij})_{1\le i,j\le f}\) be a Coxeter matrix and write

\[
 W_M=\left\langle s_1,\ldots,s_f\ \middle|\
 s_i^2=1,\ (s_is_j)^{m_{ij}}=1\text{ when }m_{ij}<\infty
 \right\rangle .                                                     \tag{3.1}
\]

### Theorem 3.1 (marked finite-volume rigidity)

Fix \(n\ge3\).  Any two faithful representations

\[
 \rho_a:W_M\longrightarrow\operatorname {Isom}(\mathbf H^n),
 \qquad a=1,2,                                                     \tag{3.2}
\]

whose images are discrete of finite covolume are conjugate in
\(\operatorname {Isom}(\mathbf H^n)\).  This applies in particular when
the \(\rho_a(s_i)\) are the side reflections of finite-volume Coxeter
chambers.

#### Proof

Let

\[
 \chi_a=\operatorname {or}\circ\rho_a:W_M\longrightarrow\{\pm1\}
 \qquad(a=1,2)
\]

be the two orientation characters.  Choose a torsion-free finite-index
subgroup, take its normal core in \(W_M\), and intersect that core with
\(\ker\chi_1\cap\ker\chi_2\); call the result \(H\).  Coxeter groups are
finitely generated linear, so the first subgroup exists by Selberg's lemma.
The resulting \(H\) is torsion-free, normal and finite index in \(W_M\), and
faithfulness makes both \(\rho_a(H)\) torsion-free finite-covolume lattices in
\(\operatorname {Isom}^+(\mathbf H^n)\).

For the reflection realizations used below, both \(\chi_a\) are just Coxeter
word parity, so one may equivalently start with the canonical parity kernel
\(W_M^+\).  The common-orientation-kernel formulation is needed for the
slightly broader statement of Theorem 3.1 concerning arbitrary faithful
lattice representations.

Mostow--Prasad rigidity applied to the abstractly specified isomorphism

\[
 \rho_1(h)\longmapsto\rho_2(h)\qquad(h\in H)                          \tag{3.3}
\]

gives \(g\in\operatorname {Isom}(\mathbf H^n)\) with

\[
                  g\rho_1(h)g^{-1}=\rho_2(h)\quad(h\in H).           \tag{3.4}
\]

For \(s\in W_M\), normality of \(H\) shows that

\[
 \rho_2(s)^{-1}g\rho_1(s)g^{-1}
       \in C_{\operatorname {Isom}(\mathbf H^n)}(\rho_2(H)).          \tag{3.5}
\]

The centralizer of a finite-covolume hyperbolic lattice in the projective
isometry group is trivial.  For example, this follows from its full limit
set, or from Borel density and the trivial center of
\(\operatorname {PO}(n,1)\).  Hence (3.5) is the identity for every \(s\),
and \(g\rho_1g^{-1}=\rho_2\) on all of \(W_M\).  \(\square\)

For a Coxeter polytope, the Poincare polyhedron theorem gives exactly the
presentation (3.1); there are no hidden extra relations in its side group.
After Theorem 3.1 the **marked homomorphism**, hence every labelled
reflection hyperplane, agrees up to conjugacy.  There are at most \(2^f\)
choices of oriented halfspaces bounded by those hyperplanes, and the
nonempty irredundant Coxeter fundamental-chamber condition retains only a
subset of them.  Thus a fixed \((n,M)\) contributes at most one conjugacy
class of marked reflection homomorphisms and finitely many marked oriented
chambers.  Forgetting the labels and/or the side orientations only takes a
quotient of that finite set; no claim is made that all \(2^f\) choices are
valid chambers.

### Caveat: do not replace global rigidity by a character-variety slogan

The scheme \(\operatorname {Hom}(W_M,\operatorname {PO}(n,1))\) is of
finite type, but the complete finite-covolume locus is not a Zariski-closed
finite-type component.  In particular, a cusped hyperbolic three-manifold's
complete holonomy can lie on a positive-dimensional character-variety
component whose nearby points are incomplete Dehn-surgery deformations.
Thus the argument

```text
finite-type character variety + "local rigidity" = finitely many points
```

is not valid uniformly in the noncompact \(n=3\) case.  Theorem 3.1 instead
uses the global Mostow--Prasad statement directly and only compares faithful,
discrete, **complete finite-covolume** representations.  This also handles
cusps and orbifold torsion.

## 4. The finite-presentation chamber theorem

### Theorem 4.1

Let \(\mathfrak C\) be any class of finite-sided finite-volume arithmetic
Coxeter chambers with

\[
 n(P)\in\mathcal N\subset\{3,4,\ldots\},\qquad
 k(P)\in\mathcal K,\qquad f(P)\le F,                                 \tag{4.1}
\]

where \(\mathcal N\) and \(\mathcal K\) are finite.  Then
\(\mathfrak C\) contains only finitely many unmarked normalized chambers and
only finitely many conjugacy classes of side-reflection groups.

#### Proof

Put \(D=\max_{k\in\mathcal K}[k:\mathbf Q]\).  Section 2 gives a finite
edge-label alphabet.  For each \(f\le F\), only finitely many symmetric
Coxeter matrices use that alphabet.  The side-reflection group of every
chamber with one such matrix is the marked group \(W_M\) in (3.1).  Apply
Theorem 3.1 for each of the finitely many pairs \((n,M)\).  \(\square\)

### Corollary 4.2 (all remaining escapes use presentation complexity)

Fix finitely many fields and dimensions \(n\ge3\).  Any infinite sequence
of pairwise nonisometric finite-volume arithmetic Coxeter chambers has
unbounded facet count.  After passing to a subsequence, the Coxeter
presentations, rather than merely the dotted entries inside one
presentation, must escape.

Consequently a facet bound automatically bounds:

* every normalized dotted Gram entry and its Weil height;
* the intrinsic finite bad-prime set of the normalized marked arrangement;
* chamber covolume; and
* every other invariant determined by the exact normalized chamber.

These are consequences of finite enumeration, not a priori inequalities.

### Corollary 4.3 (genuine bounded-AG-stratum finiteness)

Fix finitely many dimensions \(n\ge3\) and a degree bound for the Vinberg
fields.  On any finite union of exact quadric/hyperplane arrangement strata
with a bounded number of marked hyperplane components, the finite-volume
arithmetic Coxeter return image is finite.

In particular this applies to a finite list of labelled stable-pair Hilbert
types.  For the standard exact carrier the numerical identity

\[
 K_{\mathbf P^n}+Q+\frac12\sum_{i=1}^fH_i
       =\left(\frac f2-n+1\right)H                                  \tag{4.2}
\]

shows directly that the boundary Hilbert data records \(f\).  The stable
pair strata themselves may be positive-dimensional and have infinitely many
algebraic points; Theorem 4.1 says that their **complete finite-volume
Coxeter return image** is finite.  This is an actual AG-bounded-stratum
pushback theorem.  Its missing global hypothesis is that absolute reflection
maximality places every carrier in finitely many such strata.

## 5. Native integral lattices after the chamber is fixed

There is also no residual infinite native-lattice problem for a fixed
arithmetic reflection representation, once lattices are counted up to the
similarity that leaves the hyperbolic reflections unchanged.

### Proposition 5.1 (Jordan--Zassenhaus closure)

Let \(\rho:W_M\to O(q,k)\) be one of the fixed arithmetic representations
above and let \(L_0\subset V\) be a full \(\mathcal O_k\)-lattice preserved
by \(\rho(W_M)\).  Then the full \(\rho(W_M)\)-stable
\(\mathcal O_k\)-lattices in \(V\) form only finitely many homothety classes.

#### Proof

The algebra

\[
 \Lambda=\mathcal O_k[\rho(W_M)]\subseteq\operatorname {End}_{\mathcal O_k}(L_0)
                                                                         \tag{5.1}
\]

is an \(\mathcal O_k\)-order in the finite-dimensional algebra generated by
the representation: it is an \(\mathcal O_k\)-submodule of the finite module
\(\operatorname {End}_{\mathcal O_k}(L_0)\), hence is finite, and its
\(k\)-span is the full representation algebra.  Jordan--Zassenhaus gives
finitely many isomorphism classes of \(\Lambda\)-lattices whose fixed
rational module is \(V\).  A finite-covolume hyperbolic lattice is Zariski
dense.  Consequently the standard module is absolutely irreducible; after
extending scalars to an algebraic closure, Burnside's theorem, followed by
descent, gives

\[
 k[\rho(W_M)]=\operatorname {End}_k(V).                               \tag{5.2}
\]

A \(\Lambda\)-module isomorphism between two full lattices therefore extends
to a scalar in \(k^\times\), proving homothety.  \(\square\)

This is ordinary order-module finiteness **after one rational
representation has been fixed**.  It is not Nikulin's theorem that fixed
\((k,n)\) leaves finitely many reflective quadratic lattices.  Arbitrary
rescaling of an integral model can of course introduce artificial primes;
bad-prime finiteness is meant intrinsically, or after choosing one primitive
representative of each similarity class.

## 6. Why the facet bound does not follow from the proposed soft inputs

### 6.1 Exact same-field, same-lattice facet escape

The frozen redoubling families give the strongest relevant control.  For
example, the project has finite-volume noncompact Coxeter chambers in
\(\mathbf H^{10}\), all over \(\mathbf Q\), in the one fixed lattice

\[
 F=U\oplus E_8(-2)\oplus\langle-2\rangle,                             \tag{6.1}
\]

with one fixed primitive root norm/divisibility alphabet and

\[
                       f(P_N)=N+12.                                  \tag{6.2}
\]

For the archived moving roots,

\[
 q_{2m}=c+2m^2a+2m(m-1)b,qquad
 u(c,q_{2m})=(m^2-1)^2.                                               \tag{6.3}
\]

Thus dotted height is unbounded as well.  Taking \(m=p+1\) for odd primes
and \(\lambda=u/(u+1)\) gives \(v_p(\lambda)=2\), so the associated marked
four-point/Gram--Belyi packet acquires a new bad prime.  This is a statement
about that exact marked carrier; it is not being promoted to a claim about
every possible integral model of the native hyperplane arrangement.

The even larger frozen \(\mathbf H^{12}\) finite-subtree family stays in one
rank-thirteen Picard lattice and one arithmetic commensurability class while
having unbounded facet count and exponentially many chambers as the index
grows.

Every descendant in these families is a proper finite-index reflection
subgroup of its seed.  They therefore do not contradict the desired theorem
for absolutely maximal reflection groups.  They prove exactly that fixed
field, bounded dimension, finite volume, integrality, a fixed native lattice,
and even a fixed AG divisor-lattice carrier do not provide (1.1).

### 6.2 Fixed field and dimension do not bound arithmetic hosts

The failure occurs even before asking for reflections.  Fix
\(k=\mathbf Q(\sqrt5)\), select the embedding with \(\sqrt5>0\), and in rank
twelve put

\[
 q_p=\langle-p\sqrt5,1,\ldots,1\rangle                               \tag{6.4}
\]

for rational primes \(p\ne5\).  At the selected embedding \(q_p\) has
signature \((11,1)\); at the other embedding it is positive definite.
Hence it is anisotropic over \(k\), and its arithmetic orthogonal groups are
cocompact lattices in \(\operatorname {Isom}(\mathbf H^{11})\).

The rank is even, so determinant squareclass is invariant under similarity.
The squareclasses

\[
                         \det(q_p)=-p\sqrt5                           \tag{6.5}
\]

are distinct: the ratio \(p/q\) has odd valuation at a prime above \(p\).
Thus the \(q_p\) give infinitely many noncommensurable compact arithmetic
classes over one field and in one dimension, with varying finite bad primes
and discriminants.  Enlarging each arithmetic group to a maximal arithmetic
lattice does not change this conclusion.

These groups are not asserted to be reflection generated.  Their role is to
show that Mostow rigidity, fixed field, fixed dimension, compactness, and
maximality as an **arithmetic lattice** do not bound the ambient quadratic
form.  Reflection maximality must contribute genuinely new information.

### 6.3 Stable arrangements and quotient stacks

For the exact quadric/hyperplane carrier,

\[
 K_{\mathbf P^n}+Q+\frac12\sum_{i=1}^fH_i
             =\left(\frac f2-n+1\right)H.                            \tag{6.6}
\]

As \(f\) grows, the Hilbert data and, when positive, the log volume grow.
The universal target

\[
             \coprod_{f\ge n+1}[\operatorname {Arr}_{n,f}/PGL_{n+1}] \tag{6.7}
\]

is an ind-stack, not one quasi-compact finite-type stack.  Quotienting by
wall permutations does not change this.  The H10/H12 families run through
infinitely many components of (6.7).

Absolute reflection maximality means that no larger discrete reflection
group contains the side group.  In the present exact arrangement carrier,
that is not known to be a constructible, proper, or quasi-compact condition.
The side group also does not preserve the finite boundary arrangement: a
simple reflection sends most simple mirrors to new mirrors in the infinite
orbit.  Thus ordinary stable-pair automorphism finiteness does not encode
reflection-overgroup maximality.

### 6.4 Integral points, Northcott and Shafarevich

The redoubling cycles (6.3) are algebraic integers over the fixed field
\(\mathbf Q\) and have unbounded height.  Fixed field and denominator one
therefore do not activate Northcott.  The prime specialization above shows
that no fixed finite bad set activates a Shafarevich theorem.  More
importantly, Theorem 4.1 shows that these are not independent failures in
dimensions at least three: after the field and dimension are fixed, an
infinite complete Coxeter sequence can escape only by changing its Coxeter
presentation, hence by unbounded facets.

### 6.5 Mostow rigidity does not bound presentations

Mostow--Prasad says that **one** abstract lattice has one complete
hyperbolic realization.  It says nothing about how many abstract Coxeter
groups \(W_M\) occur.  The redoubling examples change \(M\) at every stage,
which is exactly how they evade Theorem 3.1.  Rigidity is the terminal
zero-dimensional gate, not the global boundedness theorem.

## 7. The minimal missing theorem

After the finite-ground-field result and any dimension ceiling, the global
higher-dimensional completion reduces to the following statement.

> **MAX-FAC.**  For every finite set \(\mathcal K\) of totally real fields
> and finite set \(\mathcal N\subset\{3,4,\ldots\}\), there exists
> \(F(\mathcal K,\mathcal N)\) such that every finite-volume arithmetic
> Coxeter chamber \(P\) with \(k(P)\in\mathcal K\),
> \(n(P)\in\mathcal N\), and side-reflection group maximal among arithmetic
> reflection groups has at most \(F(\mathcal K,\mathcal N)\) facets.

Here this maximality has the same force as absolute maximality among
discrete reflection groups.  Indeed, any discrete overgroup of a
finite-covolume lattice is again a finite-covolume lattice and contains it
with finite index; arithmeticity is invariant under commensurability.

Equivalently, only finitely many Coxeter presentations occur.  Theorem 4.1
then proves finite chambers and groups without a reflective-lattice theorem,
and Proposition 5.1 closes stable native lattices after a fixed
\(k\)-rational reflection representation has been selected.  The
proposition neither asserts finiteness of all possible rational forms nor
removes the diagonal rescaling freedom in primitive root/Cartan
decorations; neither issue is needed for chamber or reflection-group
finiteness.

This is strictly sharper than asking separately for facet, dotted-height,
bad-prime and lattice-discriminant bounds.  Any one of the following AG
statements would imply MAX-FAC:

1. the exact maximal arrangement image factors through a finite union of
   the stacks in (6.7);
2. a sheet-separating automorphic divisor has uniformly bounded total
   automorphic degree, since each marked facet consumes a positive unit of
   that degree; or
3. the exact Gram--Belyi carrier has uniformly bounded polarization degree,
   which already dominates a quadratic function of \(f\).

No current BAB, KSBA, Shafarevich, Northcott, automorphic-degree, or
quotient-stack theorem derives any of these bounds from absolute arithmetic
reflection maximality.  Importing a uniform reflection-orbifold covolume
bound followed by general lattice finiteness is another sufficient route,
but it is the classical spectral/arithmetic route rather than the requested
AG pushback.

For dimension two Theorem 3.1 is false and MAX-FAC is not enough by itself:
fixed Coxeter presentations can move in Teichmuller families.  The exact
orbifold-curve and Lambert arguments must remain the separate two-dimensional
terminal branches.

## 8. Primary references and exact controls

* G. D. Mostow, *Strong Rigidity of Locally Symmetric Spaces*, Annals of
  Mathematics Studies 78, Princeton University Press, 1973.
* G. Prasad, [*Strong rigidity of Q-rank 1
  lattices*](https://www.maths.gla.ac.uk/~mpowell/Prasad-on-Mostow-rigidity.pdf),
  Invent. Math. **21** (1973), 255--286.
* I. Reiner, *Maximal Orders*, London Mathematical Society Monographs,
  Theorem 26.4 (Jordan--Zassenhaus).
* D. Allcock,
  [*Infinitely many hyperbolic Coxeter groups through dimension
  19*](https://arxiv.org/abs/0903.0138), Geom. Topol. **10** (2006),
  737--758.

Project-native exact controls:

* [H10 Gram--Belyi redoubling lower bound](gram_belyi_coxeter_dessin_theorem.md);
* [H12 fixed-lattice finite-subtree family](jacobian_kummer_h12_maximal_disjoint_doubling_family.md);
* [bounded Schreier-degree theorem and same-field no-go](bounded_degree_schreier_coxeter_dessin_theorem.md);
* [facet-sheet automorphic-degree injection](facet_sheet_automorphic_degree_injection.md); and
* [exact global arrangement carrier](exact_agcat_functor_and_ag_finiteness_pushback_synthesis.md).
