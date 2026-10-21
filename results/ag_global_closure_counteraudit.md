# Counteraudit of the proposed global AG closure

**Date:** 2026-08-12  
**Verdict:** the exact real-framed arrangement functor survives, but the
current two-gap formulation does not.  In its stated form it omits an
independent dimension bound.  No alternative orbifold, stable-pair,
VHS/local-system, K3/derived, motivic, automorphic, or Shimura theorem found
in this audit closes that gap or the facet gap without assuming essentially
the desired boundedness.

This audit forbids the Agol--Belolipetsky--Storm--Whyte theorem, Nikulin's
reflective-lattice finiteness theorem, Vinberg's global arithmetic dimension
cutoff, and any other pre-existing global dimension/facet finiteness theorem.
The issue is whether algebraic geometry supplies the missing implication.

## 1. What an exact first arrow really means

For a marked chamber with minimal Vinberg field `k`, Lorentz space `(V,q)`,
ordered inward root rays `r_i`, and selected time component, put

\[
 F_{\rm arr}(P)=
 \bigl(\mathbf P(V),Q_q,\mathbf P(r_1^\perp),\ldots,
 \mathbf P(r_f^\perp),\sigma_0,C_P\bigr).                 \tag{1.1}
\]

The smooth quadric determines `q` up to scale, and quadric polarity returns
the root lines.  The selected real component and inward orientations return
the signs, hence the complete normalized Gram, all finite angles, parallel
relations, ultraparallel distances, and the same half-space intersection.
Thus there is an exact return on the image.

There are two qualifications.

1. This is a functor to **real-framed algebraic arrangements**, not to bare
   complex schemes or stacks.  The selected real component, uniformizing
   sheet, and half-space orientations are needed for exact return.  Forgetting
   them loses the chamber.
2. As `n` and `f` vary, the universal target is the ind-stack

   \[
     \coprod_{n,f}[\operatorname{Arr}_{n,f}/\operatorname{PGL}_{n+1}],
                                                               \tag{1.2}
   \]

   not one quasi-compact finite-type stack.  The exact first arrow therefore
   does not itself supply boundedness.

The real framing is finite after `n,f` are fixed, so it is not the principal
finiteness problem.  But it must be stated: a purported exact functor to a
bare AG category is not currently proved.

## 2. The actual downstream implication

The valid rigidity theorem is the following.

> **Finite presentation theorem.**  In dimensions `n >= 3`, a finite set of
> Vinberg fields, a finite set of dimensions, and a uniform facet bound imply
> finiteness of finite-volume Coxeter chambers and side-reflection groups.

Indeed, if `[k:Q] <= D`, every finite Coxeter label `m` satisfies

\[
  4\cos ^2(\pi/m)\in k,
  \qquad [\mathbf Q(4\cos ^2(\pi/m)): \mathbf Q]
       =\varphi(m)/2,                                       \tag{2.1}
\]

so only finitely many labels occur.  A facet bound gives finitely many
Coxeter matrices.  For a fixed matrix `M`, the Poincare polyhedron theorem
identifies the side group with the abstract Coxeter group `W_M`.
Mostow--Prasad rigidity, applied to a common torsion-free normal core, and
the trivial centralizer of a finite-covolume lattice show that two complete
faithful finite-covolume representations of `W_M` in the same `H^n` are
conjugate.  Hence a fixed `(n,M)` returns at most one marked reflection
representation and finitely many choices of sides.

No separate dotted-height, bad-prime, or native-lattice bound is needed to
count the real chambers/groups after this step.

### 2.1 The hidden third gap

The existing `MAX-FAC` statement is quantified as follows:

\[
 \text{for finite }\mathcal K\text{ and finite }\mathcal N,
 \quad f(P)\le F(\mathcal K,\mathcal N).                    \tag{2.2}
\]

The new field funnel proves a finite field set for `n >= 6`, and the
low-dimensional descent would extend this to all dimensions if compact
`H^3` field finiteness were proved.  It does **not** bound `n`.  A finite
field set does not imply a dimension bound; the published narrow-face
arguments also append a separate dimension theorem at exactly this point.

Consequently

```text
compact H3 field finiteness + MAX-FAC as stated
```

is not sufficient for the all-dimensional theorem.  One must add

\[
  \boxed{\text{DIM-AG: } n(P)\le N}                         \tag{2.3}
\]

or strengthen (2.2) to a dimension-free statement

\[
 \boxed{\text{MAX-FAC}^{\rm glob}:\quad
 f(P)\le F(\mathcal K)\text{ in every dimension}.}          \tag{2.4}
\]

The stronger statement includes the dimension bound because every
full-dimensional finite-volume polytope has at least `n+1` facets.

Thus the honest remaining ledger is either

```text
H3-FIELD + DIM-AG + conditional MAX-FAC,
```

or

```text
H3-FIELD + dimension-free MAX-FAC^glob.
```

## 3. A precise sufficient AG theorem

There is a useful way to combine the dimension and facet gaps without
hiding either one.

> **AG-DEG theorem (conditional implication).**  Assume the global field
> set is finite.  Suppose every exact carrier has a projective
> sheet-separated model `(X_P,L_P)` and an effective divisor
>
> \[
>  Z_P=\sum_{i=1}^{f(P)}q_i Z_i+Z_{\rm aux},\qquad q_i\ge q_{\min}>0,
>                                                               \tag{3.1}
> \]
>
> such that distinct facets give distinct prime divisors, `L_P` is an ample
> integral Cartier polarization (or stack denominators are uniformly
> bounded), and
>
> \[
>   {Z_P\cdot L_P^{\dim X_P-1}\over q_{\min}}\le B           \tag{3.2}
> \]
>
> with one absolute `B`.  Then all maximal arithmetic Coxeter chambers and
> reflection groups are finite in number.

On a scheme, every prime divisor has positive integral `L`-degree.  Hence

\[
  f(P)q_{\min}
  \le \sum_iq_i(Z_i\cdot L^{r-1})
  \le Z_P\cdot L^{r-1},                                    \tag{3.3}
\]

so `f(P) <= B` and `n(P) <= B-1`.  Equation (2.1) gives a finite angle
alphabet, and the finite-presentation theorem finishes dimensions at least
three.  Dimension two is the already separate orbifold-curve theorem.

The hypotheses in (3.1)--(3.2) are not bookkeeping.  Each is essential:

* without sheet separation, many real facets may be sheets of one complex
  divisor;
* on a DM stack, unbounded stabilizer denominators destroy the unit lower
  degree in (3.3);
* normalizing degree by an unbounded separating-cover index destroys the
  sheet count; and
* an effective divisor class without a uniform upper degree gives no bound.

`AG-DEG` is the exact algebraic theorem that would achieve the requested
pushback.  What is missing is the implication

\[
 \text{arithmeticity + absolute reflection maximality}
 \Longrightarrow (3.2).                                    \tag{3.4}
\]

No theorem located in this audit proves (3.4).

There is a second, genuinely different sufficient version which bypasses
the ground-field problem altogether.

> **AG-DEG-INERTIA theorem (conditional implication).**  Replace the finite
> field-set hypothesis by the assertion that the exact carrier remembers the
> two-wall Coxeter label as a distinguished dihedral inertia subgroup and
> that all such inertia groups have order at most `I`.  Together with
> (3.1)--(3.2), this already implies global Coxeter finiteness.

Indeed, (3.3) bounds `f` and `n`, while distinguished inertia
`I_2(m)` has order `2m`, so all finite labels satisfy `m <= I/2`.  There are
then finitely many Coxeter matrices, and Mostow--Prasad finishes in
dimensions at least three.  The separate orbifold-curve theorem finishes
dimension two.  In particular the fields become finite as an **output** of
the finite returned Gram list; compact `H^3` field finiteness need not be
proved first.

A quasi-compact finite-type separated DM target with finite inertia supplies
a uniform inertia-order bound.  Thus a single finite-type, sheet-separated,
degree-bounded, label-faithful DM envelope would prove the user's desired
theorem directly.  None of the current carriers has all these properties:
stable-pair schemes do not encode `m` in inertia, while the orthogonal
reflection stacks do encode it but vary in lattice, level, dimension, and
absolute degree.

## 4. Alternative AG categories

### 4.1 The hyperbolic quotient/orbifold

The real orbifold `[H^n/W^+]` is an exact carrier, and its coarse double is
topologically `S^n`.  For `n=2` it is an algebraic orbifold curve and the
genus-zero arithmetic Fuchsian theorem gives genuine finiteness.

For `n>=3`, `SO(1,n)/SO(n)` is not Hermitian.  More decisively, a compact
arithmetic reflection-maximal `H^6` example has a torsion-free finite cover
which is a closed real-hyperbolic six-manifold.  If the orbifold were a
smooth proper complex algebraic DM stack, that cover would be a compact
Kahler group, contradicting Carlson--Toledo.  Singular projective models do
not repair the method: projective SNC varieties can realize every finitely
presented group, so allowing such singularities removes the desired
finiteness leverage rather than supplying it.

Real algebraization is also too flexible: Nash--Tognoli algebraizes every
closed smooth manifold as a real locus.  It places no bound on fields,
dimension, or facets.

### 4.2 Stable pairs, KSBA, and K-moduli

The quadric/hyperplane carrier is an exact stable-pair carrier on suitable
strata.  But KSBA boundedness fixes dimension, coefficients, and numerical
data such as volume.  Here both dimension and the boundary Hilbert data vary;
for the standard carrier

\[
 K_{\mathbf P^n}+Q+\frac12\sum_{i=1}^fH_i
   =\left(\frac f2-n+1\right)H.                            \tag{4.1}
\]

Thus the numerical data already record the unbounded quantity `f`.

There is also a sharp fixed-data control: one fixed-volume KSBA curve of
exact conic/four-line pairs contains a positive-dimensional cofinite locus
of absolutely reflection-maximal Lambert chambers, and the corresponding
hyperplane pair is uniformly K-stable with constant K-moduli point.
Arithmeticity cuts that particular curve to finitely many points, but only
through a Vinberg integral-point/norm argument.  Therefore stable-pair or
K-stability by itself gives neither discreteness nor (3.4).

### 4.3 Local systems and variations of Hodge structure

On the natural complement of the complexified quadric/hyperplane carrier,
meridians about two transverse boundary hyperplanes commute.  Side
reflections with Coxeter label `m>2` do not commute.  Hence the natural
arrangement complement cannot carry the Coxeter representation as its
monodromy.

Replacing an SNC crossing by a dihedral stack chart can impose the relation,
but then the inertia has order `2m`.  A finite-type DM family with finite
inertia has bounded inertia degree.  Thus a uniform such construction would
already assume the required label/complexity bound.  Moreover, an abstract
orthogonal local system is not automatically a polarizable geometric VHS;
fixed-rank Hodge moduli require precisely the rank bound being sought.

### 4.4 K3, hyperkahler, derived, and motivic carriers

An intrinsic cohomological action carrier would be powerful.  If the full
restriction-of-scalars trace lattice acts faithfully on one fixed numerical
cohomology lattice of rank `R`, then

\[
                 [k:\mathbf Q](n+1)\le R.                    \tag{4.2}
\]

For K3 Neron--Severi or Mukai lattices this gives a genuine conditional
dimension/field bound (`R=20` or `24`, with a sharper signature capacity in
the supersingular action setup).

But there is no universal compiler from arbitrary arithmetic Coxeter data to
such a carrier.  It would have to prove a primitive trace-lattice embedding,
discriminant gluing, a common period preserved by every reflection, geometric
or Fourier--Mukai lifts, and exact return with no additional walls.  These
conditions fail automatically for many ranks/discriminants and are open in
the proposed mixed-field lane.

A quiver or dg category can be manufactured from a Coxeter presentation, but
then its numerical rank and relations grow with the presentation.  Without a
smooth-proper geometric realization in a bounded class, this is another exact
encoding, not an AG finiteness theorem.  Likewise, no theorem says that an
arbitrary orthogonal Coxeter representation is motivic; allowing unbounded
Betti rank again loses boundedness.

### 4.5 Orthogonal/unitary Shimura and automorphic carriers

The determinant-line type-IV construction and the CM unitary-ball
construction give exact real-framed algebraic carriers.  In both, the real
hyperbolic group lies in a proper `SO(1,n)` subgroup and has infinite index
in the full arithmetic group defining the Shimura variety.  The hyperbolic
space is a totally real fixed locus, not an algebraic Shimura subvariety.

Reflection maximality inside `SO(1,n)` therefore does not imply that the
ambient orthogonal/unitary modular variety is reflection generated, or that
it has a pure reflective modular form of bounded slope.  The standard
automorphic-correction implication is actually false: the complete
finite-volume maximal chamber for

\[
 U\oplus E_8(-1)\oplus E_7(-1)
\]

has no automorphic correction in Wang's classified sense.

Compensated Borcherds divisors can be constructed after adding auxiliary
components disjoint from the real slice.  On a neat sheet-separated compact
carrier they yield the valid inequality

\[
 f(P)q_{\min}\le w\,c_1(\lambda)^n.                         \tag{4.3}
\]

This is precisely (3.3), not its missing upper bound: the weight, level,
separating-cover degree, and automorphic volume can all grow with `P`.
Neither Baily--Borel nor reflective-form classification turns real
reflection maximality into a uniform bound on the right side.

Andre--Oort also does not supply the missing step.  It controls Zariski
closures of special points in a **fixed** Shimura datum; it does not make the
special points finite, and here the field, lattice, level, dimension, and
divisor data all vary.  Lawrence--Venkatesh and Shafarevich-type theorems
likewise require fixed family/numerical type, field, and bad set, and usually
give non-density rather than the global finite count required here.

## 5. The alternative field--Hodge--frame closure

There is a second route which avoids `MAX-FAC` completely.  It is logically
sound as a **conditional hybrid theorem**, but two of its proposed inputs are
not currently proved in the required form and its terminal step is a
reconstruction of Nikulin's arithmetic argument, not a new AG boundedness
theorem.

> **FIELD--DIM--FRAME theorem (conditional implication).**  Let
> `Cox_max_arith` be the class of finite-sided, finite-volume arithmetic
> Coxeter chambers in the standard admissible Vinberg lattice-form model,
> of dimension at least three, and maximal among discrete reflection groups
> (or maximal among arithmetic reflection groups).  Assume:
>
> 1. the integral-point theorem on the exact arrangement carrier produces a
>    finite set of minimal Vinberg fields, including the compact `H^3`
>    branch;
> 2. an independent theorem gives `n(P) <= N` for every finite-volume Coxeter
>    chamber; and
> 3. every complete reflective-lattice chamber admits a full-rank connected
>    root frame with
>    \[
>      c_{ij}={4(\delta_i,\delta_j)^2\over
>       (\delta_i,\delta_i)(\delta_j,\delta_j)}<62^2.       \tag{5.1}
>    \]
>
> Then the objects of `Cox_max_arith` have only finitely many unmarked
> normalized chamber classes.

Dimension two is handled by the separate algebraic orbifold-curve theorem;
it is not silently included in the bounded-frame step.

Here is the complete implication, with no fixed-data finiteness theorem used
as a black box.  Arithmeticity supplies an admissible `k`-form and a full
`W(P)`-stable integral `O_k`-lattice `L`: take a normal core of the
intersection with an arithmetic lattice and sum its finitely many translates.
Since both `W(P)` and `O(L)` have finite covolume,

\[
             W(P)\le R(L)\le O(L)                              \tag{5.2}
\]

are finite-index inclusions.  Thus `L` is reflective.  Apply (5.1) to the
complete chamber of `R(L)`, with `m=rk L=n+1`.

For fixed `k` and `m`, reflection integrality makes every `c_ij` an algebraic
integer.  At the selected embedding it lies in `[0,62^2)`, while at every
other real embedding definiteness and Cauchy--Schwarz put it in `[0,4]`.
The Minkowski embedding of `O_k` is discrete, so there are finitely many
possible connected matrices `(c_ij)`.

For completeness, the remaining integral-frame argument has four steps.

1. Write `d_ij=(delta_i,delta_j)`.  Outside the primes dividing `2` or a
   nonzero `c_ij`, reflection divisibility and connectedness force all
   valuations of the `d_ij` to differ by one common similarity valuation.
2. At the finitely many bad primes, the same relations bound all valuation
   differences along a nonzero path of length at most `m-1`.
3. The ideal class group and the unit group modulo squares are finite.  Row
   rescaling by units absorbs square-unit freedom.  Hence only finitely many
   full-rank frame lattices `Delta=[delta_1,...,delta_m]` occur up to
   similarity.  Projective `O_k`-module and sign data contribute only finite
   class-group and unit-sign choices.
4. Reflection divisibility gives
   \[
      [L:\Delta]\le [\Delta^*:\Delta]\,[O_k:2O_k]^m,       \tag{5.3}
   \]
   so only finitely many integral superlattices `L` occur.

This is enough.  Reflection maximality in (5.2) gives `W(P)=R(L)`; in the
arithmetic-only version `R(L)` is itself arithmetic.  The complete locally
finite mirror arrangement of `R(L)` has one chamber orbit by the usual finite
gallery argument.  Thus one gets finitely many chambers from finitely many
similarity hosts.  Notice that neither a facet bound nor Mostow rigidity is
needed in this lane.

### 5.1 What is and is not circular

There is no logical circularity **if all three assumptions are proved
independently**:

* the field theorem must not assume bounded dimension, bounded facets, the
  bounded-frame theorem, or fixed-data reflective-lattice finiteness;
* the dimension theorem must not assume a finite arithmetic field list or
  any reflection-group finiteness theorem; and
* the connected-frame lemma must be reproved geometrically before the
  valuation argument.

The present compact `H^3` integral-point lane still fails the first test.  A
single dotted two-cycle has one conjugate above `4` and all other conjugates
in `[0,4]`, an interval of integer capacity one.  A new determinant coupling
or another strict-capacity coordinate is required; bounded basis size alone
does not provide it.

The presently written "minimality 14 + h-vector" lemma also fails the second
test.  In the compact branch it obtains a small dimension only after assuming
that every finite vertex-diagram component along the narrow facet has at most
two nodes.  Neither minimality 14 nor reflection maximality proves that
premise.  In the cusped branch the chamber need not be simple at an ideal
vertex, so the compact incidence count does not apply.  Finite fields do not
repair either defect: the spherical families `A_r`, `B_r`, and `D_r` have
unbounded rank over a fixed finite label alphabet.

There is a rigorous universal replacement.  The Khovanskii--Prokhorov
edge-simple argument proves

\[
                         n\le 995                              \tag{5.4}
\]

for **every** finite-volume hyperbolic Coxeter polytope, without arithmeticity,
maximality, a field theorem, or reflective-lattice finiteness.  A
self-contained replay needs all of the following, not merely the simple
polytope `h`-vector inequality:

* the compactified chamber is almost simple, hence edge-simple; ideal vertex
  links are products of simplices;
* the Nikulin--Khovanskii face-average inequality for edge-simple polytopes;
* the good/bad three-face incidence lemma;
* the finite classifications of elliptic and quasi-Lanner Coxeter diagrams
  (quasi-Lanner rank at most ten and the required diameter bounds); and
* Prokhorov's weighted double count with `C=29/3`, giving
  `n < 96C+68=996`.

Thus (5.4) has no arithmetic-finiteness hypothesis and creates no circularity.
But importing it is exactly importing a classical global dimension theorem,
which this audit's rules forbid.  Reproving it line by line is legitimate as
an independent proof, but the result and proof architecture must still be
attributed to Khovanskii and M. N. Prokhorov.

The same attribution issue occurs in the terminal step.  The connected
bounded-frame theorem is Theorem 1 in the Appendix of Nikulin's 1981 paper;
the normalization in his 1982 Theorem 5.2.1 gives `A=62`.  Steps 1--4 above
are, essentially verbatim, the proof of that 1982 fixed-`(k,m)` theorem.
They do not *invoke its conclusion*, so the chain is not circular.  They do
reprove the forbidden theorem's core mechanism.  Calling this a new AG
native-host theorem would therefore be false.

### 5.2 Does this answer the AG-mediated question?

At its full conditional strength, yes as an **AG-mediated hybrid**, not as an
AG-only proof.  The exact quadric-arrangement functor is nondecorative if its
integral cycle-coordinate geometry proves the global finite field theorem.
The dimension step can also be expressed through polytopal Chow/Hodge theory:
Timorin derives the Khovanskii face-average estimate from the Lefschetz and
Hodge--Riemann package for polytope algebras/toric intersection cohomology.
The remaining diagram weights are Coxeter geometry, and the host step is
algebraic number theory of integral quadratic lattices.

Two qualifications prevent overstatement.

1. For a nonrational polytope, the polytopal Chow ring is a combinatorial
   Hodge object, not necessarily the Chow ring of an actual projective toric
   variety.  One must either enlarge `AGCat` to include this Hodge package or
   supply a genuine algebraic realization; merely renaming an `h`-vector
   calculation "AG" is not enough.
2. This dimension carrier is separate from the exact real-framed arrangement
   carrier.  A clean functor should therefore return both the exact
   arrangement and its canonical face-poset/Hodge object, with morphisms and
   real framing specified.

With those disclosures, the proposed architecture is mathematically honest:

```text
exact AG arrangement
  -> genuinely new all-field integral-point theorem (including compact H3)
  + independent Khovanskii--Prokhorov/Chow--Hodge dimension replay
  -> finitely many (k,m)
  -> independently reconstructed Nikulin A=62 frame arithmetic
  -> finitely many native hosts
  -> reflection maximality and one chamber orbit
  -> global finiteness.
```

As of this audit, the first arrow at compact `H^3` and the proposed new
minimality-14 dimension proof are open.  Replacing the latter by the fully
reproved Khovanskii--Prokhorov theorem makes the logic complete conditional
only on the `H^3` field theorem and the bounded-frame replay, but it no longer
satisfies a ban on classical global dimension methods.

## 6. Exact global status after the counteraudit

The proved statements are:

1. an exact functor to a real-framed algebraic quadric arrangement, with
   exact return;
2. complete AG-mediated finiteness in `H^2`;
3. finite arithmetic locus on a nontrivial exact KSBA curve;
4. a finite global ground-field set in `n>=6`, by integral arrangement
   charts plus the polygon theorem; and
5. finite fields and dimensions plus bounded facets imply chamber/group
   finiteness by Mostow--Prasad.

The requested global theorem is not yet proved.  The common first missing
input is

\[
 \boxed{\text{compact H3 ground-field finiteness}}             \tag{6.1}
\]

and one then has three possible terminal routes:

\[
 \boxed{\text{DIM-AG + MAX-FAC for the resulting finite dimensions}},\tag{6.2}
\]


\[
 \boxed{\text{dimension-free MAX-FAC}^{\rm glob},
 \text{ preferably via AG-DEG}},                              \tag{6.3}
\]

or the alternate hybrid of Section 5,

\[
 \boxed{\text{independent universal DIM + connected bounded frame
 + integral-frame replay}.}                                  \tag{6.4}
\]

No alternative AG category inspected here supplies (6.1), (6.2), or (6.3).
Route (6.4) is classically available only by the Khovanskii--Prokhorov and
Nikulin mechanisms just identified; it is not a new AG closure.  The cleanest
strictly AG target remains (3.4): a class-uniform upper bound on an
**absolute, sheet-separated polarized facet degree**.  Such a bound would
simultaneously control facets and dimension and would let exact return plus
Mostow finish the original Coxeter problem without reflective-lattice
finiteness.

## 7. Primary controls

* J. Carlson and D. Toledo, [*Harmonic mappings of Kahler manifolds to
  locally symmetric spaces*](https://numdam.org/item/PMIHES_1989__69__173_0/),
  Publ. Math. IHES 69 (1989), 173--201.
* M. Kapovich and J. Kollar, [*Fundamental groups of links of isolated
  singularities*](https://arxiv.org/abs/1109.4047), especially the projective
  SNC realization of arbitrary finitely presented groups.
* P. Hacking, [*Compact moduli of hyperplane arrangements*](https://arxiv.org/abs/math/0310479),
  and P. Hacking--S. Keel--J. Tevelev,
  [*Compactification of the moduli space of hyperplane arrangements*](https://arxiv.org/abs/math/0501227).
* J. Bruinier, [*Regularized theta lifts for orthogonal groups over totally
  real fields*](https://arxiv.org/abs/0908.3076).
* H. Wang, [*The classification of 2-reflective modular forms*](https://arxiv.org/abs/1906.10459).
* G. Mostow, *Strong Rigidity of Locally Symmetric Spaces*, and G. Prasad,
  *Strong rigidity of Q-rank 1 lattices*, Invent. Math. 21 (1973), 255--286.
* V. V. Nikulin, [*On arithmetic groups generated by reflections in
  Lobachevsky spaces*](https://doi.org/10.1070/IM1981v016n03ABEH001321),
  Math. USSR-Izv. 16 (1981), 573--601, Appendix Theorem 1.
* V. V. Nikulin, *On the classification of arithmetic groups generated by
  reflections in Lobachevsky spaces*, Math. USSR-Izv. 18 (1982), 99--123,
  Theorem 5.2.1, pp. 119--122.
* A. G. Khovanskii, *Hyperplane sections of polyhedra, toroidal manifolds,
  and discrete groups in Lobachevsky space*, Funct. Anal. Appl. 20 (1986),
  41--50, for the edge-simple face-average inequality.
* M. N. Prokhorov, [*The absence of discrete reflection groups with
  noncompact fundamental polyhedron of finite volume in Lobachevsky space of
  large dimension*](https://doi.org/10.1070/IM1987v028n02ABEH000889), Math.
  USSR-Izv. 28 (1987), 401--411.
* V. Timorin, [*On polytopes simple in edges*](https://arxiv.org/abs/math/0010213),
  an independent polytope-algebra/Hard-Lefschetz proof of Khovanskii's
  face-average estimate.
