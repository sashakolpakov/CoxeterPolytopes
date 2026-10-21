# Log-Fano cohomological carriers: an AG-native dimension and field theorem

**Status (2026-08-10): proved class theorem and explicit weak-del-Pezzo
corollary; essential surjectivity for all maximal arithmetic reflection
chambers remains open.**  This gives a precise route to an ABSW-type
finiteness theorem from `POLYTOPE -> SURFACE`, without importing arithmetic
volume finiteness.  The decisive geometric conditions are a uniform
log-discrepancy gap, anti-log-canonical positivity, separated facet boundary,
and a bounded **cohomological** coefficient-field marking.

The theorem also identifies two insufficient substitutes.  Reflection
maximality cannot merely be translated into indecomposability of a Belyi
cover, and bounded underlying surface data cannot merely retain the Gram in
moving cross-ratios.  Both conditions admit elementary unbounded families.

## 1. Carrier class

Let `P` be a marked finite-volume arithmetic Coxeter polytope of dimension
`n`, with `f` facets and Vinberg ground field `k` of degree
`d=[k:Q]`.  Fix positive real numbers `epsilon,delta`.

An **`(epsilon,delta)` log-Fano facet carrier** consists of:

1. a smooth projective surface `S` and an effective rational boundary `B`
   such that `(S,B)` is `epsilon`-log terminal and `-(K_S+B)` is nef and
   big;
2. a distinct reduced component `D_i` of `Supp(B)` for every facet of `P`,
   whose coefficient in `B` is at least `delta`;
3. the real place, orientation, side, and finite-volume markings needed to
   return the original chamber rather than a conjugate; and
4. an exact return rule recovering the complete unit-normal Gram of `P`.

An **intrinsic** carrier also contains an isometric injection of the entire
torsion-free restriction-of-scalars quadratic module underlying `P` into
`NS(S)/torsion`.  After the declared sign/scalar convention, its trace form
has the distinguished summand of signature `(1,n)` and all other real
summands negative definite.  Its rational rank is `d(n+1)`.  A list of only
`n+1` divisor roots or an external local system is not this datum.

For the field-finiteness statement, fix two more integers `t,E`.  A
**bounded cohomological carrier** has at most `t` additional marked effective
integral divisor cycles and effective integral algebraic correspondences on
`S` and `S x S`, each of degree at most `E` in one bounded projective
presentation.  A signed cycle is allowed only when its positive and negative
parts are separately marked and degree-bounded.  Across the whole class
there is one global finite return grammar; it is not chosen separately for
each `P`.  Any rational normalizing coefficients are part of that grammar,
not free parameters.  Coefficient algebras and projectors are generated
functorially by the bounded integral correspondence actions.  Moving
periods, repeated-eigenspace choices, cross-ratios, and real scalar markings
cannot store coefficients.  The exact Gram is recovered from:

```text
the intersection classes of the D_i,
the induced actions of the correspondences on those classes,
their finite coefficient algebra and distinguished real eigenspace.
```

In particular, the Gram is not stored in a freely moving point or
cross-ratio.  The degree bound is a degree of algebraic cycles on the
surface, not a prelisted number field or lattice.

## 2. The AG inequalities

There are first completely local design restrictions.  The valuation of a
facet component itself gives `1-b_i>epsilon`, so a nonempty carrier class
requires

\[
                           \delta<1-\epsilon.
\]

Next suppose `r` reduced
facet branches pass through one smooth point of `S`, each with multiplicity
one there.  Blowing up that point gives an exceptional divisor of log
discrepancy

\[
                         2-\sum_{i=1}^r b_i.
\]

Consequently `epsilon`-log terminality and `b_i>=delta` force

\[
                 \boxed{r\delta<2-\epsilon}.             \tag{2.0}
\]

In particular, if a finite vertex of `P` is represented literally by the
common point of all its incident facet curves, then at least `n` branches
meet there and

\[
                         n\delta<2-\epsilon.
\]

Thus a high-dimensional compiler cannot use naive point concurrence.  It
must use a controlled exceptional cluster or another decoration, and every
exceptional valuation in that cluster supplies a further explicit linear
log-discrepancy inequality.  More generally, for every divisorial valuation
`E` above `S`,

\[
 \sum_i b_i\operatorname{ord}_E(D_i)
   < A(E,S,0)-\epsilon,
 \qquad
 \delta\sum_i\operatorname{ord}_E(D_i)
   < A(E,S,0)-\epsilon.                                  \tag{2.0a}
\]

For a P-derived blowup/proximity graph, all orders and the smooth-surface log
discrepancy on the right are integers computed recursively.  Equation
(2.0a) is therefore a finite rational feasibility system before the surface
moduli search begins.  This is the first direct translation from Coxeter
face combinatorics to surface restrictions.

There is also a dimensional obstruction before any discrepancy calculation.
Because an `epsilon`-klt pair has no coefficient-one boundary, its standard
dlt dual complex is empty.  If one instead remembers the incidence complex
of the full SNC support, an `s`-dimensional carrier has support complex of
dimension at most `s-1`.  Literal recovery of the boundary complex of an
`n`-polytope by boundary strata would therefore require `n<=s`.  In
particular, a surface cannot intrinsically recover the full boundary sphere
of a polytope in dimension `n>=3` by ordinary divisor intersections.  A
higher-dimensional Coxeter carrier on a surface must keep the missing face
poset in a bounded exceptional cluster, cohomological correspondence, or
other declared decoration.  It cannot silently call that external datum the
surface dual complex.

For the natural reflection-orbifold coefficient `b_i=1/2`, this becomes
especially sharp.  A literal finite vertex on a smooth surface has log
discrepancy

\[
                              2-{n\over2}.
\]

Hence such a pair is klt only for `n<=3`; at `n=4` it is log canonical but
not klt, and for `n>4` it is not log canonical.  This independently recovers
the dimension-three boundary of the direct full-nef surface theorem for the
naive mirror-incidence model.  It also explains why the compact `H8`
eigen-K3 carrier and every prospective higher-dimensional surface carrier
must retain higher faces through ambient incidence, exceptional clusters, or
another decoration rather than ordinary concurrence of all facet curves.

> **Theorem 2.1 (log-Fano facet and trace bounds).**  For fixed `epsilon`
> there are constants `C_epsilon` and `R_epsilon`, depending only on the
> two-dimensional log-Fano bounded family, such that every
> `(epsilon,delta)` carrier satisfies
>
> \[
>       f\le {C_\epsilon\over\delta},\qquad
>       n+1\le f\le {C_\epsilon\over\delta}.             \tag{2.1}
> \]
>
> If the carrier is intrinsic, then also
>
> \[
>       \boxed{d(n+1)\le\rho(S)\le R_\epsilon}.           \tag{2.2}
> \]

**Proof.**  Alexeev's boundedness theorem for log surfaces says that the
underlying surfaces admitting an `epsilon`-log-terminal boundary with
`-(K+B)` nef are bounded, apart from the numerically Calabi--Yau case.  The
present nef-and-big hypothesis excludes that exception.  Choose a finite
type bounding family and one relative integral very ample divisor `H` after
finite stratification.  There are uniform bounds

\[
                -K_S\mathbin{\cdot}H\le C_\epsilon,
       \qquad \rho(S)\le b_2(S)\le R_\epsilon.            \tag{2.3}
\]

Intersecting `-(K_S+B)` with `H` gives

\[
  0\le -(K_S+B)\mathbin{\cdot}H,
  \qquad B\mathbin{\cdot}H\le -K_S\mathbin{\cdot}H.
\]

Every distinct facet component has positive integral `H`-degree and
coefficient at least `delta`.  Hence

\[
 \delta f\le
 \delta\sum_i H\mathbin{\cdot}D_i
 \le B\mathbin{\cdot}H\le C_\epsilon.
\]

An `n`-dimensional full-dimensional polytope has at least `n+1` facets,
proving (2.1).
For an intrinsic carrier, the injected trace module has rational rank
`d(n+1)`, which cannot exceed `rho(S)`.  This proves (2.2).  `square`

One part of the geometry is already numerically effective.  If
`N=floor(2/epsilon)`, Jiang's optimal surface bound gives

\[
  0<\bigl(-(K_S+B)\bigr)^2
  \le \max\left\{9,N+4+{4\over N}\right\}.               \tag{2.3a}
\]

This volume bound does **not** by itself bound the facet count: a nef
anti-log-canonical class may have degree zero on boundary curves.  That is
why (2.1) uses an integral very ample `H` from the bounded family.  An
explicit closed formula for `C_epsilon` or `R_epsilon` requires effective
very-ampleness/Betti bounds in addition to (2.3a); none is fabricated here.

There is an immediate finite-angle consequence.  A finite Coxeter label
`m>=3` contributes the cyclic Gram product

\[
  4\cos^2(\pi/m)\in k,
  \qquad [\mathbb Q(4\cos^2(\pi/m)):\mathbb Q]=\varphi(m)/2.
\]

Therefore every intrinsic carrier satisfies

\[
       \varphi(m)\le 2d\le {2R_\epsilon\over n+1}.        \tag{2.4}
\]

Thus the same surface geometry simultaneously bounds dimension, ground-field
degree, and the finite Coxeter-label alphabet.  It does not yet bound dotted
entries; that is the purpose of the cohomological finiteness theorem below.

## 3. Positive-dimensional surfaces can still give a finite Gram image

> **Theorem 3.1 (bounded cohomological image).**  Fix
> `epsilon,delta,t,E`.  Arithmetic Coxeter polytopes admitting bounded
> cohomological `(epsilon,delta)` carriers have only finitely many returned
> unit-normal Grams and Vinberg fields.  Consequently they form a finite set
> up to hyperbolic isometry.  The carrier surfaces themselves may move in
> positive-dimensional families.

**Proof.**  The proof of Theorem 2.1 puts the surfaces in a bounded
projective family and bounds the number and projective degrees of all facet
cycles.  Ordered facet cycles are therefore parameterized by a finite union
of relative Chow spaces of bounded degree.  The at-most-`t` correspondence
cycles of degree at most `E` are parameterized by another finite union of
relative Chow spaces on the fibre square.  These spaces are of finite type.

After a finite stratification, the surface family and every effective cycle
(including both separately bounded parts of a signed cycle) are flat.  On
each connected stratum, intersection numbers of the marked cycles and the
integral cohomological actions induced by the correspondences are locally
constant.  Hence their characteristic polynomials and the integral
coefficient algebras generated by those actions range over a finite set.
The single global return grammar then permits only its fixed finite rational
normalization alphabet and functorially determined projectors.  It does not
choose a new projector, repeated-eigenspace vector, period, or scalar for an
individual polytope.  Thus all pairings of the returned facet classes, and
therefore all complete unit-normal Grams, range over a finite set.  The real
marking only selects among finitely many embeddings.  Cyclic Gram products
compute the Vinberg field, and a spanning Gram determines its marked
hyperplane configuration up to ambient isometry.  Exact irredundancy and
literal finite-volume tests select a subset of this finite list.  `square`

This argument needs neither rigidity of the surface nor zero-dimensionality
of its moduli space.  It proves finiteness of the **cohomological Gram
image**, not of carrier fields or carrier isomorphism classes.  It becomes
an effective enumeration once an explicit bounded projective atlas and the
constants for the cycle Chow spaces are supplied.  Alexeev boundedness alone
gives existential constants, so no fake numerical value for
`C_epsilon,R_epsilon` is asserted here.

## 4. Explicit weak-del-Pezzo endpoint

The intrinsic rank bound is completely explicit when the carrier surface is
a smooth weak del Pezzo surface.

> **Corollary 4.1.**  If `-K_S` is nef and big and the trace module of an
> arithmetic `H^n` chamber over `k` embeds in `NS(S)`, then
>
> \[
>                    \boxed{d(n+1)\le9}.                 \tag{4.1}
> \]
>
> In particular, a rational-field carrier has `n<=8`, while a real-quadratic
> carrier has `n<=3`.

Indeed a smooth weak del Pezzo surface is rational, `K_S^2` is a positive
integer, and Noether's formula gives

\[
                    \rho(S)=10-K_S^2\le9.
\]

This is a genuine explicit `P -> surface -> dimension bound`.  It is a
carrier theorem, not a claim that every arithmetic Coxeter polytope has a
weak-del-Pezzo carrier.

## 5. Conditional ABSW replacement

The class theorem has the following precise completion target.

> **Corollary 5.1.**  Suppose every reflection-maximal arithmetic
> finite-volume Coxeter chamber has an exact bounded cohomological
> `(epsilon,delta)` log-Fano facet carrier for one uniform tuple
> `(epsilon,delta,t,E)`.  Then there are only finitely many such maximal
> chambers, in all dimensions and over all ground fields.  Their dimensions,
> field degrees, finite labels, exact Grams, and fields are obtained from the
> bounded surface/cycle atlas.  No ABSW or Prasad volume theorem is used.

The unproved premise is now geometric and falsifiable.  A P-first compiler
must construct the facet boundary and coefficient correspondences, prove a
uniform discrepancy gap and anti-log-canonical nefness/bigness, and prove
exact return.  Failure of log canonicity is witnessed by a valuation or
boundary cluster; failure of nefness is witnessed by a curve.  Either is a
finite algebraic obstruction that can be translated back into a Coxeter
cut.  Surviving bounded cohomological packets are a finite positive search.

Reflection maximality should therefore be translated, if possible, into the
uniform log-Fano and correspondence-degree assertions.  Merely declaring a
carrier “minimal” or “primitive” is insufficient.

## 6. Two sharp no-go checks

### Primitivity does not bound dessin degree

For every prime `p`,

\[
                       x\longmapsto x^p
\]

is a degree-`p` Belyi map on `P1`, branched only over `0` and `infinity`.
Its monodromy is the regular action of the cyclic group of prime order, hence
primitive; equivalently the prime-degree map is indecomposable.  Its source
curve and all ordinary birational invariants are constant while its degree
is unbounded.  Thus

```text
reflection maximality -> primitive Belyi carrier
```

could not by itself bound the Gram--Belyi or Schreier degree.

### Bounded surfaces do not bound moving residue fields

The four-pointed curves

\[
                  (\mathbf P^1;0,1,\infty,\lambda)
\]

have fixed Hilbert data while algebraic `lambda` ranges over fields of
unbounded degree.  Hence log boundedness alone does not give a field list
when the Gram is stored in cross-ratios.  Theorem 3.1 avoids this defect by
requiring a bounded cohomological return.  Alternatively one must prove the
scheme-theoretic Gram image zero-dimensional, as in the separate
dessin-native field theorem.

The universal Gram--Belyi carrier lies on the cross-ratio side of this
boundary, and its frozen `H10` redoubling family has unbounded absolute
Belyi degree.  K3 carriers lie outside the log-Fano hypothesis because
`K_S` is numerically trivial.  Neither family contradicts the theorem.

## 7. Concrete completion programme

1. Starting only from the arithmetic Gram/lattice of `P`, build every
   candidate log-Fano surface, separated facet boundary, and coefficient
   correspondence without prelisting a field or known surface.
2. Compute log discrepancies and the cone of curves exactly.  Prove uniform
   `epsilon,delta`, and `-(K+B)` nef and big, or return the valuation/curve
   as a P-side obstruction.
3. Prove that the trace module embeds and that bounded-degree
   correspondences recover the same unit-normal Gram at the distinguished
   real place.
4. Use (2.1)--(2.4) as the theorem boundary.  Enumerate the finite Chow
   packet of Theorem 3.1; exact Vinberg tests turn its survivors into
   examples and emptiness into nonexistence.
5. In parallel, test whether reflection overgroups correspond to explicit
   log-MMP operations on the carrier.  This is the correct place to translate
   maximality; it is not assumed in Theorems 2.1 or 3.1.

## 8. Other algebraic objects: fixed-dimensional Fano carriers

The boundedness argument is not peculiar to surfaces.  Fix a carrier
dimension `s`.  Let `X` be a smooth Fano `s`-fold in a bounded Fano class,
let `B` contain one distinct facet divisor of coefficient at least `delta`,
assume `(X,B)` is `epsilon`-log terminal and `-(K_X+B)` is nef, and retain a
bounded integral cohomological correspondence grammar with exact Gram
return.  Birkar's Fano boundedness gives a bounded projective family.  For a
uniform integral very ample `H`,

\[
 \delta f
 \le B\mathbin{\cdot}H^{s-1}
 \le -K_X\mathbin{\cdot}H^{s-1}
 \le C_{s,\epsilon}.                                    \tag{8.1}
\]

Consequently
`n+1<=f<=C_{s,epsilon}/delta`.  If the trace module embeds in
`N^1(X)` with the Hodge-index pairing

\[
                (D_1,D_2)_H=D_1D_2H^{s-2},
\]

then

\[
             [k:\mathbb Q](n+1)\le\rho(X)\le R_{s,\epsilon}. \tag{8.2}
\]

Relative Chow spaces for the bounded divisor/correspondence cycles prove the
same finite cohomological-Gram and Vinberg-field conclusion as Theorem 3.1.
Thus higher-dimensional Fano or Fano-type carriers are a legitimate parallel
route.  They do not evade the hard gate: `epsilon,delta` and correspondence
degree must be uniform, and exact essential surjectivity from raw maximal
Coxeter data remains to be proved.

The surface case remains special because its discrepancy inequalities are
read directly from a blowup/proximity graph and the weak-del-Pezzo endpoint
has the explicit constant nine.

Primary AG inputs:

- V. Alexeev,
  [*Boundedness and `K^2` for log surfaces*](https://arxiv.org/abs/alg-geom/9402007);
- C. Birkar,
  [*Singularities of linear systems and boundedness of Fano varieties*](https://arxiv.org/abs/1609.05543).
- C. Jiang,
  [*Bounding the volumes of singular weak log del Pezzo surfaces*](https://arxiv.org/abs/1305.6435).

The reflection/nef-cone ancestry of this bridge is also explicit in
V. Nikulin,
[*Algebraic surfaces with log-terminal singularities and nef anticanonical
class and reflection groups in Lobachevsky spaces*](https://arxiv.org/abs/alg-geom/9405011).
