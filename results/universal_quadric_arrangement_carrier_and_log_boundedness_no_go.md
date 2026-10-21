# Universal quadric-arrangement carriers and the exact limit of log/toric boundedness

**Status (2026-08-12): unconditional exact carrier theorem and sharp
boundedness no-go for the standard log-Fano, KSBA, toric, wonderful-model,
and finite-inertia-stack routes.  This is not an AG proof of global
reflection-group finiteness.**

## 1. Outcome

There is a simple genuinely functorial algebraic-geometric carrier for every
marked arithmetic Coxeter polytope.  It is more canonical than the universal
blow-up/Chow carrier and it remembers the hyperbolic metric, not just a copy
of a lattice:

\[
 P\longmapsto
 \bigl(\mathbf P(V),\;Q_q,\;H_1,\ldots,H_f,\;\sigma,\;C_P\bigr).
 \tag{1.1}
\]

Here `q` is the native quadratic form, `Q_q=(q=0)` is its projective
quadric, `H_i=P(r_i^perp)` is the polar hyperplane of the `i`-th
facet root, `sigma` is the distinguished real place, and `C_P` is the
chosen component/side marking.  The marked quadric arrangement returns every
facet, every finite angle, every parallel or ultraparallel distance, and the
original chamber exactly.

This universal arrow does **not** yield finiteness from the standard AG
boundedness theorems.  Its natural reflection boundary

\[
                  B_{1/2}=\frac12\sum_{i=1}^f H_i
\]

satisfies the exact identity

\[
 K_{\mathbf P^n}+B_{1/2}
       =\left(\frac f2-(n+1)\right)\mathcal O_{\mathbf P^n}(1).
 \tag{1.2}
\]

The absolutely reflection-maximal Borcherds chamber in (H^{21}), with 210
facets, therefore gives

\[
                    K+B_{1/2}=83H,                         \tag{1.3}
\]

which is ample rather than anti-ample.  Thus absolute reflection maximality
does not imply that the universal arrangement carrier is log Fano.  In the
stable/general-type direction its natural log volume is

\[
                  \left(\frac f2-(n+1)\right)^n,           \tag{1.4}
\]

which is input complexity, not a bounded output.  KSBA assumes fixed
dimension and fixed log volume; it does not manufacture either bound.

Changing the coefficients does not repair this:

* forcing log Calabi--Yau uses coefficient `(n+1)/f`, so a uniform positive
  coefficient floor is exactly a bound on `f/(n+1)`;
* forcing a fixed stable volume still leaves the carrier dimension variable
  and makes the coefficients tend to zero if `f/(n+1)` is unbounded; and
* bounded finite-type moduli do not make the returned Gram or its field
  locally constant.  A fixed-volume stable-pair family with a moving marked
  cross-ratio has fields of unbounded degree.

The exact missing bridge is therefore not another compactification.  One
must prove from reflection maximality a uniform bound on carrier dimension
and log/Hilbert/stabilizer complexity, or prove that the exact Gram map is
zero-dimensional (equivalently, cohomologically discrete on a bounded
family).  None of BAB, KSBA, toric boundedness, wonderful compactification,
or arrangement Miyaoka--Yau supplies that bridge.

## 2. The exact universal arrangement functor

### 2.1 Input category

Let `k` be a totally real number field and let `sigma:k -> R` be a
distinguished embedding.  A **marked quadratic Coxeter
datum** consists of

* a nondegenerate `k`-quadratic space `(V,q)` of dimension `n+1`, with
  signature `(1,n)` at `sigma` and definite at every other real place;
* negative root lines `k r_1,...,k r_f` in `V`, together with their
  inward real rays at `sigma`; and
* the nonempty finite-volume chamber

  \[
  P=\{[x]\in\mathbb H_q:\ q(x,r_i)\geq0\text{ for all }i\}.
  \tag{2.1}
  \]

Arithmetic integrality is additional structure on a lattice in `V`, but
is not needed for the exact-return statement.  Morphisms are field-compatible
projective similitudes preserving the ordered root rays and the selected
positive component.  Passing from ordered to unmarked facets simply means
quotienting by the finite permutation action.

The inward-ray or component marking is indispensable.  Algebraically,
`H_i` does not distinguish `r_i` from `-r_i`, while those choices select
opposite half-spaces.  This is a finite real enhancement of an algebraic
object, not storage for any Gram coefficient.

### 2.2 Carrier category

Define a marked quadric-arrangement carrier over `(k,sigma)` to be

\[
 \mathcal A=(\mathbf P(V),Q,H_1,\ldots,H_f,C),             \tag{2.2}
\]

where `Q` in `P(V)` is a smooth quadric, the `H_i` are ordered
hyperplanes, and `C` is a chosen component of the appropriate real
hyperbolic locus cut out by them.  Morphisms are ambient projective
isomorphisms preserving all markings.

> **Theorem 2.1 (universal exact quadric-arrangement carrier).**
> For every marked quadratic Coxeter datum of dimension `n>=2`, the rule
>
> \[
> (V,q,(r_i),P)\longmapsto
> (\mathbf P(V),Q_q,(\mathbf P(r_i^\perp)),C_P)            \tag{2.3}
> \]
>
> is functorial and faithful.  The carrier recovers the same oriented
> hyperbolic chamber, including the complete normalized Gram matrix, up to
> the expected common quadratic similarity and positive rescaling of
> individual facet normals.

#### Proof

The homogeneous ideal of the nondegenerate projective quadric `Q_q` is
generated in degree two by `q`.  Hence `Q_q` determines `q` up to a
nonzero scalar.  Its polarity sends a hyperplane `H_i` to the root line
`k r_i`.  At the distinguished real place the quadric determines the
projective hyperbolic ball, and the real ray/component marking chooses the
correct side of each polar hyperplane.  Intersecting those sides gives
exactly (2.1).

Choose any representatives `r_i` of the marked rays.  The quantities

\[
       \frac{q(r_i,r_j)^2}{q(r_i,r_i)q(r_j,r_j)}           \tag{2.4}
\]

are invariant under the undetermined common scaling of `q` and under
individual rescaling of the `r_i`.  Together with the inward-ray signs,
they recover the signed unit-normal Gram.  Thus (2.4) distinguishes finite
angles, parallel facets, and ultraparallel distances.  An ambient projective
isomorphism preserving `Q_q` pulls its defining equation back to a scalar
multiple, so it is a projective similitude and preserves the recovered
chamber.  Conversely every marked projective similitude induces such a
carrier isomorphism.  This proves functoriality, faithfulness, and exact
return.  \(\square\)

### 2.3 What has and has not been algebraized

The algebraic part `(P(V),Q_q,(H_i))` contains the complete metric
arrangement.  The distinguished real place and component are required to
select the hyperbolic real form and one chamber.  Neither can be recovered
from the complexification alone: conjugate real places and different
components have the same complex variety.

This is already a useful answer to the first half of the programme:

```text
every marked arithmetic Coxeter chamber
  -> exact algebraic quadric/hyperplane arrangement
  -> the same marked chamber.
```

It is not yet the desired second half, because the dimensions, number of
marked hyperplanes, base fields, and arrangement coefficients all vary.

## 3. Exact log-class calculation

Let `H=O_{P(V)}(1)`.  Each arrangement hyperplane is
linearly equivalent to `H`, and `K_{P^n}=-(n+1)H`.  For an equal
coefficient `a` in `(0,1]`, put

\[
                         B_a=a\sum_{i=1}^fH_i.
\]

> **Proposition 3.1 (log-type trichotomy).**
> Numerically,
>
> \[
> K_{\mathbf P^n}+B_a=(af-(n+1))H.                        \tag{3.1}
> \]
>
> Consequently the pair is of log-Fano numerical type when `af<n+1`,
> log-Calabi--Yau numerical type when `af=n+1`, and log-general type when
> `af>n+1`.  In the last case,
>
> \[
>                    (K+B_a)^n=(af-(n+1))^n.              \tag{3.2}
> \]

The assertion is numerical: applying a log-Fano or KSBA theorem additionally
requires the relevant klt/slc hypotheses.  Those singularity hypotheses are
not automatic for an arbitrary facet arrangement and are another, rather
than a compensating, gate.

For the reflection-orbifold coefficient `a=1/2`, (3.1) gives (1.2).
The Borcherds chamber on `U + D_20(-1)` is an unconditional required
countercontrol: it is finite-volume, arithmetic, absolutely
reflection-maximal, has `n=21`, and has `f=210`.  Therefore

\[
 af-(n+1)=105-22=83.                                      \tag{3.3}
\]

More generally it fails the log-Fano numerical inequality for every fixed
`a>=11/105`.  In particular, absolute reflection maximality cannot be
translated into log-Fano positivity of the natural half-boundary.

If one instead forces `K+B=0`, the equal coefficient must be

\[
                             a_P=\frac{n+1}{f}.            \tag{3.4}
\]

A class-uniform lower bound `a_P>=delta>0` is precisely

\[
                             f\le\frac{n+1}{\delta}.       \tag{3.5}
\]

Thus the coefficient floor needed to count components in a bounded family
is not produced by the construction; it is the desired facet-to-dimension
bound in another notation.

## 4. Absolute reflection maximality does not supply BAB/KSBA inputs

It is tempting to interpret reflection maximality as stability or minimality
of the arrangement pair.  The notions point in different directions.
Reflection maximality forbids a larger **discrete real reflection group**
containing the side group.  Log-Fano positivity, klt discrepancies, KSBA
stability, Hilbert polynomials, and log volume are properties of the
projective arrangement and do not test the existence of such a real
reflection overgroup.

> **Proposition 4.1 (standard boundedness nonclosure).**
> Apply the exact functor of Theorem 2.1 to all absolutely
> reflection-maximal finite-volume arithmetic Coxeter chambers.  Standard
> BAB or KSBA boundedness cannot, by itself, prove that the image is finite.
> More precisely:
>
> 1. BAB fixes the dimension and a positive singularity gap before it
>    concludes boundedness.  The carrier dimension is the unknown `n`.
> 2. The natural half-boundary is not universally log Fano, by (3.3).
> 3. KSBA boundedness fixes the dimension, coefficient data, and log volume.
>    On the natural carrier the dimension is `n` and the volume is (3.2),
>    both varying input invariants.
> 4. Choosing coefficients to normalize the log class replaces the missing
>    facet bound by the missing uniform coefficient floor (3.5).
> 5. Even a bounded finite-type family does not give a finite set of Grams
>    unless the Gram-return map is locally constant on a finite
>    stratification or has zero-dimensional image.

#### Proof

Items 1 and 3 are hypotheses of the boundedness theorems, not conclusions
of them.  Items 2 and 4 follow from Proposition 3.1 and the `H^21`
countercontrol.  For item 5, finite type permits positive-dimensional moduli;
intersection-cohomological data are locally constant in a suitable bounded
cycle family, but coordinates, periods, and cross-ratios are not.  The
explicit stable-pair family in Section 5 proves that even fixed dimension,
boundary combinatorics, coefficients, and log volume do not finite-ize a
marked coefficient field.  \(\square\)

The proposition does not say that a new theorem can never relate maximality
to AG boundedness.  It identifies exactly what that new theorem must prove.
Any implication

\[
 \text{absolute reflection maximality}
 \Longrightarrow
 \begin{cases}
 n\le N,\\
 \text{uniform log/Hilbert or cycle complexity},\\
 \text{discrete exact Gram image}
 \end{cases}                                               \tag{4.1}
\]

is the missing group-to-geometry bridge.  It is not a formal consequence of
BAB or KSBA, and the first line is already a global dimension theorem of the
kind the programme is trying to derive.

There is a simple reason that variable-dimensional BAB cannot help even in
the absence of arrangements.  For every `n`, let `L_0,...,L_n` be
the coordinate hyperplanes in `P^n`.  Then

\[
 \left(\mathbf P^n,\frac12\sum_{i=0}^nL_i\right)
\]

is uniformly `epsilon`-klt for every fixed `epsilon<1/2`, and

\[
 -\left(K+\frac12\sum L_i\right)=\frac{n+1}{2}H
\]

is ample.  Uniformly mild log-Fano pairs therefore exist in every dimension.
The fixed-dimension hypothesis in BAB is essential.

## 5. Fixed-volume stable pairs do not finite-ize the Gram field

Let `X=P^1 x P^1`.  Denote vertical and horizontal
fibres by `V_t` and `H_t`.  For an algebraic
`lambda` not in `{0,1}`, take the ordered reduced boundary

\[
 \Delta_\lambda=
 V_0+V_1+V_\infty+V_\lambda+H_0+H_1+H_\infty.             \tag{5.1}
\]

It is SNC and

\[
 K_X+\Delta_\lambda=2V+H,
 \qquad (K_X+\Delta_\lambda)^2=4.                         \tag{5.2}
\]

Thus these are marked stable log surfaces with fixed dimension, coefficient
set, ordered boundary combinatorics, Hilbert type, and log volume.  The three
marked points `0,1,infinity` rigidify both projective-line factors, while the
fourth vertical point recovers `lambda`.

For odd primes `p`, put

\[
                  \lambda_p=\frac{1+\cos(\pi/p)}2.
\]

Then `[Q(lambda_p):Q]=(p-1)/2`.  Hence the fields of the
marked stable pairs have unbounded degree inside this fixed-volume bounded
family.  If the ordering is forgotten, the anharmonic orbit changes the
field degree by at most a fixed factor and the degrees remain unbounded.

This example is directly relevant to the quadric-arrangement carrier.  A
moving arrangement coefficient or cross-ratio can store a dotted Gram entry.
KSBA boundedness places that entry in a finite-type parameter space; it does
not make it take finitely many values.  To obtain Coxeter finiteness one needs
an additional theorem that the exact Gram image is zero-dimensional or is
generated by bounded integral cohomological actions and therefore locally
constant.

## 6. Toric, wonderful, and stack enhancements

### 6.1 Toric incidence loses the metric

The normal fan of a rationalized Klein polytope remembers its facet nerve.
It does not remember the Lorentz form, affine support constants, normalized
Gram, ground field, or dotted distances.  Rational perturbation is also a
choice rather than a functor on marked quadratic data.  Therefore the
automatic toric construction is an incidence carrier, not an exact carrier.

Adding the marked quadric `Q_q` and the original polar hyperplanes restores
exact return, but then the toric variety has not reduced the varying data: it
has merely been appended to the universal carrier of Theorem 2.1.  Its fan
boundedness controls the chosen rational polytope, while the metric remains
in the moving quadric/arrangement parameters.

The complete-intersection-surface escape has two further exact failures.  A
general surface section misses ambient strata of codimension at least three,
so it retains only pairwise facet incidence, and its polarization may be
replaced freely by `H -> kH`.  The logarithmic BMY expression for a
toric boundary section is a positive identity in the complete-intersection
classes; its leading term swamps any pairwise label correction under that
scaling.

### 6.2 Wonderful models preserve, but do not bound, arrangement complexity

The De Concini--Procesi construction functorially resolves a projective
hyperplane arrangement after a building set is chosen.  With the blowdown,
quadric, and labels retained, it still recovers the exact carrier.  Without
those markings its boundary incidence and standard cohomological invariants
are controlled largely by the intersection poset and forget the Lorentzian
metric.

The number of hyperplanes and arrangement flats controls the number of
blowups and boundary components.  Neither is uniformly bounded here.  Thus
the wonderful model moves the input complexity into its Hilbert/log-volume
data; it does not create a fixed bounded family containing all inputs.

### 6.3 Faithful inertia stores labels only at unbounded stack complexity

Root stacks along an SNC boundary have diagonalizable, hence abelian,
inertia.  At a transverse crossing, two order-two roots produce
\(\mu_2\times\mu_2\), which models only `I_2(2)`.  A faithful local
`I_2(m)` reflection quotient has coarse discriminant

\[
                         v^2-4u^m=0,                      \tag{6.1}
\]

which is irreducible for odd `m` and has two tangent branches for even
`m>2`; it is not the original transverse toric crossing.

Even if these quotient charts are globalized, a bounded finite-type family
of Deligne--Mumford stacks with quasi-finite inertia has uniformly bounded
geometric stabilizer length: the inertia morphism is quasi-finite of finite
type over a quasi-compact base, so its fibre lengths are bounded.  Faithfully
encoding `I_2(m)` would therefore make bounded stack complexity imply a
bound on `m`.  It does not prove that the carrier stacks form such a
bounded family.  That missing assertion is exactly the finite-label part of
the desired theorem.

## 7. Arrangement Miyaoka--Yau does not close the gap

The 2026 de Borbon--Panov theorem associates a quadratic form, determined by
the intersection poset, to a stable weighted hyperplane arrangement and
proves a Miyaoka--Yau-type inequality.  This is genuine AG information, but
it is intentionally poset-level information.  It cannot distinguish two
arrangements with the same intersection lattice and different marked
quadrics or cross-ratios, so it cannot by itself recover the dotted Gram.

Nor does the inequality give an absolute dimension bound.  For every
`n>=1`, take `N=n+2` hyperplanes in general position in
\(\mathbf P^n\), each with weight

\[
                            a=\frac{n+1}{n+2}.
\]

Every nonempty intersection of `r<=n` hyperplanes has codimension `r`
and total weight `ra<r`, while `Na=n+1`.  Hence the arrangement is klt
and Calabi--Yau in every dimension and falls within the theorem's basic
hypotheses.  The inequality constrains its intersection combinatorics but
does not exclude arbitrarily large `n`.

To turn this result into a Coxeter theorem one would need an additional
reflection-maximality condition on the complex arrangement that both
controls its Lorentzian metric parameters and forces a dimension/facet
bound.  No such translation is supplied by arrangement stability: matroid
stability is a condition on weighted flats, whereas absolute reflection
maximality forbids discrete real reflection overgroups.

## 8. Exact surviving positive theorem

The valid boundedness statement is conditional and sharp.

> **Theorem 8.1 (bounded exact arrangement image).**
> Let \(\mathcal C\) be a class of marked Coxeter chambers whose exact
> quadric-arrangement carriers satisfy all of the following uniform
> conditions:
>
> 1. their projective dimensions belong to a fixed finite set;
> 2. the number and projective degrees of all marked divisors and
>    correspondences are bounded;
> 3. they lie in a bounded log/stable family with the required coefficient,
>    singularity, and volume/Hilbert data;
> 4. the complete Gram return is generated by bounded integral
>    cohomological cycle actions under one fixed finite grammar, or its
>    scheme-theoretic image is zero-dimensional; and
> 5. the real place/component marking is only a finite selector.
>
> Then only finitely many normalized Grams, fields, and marked hyperbolic
> chambers occur in \(\mathcal C\).

#### Proof

The bounded family and bounded cycle degrees are parameterized by a finite
union of finite-type Hilbert/Chow spaces.  After finite stratification, the
integral cohomological actions and intersection matrices are locally
constant.  Under the first alternative in item 4, the fixed return grammar
therefore yields finitely many Grams.  Under the second alternative, a
zero-dimensional finite-type image has finitely many geometric points.
There are only finitely many real selectors on each returned datum.  A
spanning normalized Gram determines the marked hyperplane configuration up
to hyperbolic isometry, and its cyclic products determine the Vinberg field.
\(\square\)

The theorem is AG-native after its hypotheses are established.  For the
global maximal arithmetic reflection class, however, conditions 1--4 are
precisely the missing leverage bridge.  In particular:

```text
absolute reflection maximality
  -/-> log Fano for the universal exact carrier (H21/210 disproves it);
BAB does not bound the variable carrier dimension;
KSBA does not bound the variable log volume;
finite type does not discretize moving Gram coefficients;
bounded stack inertia does not prove the stacks are bounded.
```

## 9. Scope-safe conclusion

The requested first arrow exists unconditionally:

\[
 \mathsf{Cox}^{\mathrm{marked}}_{\mathrm{arith}}
   \longrightarrow
 \mathsf{QuadricArrangement}^{\mathrm{marked},\mathbf R}
\]

is functorial, faithful, and has exact same-chamber return.  The standard
log/toric/stable-pair upgrades do not provide the second arrow to global
finiteness.  The natural log-Fano proposal is directly contradicted by the
absolutely maximal (H^{21}/210) chamber; the KSBA and wonderful-model
proposals retain unbounded input complexity; and finite-type moduli alone
permit moving fields and Grams.

Therefore this lane isolates a precise research theorem rather than closing
it:

> Prove that absolute reflection maximality forces a uniformly bounded
> algebraic complexity package for the exact quadric arrangement, together
> with a zero-dimensional or cohomologically discrete Gram image.

Without that new bridge, invoking BAB or KSBA assumes the required bounds in
their hypotheses and is not an AG derivation of arithmetic reflection
finiteness.

## Primary sources and internal controls

Primary AG inputs:

* C. Birkar, [*Singularities of linear systems and boundedness of Fano
  varieties*](https://arxiv.org/abs/1609.05543).
* V. Alexeev, [*Boundedness and `K^2` for log
  surfaces*](https://arxiv.org/abs/alg-geom/9402007).
* C. De Concini and C. Procesi, [*Wonderful models of subspace
  arrangements*](https://doi.org/10.1007/BF01589496).
* M. de Borbon and D. Panov, [*A Miyaoka--Yau inequality for hyperplane
  arrangements in CP^n*](https://arxiv.org/abs/2411.09573),
  published in 2026.

Internal exact controls:

* [automatic toric complete-intersection audit](automatic_toric_ci_polytope_to_surface_audit.md);
* [log-Fano cohomological carrier theorem](log_fano_cohomological_coxeter_carrier_theorem.md)
  and its [hostile audit](log_fano_cohomological_coxeter_carrier_hostile_audit.md);
* [stable-pair/field obstruction](gram_belyi_carrier_dimension_field_bound_obstruction.md);
* [Borcherds H^21 outer-routing theorem](borcherds_h21_outer_routing_counterexample.md)
  and its [independent audit](borcherds_h21_outer_routing_counterexample_independent_hostile_audit.md);
* [global AG universality--leverage no-go](unconditional_ag_absw_universality_leverage_no_go.md).
