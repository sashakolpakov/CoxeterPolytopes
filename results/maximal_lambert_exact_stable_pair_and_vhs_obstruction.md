# Maximal Lambert chambers move in one fixed exact stable-pair family

**Status (2026-08-12): proved sharp obstruction and one surviving local
rigidity statement.  This is not an arithmetic counterexample and not a
proof of global arithmetic-reflection finiteness.**

## 1. Result

Absolute reflection maximality does not make the exact quadric-arrangement
carrier discrete, even after imposing substantially more algebraic geometry
than BAB or KSBA requires.

There is a one-real-parameter family of compact Coxeter quadrilaterals

\[
                         P_x\subset \mathbb H^2,\qquad x>1,
\tag{1.1}
\]

with angle sequence

\[
                  (\pi/2,\pi/2,\pi/2,\pi/3)
\tag{1.2}
\]

such that all but finitely many `P_x` are maximal among **all** discrete
hyperbolic reflection groups.  Nevertheless their normalized Grams vary
nontrivially with `x`.

For the exact universal carrier

\[
 (\mathbf P^2,Q_x,H_{1,x},\ldots,H_{4,x}),
\tag{1.3}
\]

put

\[
 \Delta_x=Q_x+\frac12\sum_{i=1}^4H_{i,x}.
\tag{1.4}
\]

Then every pair `(P^2,Delta_x)` is an SNC stable pair with the same
coefficient set, Hilbert data, boundary incidence, and

\[
 K_{\mathbf P^2}+\Delta_x=H,\qquad
 (K_{\mathbf P^2}+\Delta_x)^2=1.
\tag{1.5}
\]

It still recovers `P_x` exactly: the coefficient-one conic is the marked
quadric and the four coefficient-one-half lines are its polar facet
hyperplanes.  Hence the exact Gram map has a positive-dimensional image on
the absolutely maximal locus of this single fixed-volume KSBA family.

There is an even sharper K-stability collapse.  If

\[
                         B_x=\frac12\sum_iH_{i,x},
\tag{1.6}
\]

then `(P^2,B_x)` is uniformly K-stable for every `x`, and all the
hyperplane pairs are projectively isomorphic.  Thus their K-moduli point is
literally constant while the conic-relative Gram varies.

This is not an arithmetic counterexample.  It proves the more precise and
useful statement

```text
absolute reflection maximality
  + exact KSBA boundedness with fixed volume
  + uniform K-stability of the reflection boundary
  + fixed incidence, Chow classes, and ground field Q along an infinite sequence
  -/-> discrete Gram image.
```

Arithmetic integrality and the nongeometric-conjugate conditions must do
the missing work.  None of maximality, Hodge type, stable-pair boundedness,
K-stability, GIT stability, or finite automorphisms can be used as a synonym
for that arithmetic step.

## 2. The exact Lambert family

Use unit spacelike normals and the Lorentz-sign convention `(2,1)`.  For
`x>1` set

\[
 y=y(x)=\sqrt{\frac{4x^2-3}{4(x^2-1)}}>1
\tag{2.1}
\]

and consider

\[
G_x=
\begin{pmatrix}
1&0&-x&-1/2\\
0&1&0&-y\\
-x&0&1&0\\
-1/2&-y&0&1
\end{pmatrix}.
\tag{2.2}
\]

Direct expansion gives

\[
 \det G_x=\frac{4x^2y^2-4x^2-4y^2+3}{4}=0.
\tag{2.3}
\]

Its four principal three-by-three minors are

\[
 1-y^2,\quad \frac{3-4x^2}{4},\quad
 \frac{3-4y^2}{4},\quad 1-x^2,
\tag{2.4}
\]

so they are all negative.  The matrix has rank three and inertia
`(2,1,1_zero)`.  Adjacent entries are

\[
 0,0,0,-\cos(\pi/3),
\tag{2.5}
\]

while the two opposite entries are `-x` and `-y`, both strictly below
`-1`.  The standard Gram criterion therefore gives a compact hyperbolic
quadrilateral with (1.2), whose two opposite-facet distances have hyperbolic
cosines `x` and `y`.

For a direct realization, put `a=sqrt(x^2-1)` and use
`q=-t^2+u^2+v^2` with

\[
 e_1=(0,1,0),\quad e_2=(0,0,1),\quad
 e_3=(a,-x,0),\quad
 e_4=\left(\frac{x}{2a},-\frac12,-y\right).
\tag{2.5a}
\]

These vectors have Gram (2.2).  Choosing the four inequalities
`q(z,e_i)<=0`, the adjacent intersections are future timelike and satisfy
all the other inequalities strictly.  For example
`v_12=(1,0,0)` and `v_23=(x,-a,0)`, with
`q(v_23,e_1)=-a` and `q(v_23,e_4)=-1/(2a)`.  The other two vertices follow
cyclically.  This also verifies compactness without an abstract existence
appeal.

Conversely, the rank-three determinant relation for a quadrilateral with
these four fixed angles is exactly (2.3).  Thus (2.1) is its full
one-dimensional deformation space.  Notice that the involution exchanging
the two opposite pairs sends `x` to `y`, since for

\[
                  T(z)=1+\frac1{4(z-1)}
\tag{2.6}
\]

one has `y^2=T(x^2)` and `T(T(z))=z`.  Quotienting by all finite facet
permutations therefore still leaves a curve.

The hyperbolic area is constant:

\[
 \operatorname{area}(P_x)
 =2\pi-\left(\frac\pi2+\frac\pi2+\frac\pi2+\frac\pi3\right)
 =\frac\pi6.
\tag{2.7}
\]

Reflection in its sides gives a compact finite-covolume Coxeter group for
every real `x>1`.

## 3. Absolute maximality is Zariski-generic on this curve

> **Theorem 3.1 (cofinite absolute-maximal locus).**  Only finitely many
> isometry classes in the family `(P_x)_{x>1}` admit a proper discrete
> reflection overgroup.  Equivalently, outside a finite set modulo the
> involution `x<->y(x)`, the side group of `P_x` is absolutely
> reflection-maximal.

### Proof

Suppose `W(P_x)` is properly contained in a discrete reflection group
`W(F)`.  Since `W(P_x)` is a cocompact lattice, the index

\[
                          d=[W(F):W(P_x)]
\tag{3.1}
\]

is finite, `d>=2`, and a chamber `F` is compact.  Its reflected copies tile
`P_x`, and

\[
                     \operatorname{area}(F)=\frac\pi{6d}.
\tag{3.2}
\]

A compact Coxeter `r`-gon has

\[
 \operatorname{area}(F)=\pi\left(r-2-\sum_{j=1}^r\frac1{m_j}\right).
\tag{3.3}
\]

If `r>=5`, then `sum 1/m_j<=r/2`, so its area is at least `pi/2`.
For `r=4`, the all-right tuple is Euclidean; the largest reciprocal sum
strictly below two is

\[
                       \frac12+\frac12+\frac12+\frac13=\frac{11}{6}.
\tag{3.4}
\]

Thus every compact Coxeter quadrilateral has area at least `pi/6`.  Since
(3.2) is at most `pi/12`, `F` must be a triangle.

The smallest compact Coxeter-triangle area is `pi/42`, attained by
`(2,3,7)`.  This follows directly by ordering `2<=a<=b<=c` and maximizing
`1/a+1/b+1/c<1`: the extremal hyperbolic tuple is `(2,3,7)`.  Therefore

\[
                         \frac\pi{6d}\ge\frac\pi{42},
 \qquad 2\le d\le7.
\tag{3.5}
\]

For any fixed `d` in this finite interval, the triangle orders obey

\[
                  \frac1a+\frac1b+\frac1c=1-\frac1{6d}.
\tag{3.6}
\]

There are only finitely many ordered integer solutions: `a` is bounded by
three divided by the right side, then `b` is bounded after `a` is fixed,
and `c` is determined.

The complete elementary reciprocal census is

\[
\begin{array}{c|c}
d& (a,b,c)\\ \hline
2&(2,3,12),(2,4,6),(3,3,4)\\
3&(2,3,9)\\
4&(2,3,8)\\
5&\varnothing\\
6&\varnothing\\
7&(2,3,7).
\end{array}
\tag{3.7}
\]

Each resulting Coxeter triangle is rigid.  In its locally finite chamber
tessellation, there are only finitely many rooted connected unions of at
most seven chambers, up to its reflection group.  Consequently only
finitely many of those unions can be convex quadrilaterals, and they give
only finitely many isometry classes of `P_x`.  The unordered opposite
distance pair `{x,y(x)}` recovers `x` up to the involution (2.6), so each
isometry class has finite preimage in the parameter curve.  This proves the
claim.  \(\square\)

In particular, maximality is not a special or zero-dimensional condition
on this Coxeter moduli curve.  It is the complement of finitely many points.
Any proposal that treats reflection maximality as a Hodge locus, a GIT
polystable closed orbit, or a rigidity locus has the geometric direction
backwards in this basic family.

## 4. The exact carrier is one fixed-volume KSBA family

Let `V_x` be the rank-three quotient of the four-dimensional Gram space by
the radical of `G_x`, let `q_x` be its induced Lorentz form, and let `e_i`
denote the four root lines.  The exact carrier is

\[
 Q_x=(q_x=0)\subset\mathbf P(V_x),\qquad
 H_{i,x}=\mathbf P(e_i^\perp).
\tag{4.1}
\]

The conic is smooth.  No three of the four lines are concurrent, because
the four principal minors in (2.4) are nonzero.  A polar line `H_i` is
tangent to `Q_x` exactly when its pole `e_i` is isotropic, whereas
`q_x(e_i)=1`; hence every `H_i` meets `Q_x` transversely.  A line-line
intersection lies on `Q_x` exactly when the corresponding two-root Gram
block is degenerate.  The adjacent determinants are `1` or `3/4` and the
opposite determinants are `1-x^2` and `1-y^2`, all nonzero.  Thus there is
no conic-line-line triple point.

It follows that

\[
 \left(\mathbf P^2,
 Q_x+\frac12\sum_{i=1}^4H_{i,x}\right)
\tag{4.2}
\]

is SNC and log canonical.  Since `[Q_x]=2H` and `[H_i]=H`,

\[
 K+Q_x+\frac12\sum_iH_i
 =-3H+2H+2H=H.
\tag{4.3}
\]

The pair is therefore KSBA stable and has log volume one.  Its dimension,
coefficient set, Hilbert polynomial, number and degrees of boundary
components, dual incidence data, and all ambient intersection numbers are
fixed.  The coefficient-one component distinguishes the conic, so the
stable pair retains the full exact carrier.  By quadric polarity it returns
the varying Gram (2.2).

Theorem 3.1 now gives the promised sharp no-go:

> **Corollary 4.1.**  A single bounded fixed-volume family of exact stable
> pairs contains a positive-dimensional, Zariski-open real locus returning
> absolutely reflection-maximal compact Coxeter chambers.  Hence absolute
> reflection maximality plus all standard KSBA input data does not make the
> exact Gram image zero-dimensional.

This improves the earlier moving-cross-ratio control in two ways: every
point is itself an exact Coxeter carrier, and maximality holds off a finite
set.

## 5. Uniform K-stability and GIT forget the Gram entirely

Discard the conic temporarily and set

\[
                         B_x=\frac12\sum_{i=1}^4H_{i,x}.
\tag{5.1}
\]

The pair is klt log Fano and

\[
                      -(K_{\mathbf P^2}+B_x)=H.
\tag{5.2}
\]

Fujita's exact criterion for an SNC log-Fano hyperplane arrangement
`(P^n,sum d_i H_i)` says that uniform K-stability is equivalent to

\[
 k\sum_i d_i>(n+1)\sum_{j=1}^k d_{i_j}
\tag{5.3}
\]

for every `1<=k<=n` and every choice of `k` components.  Here `n=2`,
`d_i=1/2`, and `sum d_i=2`.  The two inequalities are

\[
                         2>\frac32,\qquad 4>3.
\tag{5.4}
\]

Thus every `(P^2,B_x)` is uniformly K-stable.  Moreover, four ordered lines
in general position are a projective frame in the dual plane, so all these
hyperplane pairs are projectively isomorphic.  Their K-moduli and GIT point
is constant.

The normalized Gram moves only in the relative position of `Q_x` and this
projective frame.  Consequently K-stability, GIT stability, log canonical
thresholds of the hyperplane arrangement, and finiteness of its
automorphism group cannot see the missing parameter.  Adding the conic
restores exact return, but Section 4 shows that it also restores a genuine
positive-dimensional stable-pair moduli curve.

## 6. Rational points of unbounded height and the exact arithmetic cut

Take an integer `N>=2` and put `x=N`.  Then

\[
 y_N^2=\frac{4N^2-3}{4(N^2-1)}\in\mathbf Q.
\tag{6.1}
\]

Although the unit Gram uses `y_N`, rescale the second root by `y_N`.
The resulting Gram has entries `y_N^2` on its second diagonal and
`-y_N^2` in position `(2,4)`, with every other entry rational.  Its radical,
rank-three quotient, quadric, and four root lines are therefore defined over
`Q`.  The normalized cyclic Gram field is also `Q`.

Only finitely many integers can lie in the exceptional set of Theorem 3.1.
Hence there are infinitely many absolutely reflection-maximal members of
this exact stable-pair family defined over `Q`, while the invariant

\[
                              4x^2=4N^2
\tag{6.2}
\]

has unbounded height.

They are deliberately **not arithmetic reflection groups**.  Vinberg's
integrality criterion requires every cyclic product of the unit Gram to be
an algebraic integer.  The length-two cycle for the opposite pair `(2,4)`
is

\[
 4y_N^2=\frac{4N^2-3}{N^2-1}
       =4+\frac1{N^2-1},
\tag{6.3}
\]

which is rational but not integral for `N>=2`.

This separates three notions that a height argument must not conflate:

```text
defined over Q + absolute reflection maximality + fixed AG boundedness
  does not bound height;

arithmeticity additionally imposes integral cyclic products and
definiteness at every nongeometric embedding.
```

The arithmetic condition has an instructive exact form on this curve.  Set

\[
                         A=4x^2,\qquad B=4y^2.
\tag{6.4}
\]

The rank equation becomes

\[
                         (A-4)(B-4)=4.
\tag{6.5}
\]

> **Proposition 6.1 (exact rational arithmetic locus).**  Among the members
> of this Lambert family whose normalized cyclic Gram field is `Q`, the
> arithmetic locus consists of the three labeled parameters
>
> \[
> (A,B)=(5,8),(6,6),(8,5),
> \tag{6.6}
> \]
>
> and hence two unlabelled isometry classes.

Indeed, Vinberg arithmeticity makes the length-two cycles `A` and `B`
algebraic integers.  Over `Q` they are integers, and the real Coxeter locus
has `A,B>4`.  Thus `A-4` and `B-4` are positive integer divisors whose
product is four, giving exactly (6.6).  Conversely the nonzero cycle graph
of (2.2) is a tree, so its cyclic products are generated by `A`, `B`, and
the finite-edge value one.  At the three points (6.6) all are integral;
there is no nonidentity real embedding of `Q`.  Vinberg's criterion therefore
makes these points arithmetic.  The first and third are exchanged by
`x<->y`, while the middle point is fixed.

This is a genuine tiny arithmetic-finiteness theorem inside one exact AG
moduli curve.  It comes from the integral-point equation (6.5), not from
KSBA or maximality.

Over varying totally real fields, the nongeometric real-place inequalities
in fact finite-ize this particular curve without an externally supplied
degree bound.  Put `u=A-4`, `v=B-4` and exchange the two opposite pairs so
that `0<u<=2` at the geometric place.  At every other embedding,
positive definiteness and `uv=4` force

\[
                         -4<\sigma(u)<-1.
\tag{6.7}
\]

For

\[
 p(T)=(T+1)(T+2)(T+3)(T+4),\qquad F(T)=p(T)(p(T)+1),
\tag{6.8}
\]

one has

\[
 \|F\|_{[-4,-1]}=225/256<1,
 \qquad \|F\|_{[0,2]}=129960.
\tag{6.9}
\]

The integer norm of `F(u)` then gives `[k:Q]<=92`; all conjugates of `u`
are bounded by four, so its possible monic integral minimal polynomials are
finite.  The disjoint sign intervals also prove `k=Q(u)`.  Therefore the
full all-field arithmetic locus of this Lambert curve is finite.  See the
[integral-point theorem and global height audit](vinberg_integral_points_on_quadric_arrangement_moduli.md)
for the proof and exact certificate.

Thus (6.5) is a positive miniature of the desired arithmetic geometry:
maximality leaves a moving `G_m`-curve, while Vinberg integrality plus a
strict conjugate-interval contraction cuts it to finitely many points.  The
global problem is to produce such a contracting, Gram-generating coordinate
packet uniformly for every maximal arithmetic chamber.

## 7. Natural integral VHS/local-system monodromy is impossible

There is a separate structural obstruction to turning the finite carrier
complement into an orthogonal local system.

> **Proposition 7.1 (commuting-meridian obstruction).**  Let
>
> \[
> U=\mathbf P(V)_\mathbf C\setminus
>       \left(Q\cup\bigcup_iH_i\right).
> \]
>
> If two carrier hyperplanes meet transversely away from the other boundary
> components, their small meridians commute in `pi_1(U)`.  Therefore no
> local system on `U` can send those meridians to two noncommuting side
> reflections.

### Proof

Near a generic point of `H_i cap H_j`, the boundary has equation `z_1z_2=0`.
The local complement is `(Delta^*)^2 times Delta^{n-2}`, whose two standard
meridians commute.  Monodromy preserves this relation.  Two distinct
orthogonal reflections commute only when their mirrors meet at angle
`pi/2`; in particular side reflections with finite Coxeter label `m>2`, or
with an ultraparallel relation, do not commute.  \(\square\)

In the Lambert family, the `(1,4)` pair has label three and already gives a
contradiction.  In the Borcherds `H^21/210` chamber, the 42 short-root odd
graph has label-three edges, so the same obstruction applies to the
mandatory high-rank control.

Root stacks do not repair this.  The inertia at an SNC crossing is abelian;
two order-two roots give `mu_2 times mu_2`.  Faithful `I_2(m)` monodromy
requires replacing the crossing by a non-SNC dihedral discriminant or a
nonabelian stack chart of order `2m`.  A bounded finite-type DM family with
finite inertia has bounded stabilizer length, so assuming such a family
would already assume a bound on `m`.

Consequently the natural complement of the finite exact carrier does not
carry the Coxeter representation as an integral VHS.  A new construction
would have to alter the local discriminant and then separately prove its
global boundedness.

## 8. Hodge loci and automorphisms have only the finite shadow

The ambient Hodge and Chow data of the exact carrier are constant:

\[
 \operatorname{CH}^*(\mathbf P^n)=\mathbf Z[H]/(H^{n+1}),
 \qquad [Q]=2H,\qquad [H_i]=H.
\tag{8.1}
\]

All of these classes remain algebraic of type `(p,p)` throughout the
arrangement moduli.  In the Lambert family, even the complete boundary
incidence and all log Chern/intersection numbers are constant while the
Gram moves.  Therefore a Hodge locus formed only from the ambient quadric
and divisor cycle classes is the whole family, not the arithmetic or
maximal locus.

This does not rule out manufacturing an auxiliary cover whose periods see
`x`.  It identifies the two new lemmas that such a proposal would require:

1. its period must recover the complete normalized Gram; and
2. arithmetic reflection maximality must force those periods into a
   bounded zero-dimensional special locus.

Neither follows from the divisor classes of the universal carrier.

Likewise, the side group is not the automorphism group of the finite
carrier.  A simple reflection sends most simple-facet hyperplanes to new
mirrors in the infinite `W`-orbit; it does not preserve the finite list
`{H_1,...,H_f}`.  The finite carrier automorphism group records diagram
symmetries, not reflection-overgroup maximality.  Adding the full mirror
orbit creates an infinite arrangement rather than a finite-type divisor.

The Lambert family makes the distinction exact: generic stable pairs have
finite automorphism group, and their reflection side groups are absolutely
maximal, but the moduli curve remains positive-dimensional.

For an arithmetic lattice, adding its full commensurator is worse: the
arithmetic commensurator is dense in the real isometry group.  It cannot act
properly discontinuously to produce the sought finite-type DM quotient.

## 9. What representation rigidity does salvage

There is one valid local discreteness statement in dimensions at least
three.  Mostow--Prasad rigidity says that two discrete faithful
finite-covolume representations of the same lattice in
`Isom(H^n)`, `n>=3`, are conjugate.  Thus, for a **fixed abstract group and
inside the complete lattice locus**, the normalized Gram has no continuous
deformation.  The `H^21/210` control is compatible with this statement.

This is not global boundedness:

* the complete-lattice locus is not a Zariski-defined finite-type component
  of the representation scheme (cusped representations can deform out of
  completeness);
* the theorem fixes the abstract group, hence its presentation, rank, and
  dimension; and
* it gives no finite list of Coxeter matrices, fields, or components.

Dimension two is precisely the exceptional deformation case, and Sections
2--6 show that absolute maximality does not remove it.  Hence local rigidity
can serve as the zero-dimensional return gate only **after** a genuinely
AG-native theorem has bounded the possible dimensions and presentation
types.

## 10. Consequence for the global programme

The proposed implication

\[
 \text{absolute arithmetic reflection maximality}
 \Longrightarrow
 \text{bounded/discrete exact AG image}
\tag{10.1}
\]

has not been proved.  This report determines which part cannot possibly be
formal:

```text
reflection maximality contributes no moduli discreteness in dimension 2;
it is Zariski-open dense on an exact fixed-volume stable-pair curve.

The finite carrier's natural Hodge classes, K/GIT stability, lct data,
automorphisms, and local-system complement do not encode the missing
arithmetic integrality.
```

The remaining viable statement must use the arithmetic condition as a
Diophantine constraint on the exact Gram map.  At minimum it must prove,
without ABSW/Nikulin finiteness conclusions,

1. a uniform bound on the degree/height or integral-model complexity of the
   arithmetic points of the arrangement moduli;
2. a dimension/presentation bound that is not a BAB hypothesis in disguise;
3. finite arithmetic points on each surviving exact Gram stratum; and
4. compatibility with the absolutely maximal `H^21/210` chamber, for which
   the natural half-boundary is of general type rather than log Fano.

Equation (6.5) exhibits the desired arithmetic geometry in its smallest
nontrivial form.  A successful global theorem needs a higher-dimensional
replacement for that integral `S`-unit cut, together with a uniform degree
bound.  Stable-pair boundedness and maximality alone provably cannot supply
either one.

## Primary and internal references

* K. Fujita, *K-stability of log Fano hyperplane arrangements*, J. Algebraic
  Geom. 30 (2021), 603--638, [arXiv:1709.08213](https://arxiv.org/abs/1709.08213),
  especially Theorem 1.5 and Corollary 1.6.
* A. Felikson and P. Tumarkin, *Reflection subgroups of Coxeter groups*,
  Trans. Amer. Math. Soc. 362 (2010), 847--858.  (The proof of Theorem 3.1
  above does not need their facet-rank theorem; elementary area bounds are
  enough.)
* E. B. Vinberg, *Discrete linear groups generated by reflections*, Math.
  USSR-Izv. 5 (1971), 1083--1119, especially the cyclic-product arithmetic
  criterion in Theorem 5.
* [Universal exact quadric-arrangement theorem and standard boundedness
  no-go](universal_quadric_arrangement_carrier_and_log_boundedness_no_go.md).
* [Borcherds `H^21/210` outer-routing and absolute-maximality
  control](borcherds_h21_outer_routing_counterexample.md).
* [Exact symbolic and reciprocal-sum
  certificate](../scripts/certify_maximal_lambert_stable_pair.py).
