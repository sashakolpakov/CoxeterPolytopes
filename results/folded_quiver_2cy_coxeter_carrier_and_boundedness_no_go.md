# Folded quiver and 2-Calabi--Yau Coxeter carriers: exact return and the boundedness barrier

**Status (2026-08-12).**  This route gives a genuine exact categorical
carrier for a nontrivial class, including noncrystallographic finite labels
such as $m=5$.  It does **not** prove the requested global finiteness theorem.
The direct Kac--Moody/preprojective construction is total only on integral
crystallographic wall data; the fusion-quiver repair handles arbitrary finite
Coxeter labels but not general ultraparallel distances without adding the
distance scalars as new decorations.  In every version presently available,
the quiver, unfolding, spherical-functor source, or marked kernel complexity
is unbounded.  Smoothness, properness, the Calabi--Yau condition, and
reflection maximality do not supply the missing uniform bound.

No ABSW, Nikulin, Vinberg dimension bound, field-degree bound, or reflective
lattice finiteness theorem is used below.  The named `X3`, `H21`, and
$(2,5,6)$ objects are controls, not inputs to a finiteness proof.

The strongest conclusions are:

1. a folded 2-Calabi--Yau construction exactly realizes every marked
   Coxeter datum whose normalized Gram is determined by finite labels and
   parallel $m=\infty$ pairs; in particular it covers hyperbolic Coxeter
   simplices of the usual finite-volume type;
2. an individual 2-spherical twist can only produce finite rank-two labels
   $2$ and $3$, while symmetrizable Kac--Moody theory only adds $4$ and
   $6$; the absolutely maximal arithmetic $(2,5,6)$ triangle therefore
   forces a genuine folding/fusion or higher-rank spherical functor;
3. restriction of scalars makes the cost precise: a native reflection over a
   degree-$d$ field has a $d$-dimensional $(-1)$-eigenspace, so it cannot
   be one spherical twist on an ordinary integral numerical group;
4. the direct facet quivers for `H21` and `X3` have respectively
   $210$ vertices/corank $188$ and $36$ vertices/corank $16$.  They
   are exact, but they store the wall list rather than compress it; and
5. a quasi-compact finite-type stack of **fully decorated** exact carriers
   would imply finiteness, but no theorem puts all maximal arithmetic
   reflection chambers in one such stack.

## 1. The categorical acceptance gates

Let $P\subset\mathbb H^n$ have ordered facets with outward roots
$r_1,\ldots,r_f$.  Write

\[
 b_{ii}=2,\qquad
 b_{ij}=-\frac{2(r_i,r_j)}{\sqrt{(r_i,r_i)(r_j,r_j)}}
 \quad(i\ne j),                                      \tag{1.1}
\]

with the harmless simultaneous sign convention chosen so that off-diagonal
simple-root entries are nonpositive.  Thus

\[
 b_{ij}=-2\cos(\pi/m_{ij})
\]

for intersecting walls, $b_{ij}=-2$ for a parallel pair, and
$|b_{ij}|>2$ for an ultraparallel pair.  The normalized Gram $B_P$, its
Lorentzian quotient, and a positive-component marking recover the marked
polytope.

A categorical carrier is useful for the requested theorem only if it passes
all five gates:

* `FORM`: a categorical numerical pairing recovers the complete $B_P$,
  including every dotted distance rather than only the abstract Coxeter
  exponent $m_{ij}=\infty$;
* `LIFT`: every facet reflection is an actual reflection functor, spherical
  twist product, spherical-functor twist, or Fourier--Mukai action;
* `RETURN`: the selected categorical/stability chamber is exactly $P$, not
  a supercone cut by hidden categorical walls;
* `BOUND`: the category **and** the marked objects, foldings, functors,
  kernels, stability data, and exact-return certificate lie in one
  quasi-compact finite-type family; and
* `CONSTR`: extraction of the discrete marked Gram is constant on a finite
  constructible stratification.

The first three gates give an encoding.  Only the last two give an
AG-native finiteness argument.

## 2. The rank-two spherical obstruction

Let \(\mathcal D\) be a Hom-finite 2-Calabi--Yau category and let $E,F$
be 2-spherical objects.  On

\[
 \Lambda=K_{\rm num}(\mathcal D)
\]

the Euler form is symmetric,
$\chi(E,E)=\chi(F,F)=2$, and the spherical twist in $E$ acts by

\[
 s_E(x)=x-\chi(E,x)[E].                               \tag{2.1}
\]

Put $a=\chi(E,F)\in\mathbf Z$.  In the span of their two numerical
classes,

\[
 \operatorname {tr}(s_Es_F)=a^2-2,
 \qquad \det(s_Es_F)=1.                               \tag{2.2}
\]

### Proposition 2.1 (individual spherical-twist angle obstruction)

If $s_Es_F$ has finite order, that order is $2$ or $3$.  More
precisely,

\[
 a=0\Longrightarrow m=2,\qquad |a|=1\Longrightarrow m=3,
\]

whereas $|a|\ge2$ gives infinite order.

#### Proof

For finite order the two eigenvalues of (2.2) lie on the unit circle, so
$|a^2-2|<2$, apart from the scalar endpoints.  Since $a\in\mathbf Z$,
the only possibilities are $a=0,\pm1$.  Direct substitution gives orders
two and three.  At $|a|=2$ the product is a nontrivial unipotent, and for
$|a|>2$ it is hyperbolic.  \(\square\)

Thus the slogan “facets are spherical objects” is already false for a
single $m=4,5,$ or $6$ edge.  `H21` contains $m=4$ pairs and `X3`
contains $m=6$ pairs, so neither control is an individual-spherical-object
carrier.

For a symmetrizable generalized Cartan matrix, the rank-two block is

\[
 A_{ij}=\begin{pmatrix}2&-p\\-q&2\end{pmatrix},
 \qquad p,q\in\mathbf Z_{\ge0}.                       \tag{2.3}
\]

The same calculation gives

\[
 \operatorname {tr}(s_is_j)=pq-2.                    \tag{2.4}
\]

Hence

\[
\begin{array}{c|ccccc}
 pq&0&1&2&3&\ge4\\ \hline
 m_{ij}&2&3&4&6&\infty.
\end{array}                                           \tag{2.5}
\]

This is the exact crystallographic barrier.  For $m=5$,

\[
 pq=4\cos^2(\pi/5)=\frac{3+\sqrt5}{2},               \tag{2.6}
\]

which is not an integer.  Ordinary symmetrizable Kac--Moody,
preprojective, and Nakajima-quiver constructions therefore cannot take an
$m=5$ facet pair to two simple Kac--Moody reflections.

The compact arithmetic $(2,5,6)$ triangle in
[`k3_compact_maximal_256_phys_empty.md`](k3_compact_maximal_256_phys_empty.md)
is an especially useful calibration.  It is absolutely reflection-maximal,
has ground field \(\mathbf Q(\sqrt5)\), and contains the forbidden edge.
Its maximality is used only to refute “maximality forces crystallographic
simple objects”; it is not used to infer any global finiteness.

## 3. Fusion quivers repair the $m=5$ defect

Heng's Coxeter-quiver construction attaches to an acyclic quiver whose
arrows have labels $m\ge3$ a fusion category

\[
 \mathcal C(Q)=\boxtimes_{m\in\operatorname {im}(v)}\mathcal C_m, \tag{3.1}
\]

where \(\mathcal C_m\) is a Temperley--Lieb--Jones fusion category (or its
even part), and an arrow labeled $m$ is weighted by an object of
Frobenius--Perron dimension $2\cos(\pi/m)$.  On

\[
 L_Q=\bigoplus_{i\in Q_0}K_0(\mathcal C(Q))e_i        \tag{3.2}
\]

there is a fusion-ring-valued symmetric form

\[
 \mathsf B(e_i,e_i)=2[\mathbf1],\qquad
 \mathsf B(e_i,e_j)=-\sum_{\alpha:i-j}[X_\alpha].     \tag{3.3}
\]

Applying Frobenius--Perron dimension to (3.3) gives exactly the real Tits
form with off-diagonal entry $-2\cos(\pi/m)$.  Two ordinary unlabeled
arrows give the limiting parallel weight $-2$.

The crucial point is that this is not merely a real-number trick.  Heng
constructs a finite ordinary **unfolded quiver** \(\check Q\) with vertices

\[
 \check Q_0=\operatorname {Irr}(\mathcal C(Q))\times Q_0              \tag{3.4}
\]

and proves an abelian equivalence

\[
 \operatorname {rep}_{\mathcal C(Q)}(Q)
       \simeq \operatorname {rep}_{\mathbf C}(\check Q).             \tag{3.5}
\]

The reflection functors act on (3.2) by the Coxeter reflections of
(3.3).  For $I_2(5)$, this is the familiar folding into type $A_4$;
the Fibonacci fusion category supplies the golden coefficient.

### Theorem 3.1 (exact folded 2-CY carrier)

Let $P$ be a marked hyperbolic Coxeter datum for which every off-diagonal
normalized Gram entry is determined by a finite Coxeter label, or is the
parallel value $-2$.  Then $P$ has a finite decorated 2-Calabi--Yau
quiver carrier with exact normalized-Gram return.

This applies in particular to the usual compact or finite-volume hyperbolic
Coxeter simplices (Lannér/quasi-Lannér type), where no independent
ultraparallel distance parameter occurs.

#### Construction and proof

Order the facets and orient every non-right edge from the lower to the
higher index.  This gives an acyclic Coxeter quiver $Q_P$; use two
unlabeled arrows for a parallel pair.  Form (3.1)--(3.5).  Let
\(\mathcal D_2(\check Q_P)\) be the standard 2-Calabi--Yau quiver category.
Its standard simple objects

\[
 S_{A,i},\qquad (A,i)\in
 \operatorname {Irr}(\mathcal C(Q_P))\times(Q_P)_0,
\]

are 2-spherical, their Euler form is the simply-laced unfolded root form,
and their twists give the standard braid/Weyl action on $K_0$.

No two vertices in a fibre

\[
 I_i=\{(A,i):A\in\operatorname {Irr}(\mathcal C(Q_P))\}              \tag{3.6}
\]

are adjacent.  The corresponding twists commute, so

\[
 T_i=\prod_{(A,i)\in I_i}T_{S_{A,i}}                                \tag{3.7}
\]

is defined independently of order.  The unfolding/folding identity says
that the **numerical action** of (3.7) on the folded module (3.2) is precisely
the reflection \(\sigma_i\) in $e_i$.  Concretely, after expanding every
fusion-ring coefficient in the simple-object basis, the product over the
fibre fixes all $v_j$ with $j\ne i$ and sends

\[
 v_i\longmapsto -v_i+
   \sum_j\sum_{\alpha:i-j}[X_\alpha]v_j
   =v_i-\mathsf B(e_i,v),                              \tag{3.7a}
\]

which is exactly the formula in Definition (3.3).  Consequently the numerical actions
satisfy the original Coxeter relations, including noncrystallographic
relations such as $m=5$.  As usual for spherical twists, the categorical
autoequivalences are an Artin/braid lift of the involutive numerical action;
one should not assert $T_i^2=1$ as an autoequivalence.

Finally, (3.3) followed by Frobenius--Perron dimension is $B_P$.  Quotient
by its radical, select the marked Lorentzian component, and impose the
inequalities against the $e_i$.  This returns the original marked chamber.
Thus, on this input class,

\[
 P\longmapsto
 (\mathcal C(Q_P),\check Q_P,\mathcal D_2(\check Q_P),
   \{I_i,T_i\},B_P^+)
 \longmapsto P                                             \tag{3.8}
\]

is the identity.  The facet ordering makes the first arrow functorial for
order-preserving maps; alternatively one retains the finite orientation
groupoid.  \(\square\)

The theorem is a real categorical improvement over an ordinary generalized
Cartan matrix: it handles $m=5$, $H_3$, $H_4$, and $I_2(m)$ without
replacing the golden coefficient by a rational approximation.

It is not yet an AG finiteness theorem.  Its cost is visible in (3.4): the
number of unfolded vertices is

\[
 |Q_0|\,\operatorname {rk}_{\mathbf Z}K_0(\mathcal C(Q)).            \tag{3.9}
\]

For a label $5$, the coefficient category is Fibonacci and already has
Grothendieck rank two.  More generally TLJ at level $m$ has $m-1$
simple objects before passage to the even part, and Deligne products multiply
ranks.  The field information has been categorified, not erased.

## 4. General polytopes: the dotted-distance obstruction

The abstract Coxeter matrix records only $m_{ij}=\infty$ when two facet
hyperplanes do not meet.  It does not record

\[
 |b_{ij}|=2\cosh d(H_i,H_j)>2.                         \tag{4.1}
\]

Heng's construction formally treats the limit $b_{ij}=-2$ by a Kronecker
pair, but (4.1) is additional continuous/algebraic data.  Therefore the
unadorned fusion quiver does not pass `FORM` for a general finite-volume
polytope.

There are two possible repairs, neither bounded.

1. In the rational crystallographic case, choose the integral Cartan
   coefficients
   \(a_{ij}=2(r_i,r_j)/(r_i,r_i)\).  A valued quiver or species remembers
   the integer product $a_{ij}a_{ji}\ge4$, hence the dotted distance.
2. Over a general ground field, mark a self-dual coefficient object, a
   scalar in a relation/potential, or a ground-order action whose numerical
   value is (4.1).

The first repair excludes noncrystallographic coefficients such as (2.6).
The second simply moves every missing Gram entry into categorical structure
constants or marked endofunctors.  Without a uniform bound on their ground
fields, heights, fusion ranks, or relation degrees, it is an exact but
tautological encoding.

## 5. Restriction of scalars measures the unavoidable cost

Let $k/\mathbf Q$ be totally real of degree $d$, let $V$ be the native
$k$-quadratic space, and let $s_r$ be a facet reflection.  On the
restriction of scalars $V_{\mathbf Q}$,

\[
 \ker(s_r+1)=kr,
 \qquad \dim_{\mathbf Q}\ker(s_r+1)=d.                \tag{5.1}
\]

### Proposition 5.1 (spherical rank cost)

A native degree-$d$ facet reflection cannot be induced on an ordinary
integral numerical Grothendieck group by one 2-spherical twist when $d>1$.
If it is a product of twists in pairwise Euler-orthogonal, numerically
independent spherical objects, at least $d$ objects are required.  More
generally, a spherical-functor realization must have numerical source rank
at least $d$.

#### Proof

Equation (2.1) has a one-dimensional $(-1)$-eigenspace.  A product of
commuting reflections in $t$ independent orthogonal vectors has
$(-1)$-eigenspace of dimension $t$.  Comparison with (5.1) gives
$t\ge d$, with equality for an exact reflection and no redundant negative
directions.  For a spherical functor, the nontrivial part of its twist lies
in the numerical image of its source, so that source has rank at least the
rank $d$ of (5.1).  \(\square\)

Fusion unfolding pays the same bill in a different currency.  For the
$(2,5,6)$ triangle, $d=2$ and the Fibonacci coefficient ring has
$\mathbf Z$-rank two.  A bounded saturated category with numerical rank
$R$ carrying the complete trace module would force

\[
                         d(n+1)\le R.                 \tag{5.2}
\]

But proving that all maximal arithmetic chambers admit such a carrier with
one uniform $R$ is already the desired dimension/field breakthrough.  It
cannot be assumed as a harmless categorical realization lemma.

## 6. Exact `H21` and `X3` replay

For a wall pair put

\[
 p_{ij}=\frac{4(r_i,r_j)^2}{(r_i,r_i)(r_j,r_j)}.       \tag{6.1}
\]

The frozen Gram packets give the following complete pair histograms.

| control | facets | Gram rank | corank | $p_{ij}$ histogram |
|---|---:|---:|---:|---|
| Borcherds `H21` | 210 | 22 | 188 | $0^{5796},1^{3465},2^{2016},4^{3780},9^{6048},16^{840}$ |
| Vinberg `X3` | 36 | 20 | 16 | $0^{501},1^{27},3^{36},4^{24},16^{30},25^{12}$ |

Here $p=0,1,2,3$ mean respectively $m=2,3,4,6$, $p=4$ is the
parallel threshold, and $p>4$ is ultraparallel.

Consequences:

* `H21` has $2,016$ $m=4$ pairs, so Proposition 2.1 rules out one
  spherical object per wall.  Its rational integral Cartan data do admit a
  valued-quiver carrier, but the direct facet quiver has 210 vertices while
  the actual Lorentzian quotient has rank 22.
* `X3` has $36$ $m=6$ pairs and the same individual-spherical failure.
  Its direct facet quiver has 36 vertices and a 16-dimensional radical.
* Compressing either quiver to an ambient lattice basis does not remove the
  cost: exact chamber extraction must then retain 210 or 36 marked spherical
  complexes, functors, Fourier--Mukai kernels, or an independently proved
  algorithm that reconstructs precisely those walls and no hidden ones.

The replay is executable:

```text
python3 scripts/certify_quiver_2cy_coxeter_route.py --check
python3 -m pytest -q tests/test_quiver_2cy_coxeter_route.py
```

The terminal report/script/test digests are pinned in
[`quiver_2cy_coxeter_route_freeze.sha256`](../data/enumerations/quiver_2cy_coxeter_route_freeze.sha256).

This is also why the controls do not prove an absolute categorical no-go.
As shown in
[`bounded_universal_mixed_ag_carrier_no_go_and_next_species.md`](bounded_universal_mixed_ag_carrier_no_go_and_next_species.md),
both controls have mixed K3/Fourier--Mukai action carriers.  The 42 short
`H21` reflections and 24 short `X3` reflections are spherical twists; the
remaining reflections are actual automorphism or Fourier--Mukai actions.
What the controls refute is the **uniform simple-spherical quiver** and the
claim that the direct quiver is a bounded compression.

## 7. Properness and Calabi--Yau structure do not imply boundedness

There are three distinct categories in the standard quiver package, and
conflating them creates a false finiteness argument.

1. For a finite acyclic ordinary quiver $Q$, the path algebra $\mathbf C Q$
   is finite-dimensional and hereditary.  Hence
   \(\operatorname {Perf}(\mathbf C Q)\) is smooth and proper.  Its
   symmetrized Euler form remembers the graph, but it is not 2-Calabi--Yau
   and vertex reflections are equivalences between different orientations,
   not automatically a fixed-category wall action.
2. The preprojective/2-Calabi--Yau completion supplies the desired 2-CY
   duality and spherical-twist action.  For a non-Dynkin (in particular
   hyperbolic/indefinite) graph, the preprojective algebra is
   infinite-dimensional.  One obtains useful Hom-finite or finite-dimensional
   subcategories, but not a single proper finite-dimensional algebra whose
   boundedness follows from the input's maximality.
3. A Ginzburg 3-CY completion is homologically smooth, and a Jacobi-finite
   potential can yield a Hom-finite 2-CY cluster quotient.  The 3-CY Euler
   form is skew-symmetric, so its individual spherical twists act as
   transvections rather than the orthogonal reflections of $B_P$.  Adding
   the symmetric Gram as separate quiver/potential decoration returns to the
   `BOUND` gate.

Two elementary families dispose of “proper/smooth should be enough.”

### Example 7.1 (fixed rank, unbounded Euler entries)

Let $Q_a$ have two vertices and $a$ parallel arrows from the first to the
second.  The algebras $A_a=\mathbf C Q_a$ are finite-dimensional,
hereditary, smooth, and proper for every $a\ge0$, while

\[
 \chi_{A_a}+\chi_{A_a}^{t}
   =\begin{pmatrix}2&-a\\-a&2\end{pmatrix}.           \tag{7.1}
\]

Thus even numerical rank two plus smoothness and properness does not bound
the Gram entry.  The endomorphism-algebra dimensions $a+2$ escape every
fixed bounded family.

### Example 7.2 (2-CY does not bound categorical size)

The Dynkin cluster/2-CY quiver categories of type $A_r$ are Hom-finite
2-Calabi--Yau categories for every $r$.  Their numbers of simple/cluster
directions grow with $r$.  If one insists on saturated geometric examples,
finite orthogonal sums of K3 categories are smooth proper 2-CY and have
unbounded numerical rank.  Connectedness or indecomposability would have to
be imposed and then proved from $P$; it is not a consequence of the words
“2-Calabi--Yau.”

Reflection maximality does not repair this.  A finite-volume hyperbolic
reflection group is infinite, hence its direct Cartan graph is not of finite
Dynkin type.  The absolutely maximal `H21` chamber is an explicit maximal
input whose 210-node direct Cartan/preprojective construction is indefinite
and nonproper at the algebra level.  Maximality says that no new hyperbolic
mirror enlarges the group; it says nothing about the dimension of a path
algebra, Jacobi algebra, unfolding, or Fourier--Mukai kernel.

## 8. Quiver varieties and stability do not close the gate

Nakajima quiver varieties geometrize representations of a fixed symmetric
Kac--Moody algebra and carry its Weyl action on homology.  For fixed

\[
                        (Q,\mathbf v,\mathbf w,\theta) \tag{8.1}
\]

they are finite-type quasi-projective symplectic varieties, with projective
Lagrangian fibres in the standard setup.  This does not provide the needed
global theorem because:

* the quiver $Q$, dimension and framing vectors \(\mathbf v,\mathbf w\),
  and stability parameter \(\theta\) are unbounded;
* the Kac--Moody Weyl group retains only crystallographic Cartan data;
* a general dotted distance or noncrystallographic finite label still needs
  folding/fusion/scalar decoration; and
* fixed-parameter GIT boundedness does not prove that reflection maximality
  places every $P$ in finitely many parameter packets (8.1).

Likewise, Bridgeland stability is an exact-return language, not by itself a
boundedness theorem.  For an indefinite quiver there can be infinitely many
real roots and spherical walls.  One must bound a chamber and prove that its
finite selected wall set is complete.  That is exactly `RETURN+BOUND`, not a
formal consequence of the existence of a stability manifold.

## 9. The conditional AG-native finiteness theorem

The part that would genuinely prove finiteness is short once the missing
boundedness theorem is supplied.

### Theorem 9.1 (bounded decorated categorical extraction)

Let \(\mathfrak M\) be a quasi-compact finite-type algebraic stack of
decorated smooth/proper dg categories or bounded quiver/fusion data.  Suppose
its objects include uniformly finite marked perfect objects/functors/kernels,
a stability or positive-cone marking, and an exact-return certificate.
Assume the complete marked Coxeter Gram extracted from an object is constant
on each stratum of a finite constructible stratification of
\(\mathfrak M\).  Then only finitely many marked Coxeter polytopes occur.

#### Proof

The extraction map has discrete target.  Each constructible stratum
contributes one value, and there are finitely many strata.  Equivalently, if
the map is locally constant, its point fibres form an open cover of the
quasi-compact underlying space and a finite subcover suffices.  No arithmetic
reflection finiteness enters.  \(\square\)

For quiver carriers, a usable version of the hypothesis must uniformly bound
at least

\[
\begin{gathered}
 |Q_0|,\ |Q_1|,\ \operatorname {rk}K_0(\mathcal C(Q)),
 \text{ all relation/potential degrees and coefficient fields},\\
 \text{dimension/framing vectors, marked spherical objects or functor
 sources, kernel classes and amplitudes},\\
 \text{and the complexity of the exact stability-wall certificate}.
\end{gathered}                                                        \tag{9.1}
\]

Under (9.1), standard representation schemes, GIT quotients, and moduli of
perfect complexes do produce finite-type parameter spaces.  Without (9.1),
the union is only an unbounded disjoint union of finite-type pieces.

## 10. Exact missing theorem and verdict

The quiver/2-CY route would achieve the user's global target if one proved:

> **Uniform folded-CY essential-surjectivity.**  Every absolutely
> reflection-maximal finite-volume arithmetic Coxeter chamber admits an exact
> folded 2-CY or smooth-proper categorical carrier satisfying one tuple of
> bounds (9.1), independent of the dimension, field, facet count, Coxeter
> labels, and dotted distances.

No known Kac--Moody, preprojective, Ginzburg, fusion-quiver, Nakajima-variety,
or stability theorem supplies this assertion.  In fact:

* ordinary symmetrizable theory is not total because of $m=5$;
* fusion unfolding is total for the finite-label part but its rank grows and
  general dotted distances require new decorations;
* the direct exact quiver has one vertex per facet and is visibly
  noncompressive on `H21` and `X3`;
* 2-CY completion of the hyperbolic quiver loses algebra-level properness;
* smooth proper path categories retain unbounded quiver/Euler complexity;
  and
* a hypothetical uniform numerical-rank bound would already imply
  $d(n+1)\le R$, the desired arithmetic cutoff in another form.

Therefore this route advances the real target by producing an exact
noncrystallographic folded 2-CY carrier and by identifying the precise rank
cost.  It does **not** currently provide the universal AG-native finiteness
pushback.

## 11. Primary sources and local dependencies

* Edmund Heng,
  [*Coxeter quiver representations in fusion categories and Gabriel's
  theorem*](https://arxiv.org/abs/2302.01866): Coxeter quivers, TLJ/Fibonacci
  coefficient categories, the unfolded ordinary quiver, abelian equivalence,
  reflection functors, and noncrystallographic root systems.
* Asilata Bapat, Anand Deopurkar, Anthony Licata,
  [*Spherical objects and stability conditions on 2-Calabi--Yau quiver
  categories*](https://arxiv.org/abs/2108.09155): the arbitrary-quiver 2-CY
  category, root-lattice Grothendieck group, and spherical-twist braid action.
* Claire Amiot,
  [*Preprojective algebras and Calabi--Yau duality*](https://arxiv.org/abs/1404.4764):
  the Dynkin/non-Dynkin preprojective distinction and 2-CY categories.
* Bernhard Keller,
  [*Deformed Calabi--Yau completions*](https://arxiv.org/abs/0908.3499):
  preprojective and Ginzburg Calabi--Yau completions and homological
  smoothness.
* Claire Amiot,
  [*Cluster categories for algebras of global dimension 2 and quivers with
  potential*](https://arxiv.org/abs/0805.1035): Hom-finite 2-CY cluster
  categories in the Jacobi-finite setting.
* Hiraku Nakajima,
  [*Quiver varieties and tensor products*](https://arxiv.org/abs/math/0103008):
  the geometric Kac--Moody action from quiver varieties.

Frozen local control packets:

```text
data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json
data/enumerations/vinberg_x3_h19.json
results/borcherds_h21_outer_routing_counterexample.md
results/projected_k3_universality_rank_saturation_no_go.md
results/k3_compact_maximal_256_phys_empty.md
results/bounded_universal_mixed_ag_carrier_no_go_and_next_species.md
```
