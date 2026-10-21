# A Divisor-Lattice Factory for Coxeter Polytopes

## 1. Purpose

The example track should search algebraic varieties, not isolated Gram
matrices.  The desired object is a projective variety whose intrinsic divisor
geometry produces a hyperbolic Coxeter polyhedron.  Its walls then come with
geometric completeness: they are the walls of a nef, ample, or movable cone,
not a guessed subset of roots.

The factory has now delivered a closed positive result in rank 12: effective
order-four K3 actions realize three archived complete-nef folds as genuine
finite-volume noncompact Coxeter polytopes in `H10`, with exact lattice-action,
glue, chamber-transport, and Vinberg certificates.  See the
[rank-12 Brandhorst--Hofmann/K3 lift theorem](../results/roulleau_rank12_bh_k3_lifts.md).

For surfaces this proposal is completely natural.  If \(X\) is a smooth
projective surface and

\[
N=\operatorname{NS}(X)/\text{torsion},
\]

then the intersection form has signature \((1,\rho-1)\) by the Hodge index
theorem.  A component of

\[
\{x\in N_{\mathbb R}:x^2>0\}/\mathbb R_{>0}
\]

is therefore hyperbolic space \(\mathbb H^{\rho-1}\).

This document distinguishes a rigorous integral construction, which is best
suited to noncompact finite-volume examples, from a totally-real enhancement
that could produce compact examples.

## 2. The integral construction

Let \(\mathcal C_X^+\) be the component of the positive cone containing an
ample class.  Suppose a finite collection of irreducible negative curves

\[
\Delta=\{C_1,\ldots,C_N\}
\]

has the following properties.

1. The nef cone is exactly
   \[
   \operatorname{Nef}(X)=
   \{x\in\overline{\mathcal C_X^+}:x\mathbin{\cdot}C_i\geq0
   \text{ for every }i\}.
   \]
2. Every wall reflection
   \[
   s_{C_i}(x)=x-\frac{2(x\mathbin{\cdot}C_i)}{C_i^2}C_i
   \]
   preserves \(N\).
3. The pairwise wall angles are Coxeter angles.

Then

\[
P_X=\mathbb P\bigl(\operatorname{Nef}(X)\cap\mathcal C_X^+\bigr)
\subset \mathbb H^{\rho-1}
\]

is a hyperbolic Coxeter polyhedron.  It is compact precisely when its closure
does not meet the projectivized isotropic cone.  Finite volume with cusps can be
certified from its parabolic vertex subdiagrams.

This is the first **factory proposition** to formalize in the verifier: an
input certificate for \((X,N,\Delta)\) proves the Néron--Severi lattice, proves
that \(\Delta\) is the complete wall set, and checks the displayed reflection
and cone conditions.  The output chamber then follows formally.

### Weak algebraic-carrier theorem

The newer examples show that full-nef equality should be the final promotion
gate, not the entry requirement.  Let an already certified rational Coxeter
cone have an even integral full-rank realization `N` of signature `(1,n)`,
with integral facet reflections and exact oriented halfspaces.  If `N`
embeds primitively in the K3 lattice, period surjectivity gives a projective
K3 surface with very-general `NS=N`.  The same cone is then automatically a
full-Picard algebraic divisor-reflection carrier.  Reflection-stable
overclosures and scaled duals are allowed provided their primitive covector
map rescales every oriented inequality by a positive number.

For this *fixed* oriented realization the statement is an equivalence: a
full-Picard weak complex-K3 carrier exists exactly when `N` embeds primitively
in the K3 lattice.  It is not yet an equivalence stated solely in terms of an
abstract Coxeter diagram.  Before claiming that stronger version one must
enumerate all positive crystallographic normal scalings.  Once a fixed
integral root span `L` is chosen, its same-form carrier lattices are the
reflection-stable even intermediate lattices `L <= N <= L^vee`; scaled duals
are useful additional realizations but are not an exhaustive substitute for
the normalization census.

Equivalently, if `r=rank(N)<=20`, the exact primitive-embedding gate is the
existence of an even lattice `T` of signature `(2,20-r)` with discriminant
form `q_T=-q_N`.  Gluing `N+T` along an anti-isometry gives the K3 lattice.
The stable sufficient shortcut is

```text
r <= 19 and 22-r >= length(A_N)+2.
```

Outside it one must run Nikulin's finite odd-prime and 2-adic determinant
tests; failure of the shortcut alone proves nothing.  Rank 20 is allowed but
has a positive-definite rank-two complement and zero-dimensional period
locus.  See the
[exact theorem audit](../results/k3_weak_algebraic_carrier_theorem_audit.md)
for the proof, local algorithm, and correct abstract-Coxeter quantifiers.

Five August 2026 calibration packages show the criterion operating without
any surface-automorphism assumption.  The `X2 H10/37` walls span
`U+E8(-2)+<-2>` and transport to an explicit nine-parameter family with that
very-general full Picard lattice.  The full `K3-S3-C2-H10/13` chamber now
has an independent determinant-one transport to that same geometric basis,
putting the 13-, 21-, and 37-wall systems on one explicit family.  The
`K3-D4 H12/17` walls span
`U+E8(-1)+A2(-1)+A1(-1)`; a visible-unit-minor primitive embedding gives a
seven-dimensional period family with that full Picard lattice.  In the
second case the order-three triality only discovers the lattice and walls;
it is not an automorphism of the returned generic surface.

The Hessian `H10/27` package proves a primitive embedding by a full
discriminant anti-isometry and gives a nine-dimensional full-Picard family.
The paired `A20 H12/19` and `H11/20` package reconstructs primitive embeddings
from one special K3 and gives seven- and eight-dimensional generic families.
In the latter pair, the reflection action makes the carrier denominator
canonical: exhaustive stable integral orders are `[1,2]`, but the only even
order is `[1]`, because the unique index-two overclosure is odd.  This is the
desired split between continuous existence and discrete rigidity: the period
domain supplies surfaces once the embedding passes, while the Coxeter group
selects the admissible Picard lattice and the Vinberg test separately decides
finite volume.

See the
[`H10/37` equation carrier](../results/x2_h10_37_ue8m2a1_full_picard_carrier.md),
[`H10/13` equation carrier](../results/k3_s3_c2_h10_13_ue8m2a1_full_picard_carrier.md),
[`H12/17` period carrier](../results/k3_d4_h12_full_picard_k3_carrier.md),
[`GHESS H10/27` period carrier](../results/ghess_h10_27_full_picard_k3_carrier.md),
and [`A20` paired period carriers](../results/a20_pair_full_picard_k3_carriers.md).

The first three rows now give a complete fixed-moduli compiler.  For
`F=U+E8(-2)+<-2>`, every primitive negative integral mirror has exactly one
of the types `(-2,div1)`, `(-2,div2)`, `(-4,div2)`, because `A_F` has
exponent two and integrality is equivalent to `|r^2|` dividing
`2*div_F(r)`.  Any genuine oriented finite-sided Coxeter polyhedron passing
this `F`-realization gate is carried weakly by a marked very-general member of
the explicit family.  The full reflection group has the nonuniform 13-wall
fundamental chamber, so finite volume is equivalent to finite reflection-
subgroup index and every finite-volume chamber in this carrier class is
noncompact.  The 13-, 21-, and 37-wall indices are exactly `1`, `6720`, and
`241920`.  See the
[common-moduli criterion](../results/ue8m2a1_common_k3_moduli_carrier.md).
The ambient surface equation does not provide equations or effective curve
representatives for arbitrary transported wall classes.

### Redoubling as a surface-to-polytope amplifier

The full 13-wall chamber on `F` contains two parallel doubling walls.  Their
primitive roots have square `-2`, divisibility two, inner product two, raw
product square four, and normalized squared ratio one.  Allcock's
finite-subtree theorem therefore applies: the two reflections generate an
infinite-dihedral Cayley line, and every `N`-vertex path gives a genuine
finite-volume noncompact `H10` Coxeter polytope with `N+12` facets,
reflection index `N`, and volume multiple `N`.  Closed root recurrences keep
every boundary reflection integral on the same full Picard lattice.  Thus
one nine-parameter surface family weakly carries an infinite tower of exact
polytopes.  See the
[redoubling-family theorem](../results/k3_s3_h10_redoubling_family.md).

This adds a reusable factory stage.  After finding one complete reflective
seed on a K3-admissible lattice, compute its graph of pairwise disjoint
doubling walls.  Every finite subtree of the corresponding reflection Cayley
tree is then a certified finite-index output, with no fresh unconstrained
Gram search.  Conversely, a proposed crystallographic Coxeter polytope has a
weak full-Picard K3 carrier exactly when its fixed oriented reflection
lattice embeds primitively in the K3 lattice.  The remaining inverse theorem
must classify which abstract Coxeter/Gram data force membership in one of
these algebraic seed classes; only then can empty surface strata imply an
abstract polytope no-go.

This yields four deliberately distinct output levels:

1. all mirrors are algebraic divisor classes in the full Picard lattice;
2. selected square-minus-two orientations are proved effective, preferably
   by exact decompositions in a complete curve stock;
3. the imported chamber contains an ample class; and
4. the imported chamber is the complete nef cone.

For level 3 it is enough to find a positive nef class `H` strict on the
imported walls.  If `A` is ample, then `mH+A` is ample and imported-strict for
all sufficiently large `m`.  For level 4, every primitive facet must be a
genuine `(-2)` root wall and no additional Picard-root hyperplane may cut the
interior.

There is now an exact unbounded-rank rational-surface theorem at level 1.
Every integral lattice `G` of signature `(1,m)` embeds after a common square
rescaling as

```text
D^2 G -> I_(1,1+4m)=NS(Bl_(1+4m)(P^2)).
```

The construction rationally diagonalizes `G`, writes the positive diagonal
coefficient as a difference of two squares, writes every negative coefficient
as minus four squares, and clears denominators.  A self-adjoint rational
operator on the embedded space extends after a scalar to an integral
self-adjoint Picard operator, represented on the rational surface by an
integral transpose-invariant sum of divisor-product correspondences.  Before
scaling, the rational orthogonal projector is itself a divisor-product Chow
idempotent, and the whole coefficient field acts on its Tate-motive image
under literal correspondence composition.  Consequently every
admissible totally-real arithmetic Coxeter datum has a semialgebraic divisor-
eigenslice carrier on a rational blowup.  For a carrier containing the full
trace form isometrically up to scale, Hodge index gives the converse signature
gate: the distinguished Lorentzian place must be the only indefinite real
place, and the Néron--Severi complement is negative definite.  This converse
does not cover an unrelated carrier that discards the conjugate summands.  The
general theorem and the explicit compact
Bugaenko carrier on `Bl_19(P^2_Q)` are in
`results/rational_surface_eigencarrier_universality.md`.

Canonical-class, effectivity, ample, Mori, and nef compatibility remain
separate requirements.  Indeed universality shows that weak surface
existence is too soft to constrain the polytope; those promotion conditions
are exactly where classification can acquire force.

### K3 specialization

For a projective K3 surface, every smooth rational curve has square \(-2\), and
the reflection is

\[
s_C(x)=x+(x\mathbin{\cdot}C)C.
\]

The nef cone is a fundamental chamber for the Weyl group generated by
\((-2)\)-reflections.  With the repository's positive-root convention, put

\[
B_{ij}=-(C_i\mathbin{\cdot}C_j).
\]

Distinct irreducible curves have nonnegative integral intersection.  Hence a
K3 wall system is strongly restricted:

- \(C_i\cdot C_j=0\) gives \(m=2\);
- \(C_i\cdot C_j=1\) gives \(m=3\);
- \(C_i\cdot C_j=2\) gives parallel walls; and
- \(C_i\cdot C_j>2\) gives ultraparallel walls.

Thus the finite Coxeter part is automatically simply laced.  This is a useful
geometric pruning theorem, not an empirical label cutoff.

### Immediate rank-10 model

For an Enriques surface,

\[
\operatorname{Num}(X)\simeq U\oplus E_8(-1)
\]

has signature \((1,9)\).  Configurations of smooth rational curves therefore
produce Coxeter polyhedra directly in \(\mathbb H^9\).  This is an excellent
calibration for the geometric factory and a source of rigid noncompact
finite-volume examples.  Isotropic genus-one classes force cusps, so this
model is not a compact dimension-9 route.

### Equivariant K3 fixed-section lift criterion

A symmetry of a divisor chamber is not automatically a surface
automorphism.  The successful fixed-section examples now give the following
reusable sufficient criterion.  Let $N\hookrightarrow\Lambda_{K3}$ be a
primitive hyperbolic lattice with orthogonal complement $T$, let
$a\in O(N)$ preserve the chosen complete mixed reflection chamber, and let
$\gamma:A_N\to A_T$ be the anti-isometry defining the primitive glue.
It is enough to find $b\in O(T)$ such that

\[
 \bar a\,\gamma=\gamma\,\bar b
\]

on the discriminant groups, together with a positive isotropic eigenperiod
$\omega\in T_{\mathbb C}$ for $b$.  The graph glue then extends
$a\oplus b$ integrally to the K3 lattice.  Compute

\[
 NS(X)=\Lambda_{K3}\cap\omega^\perp
\]

exactly, rather than assuming that it remains $N$, and require its fixed
lattice to be the intended divisor-section lattice.  Finally either prove
that the full anti-invariant lattice has no $(-2)$-vectors or exhibit a
generic invariant strict chamber point avoiding every root hyperplane.  The
resulting Hodge isometry preserves a Kähler chamber, and Strong Torelli gives
a genuine K3 automorphism.

The order of $b$ need not equal the order of the visible Picard action
$a$.  This is now a constructive search rule rather than a caveat: the
`A18-C06` and `A20-C01` Picard involutions fail the naive involution glue but
lift after using order-four transcendental actions.  The other certified
cases use special periods whose additional algebraic classes are entirely
anti-invariant.  Thus a failed order-two discriminant match is an instruction
to enumerate bounded higher-order complement isometries, not a carrier
no-go.

The dual `(16,6,1)` test gives the complementary rejection gate.  For a
generic Picard lift, the period forces the rank-six complement to have one
cyclotomic type.  Comparing its mod-two characteristic polynomial with the
required discriminant action leaves order four, but the Smith valuation of
$B-I$ would require $16\mid\det(B-I)$ while cyclotomy gives determinant
eight.  A special order-four glue exists, yet its added algebraic coinvariant
contains two negated $(-2)$-roots and therefore fails the Kähler gate.  This
is a generic-Picard boundary and a specific special-period obstruction, not a
global special-period no-go; see
[`a16_dual_c2_h11_29_k3_lift_obstruction.md`](../results/a16_dual_c2_h11_29_k3_lift_obstruction.md).

### Effective-action database intersection criterion

The successful rank-14 lifts turn the preceding construction into a second,
finite positive search criterion.  Let `(N,Delta)` be a hyperbolic lattice
with a **complete** K3 nef chamber and let `A` be a finite subgroup of its
weighted wall group.  An effective finite K3-action row supplies a K3 lattice
`L`, a symplectic kernel `H_0`, and a distinguished non-symplectic generator
`g` of transcendental value `n`.  Compute, with saturated integral kernels,

\[
 T=L^{H_0}\cap\ker\Phi_n(g),\qquad N_{\rm row}=T^\perp_L.
\]

The row realizes the requested fixed-section construction if all of the
following exact gates pass:

1. `N_row` is integrally isometric to `N`.  Signature and determinant are
   only prefilters; compare the Smith invariants and discriminant quadratic
   form, and then use either an explicit isometry or an applicable indefinite
   uniqueness theorem.
2. A subgroup of the induced action on `N_row` is conjugate to `A` as a
   lattice action.  A fixed-lattice match alone is insufficient: compare the
   fixed and coinvariant lattices, their actions, and the primitive glue
   graph, or exhibit a full integral intertwiner.
3. The action preserves a Kähler chamber.  Effectiveness of the database row
   supplies one such chamber; for transport to the named chamber, use an
   invariant strict interior class and the refined period map.  A rootless
   algebraic coinvariant is a particularly clean sufficient gate.
4. The literal fixed chamber passes the independent Vinberg endpoint and
   maximal-parabolic tests.  Completeness of `Delta`, rather than the database
   classification, is what identifies the resulting polytope and proves its
   volume type.

Under these hypotheses, chamber transport followed by Strong Torelli gives a
projective K3 surface with a genuine automorphism subgroup inducing `A`, and
the fixed nef cone is exactly the archived Coxeter polytope.  This is both a
positive theorem template and a practical search algorithm: cyclotomic rank,
discriminant form, and fixed/coinvariant glue reduce an unstructured period
search to finitely many effective rows.

The first exact intersections give nine K3 realizations of fixed sections
of complete nef chambers: three
rank-12 cases, five rank-14 cases, and the rank-16 `H14/17` case.  They also
give effective `C4 x C2`, `C8 x C2`, and `C16 x C2` enrichments, including
one K3 surface simultaneously realizing the `H11/16` and `H10/14` sections.
In each case the abstract chamber or arithmetic denominator was already known
or already certified; the new output is a genuine algebraic realization.  See
the [frozen positive status](../results/positive_breakthroughs_2026_08_07.md).

There is a further positive lesson in the rank-16 case.  A no-go for an
order-two automorphism inducing a given Picard involution does not exclude a
higher-order automorphism with that same Picard restriction.  The `H14/17`
action is induced by an order-four automorphism whose square is the identity
on `NS` and `-1` on the transcendental lattice.  Therefore a failed
same-order glue test should route the candidate into a higher-cyclotomic-row
search, not directly into the rejection ledger.

### Elliptic fiber-root inverse-design criterion

The order-six `E8+D4` realization makes the original algebro-geometric idea
constructive at equation level.  To design a target

\[
N=U+R(-1)
\]

with a prescribed finite action on the root system `R`, proceed as follows.

1. Realize `U` by a Jacobian elliptic fibration and split `R` among reducible
   Kodaira fibers.
2. Choose a CM fiber automorphism, possibly combined with a base action, whose
   permutation of fiber components is the desired Dynkin action.
3. Check minimality and total Euler number 24, then exclude unintended
   reducible fibers.
4. Prove that the very-general trivial lattice is the full Neron--Severi
   lattice.  This requires a period/moduli-dimension bound, Mordell--Weil
   rank and torsion control, and saturation; Shioda--Tate rank alone is only a
   lower bound.
5. Identify the complete nef chamber and its literal wall action, then apply
   the fixed/coinvariant/glue, invariant-Kahler, and Vinberg volume gates from
   the effective-action criterion.

For pairwise-distinct nonzero `1,a,b,c`, the family

\[
y^2=x^3+t^3(t-1)(t-a)(t-b)(t-c)
\]

has fibers `II*+I0*+4II`.  The order-six map
`(x,y,t)->(omega*x,-y,t)` fixes `U+E8(-1)` and cyclically permutes the three
outer `I0*` components, giving literal `D4` triality.  Very generally
`NS=U+E8(-1)+D4(-1)`, and the fixed complete-nef section is the certified
finite-volume noncompact `H11/13` chamber.  See the
[explicit-family theorem](../results/ue8d4_order6_weierstrass_family.md).

The component permutation is necessary but not sufficient: without generic
`NS` exactness, primitive glue, an invariant Kahler chamber, and the complete
Vinberg tests it is only a diagram coincidence, not a hyperbolic-polytope
construction.

For finite group actions these requirements now admit an exact inverse
formulation.  A specified K3-lattice extension and period character determine
an equivariant period stratum; complete separation from every ambient
square-minus-two root identifies a proposed fixed cone with the invariant
nef section.  The section is the pullback of the actual nef cone of the
normal quotient surface.  This does not by itself make it Coxeter or
finite-volume: geometric angles, integral folded reflections, and Vinberg's
endpoint/parabolic test remain separate.  See the
[finite-action inverse-carrier theorem](../results/k3_finite_action_invariant_nef_inverse_carrier_theorem.md).

Imposing the symmetric fiber positions `+/-r,+/-s` adds an order-four
base/fiber symmetry to the order-six model.  The resulting one-parameter
subfamily has a cyclic `C12` action.  At `(r,s)=(1,2)`, exact reduction at
`p=7` proves that the named rational member still has geometric Picard rank
14: four direct point counts reconstruct the rank-eight Frobenius factor,
whose normalized roots are proved noncyclotomic.  This supplies a useful
arithmetic terminal gate for the inverse-design method: after proving a
very-general statement by period dimension, a good-reduction Frobenius
certificate can pin an individual number-field carrier.  See the
[named-member theorem](../results/ue8d4_named_c12_member.md).

A second instance reaches dimension 14.  For very general parameters on an
explicit nonempty open,

\[
y^2=x^3+t^2(t^2-3)x+t^3(2+bt^2+ct^4)
\]

has fibers `I2*+II*+6I1` and exact
`NS=U+E8(-1)+D6(-1)`.  The map `(x,y,t)->(-x,i*y,-t)` exchanges the two
spinor components of the `I2*` fiber.  The local tangent-branch resolution
and the order-four Lefschetz congruence independently force this outer
`D6` action.  Fixed/coinvariant/discriminant/glue comparison then identifies
Kondo's involution, and the fixed complete-nef section is Vinberg `P14`, the
finite-volume noncompact `H14/17` chamber.  See the
[explicit rank-16 theorem](../results/ue8d6_order4_weierstrass_family.md).

The same recipe has now produced an explicit rank-12 fixed chamber.  With
`s=t^2`, the normalized four-parameter model

\[
y^2=x^3+A(t^2)x+tC(t^2)
\]

where

\[
\begin{aligned}
A(s)&=-3+\alpha(s-1)+\beta(s-1)^2,\\
C(s)&=2-(\alpha+1)(s-1)+\gamma(s-1)^2+\delta(s-1)^3,
\end{aligned}
\]

has fibers `II*+2I2+10I1` and very-general
`NS=U+E8(-1)+A1(-1)^2`.  The order-four action
`(x,y,t)->(-x,i*y,-t)` exchanges the two `I2` root components, and the exact
BH/glue handoff produces the complete finite-volume noncompact `H10/12`
fixed section.  See the
[explicit rank-12 theorem](../results/ue8a1x2_order4_weierstrass_family.md).

A second four-parameter solution replaces the `II*` root stock by `I4*`.
The explicit polynomials `P,Q` in

\[
y^2=x^3+P(t^2)x+tQ(t^2)
\]

force `I4*+2I2+10I1` and very-general
`NS=U+D8(-1)+A1(-1)^2`.  Equivariant branch resolution proves that the
order-four map fixes the `D8` components individually, while swapping the
two `A1` roots.  Shioda's height formula excludes torsion and therefore
closes the primitive-lattice gate.  The result is the complete `H10/15`
section; see the
[explicit `D8+A1^2` theorem](../results/ud8a1x2_order4_weierstrass_family.md).

The third rank-12 stock is

\[
y^2=x^3+t^2P(t^2)x+t^3Q(t^2),
\]

where `P,Q` are the four-parameter normalized polynomials in the dedicated
theorem.  It forces `2I0*+2I2+8I1` and
`NS=U+D4(-1)^2+A1(-1)^2`.  The invariant Tate-root coordinate at each
`I0*` point proves both `D4` diagrams are fixed vertexwise, while the two
`A1` roots are exchanged.  This gives the complete `H10/38` section; see
the [explicit `D4^2+A1^2` theorem](../results/ud4x2a1x2_order4_weierstrass_family.md).

Thus the complete rank-12 root-stock ladder

```text
E8+A1^2  -> H10/12,
D8+A1^2  -> H10/15,
D4^2+A1^2 -> H10/38
```

is now equation-level.  In each case the same five terminal gates—fiber
open, generic Picard exactness, Mordell--Weil/saturation, literal component
action and glue, complete Vinberg volume—are independently checkable.

There is an especially strong shortcut when the fiber-root lattice is
unimodular.  The classical Inose pencil

\[
y^2=x^3-3At^4x+t^5(t^2-2Bt+1)
\]

forces two `II*` fibers and hence `U+E8(-1)^2`.  On the explicit open where
the residual quartic has four simple roots, its two independent coefficients
give the full two-dimensional lattice-polarized moduli space.  Therefore the
very-general Picard lattice is exactly the forced unimodular lattice, with
no primitive-closure, torsion, or glue alternatives.  The full nef cone is
the certified 19-wall finite-volume noncompact chamber in `H17`.  See the
[Inose-family theorem](../results/ue8x2_inose_h17_family.md).

One controlled discriminant collision raises the Picard rank by one without
leaving the same algebraic reservoir.  Setting `A=s^2`, `B=1-s^3` gives

\[
y^2=x^3-3s^2t^4x+t^5(t^2-2(1-s^3)t+1),
\]

with fibers `II*+II*+I2+2I1` and very-general
`NS=U+E8(-1)^2+A1(-1)`.  The resulting full nef cone is the complete
24-wall finite-volume noncompact chamber in `H18`.  See the
[Inose-collision theorem](../results/ue8x2a1_inose_h18_family.md).

At Picard rank 20 the same philosophy admits a mixed-root extension.  The
rigid Inose model `y^2=x^3+t^5(t-1)^2` forces
`U+2E8(-1)+A2(-1)` directly.  The order-three quotient model
`y^2=x^3+(t^3-1)^4` reaches the same lattice from `E6(-1)^3`; its visible
`3`-torsion section constructs the diagonal index-three glue.  The complete
Vinberg chamber uses 24 square-minus-two curve walls and 12 additional
square-minus-six divisibility-three divisor cuts.  This gives an explicit
`H19/36` algebraic carrier while preserving the distinction between a
mixed reflection chamber and the full Mori-wall cone.  See the
[two-carrier theorem](../results/vinberg_x3_h19_explicit_carriers.md).

This gives a reusable **collision ladder**: start with a full-dimensional
unimodular fiber-root family, then impose a simple discriminant collision
that inserts one prescribed ADE root while leaving the specialized family
full-dimensional in the new lattice-polarized moduli space.

This suggests a priority rule for the positive search: among fiber-root
designs, test unimodular trivial lattices first.  They convert the usual
Picard/saturation/Mordell--Weil cluster into a single moduli-dimension gate,
after which only completeness and Vinberg volume remain.

## 3. Why an ordinary divisor lattice cannot give the compact target

There is a clean arithmetic obstruction.  A divisor lattice is an integral
quadratic lattice over \(\mathbb Q\).  In rank at least 5, every indefinite
rational quadratic form is isotropic by Meyer's theorem.  Consequently its
arithmetic orthogonal group is nonuniform.

If a finite-sided reflection chamber for a subgroup
\(W\leq O^+(N)\) were compact, then \(W\backslash\mathbb H^{\rho-1}\) would have
finite volume and \(W\) would have finite index in the arithmetic lattice
\(O^+(N)\).  Nonuniformity is invariant under finite index, a contradiction.
Therefore:

> A Coxeter chamber obtained from an ordinary integral Néron--Severi lattice
> cannot be compact in hyperbolic dimension at least 4.

The divisor-cone idea remains highly useful for noncompact examples, but a
compact dimension-9 construction needs extra structure.

## 4. Totally-real eigenlattice enhancement

Let \(F\) be totally real and suppose the divisor lattice \(N\) carries a
self-adjoint ring embedding

\[
\iota:\mathcal O_F\hookrightarrow
\operatorname{End}_{\mathbb Z}(N).
\]

This is an endomorphism structure on the divisor lattice, not an action on the
variety by automorphisms and not necessarily an action by isometries of the
intersection form.  A concrete construction begins with an
\(\mathcal O_F\)-module \(M\), an \(F\)-bilinear form \(h\), and a scalar
\(\alpha\) in the codifferent, and sets

\[
(x,y)_N=\operatorname{Tr}_{F/\mathbb Q}\bigl(\alpha h(x,y)\bigr).
\]

The sign condition belongs to the scaled forms

\[
q_\sigma=\sigma(\alpha)h_\sigma.
\]

Assume that \(q_{\sigma_0}\) has signature \((1,d)\), while every other
\(q_\sigma\) is negative definite.  Then the trace form has the single
positive direction required by the Hodge index theorem.  Under restriction of
scalars,

\[
N_{\mathbb R}\simeq\bigoplus_{\sigma:F\hookrightarrow\mathbb R}M_\sigma,
\]

and the distinguished summand
\(E=M_{\sigma_0}\), equipped with \(q_{\sigma_0}\), contains a copy of
\(\mathbb H^d\).  In the concrete quadratic-field examples below,
\(\alpha=-\phi/\sqrt5\) is negative at both real places; equivalently,
\(h_{\sigma_0}\) has signature \((d,1)\) and the conjugate form is positive
definite.  This is the convention used by the Bugaenko, discriminant-11, and
discriminant-19 certificates.

Define the eigen-nef section

\[
P_{X,F}=\mathbb P\bigl(
\operatorname{Nef}(X)\cap E^+
\bigr).
\]

A trace-lattice \((-2)\)-class defines a mirror of the restricted arrangement
when its hyperplane meets \(E^+\).  Such a mirror is a facet of
\(P_{X,F}\) only if it actually supports the chosen section.  If the
associated \(F\)-reflection preserves \(M\), it preserves the complete
restricted root arrangement; it need not preserve one chosen nef section.
Unlike the full rational
orthogonal group, the \(F\)-linear group is compact at all conjugate real
places and can be a **uniform** arithmetic lattice.  This is the arithmetic
mechanism used by compact hyperbolic reflection groups over totally real
fields, now placed inside a divisor lattice.

For a general totally-real field the **weak full-Picard** construction is now
an exact theorem: it exists precisely when the integral trace lattice admits
an orientation-compatible primitive embedding in the chosen K3/Beauville
lattice.  Visible divisor walls of each allowed square also lie in a finite,
explicitly bounded field-norm alphabet.  See
`results/totally_real_eigencarrier_criterion.md`.  The stronger eigen-nef
promotion still requires a complete wall theorem.  For the fixed golden
scalar and an **even** \(\mathbb Z[\phi]\)-form, that promotion is available
after allowing a K3 Weyl-conjugate marking.
The five possible visible trace-root norms collapse to the two even shells
\(2\) and \(4-2\phi\), which form one unit-square orbit.  Therefore the
restriction of the complete K3 \((-2)\)-arrangement is exactly the norm-two
field-reflection arrangement for every such lattice.  A regular time-like
vector selects one of its chambers, and period theory plus Weyl transitivity
realizes that chamber exactly as an eigen-nef section.  See
`results/k3_eigen_nef_factory_theorem.md`.

This existential Weyl framing does not say that an arbitrary fixed marking
has a nonempty section, nor does it decide whether the selected chamber is a
finite polytope.  The remaining conditions are explicit:

1. the trace lattice must be even and admit a primitive embedding in the K3
   lattice;
2. the lattice-polarized K3 moduli stratum must be nonempty;
3. outside the even golden case, all divisor walls meeting the
   \(\sigma_0\)-eigenspace must be determined;
4. outside that case, their projected normals must be integral
   \(F\)-reflection roots;
5. the selected arrangement chamber must have finitely many walls and the
   required compactness; and
6. the facet reflections must preserve the full restricted arrangement and
   carry the chamber to its adjacent chambers.  They should not be required
   to preserve one nef section setwise: a reflection fundamental chamber is
   not invariant under its own facet reflections.

These are algebraic, lattice-theoretic, and real-semialgebraic tests.  A
failure produces a precise obstruction rather than an inconclusive Vinberg
timeout.

In dimension nine the surface carrier is itself a finite arithmetic moduli
problem.  A rank-20 trace lattice \(N\) is realized precisely when an even
positive-definite binary lattice \(T\) has discriminant form
\(q_T\simeq-q_N\); gluing along an anti-isometry gives the K3 lattice.  There
are only finitely many reduced binary forms of determinant \(|\det N|\).
Each oriented survivor determines a singular K3 surface by Shioda--Inose.

### The algebraic meaning of chamber volume

For the full nef chamber of a complex projective K3 surface of Picard rank at
least three, there is an exact intrinsic criterion:

\[
\operatorname{vol}_{\rm hyp}\mathbb P\operatorname{Nef}(X)<\infty
\Longleftrightarrow
 [O^+(\operatorname{NS}(X)):W(X)]<\infty
\Longleftrightarrow
|\operatorname{Aut}(X)|<\infty
\Longleftrightarrow
\operatorname{Eff}(X)\text{ is rational polyhedral}
\Longleftrightarrow
\operatorname{Cox}(X)\text{ is finitely generated}.
\]

Indeed, the nef cone is a fundamental chamber for the \((-2)\)-Weyl group,
the integral orthogonal group is an arithmetic lattice, Sterk identifies
finiteness of its quotient by the Weyl group with finiteness of
\(\operatorname{Aut}(X)\), and Artebani--Hausen--Laface identify the latter
with rational polyhedrality of the effective cone and finite Cox-ring
generation.  Moreover, for nef divisors on a surface
\(\operatorname{vol}_X(D)=D^2\), and on the slice \(D^2=1\) one has

\[
\frac12\operatorname{Hess}\bigl(-\log\operatorname{vol}_X\bigr)(u,v)
=-u\cdot v.
\]

Thus the hyperbolic measure is itself the Hessian measure of asymptotic
section growth.  Cusps are nef isotropic rays, hence genus-one fibrations.
This explains both finite volume and the Meyer noncompactness obstruction in
algebro-geometric terms.  The full proof and precise quotient distinction are
recorded in `results/k3_nef_covolume_dictionary.md`.

This already supplies high-dimensional calibration objects. The published
finite-automorphism K3 configurations of Picard ranks 15 through 19 have
complete full-nef Coxeter chambers in hyperbolic dimensions 14 through 18;
they are noncompact of finite volume, with cusps supplied by genus-one fibre
classes. The classification used by Artebani--Hausen--Laface has no rank-20
case, so dimension 18 is sharp for a **full** finite-volume K3 nef chamber.
This does not bound eigen-nef slices or arbitrary Coxeter polytopes. Exact
wall data and source locators are in
`results/k3_high_rank_nef_chambers.md`.

The real-multiplication construction uses a weaker invariant.  If
\(E\subset\operatorname{NS}(X)_{\mathbb R}\) is the distinguished eigenspace,
restrict every effective-degree intersection functional to \(E\):

\[
\operatorname{VisEff}_{E}(X)=
\overline{\operatorname{cone}}
\{(D,\cdot)|_E:H^0(X,\mathcal O_X(D))\ne0\}\subset E^*.
\]

Then \(\operatorname{Nef}(X)\cap E\) is the dual visible cone.  The whole Cox
ring may be infinitely generated while this projected cone is polyhedral and
its unit-volume dual section is compact.  Accordingly, G4 tests projected
wall completeness, Hessian covolume, and visible isotropic rays rather than
finite generation of the ordinary Cox ring.  The self-adjoint action can be
packaged over \(\mathbb Q\) as an algebraic correspondence in
\(CH^2(X\times X)_{\mathbb Q}\), so this is an invariant of an enhanced
algebro-geometric pair, not merely of an abstract real slice.

## 5. The dimension-9 rank coincidence

For a complex K3 surface, \(\rho\leq20\).  If \(F\) is real quadratic and
\(M\) has \(\mathcal O_F\)-rank \(d+1\), then the underlying divisor lattice has
rank \(2(d+1)\).  At \(d=9\),

\[
2(d+1)=20.
\]

Thus a compact \(\mathbb H^9\) attempt over
\(F=\mathbb Q(\sqrt5)\) uses a Picard-rank-20 K3 surface exactly, with no spare
divisor directions.  Its transcendental lattice has rank 2, so the K3 is a
singular K3 surface and can in principle be constructed through its period or
Shioda--Inose data.  In dimension 8 the trace lattice has rank 18, leaving a
2-dimensional K3 moduli stratum; this makes Bugaenko's chamber the natural
calibration.

The first trace calculation is positive.  For the now-audited Bugaenko
normalization, the codifferent scalar

\[
\alpha=-\phi/\sqrt5
\]

gives the exact trace lattice

\[
\langle2\rangle\oplus\langle-2\rangle
\oplus E_8(-1)^{\oplus2}.
\]

It is even of signature \((1,17)\), determinant \(-4\), and embeds primitively
in the K3 lattice with complement
\(U\oplus\langle2\rangle\oplus\langle-2\rangle\).  Thus the associated
2-dimensional lattice-polarized K3 stratum is nonempty.  The exact basis-change
and embedding certificates are in `results/bugaenko_k3_trace.md`.  The
five-shell parity reduction and exact Vinberg completion now recover the
11-wall chamber as an eigen-nef section, existentially up to a K3
Weyl-conjugate marking.  The complete calculation, including all 320 mirrors
meeting the closed section, is in `results/bugaenko_k3_eigen_nef.md`.

The coincidence gives a bounded first search:

1. use the completed Bugaenko rank-9 construction as the mandatory
   eigen-nef calibration;
2. enumerate rank-10 forms whose rank-20 trace lattices have a
   positive-definite rank-2 K3-lattice complement; and
3. search the resulting zero-dimensional K3 strata for compact eigen-nef
   sections.

The discriminant-form embedding test is expected to be extremely selective.
It replaces a broad rank-10 lattice census by a moduli problem with a fixed
rank budget and exact local obstructions.

It already excludes the first natural seed.  For

\[
L_9^{(0)}=\langle-2\phi\rangle\oplus E_8\oplus A_1,
\]

the same scalar produces

\[
N_9\simeq\langle2\rangle\oplus\langle-2\rangle^{\oplus3}
\oplus E_8(-1)^{\oplus2},
\qquad A_{N_9}\simeq(\mathbb Z/2)^4.
\]

A primitive embedding in the rank-22 unimodular K3 lattice would have a
rank-2 complement with the same discriminant group.  A rank-2 lattice has
discriminant-group generator length at most 2, whereas \(A_{N_9}\) has length
4.  Hence this trace lattice is not K3-embeddable.  See
`results/l9_seed_k3_obstruction.md` for the exact certificate.  This rules out
the seed only from the K3 factory, not from reflectivity in its own right.

More generally, the same argument excludes every orthogonal rank-one extension
\(L_8\oplus\langle c\rangle\) with an even trace form.  Its rank-two trace
block in the basis \(e,\phi e\) has shape

\[
\begin{pmatrix}t&s\\s&t+s\end{pmatrix}.
\]

Evenness forces both \(t\) and \(s\) to be even, so this block adds at least two
generators to the 2-primary discriminant group.  Together with the two already
present in the Bugaenko trace lattice, the generator length is at least 4.
The full proof is `results/orthogonal_extension_k3_obstruction.md`.
Consequently the next rank-10 candidates must use nonorthogonal coupling or
overlattice gluing; appending a positive rank-one factor can no longer consume
search time.

The smallest coupled search already produces a K3-realizable survivor:

\[
H_{11}=
\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}\oplus E_8.
\]

It has signatures \((9,1)\) and \((10,0)\) at the two real places.  Its trace
lattice is

\[
N_{11}\simeq U\oplus
\begin{pmatrix}-2&-1\\-1&-6\end{pmatrix}
\oplus E_8(-1)^{\oplus2},
\]

of determinant \(-11\), with an explicit primitive K3 embedding.  The
corresponding singular K3 has transcendental lattice
\(\left(\begin{smallmatrix}2&1\\1&6\end{smallmatrix}\right)\) and an explicit
Inose model over \(\mathbb Q\).  This is an algebraic-geometric carrier, not
yet a Coxeter example; see `results/coupled_discriminant11_candidate.md`.

### Five-shell lemma for restricted K3 mirrors

There is a finite arithmetic filter before any Vinberg enumeration.  Write the
field norm of an integral vector as

\[
x=h(\delta,\delta)=a+b\phi,\qquad a,b\in\mathbb Z.
\]

For \(\delta\) to be a \((-2)\)-class in the trace Néron--Severi lattice,

\[
\operatorname{Tr}\!\left(-\frac{\phi}{\sqrt5}x\right)=-(a+b)=-2,
\]

so \(a=2-b\).  If its trace hyperplane meets the positive cone in \(E\), its
normal is negative for \(q_{\sigma_0}\), hence positive for
\(h_{\sigma_0}\).  The conjugate form \(h_{\sigma_1}\) is positive definite.
Thus \(x\) is positive at both real places, and

\[
2+b(\phi-1)>0,\qquad 2-b\phi>0.
\]

Exact comparison gives \(-3\leq b\leq1\).  Therefore every relevant
restricted K3 root belongs to exactly one of the five norm shells

\[
5-3\phi,\quad4-2\phi,\quad3-\phi,\quad2,\quad1+\phi.
\]

This is only a necessary condition for a restricted trace-root mirror.  A
vector in one of these shells is an integral field-reflection root precisely
when

\[
\frac{2h(M,\delta)}{h(\delta,\delta)}\subseteq\mathcal O_F.
\]

Both filters are implemented exactly in
`src/coxeter_polytopes/geometric/k3.py`.  The remaining search is a five-shell
root problem followed by chamber identification, reflection, Coxeter-angle,
and compactness tests---not an unrestricted scan of divisor classes.

For an **even** \(\mathcal O_F\)-form, meaning
\(h(v,v)\in2\mathcal O_F\) for every \(v\in M\), the mirror set simplifies
decisively.  Only the shells

\[
2\qquad\text{and}\qquad4-2\phi=2\phi^{-2}
\]

can occur.  If \(h(\delta,\delta)=4-2\phi\), then

\[
h(\phi\delta,\phi\delta)=\phi^2(4-2\phi)=2.
\]

Multiplication by \(\phi\) changes the ambient trace-lattice vector, but on
\(E=M_{\sigma_0}\) it scales the normal by the positive number
\(\sigma_0(\phi)\).  It therefore leaves the **restricted** mirror unchanged.
Conversely, a vector \(r\in M\) of field norm 2 has trace square

\[
\operatorname{Tr}_{F/\mathbb Q}(2\alpha)=-2
\]

and its restricted mirror meets \(E^+\).  Its field reflection is

\[
s_r(v)=v-h(v,r)r,
\]

which preserves \(M\) because \(h\) is integral.  Consequently, under the
stated trace-lattice and K3-realization hypotheses,

> the complete restricted K3 root arrangement of an even field form is
> exactly the arrangement of mirrors of norm-2
> \(\mathcal O_F\)-reflection roots.

This is an equality of arrangements in \(E\), not an equality of the ambient
trace hyperplanes.  It proves mirror completeness, but it does not prove that
a chamber has finitely many facets or compact closure.  The forms \(L_8\),
\(H_{11}\), and the discriminant-19 carrier are all even.

For \(L_8\), the separately audited Vinberg completion supplies exactly the
missing chamber theorem.  The restricted K3 arrangement has Bugaenko's
compact 11-facet chamber as a cell, with no additional divisor facet.  Exact
finite-parabolic closure at its 41 vertices finds all 320 mirrors meeting the
closed cell: 11 facets and 309 lower-face contacts.  A generic period has
\(\operatorname{NS}(X)=N\); after conjugating the self-adjoint embedding
\(\iota\) by a K3 Weyl element, the chosen ample chamber cuts the eigenspace
in precisely this cell.  Root reflections act trivially on \(A_N\), so the
Weyl conjugation extends across the primitive K3 embedding by the identity on
\(N^\perp\).  This is the first completed compact output of the factory; see
`results/bugaenko_k3_eigen_nef.md`.

For \(H_{11}\), this filter can be applied to existing algebraic-geometric
data.  Shimada's discriminant-11 K3 has exactly the reduced transcendental
lattice \(T_{11}^{\mathrm{red}}\).  Of his 347 published smooth-rational-curve
representatives, 17 pass both the two-place positivity and field-reflection
tests.  They span only rank 9 over \(F\), and one step under all 767 published
automorphism generators adds no survivor.

There is a stronger cone gate.  Two oriented effective \((-2)\)-classes in
the published data have timelike projections to the identity eigenspace in
opposite time-cone components.  Reverse Lorentzian Cauchy--Schwarz implies
that every positive eigenvector violates one of their two nef inequalities.
Thus the eigen-nef section of this fixed marking is empty.  This
**opposite-timecone divisor test** is finite and must be applied before root
enumeration for every future marking.  A Weyl-conjugate marking changes the
eigenspace relative to the nef cone and therefore requires a separate
chamber-walk problem.  That problem has now passed its first gate: an exact
integral conjugate has a positive identity-eigenvector strictly inside
Shimada chamber 68.  Its 17 local curve walls pass the field-reflection filter,
and the 347 published curve representatives yield 18 conjugated survivors of
field span rank 10.  This proves nonemptiness, not completeness; see
`results/shimada_disc11_conjugate_marking.md`.

### Restricted-root chambers and the finite quotient theorem

The preceding fixed-chamber view has an intrinsic formulation.  Let
\(E=M_{\sigma_0}\) be the Lorentzian summand defined above, and restrict the
complete trace-lattice root arrangement to it:

\[
\mathcal H_E=
\bigcup_{\substack{\delta\in N,\;\delta^2=-2\\
\delta^\perp\cap E^+\ne\varnothing}}
(\delta^\perp\cap E^+),
\]

where orthogonality is for the trace form on \(N\).  A **restricted cell** is
a connected component of \(E^+\setminus\mathcal H_E\).

#### K3 chamber restriction and Weyl conjugation

Let \(\mathcal H_N\) be the full \((-2)\)-root arrangement in the chosen
component of the positive cone of \(N\).  Every restricted cell \(C\) lies in
a unique ambient root chamber \(D\), and in fact

\[
C=D\cap E^+.
\]

Indeed, an ambient root chamber is convex, so its intersection with \(E^+\)
is connected; conversely, a path in a restricted cell crosses no ambient
root hyperplane.  For a projective K3 surface, the Weyl group \(W(N)\)
generated by trace-lattice \((-2)\)-reflections is transitive on the ambient
root chambers, and the ample cone is one of them.  Hence some \(w\in W(N)\)
satisfies

\[
wD=\operatorname{Amp}(X),
\qquad
wC=\operatorname{Amp}(X)\cap wE^+.
\]

Here \(wE\) is the distinguished eigenspace for the conjugated
\(\mathcal O_F\)-action \(a\mapsto w a w^{-1}\).  Conversely, every nonempty
ample section arising from such a Weyl-conjugate marking is obtained from a
restricted cell in this way.  This is an existence correspondence; it does
not say that distinct cells give inequivalent markings.

#### What arithmetic cocompactness makes finite

Work projectively and put

\[
\Gamma_h=\operatorname{PO}^+(M,h).
\]

Because \(h\) is definite at the conjugate real place, it has no nonzero
isotropic vector over \(F\).  The associated orthogonal group is therefore
\(F\)-anisotropic, and the Borel--Harish-Chandra--Godement compactness
criterion gives

\[
\Gamma_h\backslash\mathbb H(E)\quad\text{compact}.
\]

For any fixed finite set of field-norm shells, the corresponding integral
mirror arrangement is locally finite.  To see this, a mirror meeting a fixed
compact subset of \(\mathbb H(E)\), together with the fixed norm at
\(\sigma_0\), bounds the \(\sigma_0\)-coordinates of its normal.  Definiteness
and the fixed conjugate norm bound its other real coordinates.  The Minkowski
embedding of \(M\) is discrete, so only finitely many normals occur.

The following consequences are finite and no stronger:

1. there are finitely many mirror orbits and restricted-cell orbits under
   \(\Gamma_h\);
2. for a cell \(C\), its stabilizer
   \(\Gamma_C=\{\gamma\in\Gamma_h:\gamma C=C\}\) has finitely many orbits on
   the facets of \(C\); and
3. \(\overline C/\Gamma_C\) is compact, where the closure is taken inside
   \(\mathbb H(E)\).

The third statement does **not** imply that \(C\) is compact or finite-sided.
A lifted chamber can be periodic under an infinite stabilizer.  Properness of
the arithmetic action gives the precise criterion

\[
\boxed{\quad
\overline C\text{ is compact}
\quad\Longleftrightarrow\quad
\Gamma_C\text{ is finite}.
\quad}
\]

If \(\Gamma_C\) is finite, local finiteness also implies that \(C\) has only
finitely many facets.  If \(\Gamma_C\) is infinite, compactness of the finite
quotient records only periodicity; it gives no compact lifted polytope.

#### When a restricted cell is Coxeter

Suppose each facet of a restricted cell is supported by a vector \(r\in M\)
whose field reflection

\[
s_r(v)=v-\frac{2h(v,r)}{h(r,r)}r
\]

preserves \(M\).  It then preserves the trace form and the complete
trace-root arrangement, and it exchanges the two cells adjacent to its
mirror.  If the cell is finite-sided, Vinberg's reflection-chamber theorem
implies that it is a Coxeter polyhedron; the dihedral labels should still be
recomputed exactly as an independent audit.  In particular, integral facet
reflections together with finite \(\Gamma_C\) produce a compact Coxeter
polytope.

For an even field form, the preceding norm-2 mirror theorem makes every
restricted mirror reflective.  Let \(W_2\) be the projective group generated
by all norm-2 field reflections.  It is normal in \(\Gamma_h\) and acts
transitively on the restricted cells.  After projectivizing away the central
sign, the chamber action gives

\[
\Gamma_h/W_2\simeq\Gamma_C.
\]

Consequently the following conditions are equivalent in the even case:

\[
\begin{aligned}
&\text{a restricted norm-2 chamber is compact};\\
&\text{its stabilizer \(\Gamma_C\) is finite};\\
&[\Gamma_h:W_2]<\infty.
\end{aligned}
\]

Thus arithmetic cocompactness does not prove reflectivity.  It turns
compactness of the desired chamber into the still-substantive finite-index
problem for \(W_2\).

> **Restricted K3 chamber theorem (even quadratic case).**  Let
> \(F=\mathbb Q(\sqrt5)\), \(\mathcal O_F=\mathbb Z[\phi]\), and
> \(\alpha=-\phi/\sqrt5\).  Let \((M,h)\) be an integral even
> \(\mathcal O_F\)-lattice such that \(h_{\sigma_0}\) has signature
> \((d,1)\) and \(h_{\sigma_1}\) is positive definite.  Assume that the
> trace lattice \(N=(M,\operatorname{Tr}(\alpha h))\), with its
> self-adjoint \(\mathcal O_F\)-action, is realized as
> \(\operatorname{NS}(X)\) for a projective K3 surface.  Then the complete
> restricted K3 root arrangement in \(E=M_{\sigma_0}\) is the norm-2 field
> mirror arrangement.  Every one of its cells is the ample section of a
> Weyl-conjugate marking.  If
> \(W_2\leq\Gamma_h=\operatorname{PO}^+(M,h)\) is generated by the norm-2
> field reflections, then a restricted cell is a compact Coxeter polytope if
> and only if \([\Gamma_h:W_2]<\infty\).  Finite index, or equivalently a
> finite cell stabilizer, is an additional hypothesis to be certified; it
> does not follow from cocompactness of \(\Gamma_h\).

#### Required finite-quotient certificate

Cocompactness is a mathematical finiteness reduction, not an implemented
enumeration algorithm.  A complete finite-quotient computation would have to
supply, with exact bounds:

1. a compact covering set or finite-sided arithmetic fundamental polyhedron
   for \(\Gamma_h\), including side-pairing matrices;
2. an exhaustive list of the five-shell mirrors meeting that set;
3. the induced finite cell quotient and all boundary identifications;
4. the action groupoid of lifted cells, from which each \(\Gamma_C\) is
   recovered; and
5. for every candidate cell, its complete facet orbits and the integrality of
   every facet reflection.

In a labelled quotient complex, products of side-pairing labels around loops
based at a cell, together with its recorded finite isotropy, generate
\(\Gamma_C\).  Compactness requires a proof that this matrix group is finite;
the finiteness of the quotient complex alone is insufficient.  In the even
case one may instead certify the equivalent finite-index statement
\([\Gamma_h:W_2]<\infty\).

None of these fundamental-domain, exhaustive-intersection, or stabilizer
computations has yet been completed for \(H_{11}\) or \(H_{19}\).  Exhausting
the finite quotient would exclude all Weyl-conjugate markings only after
every cell orbit has received one of these exact certificates.

## 6. Hyperkähler extension

For a projective irreducible holomorphic symplectic variety, the
Beauville--Bogomolov form restricts to signature \((1,\rho-1)\) on the divisor
lattice.  MBM classes cut out the Kähler cone, and prime exceptional divisors
can induce integral monodromy reflections.  This supplies a higher-dimensional
version of the same factory.

The caveat is important: an MBM wall is not automatically a Coxeter-reflection
wall.  The certificate must check reflection integrality and the Coxeter angle
condition.  Hyperkähler varieties are therefore a second-stage extension after
the K3/Enriques implementation, not a reason to relax the proof standard.

### 6.1 High-dimensional noncompact calibration

The extension is not merely hypothetical at the Picard-lattice level.
Alexeev's reflective 2-elementary examples give finite-volume noncompact full
root chambers in:

- \(\mathbb H^{14}\), from the rank-15 Picard lattice `(15,9,1)` of a
  `Hilb^2(K3)` fourfold, with 98 walls;
- \(\mathbb H^{19}\), from the rank-20 Picard lattice `(20,2,1)` of a
  complex K3 surface, with 50 walls; and
- \(\mathbb H^{21}\), from the rank-22 lattice `(22,2,0)` of a
  supersingular K3 surface in characteristic two, with 210 walls.

All three matrices, their exact ranks and inertias, and the geometric caveat
are audited in
[results/higher_dimensional_algebraic_reflection_examples.md](../results/higher_dimensional_algebraic_reflection_examples.md).
These are existing noncompact examples, not a new compact frontier claim.

### 6.2 The first golden OG10 carrier

The rank boundary is now crossed by an explicit construction.  Over
`O=Z[phi]`, take

\[
H_{OG}=\begin{pmatrix}-2\phi&1\\1&2\end{pmatrix}
       \oplus E_8\oplus\langle2\rangle.
\]

It has signatures `(10,1)` and `(11,0)`.  With trace scalar
`-phi/sqrt(5)`, its rank-22 trace lattice is

\[
N=U\oplus T_{11}(-1)\oplus2E_8(-1)\oplus2A_1(-1),
\qquad \det N=-44.
\]

An explicit order-44 discriminant gluing embeds `N` primitively in

\[
\Lambda_{OG10}=3U\oplus2E_8(-1)\oplus A_2(-1)
\]

with positive complement `<2>+<66>`.  Period surjectivity and the positive
class in `U` therefore produce a projective OG10-type variety with
`NS(X)=N`.  Its distinguished golden eigenspace has signature `(1,10)`.
The full calculation and hostile audit are recorded in
[results/og10_golden_h10_carrier.md](../results/og10_golden_h10_carrier.md).

There is also an automatic numerical wall lift.  If `H_OG(r,r)=2`, the
underlying divisors `D=r` and `tau D` have squares `-2` and `-4`, ambient
divisibility one, and are both OG10 wall divisors.  Self-adjointness gives

\[
(\tau D,x)=\phi(D,x)\qquad(x\in V_\phi),
\]

so their ambient hyperplanes restrict to the same norm-two field mirror.
The remaining geometric gates can in fact be closed.  In the explicit
ambient gluing, the mod-`3` radical of the OG10 Beauville lattice misses
`N/3N`; hence an `N`-class of ambient divisibility three is divisible by
three in `N`.  The `(-6,3)` and `(-24,3)` wall types are therefore absent.
The complete visible Kähler-wall arrangement reduces to the three field
norms

\[
2,\qquad4,\qquad6-2\phi.
\]

Markman's chamber theory and Mongardi--Onorati's OG10 classification realize
every component of its complement as a same-period marked Kähler slice.

The movable cone is even cleaner.  Its only stably prime-exceptional walls are the
square-`-2` classes, whose two visible shells form one unit-square orbit.
After re-marking,

\[
\operatorname{Mov}^0(X)\cap V_\phi^+
\]

is any prescribed chamber of the complete norm-two field-reflection
arrangement.  If `D=r` and `E=tau^{-1}D`, then

\[
D^2=E^2=-2,\qquad(D,E)=0,\qquad (R_D R_E)|_N=s_r.
\]

Thus every golden reflection is the restriction of a product of two
orthogonal divisorial OG10 monodromy reflections.  This establishes level 3
for the movable Coxeter cone.  Finite-sidedness, finite volume, and compactness
are logically independent of that realization; for this particular carrier,
the latter two are decided negatively below.  See
[the complete wall and movable-cone theorem](../results/og10_golden_complete_wall_chambers.md).

There is now an exact dual decision theorem.  The arithmetic centralizer of
`tau` is a uniform lattice because the conjugate field form is positive
definite, and

\[
\operatorname{vol}(C)<\infty
\iff\overline C\text{ is compact}
\iff|\operatorname{Bir}(X,\tau)|<\infty.
\]

The measure is intrinsic to divisor algebra: the OG10 Fujiki relation gives
`vol_X(D)=c_X B(D,D)^5`, and the Hessian of `-(1/10) log vol_X` on the unit
slice is the hyperbolic metric.  Moreover, a rational codimension-10 cycle on
`X x X`, built from divisors and their complete intersections, acts as `tau`
on `NS(X)_Q` and zero on the period plane.  Thus the enhancement itself is a
cohomological algebraic correspondence, while the enhanced birational group
decides total chamber volume.  See
[the OG10 enhanced-covolume theorem](../results/og10_enhanced_covolume_theorem.md).

For this first carrier the dual theorem is now evaluated, not merely
available.  The complete anti-symplectic fixed chamber \(P_9\) has Bugaenko's
compact chamber \(P_8\) as a complete right-angled facet.  The nontrivial
diagram reversal of \(P_8\) extends to an integral isometry \(A\) of the fixed
lattice and preserves the selected transverse side.  Local finiteness and the
uniqueness of the chamber adjacent to the relative interior of a complete
facet therefore give \(A(P_9)=P_9\).  With the independently certified
Pell-fold chamber involution \(g\),

\[
                     \operatorname{tr}(Ag)=8+8\phi>10.
\]

Hence \(Ag\) has infinite order.  The chamber stabilizer is infinite, the
norm-two Weyl group has infinite index in its uniform arithmetic group, and

\[
 \operatorname{vol}(P_9)=\operatorname{vol}(C)=\infty,
 \qquad |\operatorname{Bir}(X,\tau)|=\infty.
\]

The explicit fixed-lattice symmetries are not asserted to lift to the archived
H10 model; the final birational statement follows from the covolume
equivalence.  See the
[infinite-stabilizer theorem](../results/og10_golden_h9_infinite_stabilizer.md).

This yields a reusable **complete-facet symmetry obstruction** for future
inverse design.  If a symmetry of a complete facet extends integrally to the
ambient divisor lattice and preserves the chosen side, it is a symmetry of the
whole chamber.  If its product with any other chamber symmetry has infinite
order—for example, because the absolute trace exceeds the Lorentzian rank—then
the chamber stabilizer is infinite.  In the uniform arithmetic setting this
forces infinite reflection index and infinite hyperbolic volume.  Candidate
surfaces should pass this finite symmetry screen before expensive wall
enumeration.

The first inverse-designed neighbour shows that this is not the only useful
screen.  Coupling the last positive line produces

\[
M_{20}=
\begin{pmatrix}-2\phi&1&0\\1&2&1\\0&1&2\end{pmatrix}\oplus E_8.
\]

Its trace lattice has determinant \(-20\) and an explicit primitive embedding
in the actual OG10 Beauville lattice with complement
\(\langle2\rangle\oplus\langle30\rangle\).  Nevertheless, for the natural
fixed \(L_9\), every ambient norm-two mirror which meets the fixed hyperbolic
space projects to an existing fixed norm-two mirror.  If \(y\) is the
transverse coordinate, exact bounds at both real embeddings force

\[
                 y\in\{0,\pm1,\pm(\phi-1)\}.
\]

For \(y^2=1\), \(2\phi p\) is an integral fixed root of norm two; for
\(y^2=2-\phi\), \(2(\phi-1)p\) is.  Hence the projected arrangement is exactly
the already infinite-volume \(L_9\) arrangement.  Uniform arithmeticity then
rules out finite volume for every norm-two \(\mathbb H^{10}\) chamber of
\(M_{20}\).  See the
[A2-coupled projection theorem](../results/og10_golden_a2_coupled_candidate.md).

This gives a second reusable **projected-root obstruction**: before running
Vinberg on a new coupled carrier, enumerate the finitely bounded transverse
coordinates at every real embedding.  A viable design must produce a positive
projected mirror not proportional to an old fixed root and then prove that it
actually cuts the selected chamber.  Merely creating a finite Coxeter edge
between a root and its involutive image is not enough.

### 6.3 Wall-first covector compiler

The next inverse-design step makes that rule constructive.  Fix

\[
F=\langle-2\phi\rangle\oplus E_8\oplus\langle2\rangle,
\qquad
H(\ell)=\begin{pmatrix}F&\ell\\ \ell^t&2\end{pmatrix},
\]

and put \(p=F^{-1}\ell\).  The adjoined basis vector \(e\) is automatically a
norm-two root whose restricted wall is \(p^\perp\).  The map

\[
J|_F=1,\qquad J(e)=2p-e
\]

is an integral isometric involution exactly when \(2p\in F\).  The Schur
complement gives the complete two-place signature gate

\[
0<p^2<2,\qquad (p')^2<2.
\]

Thus one first prescribes a covector whose wall separates the archived
interior point from the exact hyperideal ray \(v_{61}\), and only then asks
for trace-lattice gluing and effectivity.  In coordinates
\(\ell=(A,E_8w,C)\), integrality is automatic and the signature test becomes

\[
q=2p^2=-A^2/\phi+C^2+2w^2,
\qquad 0<q<4,\quad q'<4.
\]

The first exact finite boxes are negative:

| declared family | exact admissible wall tests | walls cutting \(v_{61}\) |
|---|---:|---:|
| \(w=0\), coefficient box 30 | 66 | 0 |
| \(w=\phi^n v\), \(1\le n\le8\), \(v^2\le8\), box 30 | 2,706,240 | 0 |
| \(w=a+\phi b\), \(b^2=2\), box 30 | 15,360 (from 13,920 eligible \(w\)'s) | 0 |

Every comparison and the mixed \(E_8\) sphere enumeration are exact, and the
verifier replays every wall record against the archived interior and ray.
This is a finite-box no-go for the simplest one-root carriers, not a global
obstruction to algebro-geometric inverse design.  Two further mixed shells
have since been promoted exactly.  Through \(b^2\le6\), 569,760 eligible
mixed directions produce 1,472,640 wall tests and no separator.  A
structurally independent Weyl replay replaces the 242,965,920 raw
\((a,b)\) pairs by 79,923 representative tests and compares all eligible,
coefficient, and wall records literally.  The zero, scalar, and mixed family
counts sum to 4,178,946 evaluations but overlap in 209,280 scalar/mixed
records; the exact geometric union has 3,969,666 distinct evaluations.  See
the [family-overlap audit](../results/og10_golden_projected_wall_family_overlap.md),
the
[projected-wall box theorem](../results/og10_golden_projected_wall_box_no_go.md)
and the successive
[mixed-shell-two](../results/og10_golden_projected_wall_mixed_shell2_no_go.md)
and
[mixed-shell-three](../results/og10_golden_projected_wall_mixed_shell3_no_go.md)
theorems.

The zero-\(E_8\) subfamily is now global rather than bounded.  Put
\(X=-A/\phi\).  If \(\ell=(A,0,C)\) separated the interior and escape ray,
its two evaluations would force

\[
13<\frac C X<\frac{\phi(601+199\phi)}2.
\]

Together with \(0<q,q'<4\), this gives \(|X|<1/6\) and
\(|X'|<8/3\).  Minkowski reconstruction reduces every possible separator to
189 pairs, of which 26 pass the signature gate and none cuts.  A direct
independent 525-pair covector superset agrees.  See the
[global zero-E8 theorem](../results/og10_golden_projected_wall_zero_e8_global_no_go.md).

There is now a global reason not to respond by adjoining a second norm-two
line.  A full two-line field extension of \(F\) has trace rank 24 and trace
signature \((1,23)\), which cannot occur as an OG10 Picard lattice.  Even in a
larger abstract carrier, if \(F\) is saturated and an integral involution
fixes it pointwise with no further fixed directions, then for every norm-two
root \(r\)

\[
p=\frac{r+Jr}{2}\in F_K,\qquad 2p\in F.
\]

The root therefore generates a one-line subcarrier and satisfies exactly the
same two Schur inequalities.  Cross-coupling two such walls cannot jointly
remove a ray unless one wall removes it individually.  This reduction is
conditional on saturation and on the stated involution; it does not cover a
changed fixed lattice or a noncentral action.

The first genuinely different OG10 reflective type is a primitive class of
square \(-6\) and ambient divisibility three.  Writing

\[
H_6(k)=\begin{pmatrix}F&3k\\3k^t&6\end{pmatrix},\qquad p_0=F^{-1}k,
\]

gives the exact gate \(0<2p_0^2<4/3\) at the identity place and
\(2(p_0')^2<4/3\) at the conjugate place.  Exact replay of the first zero,
scalar-unit, and mixed boxes gives respectively 18, 696,000, and 5,760 wall
tests, all nonseparating.  A selected trace lift must still be placed on the
unique mod-three radical line by an explicit primitive OG10 gluing before it
can be called a divisibility-three MBM wall.  See the
[rank-two reduction and norm-six compiler](../results/og10_golden_rank_two_reduction_no_go.md).

### 6.4 Three levels that must not be conflated

For the geometric factory, keep the following implications separate:

1. `L` is an abstract reflective hyperbolic lattice;
2. `L` occurs as the Picard or Neron--Severi lattice of a variety; and
3. the finite root chamber is actually its nef, movable, or Kahler cone.

Alexeev's full reflection group uses both short `(-2)` and long `(-4)` roots.
On a K3 surface only the effective `(-2)` roots automatically cut out the nef
cone.  On a hyperkähler variety, Markman's prime-exceptional-divisor theorem
connects certain negative classes to monodromy reflections, but every proposed
wall still needs an MBM/effectivity certificate.  Thus algebraic realization
of the lattice establishes level 2; it does not automatically establish
level 3.

For compactness, the whole rational Picard cone is the wrong ambient space in
rank at least five: rational isotropic rays force cusps.  The additional object
to design is a self-adjoint order action whose distinguished totally-real
eigenspace is anisotropic.  The construction above shows that OG10 type is
not only a possible carrier: it realizes this quadratic eigenslice in
\(\mathbb H^{10}\).  For this OG10 carrier the theorem above also reaches
level 3: the marked eigen-movable cone is the norm-two Coxeter chamber, and
the complete Kähler subdivisions are known arithmetically.  The subsequent
complete-facet symmetry argument proves that this chamber has infinite volume
and is noncompact.  Finite-sidedness remains open, but it is no longer a
compact-example gate for this carrier.  The next algebraic candidate should be
inverse-designed from a surviving H9 combinatorial/Gram pattern and screened
by both obstructions above, not selected by an unconstrained lattice census.
The A2-coupled discriminant-20 neighbour is the first completed test of this
workflow: it passes the OG10 embedding gate but fails the projected-root gate
and is therefore discarded before a Vinberg enumeration.

## 7. Relation to the Gram--slack realization space

The divisor factory is more rigid than a bare semialgebraic realization
space.  Nevertheless, its output can be exported to one.  For a fixed
combinatorial type, vertex vectors \(v_a\), wall roots \(r_i\), and slacks

\[
S_{ai}=-(v_a,r_i)
\]

give a positive point on the slack variety, while
\(G_{ij}=(r_i,r_j)\) gives the Coxeter Gram point.  Conversely, the Gram--slack
space is retained for the theorem track and for independently verifying the
polytope produced by a divisor cone.

In the \(d=9,n=14\) target, \(G\) has corank 4.  Introducing its kernel as a
point of \(\operatorname{Gr}(4,14)\) turns the rank equations into the bilinear
incidence equations \(GK=0\).  This is the primal--dual bridge: the same
semialgebraic stratum can yield an exact point (an example) or an emptiness
certificate (a theorem), while the divisor-lattice object supplies a rigid
source of distinguished points.

## 8. Evidence states and first gates

The geometric track uses five statuses:

- `geometric_compact_reflective`;
- `geometric_cofinite_noncompact`;
- `geometric_obstruction_proved`;
- `geometric_realizable_noncoxeter`; and
- `undecided`.

The first gates are:

1. **G0:** certify one known Enriques/K3 divisor-cone Coxeter chamber;
2. **G1:** prove and implement the trace-lattice parity, discriminant form, and
   primitive K3-embedding tests;
3. **G2 (complete):** Bugaenko's rank-9
   \(\mathbb Z[\phi]\)-lattice admits the required K3 trace realization;
4. **G3 (complete):** its 11 walls form a Weyl-conjugate K3 eigen-nef
   section, with all restricted K3 walls exactly enumerated; and
5. **G4:** enumerate only the rank-10 forms passing G1 and search their
   Picard-rank-20 K3 realizations.

No broad rank-10 Vinberg sweep is authorized as the primary example search
before G1 determines how restrictive the geometric embedding conditions are.

## 9. Primary references to audit

- Discrete reflection chambers are Coxeter chambers: Vinberg 1971,
  Theorems 1--2 [Vinberg1971].
- K3 positive, ample, and nef cones and the Weyl-group reduction
  [Sterk1985, PiatetskiShapiroShafarevich1971].
- Cox-ring finite generation, finite automorphisms, and the full high-rank
  finite-nef list [ArtebaniHausenLaface2010, Kondo1989, Roulleau2022].
- Arithmetic lattices and the anisotropic compactness criterion:
  Borel--Harish-Chandra 1962, especially Theorem 11.8
  [BorelHarishChandra1962].
- Enriques surfaces with finite automorphism group and explicit curve diagrams
  [Kondo1986, BrandhorstShimada2022].
- MBM walls and the hyperkähler cone conjecture
  [AmerikVerbitsky2015, AmerikVerbitsky2017].
- Slack realization spaces [GouveiaEtAl2019].

Every cited implication must receive a theorem and page locator in the
frontier ledger before it is used to prune a search.
