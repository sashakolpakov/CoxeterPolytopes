# Polytope-to-surface moduli principle: exact theorem blueprint

**Status (2026-08-08):** hostile-audited synthesis.  The general principle
below is a theorem schema, not yet a theorem for all arithmetic Coxeter
polytopes.  Its weak rational, complete-root K3, fixed- and variable-
decoration Looijenga, and canonical determinant-line automorphic-arrangement
instances are proved at their stated strengths.  The one-boundary Coble
instance is a finite classification/recognition theorem, not a broad
universality theorem.

## The precise answer

The statement one can rigorously seek is not

```text
an arithmetic Coxeter polytope P canonically determines one surface X.
```

Periods, primitive-embedding orbits, integral wall normalizations,
overlattices, Weyl branches, and boundary decorations make such a surface
neither unique nor functorial.  The correct statement is

```text
a marked arithmetic Coxeter datum P canonically determines a carrier
moduli problem M_S(P), relative to a specified surface theory S;

P has an S-carrier of the required strength
    <=>
M_S(P) has a point passing the exact surface and cone gates.
```

The moduli problem is allowed to be empty.  Calling its formal definition an
``associated surface'' would hide the main existence question; an associated
surface is obtained only after choosing a geometric point.

Nor should `M_S(P)` always be called a finite-type algebraic stack.  In the
K3 case the exact-Picard locus is a type-IV period domain with a countable
Noether--Lefschetz arrangement removed; in the Looijenga case it is a period
torus with root hypertori removed and a nonseparated chamber lift.  The safe
universal term is **marked moduli functor/locus**.  Stack or coarse-space
structures should be asserted separately in each carrier theory.

## Input: more than an abstract Coxeter diagram

The natural input is a labelled and oriented arithmetic realization

```text
P^ar = (F,M,h,sigma_0; r_1,...,r_s; chosen positive sides),
```

not merely an unlabelled real polytope or Coxeter diagram.  A wall hyperplane
does not select a primitive integral normal.  Different positive rational
rescalings can give different crystallographic root data, and a scaled dual
does not exhaust them automatically.  Consequently an unmarked abstract
polytope has a canonical carrier problem only after taking the groupoid over
**all** of its admissible arithmetic realizations.  That groupoid need not be
known to be finite or effectively enumerable.

Fix a carrier theory `S`.  Let `A_S(P)` be the groupoid of surface-free lift
packages derived from `P`: integral root normalizations, reflection-stable
saturations, discriminant gluing, characteristic classes, boundary cycles,
or other finite decorations required by `S`.  Define

```text
M_S(P) = [ disjoint_union_(a in A_S(P)) M_S(a) / package isomorphisms ].
```

This groupoid-valued construction is canonical relative to the marked input
and the stated carrier theory.  A chosen list of representatives is only a
presentation of it.  Every claim for a bare diagram must therefore prove
that its compiler for `A_S(P)` is finite and exhaustive; otherwise an empty
computed union is only a package-relative obstruction.

## Exact gates

For a fixed package `a`, the useful theorem has six independent gates.

1. **Hodge/lattice gate `H`.**  The surface lattice has one positive
   direction and the marked wall classes embed with the required saturation,
   canonical class, discriminant form, and orientation.
2. **Moduli gate `M`.**  The corresponding marked period, deformation, or
   birational moduli locus is nonempty.  This is where the softer algebraic
   geometry constructs surfaces.
3. **Effectivity gate `E`.**  The requested signed wall classes are represented
   by the stipulated irreducible curves or divisors.  Algebraicity of a class
   alone is not enough.
4. **Cone-completeness gate `C`.**  The displayed inequalities equal the
   complete nef cone, Mori chamber, or certified invariant section.  Every
   exceptional, boundary, root, or MBM wall has been exhausted; no hidden
   wall cuts the interior.
5. **Marked-recovery gate `R`.**  Extraction returns the original labelled
   integral walls, sides, Gram matrix, and Coxeter angles, rather than only an
   isometric ambient lattice or an unlabelled combinatorial type.
6. **Hyperbolic gate `V`.**  Finite volume is proved after extraction, from a
   complete rational-polyhedral nef cone in Picard rank at least three, a
   complete no-hyperideal ray certificate, or literal Vinberg endpoint and
   full parabolic-neighbourhood tests.  It is not inserted into the definition
   of the surface locus.

Thus the non-tautological target for an exhaustively compiled surface-free
class `C` is

```text
P has a strong S-carrier
  <=>
there exist a in A_S(P) and x in M_S(a)
for which H, M, E, C, and R hold;

the extracted chamber is a finite-volume Coxeter polytope
  <=>
V also holds.
```

The surface conditions may be classical--exact Neron--Severi lattice,
irreducible negative curves, rational-polyhedral nef cone, finite Cox ring,
isotropic fibre classes--but their equivalence with polytope properties is
carrier-specific.  None of these statements requires a facet reflection to
be induced by a surface automorphism.  Automorphism finiteness may occur as
an equivalent algebraic test for cone polyhedrality; it is not an
automorphism-realization assumption.

## What is genuinely universal

For an admissible totally-real quadratic-form arithmetic datum, restriction
of scalars has Hodge-compatible signature exactly when there is one
Lorentzian real place and all other real places are negative definite.  A
common square multiple of that trace lattice embeds in the divisor lattice
of a blowup of `P2`; an algebraic correspondence cuts out the distinguished
eigenspace and integral divisor halfspaces recover `P` up to a common metric
scale.  This proves a nonempty weak rational-surface carrier for every such
`P`, with no rank ceiling.  See
[`rational_surface_eigencarrier_universality.md`](rational_surface_eigencarrier_universality.md).

This universal result supplies only `H`, weak `M`, and labelled linear
recovery.  Its divisor classes need not be effective, irreducible, extremal,
or complete nef walls.  It therefore cannot turn surface classification into
a Coxeter nonexistence theorem.  There is currently no theorem that
arithmeticity alone forces a strong K3, Looijenga, Coble, or Mori-dream
surface carrier.

## Exact scoped theorem I: K3 surfaces

Fix an oriented even integral realization `a` with lattice `N` of signature
`(1,n)`.  The weak marked K3 moduli problem is nonempty exactly when

```text
N embeds primitively in Lambda_K3 = U^3 + E8(-1)^2.
```

Equivalently there is an even complement of signature `(2,20-rank(N))` with
discriminant form `-q_N`.  This is an exact finite local discriminant-form
decision; a stable-range length inequality is only a sufficient shortcut.
The exact-Picard locus has dimension `20-rank(N)`.  See
[`k3_weak_algebraic_carrier_theorem_audit.md`](k3_weak_algebraic_carrier_theorem_audit.md).

For a proposed cone `C_P`, the strong locus is nonempty exactly when, in
addition,

```text
every irredundant facet normal is a primitive (-2)-root, and
C_P is one complete chamber of the hyperplanes of all (-2)-roots of N.
```

Period surjectivity supplies the surface and the K3 Weyl theorem moves its
nef cone to that chamber.  Conversely the complete nef cone of an exact-Picard
K3 has precisely this form.  This is the clean proved necessary-and-sufficient
`P -> surface -> P` theorem, but only for the complete primitive `(-2)`-root
class.  See
[`k3_strong_inverse_carrier_transfer_theorem.md`](k3_strong_inverse_carrier_transfer_theorem.md).

Inside this class, surface algebra gives the desired intrinsic dictionary:

```text
finite hyperbolic volume
  <=> finite index of the full (-2)-reflection group
  <=> rational-polyhedral Nef/Eff
  <=> finitely generated Cox ring
  <=> finite Aut(X),

compactness
  <=> no rational isotropic class
  <=> no genus-one fibration.
```

The automorphism clause is optional: Cox-ring finite generation or nef-cone
polyhedrality already states the surface-side condition.  No root reflection
is asserted to come from `Aut(X)`.  The dictionary and its scope are frozen in
[`k3_nef_covolume_dictionary.md`](k3_nef_covolume_dictionary.md).

The remaining algorithmic caveat is important.  Primitive embedding is a
finite decision, but ``no hidden `(-2)` root'' is finite only after a root
height bound, reflective classification, or another stopping theorem is
provided.

## Exact scoped theorem II: decorated Looijenga pairs

For a fixed rational-surface decoration extracted from `P`, let `D` be the
anticanonical cycle, `Lambda=D^perp`, `R0` the intended nodal root lattice,
and `Phi` the deformation roots computed from the generic Looijenga cone
`C_D++`.  The exact-root period locus is

```text
T_R0^o = { phi in Hom(Lambda,C*) :
           phi|R0=1 and phi(beta)!=1 for beta in Phi\R0 }.
```

Writing `Q=Lambda/R0`, this locus is nonempty exactly when `Q` has a
character nontrivial on every forbidden-root image.  A branch of the
nonseparated lifted period space then makes the intended simple roots
effective.  The final independent gate is

```text
C_P = C_D++ intersect the intended simple-root halfspaces.
```

For that fixed decoration these character, branch, and complete-exceptional-
cone conditions are necessary and sufficient for a full-nef Looijenga
carrier.  The arithmetic `H4` calibration has carrier locus `C*\{1}` and an
eight-wall full nef cone with five finite and four ideal vertices.  See
[`looijenga_h4_inverse_carrier.md`](looijenga_h4_inverse_carrier.md).

The supplied-decoration restriction is now removed in the bounded class
`C_triangle(4,8)`.  Starting with only a primitive full-span wall system in
`I_(1,4)` having at most eight square-minus-one/two walls, a terminating
compiler recovers the canonical class by adjunction and exhausts all 165
anticanonical triangles, blowdown models, period/root branches, and
`C_D++` inequalities.  Its carrier locus consists of exactly four weighted
Gram types `T6a,T6b,T7,T8`.  This is a total triangular-carrier detector, not
universal carrier existence; rejected inputs may still be abstract Coxeter
polytopes.  See
[`triangular_looijenga_h4_compiler.md`](triangular_looijenga_h4_compiler.md).

## Exact scoped theorem III: complex Coble surfaces

There is a clean surface-side necessary-and-sufficient criterion for the
whole complex Coble class.  If `S` is a complex Coble surface with boundary
`B=B_1+...+B_n`, then `(S,B/2)` is a klt Calabi--Yau pair and Kondo's action
on the Coble--Mukai lattice has finite kernel.  Totaro's theorem therefore
gives

```text
Aut(S) finite
  <=> Nef(S) rational polyhedral
  <=> Cox(S) finitely generated.
```

Under these conditions the projectivized **full** nef cone is automatically
a finite-volume Coxeter polytope in `H^(9+n)`.  Indeed, every negative curve
is a boundary curve `B_i^2=-4`, a `(-1)`-curve `E`, or a `(-2)`-curve `R`.
Adjunction and nonnegative intersections between distinct curves leave only
Coxeter angles `pi/2`, `pi/3`, `pi/4`, parallel walls, and ultraparallel
walls.  The scaled facet normals

```text
B_i, 2E, R
```

have integral Cartan pairings and generate a full-rank reflection-stable
lattice, so the polytope is arithmetic even when an individual `B_i`
reflection does not preserve the unimodular Neron--Severi lattice.  Complete
rational-polyhedral nef equality gives finite volume by the convex--Hodge
theorem.  Kondo's effective-root theorem supplies a special genus-one
fibration for every finite-automorphism type, so its nef fibre class is
isotropic and the polytope is noncompact.

Conversely, if the full nef cone is a finite-sided Coble Coxeter polytope,
then it is rational polyhedral; Totaro makes the image of `Aut(S)` on
`N^1(S)` finite, and Kondo's finite-kernel theorem makes `Aut(S)` finite.
Thus the displayed surface-algebra criterion is also necessary and
sufficient for the full-nef Coble construction.  Kondo's complex
classification has exactly three carrier surfaces: Type I with `n=1`, Type
I with `n=2`, and Type II with `n=1`.  This classifies the **surface carrier
types**, but it becomes a surface-free Gram recognition theorem only after
the complete wall matrices have been extracted.

The current Coble result is a finite classification/recognition theorem, not
a general `P -> Coble surface` theorem.  For a complex Coble surface with one
boundary component and finite automorphism group, Kondo's classification
leaves Types I and II.  In full Neron--Severi rank eleven their complete nef
cones give two noncompact finite-volume arithmetic Coxeter polytopes in
`H10`:

```text
Type I:  22 facets, 315 finite rays,  38 ideal rays;
Type II: 30 facets, 1333 finite rays, 207 ideal rays.
```

The pair `(S,B/2)` is klt Calabi--Yau.  Totaro's complex-surface theorem and
finite `Aut(S)` imply a complete rational-polyhedral nef cone and finite Cox
ring; Kondo exhausts negative curves as the boundary, `(-1)`-curves, and
`(-2)`-curves.  In the one-boundary case the primitive boundary hyperplane is
normal to `-K`, so the full wall system has norm `-1` and `-2` integral
reflections.  A Kondo genus-one fibration supplies an isotropic nef ray and
hence noncompactness.  Exact cone/Vinberg enumeration supplies the displayed
counts.

Conversely every complex finite-automorphism one-boundary Coble full-nef
chamber is one of these two marked integral types.  There is now an
independent polytope-first recognizer rather than merely a Gram lookup.  From
an oriented odd-unimodular `H10` chamber it enumerates characteristic
norm-minus-one facets `b`, proves `b^perp=E10`, recovers the orthogonal
norm-minus-two root chamber, and requires the classes `E+b` from all remaining
norm-minus-one walls to equal the complete primitive-isotropic cusp set.  The
integral saturation is included in the test.

Applied to the pre-existing `D16-D5-H10-30` walls alone, this process uniquely
reconstructs Kondō Type II: one boundary class, twelve roots, seventeen cusp
lifts, and the complete 30-wall full-nef package.  A hostile implementation
imports neither the surface-first wall list nor the main Coble certifier.  Its
frozen local denominator contains six odd-unimodular `H10` chamber types, 24
integral presentations, and 154 characteristic candidates.  Only `P10/12`
and `D16-D5-H10-30` pass raw incidence and complete cusp regeneration; Kondō's
classification excludes `P10/12` only as a finite-automorphism one-boundary
Coble carrier, leaving Type II as the unique positive in that denominator.
The two possible `E10` index-two saturations are exchanged by the order-48
root-diagram automorphism group, so the integral package is unique up to
diagram automorphism.

This is a genuine but scoped `P -> surface` theorem for one carrier species.
It produces a new strong carrier of a known polytope, not a new abstract
polytope.  It does not prove that an arbitrary arithmetic `H10` polytope has a
Coble carrier, and the five negative rows are carrier obstructions rather
than polytope nonexistence theorems.  See the
[main Coble theorem](coble_one_boundary_full_ns_h10_theorem.md) and the
[independent polytope-first hostile audit](polytope_first_coble_h10_hostile_audit.md).

There is one safe extension and one unsafe extension:

- the complex Type-I surface with two boundary components now has a frozen
  complete full-`NS` 26-wall nef chamber in `H11`: exact double description is
  `508+42+0`, literal Vinberg has 3,641 edges and 106 affine audits with zero
  failures, and every arithmetic reflection is integral on
  `2 Pic~=U+E8(-2)+<-2>^2`.  This is a project-new exact positive realization
  with worldwide abstract priority open.  It remains strictly surface-first,
  not an inverse recognition theorem; see the
  [two-boundary theorem](coble_type_i_two_boundary_full_ns_h11.md);
- the `n=5` and `n=8` Kondo surfaces occur only in positive characteristic.
  Totaro's cited Cox-ring equivalence is explicitly over the complex
  numbers, and Kondo's Coble--Mukai root exhaustion does not automatically
  exhaust every full-`NS` `(-1)` wall.  Proposed `H14` and `H17` conclusions
  therefore require an independent negative-curve and nef-cone exhaustion
  before promotion.

## The completed bounded triangular class

The first completed variable bounded target is the triangular Looijenga class

```text
C_triangle(4,8): marked rational arithmetic H4 wall systems with
  at most eight facets,
  full wall span and odd-unimodular saturation I_(1,4),
  primitive wall squares in {-1,-2},
  the Coxeter intersection alphabet forced by those curve types, and
  no supplied surface, boundary decoration, or volume verdict.
```

The last phrase must be encoded without naming a surface: enumerate
characteristic vectors `k` of square five, decompositions
`-k=d1+d2+d3` with the triangle intersection matrix, and integral blowdown
bases certifying the boundary blowups.  These are finite lattice equations.
Because `k^perp` is negative definite, both deformation-root and numerical
exceptional-class searches terminate.  Thus every obligation can be made
finite:

1. enumerate all arithmetic Gram/root realizations and all decorations;
2. prove that the compiler is exhaustive up to marked lattice isometry;
3. apply the exact character criterion and enumerate chamber lifts;
4. exhaust exceptional and root walls and test complete `C_D++` equality;
5. recover the labelled Gram matrix and replay exact rays/Vinberg data; and
6. classify the finite discrete Gram image of the surviving period strata.

The completed compiler gives the genuine theorem

```text
P in C_triangle(4,8) has a full-nef Looijenga surface carrier
  <=>
one of its finitely and exhaustively derived character/cone loci is nonempty.
```

The census is exact: 165 triangles and 885 period/chamber branches reduce to
260 marked wall systems and exactly four Gram types, with facets `6,6,7,8`
and finite/ideal profiles `5+2,4+3,5+3,5+4`.  All four are finite-volume
noncompact full-nef `H4` carriers.  `T6a/T6b` are Tumarkin pyramids and `T7`
is Im Hof's doubly truncated orthoscheme; their new content is the AG carrier,
not the abstract polytope.  All four are canonical mod-two congruence
preimages in `W(I_(1,4))`, of indices `10,20,30,60`.  `T8` has no bounded-
corpus match, but worldwide priority remains open.  See the
[`compiler`](triangular_looijenga_h4_compiler.md) and
[`lineage/congruence theorem`](triangular_looijenga_h4_worldwide_lineage.md).

This is still not `U_total`: a rejected input has no triangular Looijenga
carrier, but may have another carrier or exist with none.  The next target is
therefore a jointly universal finite carrier stack, preferably compact, not
another single-species lookup table.

## Exact scoped theorem IV: determinant-line automorphic arrangements

Let a compact arithmetic `P` be presented by an even primitive totally-real
wall lattice of signature `(1,n)` at one place and definite conjugates.  Its
coordinate-free determinant line has the forced opposite sign pattern.
Adjoining it and enumerating reflection-stable discriminant gluings gives a
finite family of `(2,n)` orthogonal Shimura stacks, projective when the field
degree is greater than one.  A canonical real slice pulls each facet Heegner
divisor back to exactly its mirror.  Negative-definite intersections recover
elliptic candidates, oriented halfspace feasibility recovers actual faces,
and the marked dihedral inertia subgroup of order `2m_ij` recovers the
finite Coxeter label.  Parallel and ultraparallel data remain in the marked
real-slice Gram.  On compact `P`, hidden sheets in every fixed divisor shell
are finitely enumerable.

Bruinier's vector-valued Hilbert cusp pairing is necessary and sufficient
for a prescribed facet principal part to give a Whittaker/Borcherds divisor.
That is an exact iff for this automorphic carrier species, not automatic
correction for every reflective lattice.  A general complete-intersection
surface gives curves determined by `P`, but they encode the divisor
arrangement rather than automatically forming its full nef cone.

For Bugaenko's compact `H8/11` datum, the unique determinant package and
hidden-sheet gate are complete, but all eleven facets collapse to one
standard Heegner column with 309 additional lower-face contacts.  The exact
weight-`(5,5)` Poincare coefficient estimate is nonzero, so the pure norm-two
Whittaker principal part is obstructed.  The complete globally minimal
norm-`1/80` compensator shell is nonempty but disjoint from the real slice,
and its target-plus-shell `2x2` Poincare matrix is invertible, so this first
support is obstructed as well.  The determinant safe-span theorem now proves
that some larger finite slice-disjoint support cancels the obstruction; over
a real/integral q-structure the correction is nonnegative/effective
(conditional here on Bruinier's Hilbert q-integrality statement).  Thus this
is a first-shell obstruction, not a polytope obstruction or a global carrier
failure.

An actual dimension bound requires an independently defined compact class
with a uniform bound on the automatic carrier's height, level, weight, and
degree, or a jointly exhaustive alternative carrier disjunction, plus a
uniform high-weight obstruction estimate on every allowed support.
A facet bound additionally requires a sheet/degree inequality preventing
many real walls from collapsing to one global Heegner orbit.  These are
explicit missing lemmas, not assumptions.  See the
[`determinant-line transfer theorem`](polytope_to_automorphic_divisor_transfer_theorem.md),
the [`Bugaenko Whittaker gate`](bugaenko_h8_whittaker_gate.md), and the
[`first compensator theorem`](bugaenko_h8_whittaker_first_compensator.md).

## What would count as the broad theorem

A broad theorem still missing would exhibit a nontrivial variable class `C`,
defined without surfaces and larger than a finite list of known Gram
matrices, for which:

```text
A_S(P) is finite and exhaustive for every P in C;
H and M always construct the associated surface moduli problem;
E and C are decidable by intrinsic surface invariants;
R recovers P exactly; and
V follows from the complete surface cone.
```

Only after that forward universality theorem may emptiness or classification
of the surface loci be transferred into an abstract Coxeter nonexistence or
finiteness theorem.  Until then, K3, Looijenga, or Coble failure is a scoped
carrier obstruction, not a theorem that the arithmetic polytope does not
exist.

## Supporting synthesis

The longer association ladder and decision interface are in
[`arithmetic_coxeter_to_surface_association_criterion.md`](arithmetic_coxeter_to_surface_association_criterion.md).
The general implication from a complete rational-polyhedral nef cone to
finite hyperbolic volume is proved in
[`projective_surface_nef_hyperbolic_volume_theorem.md`](projective_surface_nef_hyperbolic_volume_theorem.md).
