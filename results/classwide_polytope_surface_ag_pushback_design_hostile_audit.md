# Hostile design audit: a class-wide `POLYTOPE -> SURFACE -> AG -> theorem on P` route

**Date:** 2026-08-09  
**Verdict:** the closest genuinely non-tautological route is the **cusped
Coble--Vinberg incidence class**.  Its finite carrier compiler, exact
same-`P` theorem for every geometric branch, and independent AG
classification are already proved.  Two linked strict-L4 gates remain:
geometric nonemptiness on a natural surface-free, non-singleton
subdenominator, and a nonthin terminal pushback not already encoded by its
finite core grammar.  The analogous compact golden-K3 route is not one technical
lemma from completion: its missing essential-surjectivity statement already
contains the desired `H9` nonexistence theorem.

Finite-group/triality folds of simply-laced chambers are useful fishing
compilers, but they fail this audit as the main theorem route.  Their local
folding data neither force a complete K3 root chamber nor make the terminal
surface classification independent of the Coxeter classification.  Adding
root completeness to the definition repairs exact return only by putting the
hard answer into the input.

The compact and cusped tracks must be stated separately.  A strong Coble
full-nef chamber is always noncompact, whereas the completed direct K3-root
compact theorem stops in `H3`.  No presently justified single surface
species spans both regimes.

**Post-audit update.**  The two four-point cores are exact: their P-derived
rank-drop designs are realizable at `P1(F3)`, construct all root seeds, and
return the same `H14/60` and `H17/48` chambers.  The advertised extension
closure is false, however: `O(2,0)` is base-point-free with a generically
separable associated map but its general divisor is two fibres.  Correct
status is core-only `L3`, not `L3+`.  The seven-row follow-up found different
finite geometries instead.  Raw `H11/206` gives a Möbius-plane/Frobenius
carrier whose Coble identities hold iff `q=3`; raw `H10/78` gives the
characteristic-five `K5` Desargues/cut design and an explicit integral
ten-node sextic carrier.  Both return the exact original chambers.  A
normalized arbitrary-input recognizer for all five cores is now frozen; it
closes extraction but remains L3 because the grammar is rigid and terminal
pushback is thin.  The fifth explicit P-first calibration, Type-VII `H10/72`,
uses a tetrahedral `4+6` design, nodal quartics, and a boundary with four
nodes and six cusps.  The formal-L4 gate is now a nonzero extension class
with uniformly nonempty moduli and genuinely new AG pushback, not another
hardwired row.  See the [`CI4` hostile audit](coble_four_point_constructible_incidence_hostile_audit.md),
[`seven-row audit`](coble_remaining_seven_contraction_design_audit.md), and
[`generic recognizer`](coble_finite_geometry_raw_recognizer.md),
[`Desargues inverse`](coble_type_vi_char5_desargues_inverse.md), and
[`Type-VII inverse`](coble_type_vii_char5_tetrahedral_inverse.md).

## 1. The contract being audited

For a surface-free class `C` of marked arithmetic Coxeter data, a strict
inverse theorem must have five logically separate gates.

1. **`C`: independent denominator.**  Membership is computed from the
   lattice, oriented walls, Gram data, Coxeter relations, and finite
   combinatorial decorations derived from `P`.  It cannot say that a surface
   exists, that a period locus is nonempty, or that `P` is one of a known
   list of answers.
2. **`Car`: finite exhaustive carriers.**  From each `P`, enumerate finitely
   many integral normalizations, saturations, complements/gluings, boundary
   frames, contraction plans, period/action packages, or incidence strata.
   Every carrier of the declared species must occur on one branch, and at
   least one branch must be nonempty for every `P in C`.
3. **`Return`: exact marked same-`P`.**  A successful branch must recover the
   same integral wall normals, orientations, Gram matrix, and cone.  This
   requires both inclusions

   ```text
   every facet of P is an effective geometric wall;
   every effective root/negative curve/MBM wall is nonnegative on P.
   ```

   The second inclusion is the hidden-wall gate.  A lattice embedding or a
   list of visible curves is not enough.
4. **`AG`: independent algebraic-geometric input.**  A classification,
   Hodge bound, deformation theorem, Cox-ring theorem, or surface invariant
   must be applied only after the carrier and exact return have been proved.
   Reusing a Coxeter classification under surface terminology does not meet
   this requirement.
5. **`Pushback`: theorem about every a-priori `P`.**  The AG conclusion must
   return a new dimension bound, allowed-dimension set, finite
   classification, or nonexistence theorem on the original denominator—not
   merely on the subset already defined by carrier success.

The hyperbolic-volume gate remains separate.  A finite wall list can have
hyperideal ends and infinite volume.  It is enough either to run Vinberg's
complete elliptic-edge and full parabolic-neighbourhood tests, or, after
proving equality with a complete rational-polyhedral nef cone in Picard rank
at least three, to invoke the convex-Hodge finite-volume theorem.  Neither
route follows from “finite-sided” alone.  See the
[`nef-volume theorem`](projective_surface_nef_hyperbolic_volume_theorem.md).

## 2. The completed baseline, and why it is not the flagship

There is already a genuine class-wide inverse theorem for the compact class

```text
K3RootComp = even integral hyperbolic lattices N of rank at least three
             with a compact chamber of the complete (-2)-root arrangement.
```

This definition contains no surface or K3 embedding.  Compactness makes the
root reflection group cocompact; Godement plus Meyer gives
`rank(N)<=4`, hence `dim(P)<=3`.  Primitive K3 embedding is then automatic,
and K3 Weyl theory returns the same complete chamber with no hidden roots.
The Nikulin--Vinberg/Roulleau classification leaves exactly eight marked
types: six in `H2` and two in `H3`.  This theorem and its hostile audit are
in
[`compact_complete_k3_root_return_classification.md`](compact_complete_k3_root_return_classification.md)
and
[`compact_complete_k3_root_return_hostile_audit.md`](compact_complete_k3_root_return_hostile_audit.md).

The finite-volume extension `Root2FV` similarly has 118 marked complete
norm-two types, eight compact and 110 cusped, with no member above `H18`.
Every type has an exact K3 return.  See
[`global_complete_k3_root_automorphic_return.md`](global_complete_k3_root_automorphic_return.md).

These are important denominator theorems, but they do not supply the desired
new terminal AG leverage.  Root completeness makes the hidden-wall gate part
of the Coxeter/lattice input, and the dimension bounds and finite lists come
first from Godement--Meyer or the Nikulin--Vinberg Coxeter-lattice
classification.  K3 theory supplies exact realization and useful
equivalences, but not an independently stronger last arrow.  They are the
correct regression standards against which a new route must be measured.

## 3. Hostile rejection of the `D4`-triality fold as flagship

A natural surface-free candidate is the class of connected integral
`m=3,6` chambers for which, after contracting the `3`-edge components, the
`6`-edge quotient is loop-free and bipartite.  This condition is exactly
equivalent to assigning root squares `2` and `6`: equal squares across a
`3`-edge and unequal squares across a `6`-edge.  Every square-six wall then
unfolds to three mutually orthogonal square-two roots, with a cyclic or
`S3` triality action.

For a fixed `P`, the natural finite packages are:

1. the two complementary long/short colourings;
2. all integral unfolded square-two Grams and their saturated radicals;
3. all even reflection-stable overclosures in the finite discriminant
   group;
4. all primitive K3 complements and discriminant anti-isometry double
   cosets;
5. all extensions of the triality action through the glue and period;
6. all complete-root chambers containing the displayed unfolded walls; and
7. the invariant-slice comparison, including every projected hidden root.

The first five stages are genuinely finite.  The last two are the hard
problem.  Local foldability and primitive K3 embeddability do not force the
unfolded walls to be the simple walls of a complete root chamber.  Nor does
finite volume of `P` force finite volume, or even Lorentzian signature, for
the unfolded chamber.

The existing Bredon census is a direct counter-calibration.  Of 24 ADEG
diagrams, 21 admit 42 integral triality orientations and 31 pass a stable
K3-embedding test, but the precompiler certifies no same-`P` return by
itself.  More sharply:

- the `N11` orientation is exactly the `C3` fold of the complete
  `U+E8(-1)+D4(-1)` K3 nef chamber and gives an equation-level
  `P -> K3 -> same P` return;
- exhaustive comparison with a 55-stock local denominator gives eight
  further exact fixed-divisor-section returns: one `N2`, three `N3`, two
  `N5`, one `N6`, and one `N7` package.  The residual 30-package census is
  `P=8,R=15,O=7`, where `R/O` are only denominator-relative outcomes;
- the apparently analogous `N9` orientation is not the unique classified
  rank-12/13-wall K3 stock: its unfolded graph is a tree with a degree-four
  vertex and determinant `-12`, whereas the stock contains a triangle, has
  maximum degree three, and determinant `-3`.

See
[`bredon_adeg24_triality_precompiler_audit.md`](bredon_adeg24_triality_precompiler_audit.md)
and
[`bredon_n9_n11_unfolded_k3_stock_lineage.md`](bredon_n9_n11_unfolded_k3_stock_lineage.md),
as well as the
[`remaining-stock promotion theorem`](bredon_adeg24_remaining_k3_stock_promotions.md).

Thus the local P-only fold data do not determine the global nef completion.
Defining the class instead by saying that the unfolding is a complete K3
root chamber would make exact return true, but then the decisive condition
has been inserted at input.  The available finite-automorphism K3
classification is also essentially the complete-root Coxeter
classification, so it does not furnish an independent AG terminal bound.

**Audit verdict:** retain triality as a positive-example and lineage search
engine; reject it as the present class-wide theorem route.

## 4. Compact route: golden coherent K3 carriers

### 4.1 Independent class

Let `C_35^conn` consist of actual oriented crystallographic realizations over
`Q(sqrt(5))` whose `m=3,5` wall graph is connected.  No surface, K3
embedding, complete-root hypothesis, or carrier-success condition occurs in
the definition.  The actual arithmetic realization, rather than only an
angle diagram with dotted edges, is essential.

The all-scaling theorem proves a unique equal-root normalization and a finite
saturation compiler for every such realization.  Infinite unit twists are
vertex coboundaries.  See
[`golden_normal_scaling_universality.md`](golden_normal_scaling_universality.md).

### 4.2 Finite exhaustive packages

For a compact `P in C_35^conn`, the complete full-trace K3 compiler is:

1. normalize the unique square-two golden wall lattice `R`;
2. enumerate every even, field-stable, reflection-stable overclosure between
   `R` and `R^vee` (or between the trace lattice and its dual);
3. enumerate every primitive K3 complement and discriminant-glue double
   coset;
4. on the complement, enumerate rational self-adjoint square-root actions
   `B_T^2=5`, integral glue-compatible order actions, and Hodge-signature
   branches;
5. choose exact-Picard periods and Weyl framings; and
6. enumerate every trace root whose mirror meets compact `P`, rejecting a
   branch unless all projected root inequalities are redundant on `P` and
   every input facet is attained.

Every stage is finite for fixed compact `P`.  The last condition is the exact
same-`P` iff: primitive embedding plus visible roots is insufficient.
Compactness makes the hidden-root search finite, but does not make it
automatically pass.

### 4.3 Independent AG input and conditional pushback

A coherent full-`H^2` quadratic action is self-adjoint.  The divisor
eigencarrier has rational rank `2(n+1)`, and the transcendental distinguished
eigenspace must contain the two-dimensional K3 period plane.  Hence

\[
                  2(n+1)+4\le 22,
                  \qquad n\le 8.                         \tag{4.1}
\]

This is genuinely Hodge-theoretic.  Bugaenko's compact `H8` chamber is the
sharp coherent endpoint, and all twelve golden Lanner chambers have coherent
same-`P` returns.  See
[`golden_lanner_full_h2_coherent_returns.md`](golden_lanner_full_h2_coherent_returns.md).

The desired compact theorem would be:

> Every compact `P in C_35^conn` has at least one coherent full-`H^2` K3
> package passing the complete restricted-root gate.  Consequently
> `dim(P)<=8`, sharply.

That theorem is not proved.  At `n=9` the Picard trace lattice has rank 20,
the transcendental complement has `Q(sqrt(5))`-rank one, and each real
eigenspace is one-dimensional.  It cannot contain the period plane.  The
exact rank-20 audit finds packages for which rational and integral glue gates
pass, but every Hodge gate fails.  See
[`golden_rank20_full_h2_extension_gate.md`](golden_rank20_full_h2_extension_gate.md).

Therefore the missing universality lemma is already equivalent, at the top
boundary, to proving that no compact `H9` member exists.  It cannot be
presented as a routine consequence of primitive embedding.  Defining a
subclass by “admits a coherent package” would give a valid recognition
criterion but make (4.1) tautological as a polytope bound.

**Compact verdict:** the finite compiler and AG obstruction are exact, but
there is no small class-wide existence lemma presently in sight.  The honest
parallel plan is to search `H9` for an `NS`-only complete-root return while
separately seeking a genuinely P-forced reason for full-cohomology
coherence.  A positive `H9` survivor would refute coherent universality but
would itself be a valuable new example.

## 5. Cusped route: the Coble--Vinberg incidence class

### 5.1 Independent non-singleton class

`D_CobVin` is defined from an oriented integral Coxeter datum and a strict
integral interior vector.  The compiler enumerates all boundary de-scalings
to

\[
 N\simeq I_{1,9+b},\qquad B_i^2=-4,\qquad
 K=-\frac12\sum_i B_i\in N,
\]

requires `K` to be characteristic, and checks the adjunction-compatible
wall alphabet `(-4),(-2),(-1)`.  It then reconstructs the rank-ten
Coble--Mukai root lattice, applies Vinberg's finite-index criterion, exhausts
primitive cusp lifts, and reassembles the original oriented walls.  No
surface, characteristic, Kondō type, boundary count from the classification,
or incidence solution occurs in the definition.

The class is already non-singleton: the archived `H14/60` and `H17/48`
inputs are distinct members.  The full construction is in
[`coble_vinberg_incidence_inverse_transfer.md`](coble_vinberg_incidence_inverse_transfer.md).

### 5.2 Finite exhaustive carrier packages

For every boundary branch, the compiler enumerates all contraction sequences
by the numerical update

\[
 D'.F'=D.F+(D.E)(F.E),\qquad K'=K-E,
\]

until reaching `P2` or a Hirzebruch surface.  Reversing each sequence gives
a finite-type iterated-blowup incidence stack.  Its locally closed conditions
impose:

```text
smooth, geometrically irreducible, distinct wall curves;
all prescribed ordinary and infinitely-near incidences;
an ample marked interior class;
disjoint rational boundaries with sum B=-2K;
h^0(-K)=0 and h^0(-2K)=1.
```

The collection is finite because every contraction lowers Picard rank and
only the finitely many raw facet classes can occur.  It is exhaustive in the
needed direction: on any exact strong-Coble carrier, every curve contracted
by a minimal-model sequence is a negative curve and hence a nef facet, so
some full contraction sequence appears in the compiler.

The optional K3-cover compiler gives a second finite description: boundary
halves, every stable/split wall lift, all cross-pairings, even involution-
stable overclosures, discriminant glue, equivariant K3 embeddings, periods,
ramification frames, and projected roots.  This is valuable for complex
promotion and crystalline analysis but is not needed for exact return once a
rational-surface incidence branch exists.  See
[`coble_k3_cover_polytope_to_surface_bridge.md`](coble_k3_cover_polytope_to_surface_bridge.md).

### 5.3 Exact same-`P` iff, including hidden curves

For every algebraically closed field of characteristic different from two,
the proved criterion is

```text
P has an exact strong-Coble full-nef carrier
  <=>
one of its finitely compiled incidence stacks has a geometric point.
```

Every geometric point is a strong Coble surface.  Vinberg finite index plus
Kondō's effective-root theorem proves that the compiled Coble--Mukai roots
are all effective irreducible roots.  This eliminates hidden `(-2)` curves
and hidden `(-1)` curves meeting two boundary components.  Any remaining
`(-1)` curve produces a primitive nef isotropic class `f=E+B_i/2`; the
complete cusp-lift enumeration then eliminates it as well.  Thus all
negative curves are exhausted and

\[
                  \mathbf P(\operatorname{Nef}S)=P
\]

with the original integral orientation and marking.  This is a genuine
uniform exact-return theorem, not a lookup in one of the nine surface
templates.

Once full-nef equality is known, rational polyhedrality also proves finite
hyperbolic volume by the convex-Hodge theorem.  The reflection integrality,
Coxeter angles, and literal marked wall equality still come from the
P-derived arithmetic compiler.  Every resulting Coble chamber is
noncompact: pulling back a ruling from a rational minimal model supplies a
nef isotropic ray.

### 5.4 Independent AG input and the sharp conditional theorem

Exact finite-sided full-nef equality forces `Aut(S)` to be finite: the image
of `Aut(S)` permutes a finite spanning set of primitive facet normals and
Kondō's action on the Coble--Mukai lattice has finite kernel.  Kondō's
classification then gives

\[
 b\in\{1,2,5,8\},\qquad
 \dim P=9+b\in\{10,11,14,17\}.                         \tag{5.1}
\]

The `H14` and `H17` rows exist only in characteristic three.  In particular
there are no compiler-positive Coble-nef members in `H12`, `H13`, `H15`,
`H16`, or above `H17`.  This restriction is AG-native; it was not inserted
in `D_CobVin`.

All nine classified full-NS chambers have now been constructed exactly:

| type | characteristic | `b` | dimension | full-nef facets |
|---|---:|---:|---:|---:|
| I | any `!=2` | 1 | `H10` | 22 |
| II | any `!=2` | 1 | `H10` | 30 |
| I | any `!=2` | 2 | `H11` | 26 |
| V | 3 | 2 | `H11` | 36 |
| VI | 5 | 1 | `H10` | 78 |
| VII | 5 | 1 | `H10` | 72 |
| MI | 3 | 2 | `H11` | 206 |
| VI | 3 | 5 | `H14` | 60 |
| MII | 3 | 8 | `H17` | 48 |

See
[`coble_nine_class_polytope_surface_transfer.md`](coble_nine_class_polytope_surface_transfer.md).
These are positive full-nef hyperbolic examples; several complete surface
realizations are project-new, while worldwide abstract-polytope priority is
separate and often already known.

The strict target theorem is now extremely clear:

> **Target `U_Cob`.**  Every member of a natural P-only, non-singleton
> subdenominator of `D_CobVin`, defined before the Kondō list and not by
> incidence-stack nonemptiness, has at least one geometric incidence branch.

Once `U_Cob` is proved, the already completed exact-return theorem and
Kondō's classification immediately turn (5.1) into a theorem on every input
polytope in that subdenominator.  If `U_Cob` held for all of `D_CobVin`, it
would give the strongest version; the programme should not assume that this
maximal statement is true.

## 6. Smallest credible completion lemma (original proposal and correction)

The smallest non-tautological target is not “the incidence stack is
nonempty.”  That merely renames carrier existence.  It is an explicit
**free-incidence realization lemma** whose hypotheses are finite numerical
conditions derived from `P`.

For each P-derived contraction branch, form its dependency data on `P2` or
`F_e`.  Call the branch freely realizable when the following hold.

1. Every reverse blowup centre is either a free point on a previously
   constructed smooth curve or a prescribed transverse intersection of two
   previously constructed curves; the dependency graph is acyclic.
2. Every new curve class lies in a base-point-free system after the earlier
   jet conditions, with the required `H^1` vanishing and evaluation-map
   surjectivity verified from its numerical class.
3. The imposed ordinary and infinitely-near conditions have independent
   expected codimension, expressed by the ranks of explicit evaluation
   matrices on the minimal model.
4. The boundary classes satisfy a numerical exact-sequence/Hilbert-function
   certificate forcing `h^0(-K)=0` and `h^0(-2K)=1` on the resulting open
   stratum.

These predicates are computed from the intersection table and contraction
formulas supplied by `P`; they do not name a Kondō type or assume a surface
point.  Induction on the acyclic dependency order, followed by Bertini and
upper semicontinuity, should prove that every freely realizable branch has a
geometric point.  The existing hidden-curve theorem then gives exact
same-`P`, and Kondō supplies (5.1).

The route to a completed theorem is finite and falsifiable.

1. Run the all-contraction compiler on the nine frozen full-NS rows and
   canonicalize their P-derived dependency data without importing their
   surface labels.
2. Test whether the free-incidence predicate contains at least two distinct
   raw inputs.  If it does, this immediately gives a non-singleton theorem
   class.  If it does not, enlarge the lemma to a controlled determinantal
   or supersolvable-matroid predicate rather than adding answer templates.
3. Prove the realization lemma uniformly on that combinatorial predicate.
4. Freeze exact elimination only for the residual characteristic-dependent
   branches.  In particular, the `H14` conic/tangent-line branch and the
   `H17` `4x4` grid branch already show how characteristic three can be
   forced by P-derived incidence equations.
5. Apply the existing exact-return theorem verbatim, then Kondō's
   classification.  No new hidden-root argument or surface classification
   is needed.

This is materially smaller than the K3-cover alternative.  The latter still
requires, uniformly, equivariant primitive embeddings, complex or
crystalline periods, uniqueness of the ramification frame, and complete
projected-root equality.

The subsequent exact tests change the implementation order.  Conditions
2--3 above cannot require evaluation surjectivity: on `H14/60` the six
`3x3` line matrices have rank two, and on `H17/48` all twenty-four `4x4`
graph matrices have rank three.  Their two finite cores are proved, but the
attempted extension induction fails the composed-pencil test and is removed.
The replacement finite-geometry grammar is now executable on arbitrary
normalized Gram/wall/orientation input.  It recognizes the polar/chord,
grid/projectivity, Möbius `S(3,4,10)`, characteristic-five Desargues/cut,
and Type-VII tetrahedral `4+6`
cores without archive names or stored answer Grams, invariantly under
determinant-one shear and wall reversal.  This closes the recognition defect
only.  The resulting five-core class is rigid and already records its four
dimensions, so Kondō's terminal classification supplies no nonthin new
polytope bound or exclusion.  Its honest score is
`L3_PLUS_CLOSED_GENERIC_RECOGNITION_GATE`, with
`Pushback=FAIL_NONTHIN_NOVELTY`, not L4.  The completed `H10/72`
tetrahedral `4+6` construction is now part of the generic recognizer.

The relative-embedding theorem has now been corrected: after resolving the
assigned base cluster, constant-rank/base-change/Fitting, relative
length-one/length-two evaluation, smooth Bertini, Jouanolou irreducibility,
and proper Hilbert bad-locus certificates give a rigorous sufficient
geometric branch.  But the class is P-indexed by the existence of that finite
witness; raw Coxeter/arithmetic conditions do not force a witness and the
certificate language is not exhaustive for all carriers.  The simultaneous
nine-row census finds five zero-step cores and **no nonzero extension**.
Therefore this is `L3 + corrected sufficient extension theorem`, not L4, and
the next completion test remains a nonzero raw-P-forced and jointly
exhaustive extension/incidence class with a genuine AG consequence beyond
the named cores.  See the
[`corrected theorem and census`](coble_corrected_extension_pushback.md).

## 7. Audit table

| route | compact/cusped | independent class `C` | finite packages | class-wide carrier | exact return | independent AG terminal input | theorem on every raw `P` |
|---|---|---:|---:|---:|---:|---:|---:|
| `K3RootComp` | compact | pass | pass | pass | pass | no: boundary is Godement--Meyer/Nikulin--Vinberg | pass, but not flagship |
| `Root2FV` | 8 compact, 110 cusped | pass | pass | pass | pass | no: 118-type bound is imported first | pass, but not flagship |
| local `D4` triality folds | mainly cusped in current stock | pass | pass | fail | isolated positives only | no independent stronger AG bound | fail |
| compact `C_35^conn` coherent K3 | compact | pass | pass for fixed `P` | fail | conditional per branch | pass: Hodge period gives `n<=8` | fail |
| `D_CobVin` incidence | cusped | pass, non-singleton | pass | **one missing gate** | pass for every geometric branch | pass: Kondō gives (5.1) | carrier-positive locus only |
| core-only `CI4` | cusped | arbitrary-input P-only recognizer passes | pass for two cores | pass for the two rigid cores only | pass with `D_CobVin` | conditional pass | actual `L3`; the old extension theorem is false |
| finite-geometry cores (Möbius/Desargues) | cusped | arbitrary-input P-only recognizer passes | pass on both cores | pass on `H11/206,H10/78` with `D_CobVin` | pass | pass after return | `L3` plus closed recognition gate; nonthin pushback fails |
| Type-VII tetrahedral `4+6` core | cusped | arbitrary-input P-only recognizer passes | pass on all five intrinsic frames | pass on the frozen datum | pass | pass after return | actual `L3`; new explicit realization, no class-wide payoff |
| corrected `C_FG_HilbExt` | cusped | P-indexed witness class, not unqualified P-only | finite certificate language passes | pass for every witnessed branch; no raw-P-forced or exhaustive existence | pass after independent `D_CobVin` gate | Kondō gives the four dimensions and characteristic-five restriction | `L3` sufficient theorem; nine-row search has zero nonzero extensions; not L4 |

## 8. Nontrivial slices already proved by local artifacts

The current workspace does contain more than formal blueprints.

1. **Completed compact inverse class:** exactly eight complete integral
   norm-two compact K3 chambers, all in `H2` or `H3`, with literal
   `P -> K3 -> same P` return.
2. **Completed global complete-root class:** 118 finite-volume marked types
   through `H18`, all with exact K3 return; this is the broadest completed
   denominator, though not a new AG-derived boundary.
3. **Coble finite-geometry inverse positives:** `H14/60` and `H17/48` pass
   the common four-point core theorem; `H11/206` passes the Möbius/Frobenius
   compiler; `H10/78` passes the characteristic-five Desargues/sextic
   compiler; and `H10/72` passes the tetrahedral `4+6` line/conic/quartic
   compiler with a four-node/six-cusp boundary.  Every branch has exact
   hidden-root/cusp-lift return.
4. **Nine Coble full-nef outputs:** exact positive chambers occur in
   `H10`, `H11`, `H14`, and `H17`, with the facet counts in the table above.
5. **Golden coherent positives:** all twelve classical golden Lanner
   chambers and Bugaenko's compact `H8` chamber have coherent K3 returns;
   every quadratic rank-20 coherent `H9` branch is Hodge-impossible.
6. **Triality slice:** Bredon `N11` has an exact order-six K3 same-polytope
   realization; eight further packages (`N2`, three `N3` rows, two `N5`
   rows, `N6`, and `N7`) have exact complete-nef fixed-divisor-section
   returns against the 55-stock local denominator.  The residual comparison
   is `P=8,R=15,O=7`; these isolated positives still do not force a
   complete-nef carrier for every ADEG input.

Items 3--6 are valuable positive or boundary results.  Five Coble
finite-geometry constructions are now exact; one arbitrary-input recognizer
covers all five without reading frozen answers.  Their AG pushback is
still thin.  They remain
`L3` until a nonthin input class and non-tautological terminal theorem are
both frozen.

## 9. Final strategy recommendation

Run two tracks in parallel, with unequal expectations.

- **Cusped theorem track:** make `D_CobVin`'s incidence nonemptiness
  input-computable and structural.  Five exact core constructions now cover
  `H10/72`, `H10/78`, `H11/206`, `H14/60`, and `H17/48`.  A normalized
  arbitrary-input recognizer covers all five.  The next concrete goal is to
  enlarge to a
  nonzero extension denominator whose Kondō pushback narrows information not
  already encoded by the selected core.  This is the
  shortest route to a genuine

  ```text
  P -> finite incidence stacks -> exact Coble nef cone -> Kondo
    -> dimensions {10,11,14,17}, char-3 endpoint restrictions.
  ```

- **Compact primal/dual track:** continue the golden finite compiler, but do
  not label coherent essential surjectivity a small lemma.  Search `H9`
  `NS`-only packages and complete roots on the example side; on the theorem
  side seek a new P-only mechanism forcing the field action on the full
  Hodge structure.  The two outcomes are clean: a survivor is a new compact
  example and refutes coherent universality, while a genuine forcing theorem
  yields the sharp `n<=8` bound.

The triality compiler should continue in parallel as a fishing and lineage
engine, not as the main inverse theorem.  The Coble incidence compiler is the
only current route where the finite carrier package, exact same-`P` theorem,
independent AG classification, high-dimensional positive endpoints, and a
plausibly isolated final lemma all coexist.
