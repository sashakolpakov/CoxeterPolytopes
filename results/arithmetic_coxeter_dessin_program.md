# Arithmetic Coxeter dessins: general inverse programme

**Status (2026-08-10):** the marked arrangement compiler, compensated
Borcherds criterion, fixed-support finite-volume cusp theorem, repaired
decorated-surface compression, finite congruence-monodromy theorem, and
literal Gram--Belyi surface return are proved.  The explicit missing gate is
classification pushback: a uniform absolute Gram--Belyi/Hilbert-degree or
dotted-Gram-height bound from maximality or rigid Borcherds geometry.  This
document does not claim a completed unconditional classification theorem.
The independent canonical-chamber-rigidity synthesis now proves the exact
conditional AG payoff and corrects literal full-nef return to an eigen-nef
return whenever the minimal Vinberg field is non-rational.

## 1. Desired dictionary

The input is a marked compact or finite-volume arithmetic Coxeter polytope
\(P\subset\mathbb H^n\), given only by its exact oriented Gram datum.
The proposed arithmetic Coxeter dessin is

\[
 {\cal C}(P)=(X_P,D_P,C_P,\mathbb V_P,\partial_P):
\]

- \(X_P\) is the orthogonal Shimura stack supplied by the sign-correct
  determinant-line extension from signature \((1,n)\) to \((2,n)\);
- \(D_P\) is the marked Heegner mirror arrangement;
- \(C_P\) is the chosen component of the real complement;
- \(\mathbb V_P\) is the orthogonal arithmetic local system retaining the
  quadratic module and facet roots; and
- \(\partial_P\) marks the rational isotropic/Baily--Borel strata.

The required inverse theorem reconstructs the normalized Gram matrix, every
Coxeter label, the oriented face poset, and every finite-volume cusp from
this decorated algebraic object.

The Belyi analogy is structural, not decorative:

| Belyi/dessin datum | Coxeter-dessin datum |
|---|---|
| curve over a number field | orthogonal Shimura stack over its reflex field |
| branch divisor | Heegner mirror divisor |
| local monodromy | distinguished reflection/dihedral inertia |
| dessin chamber | component of the real mirror complement |
| permutation labels | marked integral root local system |
| puncture/cusp | rational isotropic Baily--Borel stratum |
| Belyi function | Borcherds automorphic discriminant, when unobstructed |

The proved three-level carrier/rigidity split is stated separately in the
[Coxeter--Borcherds dessin route](coxeter_borcherds_dessin_route.md): the
arrangement dessin is universal, safe compensation is exact on the real
slice, and only pure or uniformly bounded-slope dessins are presently rigid
enough for classification pushback.

## 2. What is already proved

The existing determinant-line transfer theorem proves, for a marked even
admissible arithmetic quadratic realization:

1. a canonical sign-correct determinant line and a finite list of even
   reflection-stable gluings;
2. the type-IV orthogonal Shimura stack;
3. exact equality between each hyperbolic facet and the real pullback of its
   Heegner mirror;
4. recovery of elliptic candidate intersections, followed by the necessary
   oriented halfspace test for actual faces;
5. recovery of a finite Coxeter label \(m_{ij}\) from the distinguished
   dihedral inertia subgroup of order \(2m_{ij}\);
6. a terminating hidden-sheet enumeration for compact \(P\); and
7. Bruinier's finite linear obstruction: a prescribed Heegner principal part
   is realized by a weakly holomorphic Whittaker/Borcherds product if and only
   if it annihilates the dual vector-valued cusp space.

These statements give a general compact arrangement-level
\(P\longrightarrow{\cal C}(P)\longrightarrow P\) theorem.  The separate
fixed-support cusp theorem extends exact marked recovery to finite-volume
noncompact inputs, and the repaired surface theorem gives a faithful
decorated carrier.  Neither an undecorated ordinary surface nor the universal
compensated product gives classification pushback by itself.

## 3. The universal field split

Let \(k\) be the Vinberg/definition field and suppose \(q\) has signature
\((1,n)\) at one real embedding and is definite at every other real
embedding.

- If \(k\ne\mathbb Q\), then \(q\) is anisotropic over \(k\): a nonzero
  \(k\)-isotropic vector would remain isotropic at a definite conjugate real
  place.  Hence the arithmetic quotient, and therefore a finite-covolume
  reflection chamber, is compact.
- A finite-volume noncompact arithmetic quadratic reflection chamber must
  therefore have \(k=\mathbb Q\).
- Over \(\mathbb Q\), Meyer isotropy for indefinite forms of rank at least
  five implies that a compact quadratic arithmetic chamber can occur only
  for \(n\le3\).

Thus the general compiler naturally splits into a totally-real compact lane
and a rational cuspidal lane.  This is the first field/dimension theorem that
the algebraic dictionary must expose and strengthen.

### Dessin-native finite fields

The finite field list need not be imported from Prasad's arithmetic-volume
formula or ABSW. The separate Mostow--Prasad rigidity theorem is used only
for uniqueness after a Coxeter matrix has been fixed. Define a
bounded universal coefficient-order and Hilbert/Borcherds grammar
\(\mathcal D\), without prelisting fields, and put exact outward unit-normal
geometric Gram coordinates on its splitting/eigencharacter cover. If the
scheme-theoretic image of that geometric Gram map (equivalently Cartan data
modulo diagonal root rescaling) is zero-dimensional, its finite
\(G_{\mathbb Q}\)-set enumerates the unit-normal geometric Gram fields and
their cyclic-product Vinberg subfields. Arbitrary primitive-root and
asymmetric-Cartan coordinate fields are not covered. Surface fibres may remain
positive-dimensional.  Only the stronger goal of enumerating and descending
carrier fields requires the fully rigidified exact-return stack itself to be
zero-dimensional with trivial stabilizers.  The distinguished Lorentzian
embedding is selected only on the real enhancement.

This is an unconditional implication with a conditional input:
zero-dimensionality of the global Gram image is not yet proved.
Bounded Hilbert polynomials give finite type but can leave positive-dimensional
families.  The full marked-incidence tangent complex, or the equivalent
exhaustive Jacobian on all charts, must vanish everywhere.  See the
[dessin-native field theorem](dessin_native_field_effectivity_theorem.md).

There is now a second, strictly scalar finite-field route that does not
require the global Gram image to be zero-dimensional.  On an exact
full-Picard K3 trace carrier, reflection-word traces make all directed cyclic
Cartan products integral.  They reduce to pairs and chordless cycles, and
the positive-definite conjugate Gram matrices force every chordless cycle to
have length at most `n+3`.  After the canonical affine shift
`theta=(-1)^ell*c+2^ell`, each coordinate is totally positive and either
uniformly small or dominant.  A bound on the single bidegree
`H.tau_theta(H)` for every selected large generator cycle bounds every
conjugate and gives finite integral minimal-polynomial boxes.  For `n=2`,
`H^2>=s0`, and bidegree at most `C`, the explicit discriminant bound is
`(32*max(64,C/s0))^30`.  The unresolved input is now geometric and scalar:
exact maximal eigen-nef return must canonically supply `H` in the trace block
and bound those large-cycle degrees.  See the
[chordless-cycle height theorem](k3_chordless_cycle_height_field_bound.md)
and its
[independent hostile audit](k3_chordless_cycle_height_field_bound_independent_hostile_audit.md).

### Global faithful-packet bridge

The dimension step is no longer heuristic.  If one integral dessin packet
contains the complete restriction-of-scalars root module and has rank at
most `B`, then

```text
[k:Q](n+1) <= B.
```

If instead its distinct facet curves have total degree at most `D` against
one integral ample class, then `n<=D-1`.  Fixed `B,D`, coefficient height,
and a finite exact-return grammar give a terminating enumeration of complete
Grams, fields, and finite-volume chambers.  Exact star-Gram and rational-
surface blowup families show that bounded four-node data, bounded individual
curve degree, and the word “surface” do not yield either global ceiling.  In
the six frozen controls the full trace module is already intrinsic; what is
missing is a ceiling derived uniformly from maximality.  See the
[packet theorem](integral_packet_global_dimension_bridge.md), its
[hostile audit](integral_packet_global_dimension_bridge_hostile_audit.md),
and the [six-carrier calibration](six_carrier_integral_packet_rank_calibration.md).

### Intrinsic minimal-carrier dimension theorem

There is now one nonempty AG lane with an absolute dimension bound.  Suppose
the full trace module of `P` embeds isometrically in the numerical lattice of
a minimal smooth projective complex surface, and every simple reflection
extends to a pairing isometry fixing the canonical class.  The canonical
class is then orthogonal to the trace module.  Since that Lorentzian module
already uses the sole positive numerical direction, minimal-surface
classification forces Kodaira dimension zero and

```text
[k:Q](n+1) <= rho(S) <= 20.
```

The positive-characteristic ceiling is `22`.  Actual automorphisms are a
strict subcase; a Wehler K3 realizes an arithmetic ideal triangle in `H2` by
three deck-reflection automorphisms.  Effective nef-wall roots cannot belong
to that strict automorphism lane, because an automorphism cannot send an
effective class to its negative; canonical-preserving monodromy or
correspondence is the relevant broader possibility.

Minimality can now be replaced by an exact operator condition.  If the
exceptional lattice of a contraction to a minimal surface has a finite
numerical operator orbit, their
negative-definite direct sum has finite integral orthogonal group.  A
finite-index finite-covolume trace group has no invariant functional, so the
simultaneous equivariant projection of the trace block vanishes.  The entire
block descends isometrically.  Conversely, orthogonal descent confines the
orbit to the finite orthogonal group of the negative complement.  Finite
numerical operator orbit of the chosen exceptional lattice,
trace-orthogonality, and isometric descent are therefore equivalent.
For actual automorphism carriers with nonnegative Kodaira dimension the minimal
contraction is automatically equivariant, so the `20` bound holds even when
the original surface is nonminimal.  No complement-rank bound is possible in
the correspondence category: the reflection-maximal `Gamma16` K3 followed
by `r` redundant blowups has complement `<-1>^r` for arbitrary `r`.  Thus the
raw compiler must force a **finite numerical operator orbit of a geometric
exceptional lattice**, not bounded
complement rank.

This is not yet global: the universal arrangement and Gram--Belyi
constructions do not presently produce an intrinsic packet with such a
finite exceptional-lattice orbit for every `P`.  See the
[theorem](project_native_minimal_surface_automorphism_dimension_bound.md) and
[hostile audit](project_native_minimal_surface_automorphism_dimension_bound_hostile_audit.md),
and the
[equivariant-minimalization theorem](equivariant_minimalization_trace_descent_theorem.md)
with its
[hostile audit](equivariant_minimalization_trace_descent_hostile_audit.md).

### Canonical-cover and rational Mori-boundary dimension theorem

The operator hypothesis can be weakened further when canonical
orthogonality is supplied geometrically.  A positive-similarity injection of
the full trace space into `K_X^perp` on one minimal surface already gives
`[k:Q](n+1)<=20` over `C` and `<=22` in positive characteristic.  In
particular, the canonical double cover of every complex strong Coble surface
is K3 and carries the full extended Picard lattice, giving full-nef
dimension at most `19` without using the finite-automorphism classification.
The exact action audit separates actual automorphisms, coherent monodromy,
generatorwise monodromy, and divisor correspondences; none is silently
promoted to another.
The rank ceilings themselves are sharp in the frozen corpus: the complex
`X3-H19/36` K3 has Picard rank `20`, and the characteristic-two
`BOR-D4-H21/210` supersingular K3 has rank `22`.  This sharpness does not
supply the missing raw-P carrier theorem.

There is also a direct rational-surface dichotomy.  If the complete
polyhedral Mori cone is generated by smooth rational `(-1)` and `(-2)`
curves, adjunction makes `-K` nef, hence `rho<=10`.  If it instead has
pairwise-disjoint `(-4)` boundary curves and every remaining generator has
type `(-2,B.C=0)` or `(-1,B.C=2)`, the generator span annihilates
`K+B/2`; therefore `B=-2K` and the canonical double cover is K3.  This gives
carrier-relative ceilings `9`, `19`, and `21` for literal full-nef chambers.
See the
[canonical-cover theorem](canonical_orthogonal_k3_cover_dimension_theorem.md),
its [hostile audit](canonical_orthogonal_k3_cover_dimension_theorem_hostile_audit.md),
and the
[rational Mori-boundary theorem](rational_mori_boundary_k3_dimension_dichotomy.md)
with its
[hostile audit](rational_mori_boundary_k3_dimension_dichotomy_hostile_audit.md).

The corrected global transfer statement is now exact.  A literal full-nef
surface chamber has minimal Vinberg field `Q`; irrational fields cannot be
stored in an ordinary rational Neron--Severi nef cone and must be returned by
the distinguished eigen-nef section of the full trace action.  Conditional
on an exhaustive exact compiler, the direct rational, weak-incidence
rational, and minimalizable complex full-trace branches give

```text
n <= 11 and k=Q;
n <= 25 and k=Q;
[k:Q](n+1) <= 20,
```

respectively.  Hence `n<=25` is an AG corollary once raw essential
surjectivity is proved.  The rank inequality does not finite-ize fields;
that requires a bounded discriminant or a finite integral Gram/
correspondence or primitive spectral packet.  Finite maximal chambers then
also require finite hosts or complete packets, host preservation,
reflection saturation, and the complete one-orbit mirror arrangement.  See
the
[canonical-chamber-rigidity theorem](canonical_chamber_rigidity_absw_transfer_theorem.md)
and its
[hostile audit](canonical_chamber_rigidity_absw_transfer_theorem_hostile_audit.md).

This carrier premise is genuinely geometric.  An absolutely maximal
rank-three `(2,4,infinity)` triangle with integral spanning reflections and
the expected numerical wall alphabet nevertheless has its marked rational
canonical class inside the proposed chamber, so it cannot be a rational
full-nef carrier.  Bare maximality therefore does not even force canonical
escape, much less simultaneous facet effectivity or cone completeness.  See
the
[P-first obstruction](rational_pfirst_effectivity_realizability_no_go.md).

Host-reflection saturation nevertheless closes one nontrivial P-first class
positively.  For `I_(1,9)`, relative reflection maximality and the primitive
characteristic anticanonical cusp `F=-K` force a unique ten-facet alcove with
one section root and nine affine fibre roots.  The impossible Hamming fibre
box becomes a positive-root subdivision, while a nine-blowup flex pencil
realizes the saturated alcove as the exact complete nef cone.  This is a
scoped essential-surjectivity theorem, not a global one: other ranks,
adapted hosts, and maximal chambers whose `-K` is not a cusp remain outside.
See the
[anticanonical-cusp closure](rational_pfirst_anticanonical_cusp_saturation_closure.md).

The adjacent timelike anticanonical range is now equally P-first.  On
`I_(1,n)`, `3<=n<=8`, every relatively reflection-maximal finite-volume
chamber with ordinary physical `(-1)/(-2)` facets is `Stab(K)`-conjugate to
the standard simplex.  Here `-K` is a finite vertex of square `9-n`, not a
cusp.  Truncating the flex base cluster constructs all facets and an
anticanonical pencil with no fixed curve, which independently certifies
nefness; reflection saturation then gives exact full-nef return.  The
remaining rational property-4 gap is concentrated in twice-primitive
`(-4)` facets, `-K` outside the chamber, and nonunimodular adapted hosts.
See the
[timelike anticanonical closure](rational_pfirst_timelike_anticanonical_vertex_closure.md)
and its
[hostile audit](rational_pfirst_timelike_anticanonical_vertex_closure_hostile_audit.md).

The wrong-log branch is bounded as well under the exact weak-incidence gate.
If the complete Mori alphabet consists of smooth rational `(-1),(-2),(-4)`
curves, the `(-4)` curves are disjoint, `B.R=0` for every `(-2)` generator,
and `B.E<=2` for every `(-1)` generator, then `D=-(K+B/2)` is nef.
Contracting the `(-4)` curves gives a projective klt surface with
`1/4(1,1)` points and nef anticanonical class.  Orbifold `c2`
semipositivity and the exact formulas

```text
D^2=10-rho+b,
c2_orb=rho+2-7b/4
```

give `rho<=26`, hence `n<=25`.  This is an AG/orbifold dimension theorem,
not a Nikulin--Vinberg or volume argument.  Its carrier hypotheses are real:
the remaining unproved lane is exactly intersecting `(-4)` curves, positive
`(-4)/(-2)` incidence, or aggregate `B.E>=3`; no klt boundary supported on
the negative walls repairs those failures.  See the
[wrong-log rank theorem](rational_wrong_log_orbifold_rank_gate.md) and its
[hostile audit](rational_wrong_log_orbifold_rank_gate_hostile_audit.md).

For reflection-maximal chambers this weak-incidence result now has a
class-wide finiteness upgrade.  Every physical wall acts integrally on
\[
 M=\langle\operatorname{NS}(S),B_1/2,\ldots,B_b/2\rangle,
 \qquad q_M=2q_{\operatorname{NS}},
\]
so maximality saturates the facet group to the full integral reflection
group \(R(M)\).  Its complete mirror chambers form one orbit.  The
algebro-geometric ceiling \(\rho\leq26\), together with the finite
exponent-two interval of adapted hosts, therefore gives only finitely many
abstract chambers, all with \(n\leq25\) and Vinberg field \(\mathbf Q\).
This is an ABSW-type theorem on the already constructed carrier class, not
raw essential surjectivity.  See the
[weak-incidence maximal-finiteness theorem](rational_wrong_log_reflection_maximal_finiteness.md)
and its
[independent hostile audit](rational_wrong_log_maximal_finiteness_hostile_audit.md).

The hostile complement has also been reduced without a small-diagram census.
The finite-vertex discrepancy theorem is uniformly `epsilon`-klt for every
strict `epsilon<1/2`, while explicit free blowup chains show that it cannot
bound rank.  Minimal H1 and H2 incidences have sharp log-minimal core ranks
`14` and `12`; H3 has an exact spherical/affine/hyperbolic star trichotomy.
Finite volume closes the affine abundance gate: the null cycle
`F=4E+B_1+...+B_4` is base-point-free because its nonreduced rational-tree
fibre has torsion-free multidegree-zero Picard group.  The resulting
genus-three fibration has relatively minimal rank at most `26`.
Finite volume forces the `W_T`-saturated mixed projections to span the erased
tail, but one `W(B_m)` orbit already spans arbitrary `m`.  An effective affine
`B_m` blowup packet makes this sharp with one mixed curve of constant
coefficient energy.  Consequently no individual-facet complexity bound can
finish the proof.  The explicit ruled packet is nevertheless not a maximal
counterexample: on its full Neron--Severi host the ruling-swap root `c-d`
has a mirror through the ample cone, so it supplies a proper discrete
reflection overgroup in every rank.  Its universal visible simplex has
finite volume exactly for `m<=4`, and every later finite-volume completion
must add a hidden wall with `sum_(1..5) mu_i>2a+b`.  This is a full-host
family obstruction, not yet the adapted-host theorem for a general
genus-three vertical tail.

There is nevertheless a uniform global width theorem.  For incomparable
exceptional valuations `i,j`, the relative ample cone contains a class
orthogonal to `e_i-e_j`.  That difference reflection preserves the adapted
host exactly when every boundary pairing is even.  Hence a
reflection-maximal exceptional forest has at most `2^b` incomparable
vertices and is covered by at most `2^b` chains; with no half-boundary it is
one chain.  This does not bound chain depth, and odd parity does not certify
maximality.

The direct full-NS minimal model is now forced to be `P2` or `F2`: every
blowup of `F0` retains the ruling-swap mirror through `Amp`, and the strict
negative section on `F_e` has an integral facet reflection only when
`e+t|2`.  On the marked standard frames, visible finite volume lasts only
through chain lengths `9` and `8`, respectively.  The first escape is
genus-one and can be cut only by a doubled `(-4)` wall; one such numerical
wall removes the old rays but produces a new mixed escape two stages later.
This is an exact no-go for the diagonal proof, not a depth theorem.

The first mixed ray supplies the missing depth theorem on those exact marked
frames.  For `m>=12` over `P2`, the spacelike vector
`7H-3e1-2e2-...-2e12` is nonnegative on every complete physical
`(-1),(-2),(-4)` wall; over `F2`, `4s+9f-2e1-...-2e11` does the same for
`m>=11`.  The proof is parity-free and uses only adjunction, proximity, and
the marked line/fibre inequalities, giving dimension at most eleven on
either frame.

Reflection maximality now forces those frames on the direct full-NS host.
The exceptional forest is one free chain (satellites violate the complete
negative alphabet or facet integrality), and the free crystallographic
Coxeter root module forces `H-e1-e2-e3` or `f-e1-e2` to be an actual facet
even when the chamber is nonsimplicial.  Therefore every direct full-NS
carrier in the declared complete negative-extremal alphabet satisfies
`n<=11`.  Square-zero tangent Mori rays are allowed.  This is an AG class
theorem, not a chamber census; adapted hosts remain outside its scope.

The genus-three carrier supplies no hidden stable bound.  Vertical blowups
preserve its relative canonical algebra/model, Hodge bundle, moduli map,
and stabilized family, and relative Noether cancels their changes.  Nor do
the affine marks and canonical colors bound the cusp: formally physical
colored `D_n` packets exist in every even rank.  The coordinate refinement
is controlled by ambient boundary parity, which is not encoded in the
factor Gram, and global surface realization remains open.

The remaining adapted/hostile completion therefore needs a global total marked degree,
finite contraction orbit/equivariant descent, or bounded canonical-effective
orbit; higher/dotted H1/H2 incidences remain open.
See the
[hostile-incidence synthesis](rational_hostile_incidence_core_tail_reduction.md),
its [hostile audit](rational_hostile_incidence_core_tail_reduction_hostile_audit.md),
the [mixed-attachment lemma](rational_finite_volume_mixed_attachment_lemma.md),
the [sharp tail no-go](rational_global_exceptional_tail_bridge_no_go.md),
the [ruled-family maximality obstruction](rational_affine_b_tail_global_completion_obstruction.md),
its [hostile audit](rational_affine_b_tail_global_completion_obstruction_hostile_audit.md),
the [exceptional-forest comparability theorem](rational_exceptional_forest_comparability_theorem.md),
its [hostile audit](rational_exceptional_forest_comparability_hostile_audit.md),
the [direct minimal-model reduction](rational_direct_full_ns_minimal_model_reduction.md),
its [hostile audit](rational_direct_full_ns_minimal_model_reduction_hostile_audit.md),
the [single-chain standard-frame gate](rational_single_chain_standard_simplex_gate.md),
its [hostile audit](rational_single_chain_standard_simplex_gate_hostile_audit.md),
the [first-mixed-ray theorem](rational_single_chain_first_mixed_ray_bound.md),
its [hostile audit](rational_single_chain_first_mixed_ray_bound_hostile_audit.md),
the [direct full-NS dimension theorem](rational_direct_full_ns_n11_dimension_theorem.md),
its [hostile audit](rational_direct_full_ns_n11_dimension_theorem_hostile_audit.md),
the [genus-three vertical-invariant no-go](rational_genus3_vertical_blowup_invariant_no_go.md),
its [hostile audit](rational_genus3_vertical_blowup_invariant_no_go_hostile_audit.md),
the [colored-cusp no-go](rational_affine_h3_colored_cusp_factor_no_go.md),
its [hostile audit](rational_affine_h3_colored_cusp_factor_no_go_hostile_audit.md),
the [affine-H3 abundance theorem](rational_affine_h3_finite_volume_abundance.md),
and its
[hostile audit](rational_affine_h3_finite_volume_abundance_hostile_audit.md).

The effective log-boundary escape is independently closed on the minimal
H1/H2 pairs.  A free tail of length \(m\) lowers
\([-(K+\Delta)]^2\) by at least \(m\), with the zero tail boundary uniquely
optimal; but nefness on either genus-one core already forces coefficient one
on both core components.  The adjoint section makes \(K+\Delta\) effective
and nonzero after a tail blowup, so its negative is not nef.  Hence
log-Fano/BAB and log-BMY cannot supply the missing tail bound.  This exact
scope excludes higher/parallel/dotted incidences.  See the
[tail-coefficient theorem](rational_log_fano_tail_coefficient_no_go.md) and
its
[hostile audit](rational_log_fano_tail_coefficient_no_go_hostile_audit.md).

The reflection input to this dichotomy is now sharp.  On a rational
unimodular `NS`, integral primitive reflection plus smooth-rational
adjunction gives only `(-1),(-2)` physical walls, with a
twice-primitive `(-4)` exception.  The Coble incidence table is equivalently
the common fixed-class condition for `K+B/2`.  It is not a consequence of
maximality: the maximal full-nef `(2,3,infinity)` triangle on a blow-up of
`F_4` fails it.  Nor can one evade the K3 ceiling with higher standard
denominators.  The tame adapted cover for disjoint `B_i^2=-2m_i` has trivial
canonical bundle and `c2>=12*lcm(m_i)`; Noether forces `lcm(m_i)=2` and all
`m_i=2`.  See the
[reflection-integrality gate](rational_nef_reflection_integrality_gate.md)
and its
[hostile audit](rational_nef_reflection_integrality_gate_hostile_audit.md).

The numerical front end no longer has the old row-specific metadata leak.
Raw Stage 0 starts with only a unimodular Lorentzian Picard Gram and primitive
oriented full-span walls, derives every adjunction-compatible physical scale
and canonical class, and selects the weak or standard log branch on its
inferred reflection module.  Its hostile audit independently checks alternate
wall bases and residue-module membership.  The Stage-0--A bridge derives the
canonical class, boundary indices, and boundary weights of every accepted
index-two packet; the four Coble controls reproduce the native Stage-A
inputs, toric `H6` routes weak, and the maximal `F4` triangle is rejected.
Numerical Stage A then exhausts all supported contraction frames and emits
all `P2/F_e` endpoint-coordinate programs.  See the
[raw theorem](raw_adjunction_log_stage0.md),
[raw hostile audit](raw_adjunction_log_stage0_hostile_audit.md),
[bridge](raw_to_pb_atlas_stage0_a_bridge.md),
[bridge hostile audit](raw_to_pb_atlas_stage0_a_bridge_hostile_audit.md), and
[Stage-A theorem](pb_atlas_numerical_contraction_stage_a.md) with its
[hostile audit](pb_atlas_numerical_contraction_stage_a_hostile_audit.md).
Numerical
[Stage B1](pb_atlas_numerical_cluster_stage_b1.md) consumes every endpoint
program and emits its reverse-blowup/proximity DAG, all wall divisor rows,
and indexed formal `-K/-2K` jet templates.  All `282` orbit representatives
compile; no rank or nonemptiness is inferred from expected dimension.

[Stage B2a](pb_atlas_cluster_chart_stage_b2a.md) now supplies literal
proper/free/satellite charts and polynomial jet matrices for all `42` plane
endpoint-orbit representatives; its
[hostile audit](pb_atlas_cluster_chart_stage_b2a_hostile_audit.md)
independently reconstructs the coefficients and validates the strict input
contract.  The first full geometric B2b component is also closed.  Selected
only by its proximity forest `[2,2,2,1,1,1,1]`, the generated `H10` packet
has one collision-free incidence component, empty `|-K|`, a unique reduced
absolutely irreducible antibicanonical sextic, nineteen exact smooth-rational
noncontracted facet curves, and ten coefficientwise matched exceptional
facets partitioning all thirty walls.  It therefore constructs a strong
Coble surface over `Q/C`, so the canonical cover gives the carrier-class
bound `n<=19`.  A helper-independent hostile replay confirms every chart,
Groebner, collision, section, factorization, multiplicity, and wall-partition
gate.  See the
[generated B2b theorem](pb_atlas_h10_geometry_stage_b2b.md) and its
[hostile audit](pb_atlas_h10_geometry_stage_b2b_hostile_audit.md).

This closes one generated component, not the variable raw class.  Other
plane and ruled packets and essential surjectivity remain open; same-`P`
return additionally requires reflection maximality.

For reflection-maximal literal full-nef outputs, the terminal finiteness step
is also closed.  A missing integral host reflection would generate a discrete
reflection overgroup, so maximality forces the facet group to equal the full
integral reflection group `R(M)`.  Its complete mirror arrangement has one
chamber orbit.  The direct rational and standard Coble AG rank theorems leave
only finitely many host lattices, hence finitely many abstract chambers and
field `Q`, without an independent facet-degree or arithmetic-volume bound.
This is still carrier-relative: the raw-to-geometric incidence arrow remains
Stage B.  See the
[maximal-carrier finiteness theorem](rational_full_nef_reflection_maximal_finiteness_theorem.md)
and its
[hostile audit](rational_full_nef_maximal_finiteness_hostile_audit.md).

For the maximal standard-Coble output, hidden-wall exhaustion is no longer a
separate compiler step.  The equation `B=-2K` and adjunction force every
negative irreducible curve into the three adapted integral reflection types.
Saturation puts all their mirrors in `W(P)`, while an ample class lies in the
proposed chamber because every displayed facet is effective; hence the
proposed chamber is the complete nef cone.  Stage B must still prove the
strong-Coble incidence open and irreducibility/effectivity of every proposed
facet.  See the
[maximal-Coble exact-return theorem](reflection_maximal_coble_exact_nef_return.md)
and its
[hostile audit](reflection_maximal_coble_exact_nef_return_hostile_audit.md).

The same terminal argument works for a nonzero weak-log class.  If Stage B
constructs a disjoint `(-4)` boundary and independently certifies
`D=-(K+B/2)` nef, adjunction classifies every hidden negative curve into the
adapted integral reflection alphabet.  Relative maximality and the ample
orientation then give literal full-nef return, with `Bl_p(F4)` as the exact
non-Coble calibration.  Testing `D` only on the proposed facets is circular;
the compiler must return a basepoint-free multiple, a morphism, or an
independently complete Mori calculation.  See the
[weak-log return theorem](reflection_maximal_weak_log_exact_nef_return.md)
and its
[hostile audit](reflection_maximal_weak_log_exact_nef_return_hostile_audit.md).

There is now an exact finite certificate for the independent nef input.
Given a certified finite affine cover and returned global sections of
`O(mD)`, a localized polynomial identity
`g^N=sum(q_j n_j)+sum(r_l I_l)` on each chart proves that the sections
generate, so `mD` is base-point-free and `D` is nef.  The toric
`Bl_p(F_4)` provider implements the whole cover and certifies
`2D=3B+12F-2E`; the general PB-Atlas B1/B2a adapter remains a specified but
unimplemented contract.  See the
[P-first semiampleness certificate](pfirst_weak_log_semiampleness_certificate.md).

For every packet that actually supplies this finite certificate, the terminal
AG transfer is now complete without a complete Mori census.  Nefness of `D`,
adjunction, adapted-host reflection maximality, and the ample orientation
return the entire nef cone; the boundary contraction and orbifold Chern
inequality give `n<=25` (`n<=9` for empty boundary), the full-nef field is
`Q`, and finite hosts plus the complete mirror orbit give carrier-relative
finiteness.  The theorem is conditional on a supplied complete carrier
certificate; its producer does not yet consume arbitrary B1/B2a packets.
See the
[complete weak-log transfer](pfirst_weak_log_complete_carrier_transfer.md)
and its
[independent hostile audit](pfirst_weak_log_complete_carrier_transfer_independent_hostile_audit.md).

This certificate is genuinely extra geometry on an isotropic ray.  For a
primitive nef smooth genus-one divisor `F` with `F^2=K.F=0`, the exact
criterion is torsion of `O_F(F)`: order `m` is the first moving index, while
a non-torsion character fixes every multiple.  Two nine-point cubic blowups
have identical marked `(NS,K,F=-K)` and nefness but opposite answers.  Hence
an NS/eigenlattice packet alone cannot provide the required MMP morphism; the
period character or an explicit unit certificate must be returned.  See the
[isotropic period gate](rational_isotropic_weak_log_period_gate.md).

Conversely, the non-rational rank bound cannot replace the missing finite
arithmetic return.  Explicit rank-four minimal K3 packages now give exact
finite-length eigen-nef intervals, complete restricted-root saturation,
reflection maximality in their declared `Z[t]` trace hosts, and a primitive
ample class of fixed square `14`, while realizing infinitely many real
quadratic fields.  The result is relative to the declared order and has
hyperbolic dimension one; overorder/commensurable maximality and `n>=2`
remain open.  It nevertheless rules out every field-finiteness argument
based only on rank, periods, exact eigen-nef return, host-relative saturation,
or fixed polarization square.  A bounded discriminant or bounded primitive
spectral/operator packet remains necessary.  See the
[exact maximal eigen-interval theorem](rank4_k3_exact_maximal_eigen_interval_infinite_fields.md)
and its
[hostile scope audit](rank4_k3_exact_maximal_eigen_interval_infinite_fields_hostile_scope_audit.md).

The orthogonal higher-dimensional escape is closed.  Duplicating the
negative field line produces a compact restricted pentagon, but its exposed
mixed wall has nonintegral cyclic Cartan product and cannot preserve any
commensurable full trace lattice; every higher orthogonal duplicate contains
the same bad two-coordinate packet.  Hence `n>=2` requires a nonorthogonal
coupling or a theorem bounding the resulting integral mixed-wall height.
See the
[orthogonal-stabilization no-go](rank6_k3_orthogonal_stabilization_nonarithmetic_no_go.md)
and its
[hostile audit](rank6_k3_orthogonal_stabilization_nonarithmetic_no_go_hostile_audit.md).

The field-generator combinatorics itself is no longer an open end.  Pair
cycles and chordless cycles generate the Vinberg field, and chordless length
is at most `n+3`, irrespective of the number of facets.  Consequently the
only surviving `n>=2` escape is an unbounded polarized correspondence degree
for at least one large signed generator cycle.  Do not launch an unbounded
facet or cycle census; attack that one scalar height or construct the
nonorthogonal family that makes it diverge.

This is genuine AG pushback once such a carrier is constructed; it is not yet
a global ABSW replacement.  Two old prime-binomial branches remain exact
geometric false positives, but the generated B2b component proves that the
repaired chart/evaluation route can reach a real strong-Coble open.  The next
executable gate is uniformity: extend the same generated avoidance and
effective-facet proof to the variable plane/ruled packet class before
invoking maximal exact return and the K3-cover bound.  Its original failure
boundary is recorded in the
[P-first pilot audit](coble_prime_binomial_k3_rank_l4_pilot.md).

### Universal literal Gram--Belyi carrier

There is now a direct Belyi construction which starts from every arithmetic
chamber, without first choosing a Shimura or divisor-lattice species.  For
each facet pair retain the scaling-invariant squared normalized Gram entry

\[
u_{ij}=\frac{(r_i,r_j)^2}{r_i^2r_j^2}
\]

and, when it is nonzero, the marked point
\(\lambda_{ij}=u_{ij}/(1+u_{ij})\).  Attach all minimal normalized Belyi
polynomials taking that point to zero or one; an extra bit treats the
orthogonal case.  A Belyi backbone joins the leaves into a connected nodal
curve, and its product with \(\mathbb P^1\) is a projective SNC Belyi surface.
The marked object retains its distinguished real place, acute convention,
time sheet, normalization and conductor.  It therefore returns the complete
geometric unit-normal Gram and the same chamber, though not primitive
integral root lengths without a separate lattice marking.

Writing \(\beta_{\rm GB}\) for the total leaf degree and
\(M=\binom f2\), the pulled-back polarization has

\[
A_P^2=2(M+\beta_{\rm GB}),\qquad \beta_{\rm GB}\ge f(f-1).
\]

Thus bounded \(A_P^2\) or bounded \(\beta_{\rm GB}\) gives a finite effective
list of marked Belyi packets, exact Grams, chambers and Vinberg fields, with
no arithmetic-volume or rigidity theorem.  This is the first literal answer
in this project to the requested “dessins for Coxeter polytopes” direction.  It is also an
honest conditional theorem: the missing global input is a natural reason
that maximal or geometrically primitive chambers have bounded absolute
degree.  An exact `H10` redoubling family proves that finite-index descendants
have unbounded degree even at fixed rational field, lattice and root alphabet.
See the [theorem](gram_belyi_coxeter_dessin_theorem.md) and
[hostile audit](gram_belyi_coxeter_dessin_theorem_hostile_audit.md).

The first exact local passport table shows why this degree is not just a
relabelled facet count.  The \(m=3,4,6\) and parallel packets have exact
degrees \(5,3,7,2\); the \(m=5\) packet has degree at least \(29\), by a
\(29\)-adic Newton-polygon obstruction, with an explicit finite upper
witness.  Hence
\[
\beta_{\rm GB}\ge
2\binom f2+3N_3+N_4+27N_5+5N_6
\]
before dotted excess.  See the
[standard-packet calibration](gram_belyi_standard_packet_calibration.md).

The initial follow-up degree-\(29\) enumeration reached a rigorous scoped
boundary.  It excludes every one-, two-, and three-vertex lambda-colored
fibre.  In particular, the last `1,372` three-black/`WW` marked records
collapse to `784` coefficient systems and all are eliminated by exact
univariate or finite-extension gcd certificates.  It also excludes all `210`
quadratic-factor ansatz trials.  Hence at that stage a witness had to have at
least four lambda-colored vertices, while the exact \(m=5\) height remained **OPEN**, with
no timeout used as evidence; see the
[degree-29 reduction](gram_belyi_m5_degree29_reduction.md) and the
[three-black theorem](gram_belyi_m5_degree29_three_black_ww.md).

At four black vertices, the exact quotient has `12,752` master coefficient
systems for `42,994` marked records.  The first bounded resultant pass excludes
every two-black-endpoint and one-black/critical-white case.  The second closes
all `WW` statuses by a sharp degree-`756` critical resultant and `10,696`
exact tower-field endpoint gcds.  The final pass checks both `BW` and `WB`
fixed endpoint values on `28,159` exact critical tower orbits and eliminates
the last `3,458` systems / `10,374` records.  Thus a degree-`29` witness must
have at least five black vertices; the exact minimum remains **OPEN**.  See
the [four-black reduction](gram_belyi_m5_degree29_four_black_reduction.md),
the [four-black `WW` theorem](gram_belyi_m5_degree29_four_black_ww_simple.md),
the [four-black completion theorem](gram_belyi_m5_degree29_four_black_one_black_simple.md),
and its [hostile audit](gram_belyi_m5_degree29_four_black_one_black_simple_hostile_audit.md).

The five-black marked-passport row has now been quotiented exactly:
`190,047` raw records become `33,079` binary-status systems.  The largest
single endpoint block, `BB`, is completely excluded by `11,011` exact
projective finite-field certificates, removing `55,055` marked records.
The endpoint quotient `G=DQ+R` makes every white endpoint a zero of `Q`,
independent of simple/critical status, so `BW/WB/WW` collapse to
`4,225/4,225/1,292` distinct necessary systems.  The exact `(5,1^24)` factorization pass
then excludes all `1,292` black records / `3,876` `SS/SC/CS` marked records.
The remaining five-black target has `9,742` status-free necessary systems
(`22,068` in the conservative binary-status ledger) / `131,116` records;
the global remainder is `42,503,404`.  In the next exact partition
`(3^2,1^23)`, the uniform infinity gate is complete and the predetermined
`a>=21` sublocus is impossible: `12` systems / `48` records are removed.
This leaves `1,280` systems / `5,120` records in that partition, `131,068`
five-black records, and `42,503,356` globally.  The fixed-factor tower
`K=x(x-1)` then closes all `1,292` `CC` systems through dimension-`24` base
algebras and exact value tests.  After overlap this removes `1,280` more
records, leaving `3,840` `SS/SC/CS` records in the partition, `129,788`
five-black records, and `42,502,076` globally.  Separate one-critical towers
then close all `1,280` unresolved `CS` rows and the exact `41`-row `SC`
tranche `17<=a_lambda<=20`.  The partition residual is `2,519`, split as
`1,280 SS + 1,239 SC`; five-black and global remainders are `128,467` and
`42,500,755`.
The broader census stops there, and degree `29` stays **OPEN**.  The carrier
surface itself gives `e(S_P)=f(f-1)+4` and a conditional dimension bridge,
but fixed-volume stable pairs with unbounded marked fields prove that a
zero-dimensional Gram image or a separate absolute degree bound is still
essential.
See the [five-black quotient](gram_belyi_m5_degree29_five_black_quotient.md),
the [five-black `BB` theorem](gram_belyi_m5_degree29_five_black_bb.md), and
the [five-black `WW` partition theorem](gram_belyi_m5_degree29_five_black_ww_partition5.md),
with its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition5_hostile_audit.md),
the [`(3^2,1^23)` high-degree theorem](gram_belyi_m5_degree29_five_black_ww_partition33_high_a.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_high_a_hostile_audit.md),
the [complete `CC` tower theorem](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cc_tower_hostile_audit.md),
the [`CS` theorem / `SC` tranche](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche.md),
its [hostile audit](gram_belyi_m5_degree29_five_black_ww_partition33_cs_full_sc_tranche_hostile_audit.md),
and the [surface-bound obstruction theorem](gram_belyi_carrier_dimension_field_bound_obstruction.md).

## 4. Finite-volume cusp theorem

For a rational cusp represented by a primitive isotropic vector \(u\), fix a
full-rank parabolic translation cell, a rational Witt partner, a finite
fixed-denominator discriminant-coset support, and one of the finitely listed
determinant extensions.  Use a Witt decomposition and a standard
horospherical fundamental domain.  A
candidate mirror root \(r\) has two cases.

1. If \((r,u)\ne0\), the completed-square mirror identity applies and the
   mirror's penetration height is bounded explicitly
   in terms of \(r^2\), the denominator of \((r,u)\), and the compact
   Euclidean cusp cross-section.  Such mirrors lie in a computable truncated
   thick part.
2. If \((r,u)=0\), use the original affine pairing equation rather than
   dividing the completed-square identity by zero.  The mirror is vertical
   through the cusp.  At fixed norm
   and discriminant class, its image in \(u^\perp/u\) lies in a
   positive-definite lattice shell and is finite modulo the parabolic
   translation lattice.

The theorem combines these two finite enumerations with Vinberg's
complete affine/parabolic-neighborhood criterion.  It must prove that the
Heegner arrangement has exactly the listed cusp branches and no hidden
interior cut.  No timeout or bounded-height scan can replace this proof.

The general fixed-support finiteness theorem and exact three-cusp toric
calibration are now frozen in the
[finite-volume cusp transfer](finite_volume_coxeter_dessin_cusp_transfer.md).
Input-specific support compilation and literal parabolic replay remain
separate obligations.
The [independent hostile audit](finite_volume_coxeter_dessin_cusp_transfer_hostile_audit.md)
passes after repairing the vertical edge case and independently replaying
the complete three-cusp toric calibration.

## 5. Why an undecorated surface is not enough

For \(n=2\), \(X_P\) is already a complex surface.  For \(n>2\), a sufficiently
general complete-intersection surface \(S_P\subset X_P\) retains restricted
Heegner curves and the orthogonal local system.  It does not, through its
ordinary curve-intersection graph alone, retain triple and higher mirror
intersections, hence cannot reconstruct the rank-\(n\) face poset.

The honest surface object is therefore

\[
 (S_P,\mathbb V_P|_{S_P},D_P|_{S_P},
   \text{ambient higher-intersection and real-chamber marking}).
\]

Proving that a smaller intrinsic surface decoration carries the same
information is an open theorem target.  Merely invoking Belyi's theorem on an
arbitrary surface section would forget the marked mirrors and is rejected.

The first one-point cusp statement was false: after SNC resolution, the cusp
walls and toroidal boundary cannot all pass through one point (the
`(2,3,infinity)` triangle is the minimal counterexample).  The hostile-audited
[repaired compression theorem](coxeter_dessin_marked_surface_compression_theorem.md)
uses the full labeled exceptional boundary cluster and contraction map, and
retains the original stacky root/reflection/descent and ambient chamber data.
It is a noncanonical faithful decorated carrier, not intrinsic reconstruction
from an ordinary surface.  Quantitative cluster bounds remain open.

## 6. Experimental calibration matrix

The SHA-pinned numerical and directional entries are maintained in the
[exact calibration ledger](arithmetic_coxeter_dessin_calibration_ledger.md).

| Input | Field/compactness lane | Existing information | Next dessin test |
|---|---|---|---|
| Bugaenko \(H^8/11\) | real quadratic, compact | determinant extension, full real shell, pure and first compensator obstructions | construct decorated surface and quantify the least safe Borcherds correction |
| toric \(H^6/9\) | rational, noncompact | exact full-nef surface and 3 cusps | compare surface curves with the orthogonal cusp marking |
| Coble \(H^{10}\) | rational, noncompact | exact inverse/full-nef carrier | compile its determinant dessin without using the surface |
| Coble \(H^{11}/26\) | rational, noncompact | raw plane inverse; accepted line--quintic chart point has finite algebra and Vinberg field \(\mathbb Q\) | promote the one-point Galois calibration to a zero-dimensional variable-class stack |
| characteristic-three \(H^{14}/H^{17}\) | rational lattice, noncompact | vertical Coble carriers | identify what the characteristic-three carrier adds beyond the characteristic-zero Shimura dessin |

The matrix must be enlarged by at least one further compact field and one
non-Coble rational cusp example before general hypotheses are frozen.

## 7. Clear route to completion

1. Completed at fixed support: prove the rational cusp theorem in Section 4.
   Next compile the finite Heegner support and replay all parabolic
   neighborhoods for each noncompact input class.
2. Prove the decorated-surface restriction/faithfulness theorem, explicitly
   stating which higher incidence data remain external.
3. Compute the vector-valued obstruction on the calibration matrix and use
   the proved finite congruence-dessin compiler to extract invariant patterns
   depending on \(n,[k:\mathbb Q]\), discriminant, and level.  The exact
   theorem and \(E_{10}\) level-five calibration are in the
   [congruence-monodromy report](coxeter_congruence_dessin_theorem.md).
4. Define a non-singleton P-only bounded class and prove finite exhaustive
   carrier alternatives.  The natural first coordinate is the Coxeter-cover
   degree \(\deg_{\rm Cox}\): the
   [bounded-degree theorem](bounded_degree_schreier_coxeter_dessin_theorem.md)
   already gives a finite exact permutation census and an explicit Gram-height
   bound at fixed degree.  What remains is to bound that degree from marked
   Borcherds/Hilbert complexity rather than assume it.
   In parallel use the absolute Gram--Belyi degree as the literal dessin
   coordinate: prove that a natural maximal/intrinsic carrier class bounds
   \(A_P^2\), then enumerate its finite passport list directly.
5. Transfer an algebraic or automorphic classification back to a new
   example, dimension/facet bound, finiteness theorem, or nonexistence
   theorem.  Without this nonzero pushback, the programme remains an
   arrangement-level L3 theorem rather than the requested general inverse.
6. In parallel, construct the full rigid carrier equations over
   \(\mathbb Q\), prove at least their unit-normal geometric Gram image
   zero-dimensional,
   and enumerate its Galois orbits and cyclic-product Vinberg subfields.
   Prove the whole carrier locus zero-dimensional only for the stronger
   surface-field/descent conclusion.  The current H11/26 calculation
   certifies only one accepted \(\mathbb Q\)-point.

## 8. Primary sources

- E. B. Vinberg, *Hyperbolic reflection groups*, Russian Math. Surveys 40
  (1985), for arithmetic reflection data and finite-volume diagram tests.
- R. Borcherds, *Automorphic forms with singularities on Grassmannians*,
  Invent. Math. 132 (1998), for automorphic products with special-divisor
  support.
- J. H. Bruinier, *Regularized theta lifts for orthogonal groups over totally
  real fields*, J. Reine Angew. Math. 672 (2012), for Whittaker lifts,
  special divisors, and the exact cusp-form obstruction.
- The project report
  [polytope-to-automorphic-divisor transfer](polytope_to_automorphic_divisor_transfer_theorem.md)
  for the proved determinant-line, real-slice, inertia, and compact
  hidden-sheet statements.
