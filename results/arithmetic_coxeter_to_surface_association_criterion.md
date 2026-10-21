# Arithmetic Coxeter polytope to algebraic surface: the exact association ladder

**Status:** proved synthesis of the frozen rational-surface, K3,
totally-real, fixed-decoration Looijenga-pair, and one-boundary Coble carrier theorems
(2026-08-08).

## Direction of the theorem

The primary map is

```text
arithmetic Coxeter polytope P
    -> polytope-derived lattice/action packages A(P)
    -> moduli spaces of marked algebraic surfaces M_a(P).
```

It is not merely a construction of a polytope from a surface.  The surface
and every divisor/action marking are constructed from the arithmetic Coxeter
data.  Surface theory then decides how strongly the original polytope is
realized and supplies necessary-and-sufficient carrier criteria.

This file concerns `P -> surface`.  Surface-first constructions are useful
positive examples, but they are not evidence for inverse universality.  A
recognizer for one carrier species is a genuine inverse theorem only in that
species; its failure cannot exclude the polytope until the list of carrier
species and packages is exhaustive for a surface-free input class.

There are three exact levels.  The first is universal but weak; the third is
restrictive but gives a complete nef cone.  The project's main open problem
is to construct and classify useful intermediate levels.

## Arithmetic input

Let `F` be a totally real field of degree `d`, let `M` be an
`O_F`-lattice of rank `n+1`, and let `h` be an integral symmetric field form.
Fix a distinguished real embedding `sigma_0` and an oriented finite-sided
Coxeter chamber `P` cut out by roots `r_1,...,r_s`.  Use the sign convention

```text
signature(h^sigma_0) = (1,n),
h^sigma is negative definite for sigma != sigma_0.
```

The reflections in the `r_i` are required to preserve `M`.  For `F=Q`, the
condition at the other real places is empty.  For general `F`, the
one-indefinite-place condition is exactly what makes a trace lattice
compatible with Hodge index.

Here “arithmetic” means the quadratic-form/lattice-form class just specified.
No universality statement is made for non-simple arithmetic constructions
not encoded by such a totally-real quadratic lattice; extending the
association to those forms is a separate problem.

Choose a totally positive codifferent normalization `beta` and form

```text
N(P,beta) = (M_Z, Tr_F/Q(beta h)).
```

After the usual common sign change this agrees with the totally-negative
trace convention used elsewhere in the project.  Its signature and rank are

```text
signature N = (1, d(n+1)-1),
rank N      = d(n+1).
```

All positive wall rescalings, even saturations, discriminant gluings, and
orientation-compatible trace normalizations admitted by a bounded input
specification belong to `A(P)`.  For a fixed exact integral realization this
set starts with one root datum; for an abstract diagram, proving that the
normalization list is finite and exhaustive is a separate arithmetic theorem.

## Level I: every admissible arithmetic polytope has a weak rational surface

### Universal association theorem

Every input above has a constructible tuple

```text
(X_P, Gamma_P, pi_P; D_1,...,D_s),
```

where

- `X_P` is a blow-up of `P^2`;
- a common square multiple of the trace lattice embeds in `NS(X_P)`;
- `Gamma_P` is an integral self-adjoint divisor-product correspondence;
- `pi_P` is the rational Chow projector onto the trace-lattice Tate summand;
- the distinguished real eigenspace of `Gamma_P` has signature `(1,n)`; and
- the divisor halfspaces defined by the `D_i` recover the labelled, oriented
  chamber `P` exactly, up to one common positive metric scale.

If the trace lattice has signature `(1,m)`, the general construction uses at
most `1+4m` blow-ups.  Thus one may take

```text
b <= 4 d(n+1) - 3.
```

Conversely, a surface carrier required to contain the **full trace form** has
positive index one only if every nondistinguished conjugate field form is
definite of the opposite sign.  Therefore

```text
one Lorentzian place + definite conjugate places
  <=> Hodge-compatible full-trace weak surface association.
```

This is already a broad arithmetic-polytope-to-surface theorem with no K3
rank ceiling.  It is deliberately weak: the `D_i` need not be effective or
irreducible and the extracted chamber need not be the nef cone.  Precisely
because this level is universal, its existence alone cannot obstruct
arithmetic Coxeter polytopes.

## Level II: exact weak K3 association

For each oriented even integral package `a in A(P)`, with lattice `N_a` of
signature `(1,r-1)`, define `M_a^K3(P)` to be the moduli space of marked
projective complex K3 surfaces with exact Picard lattice `N_a`, retaining the
labelled divisor mirrors and their oriented halfspaces.  Put

```text
M_weak^K3(P) = disjoint_union_(a in A(P)) M_a^K3(P).
```

### Necessary-and-sufficient criterion

For a fixed package `a`,

```text
M_a^K3(P) is nonempty
  <=>
N_a embeds primitively in Lambda_K3 = U^3 + E8(-1)^2.
```

Equivalently, `r<=20` and there is an even complement `T_a` of signature
`(2,20-r)` with

```text
q_(T_a) = -q_(N_a).
```

This is a finite local discriminant-form decision, including the nonstable
odd-prime and 2-adic cases.  In the stable range, the sufficient shortcut is

```text
22-r >= length(A_(N_a)) + 2.
```

When nonempty, each primitive-embedding/period component has complex
dimension

```text
20-r.
```

Thus the required association really does replace a rigid chamber by a
usually positive-dimensional surface family.  It uses no named surface,
explicit equation, or surface automorphism.

For a bounded exhaustive `A(P)`, this gives the exact polytope-derived
criterion

```text
P has a weak full-Picard K3 carrier
  <=> some a in A(P) passes the discriminant-complement test.
```

Failure excludes K3 carriers of the exhaustively listed integral
realizations; it does not exclude the polytope or its universal rational
carrier.

## Level III: exact full-nef K3 association

Let `a` be a surviving K3 package and let `C_a` be the transported oriented
cone.  Define `M_a^nef(P)` by requiring the marking to identify `C_a` with
the complete nef cone.

### Necessary-and-sufficient criterion

```text
M_a^nef(P) is nonempty
  <=>
(1) N_a embeds primitively in Lambda_K3;
(2) every irredundant facet of C_a is a primitive square-minus-two root;
(3) C_a is a component of the complement of all square-minus-two
    hyperplanes of N_a—equivalently, no hidden root cuts its interior.
```

When these conditions hold, period surjectivity constructs `X`, the K3 Weyl
theorem moves its ample chamber to `C_a`, and every facet is represented by
an irreducible smooth rational curve.  Conversely every full-nef K3
association has exactly these properties.  Therefore, after exhaustive
normalization,

```text
P has an associated full-nef K3 surface
  <=> some a in A(P) passes embedding + primitive (-2) facets
      + complete-root separation.
```

This is the cleanest existing answer to the requested direction.  It is a
theorem about the surface associated **to the polytope**, and all its gates
are exact lattice/root conditions derived from `P`.  The primitive-embedding
gate is a finite decision; the no-hidden-root gate becomes a finite decision
only after a root-height bound, reflective classification, or other certified
stopping theorem is supplied.

## The first classification-sensitive middle association: Looijenga pairs

The gap between the universal weak rational carrier and the complete-root K3
carrier is no longer empty.  Suppose a rational arithmetic Coxeter datum `P`
passes a finite, entirely lattice-combinatorial **Looijenga-decoration**
compiler: an odd unimodular rational-surface saturation, characteristic
canonical vector, anticanonical cycle, and explicit toric blowdown certificate
are recovered from its walls.  Fix one resulting decoration, put

```text
Lambda = D^perp,
T' = Hom(Lambda,C*),
```

and let `R0` be the lattice generated by the internal nodal walls requested by
`P`.  If `Phi` is the deformation-root set independently computed from the
generic Looijenga cone `C_D++`, define

```text
T_R0^o = {phi in T' :
          phi|R0=1 and phi(beta)!=1 for every beta in Phi\R0}.
```

For `Q=Lambda/R0`, this locus is nonempty exactly when `Q` has a character
nontrivial on every forbidden-root image.  This is a finite character check
on the torsion of `Q` followed, on a surviving component, by deletion of
proper root hypertori.  When `R0` is primitive it reduces to the transparent
condition that no forbidden root is forced into `R0`.

Gross--Hacking--Keel/Friedman period surjectivity then produces a marked
Looijenga pair at every surviving character.  A choice of chamber in the
nonseparated lifted period space makes the intended simple roots themselves
effective.  Finally, and independently, one tests

```text
C_P = C_D++ intersect the selected simple-root halfspaces.
```

This last equality exhausts all numerical exceptional, negative-boundary,
and hidden nodal walls.  Thus, for each fixed decoration,

```text
P has a marked full-nef Looijenga carrier
  <=> exact character stratum nonempty
      + intended Weyl branch exists
      + complete C_D++ equality.
```

After a bounded exhaustive decoration enumeration, the disjunction of these
tests is an exact `P -> surface` criterion for the whole bounded Coxeter class.
It is already realized by an eight-facet arithmetic polytope in `H4`.  Its
walls reconstruct the blowup of `P2` at four boundary points, the exact carrier
stratum is `C*\{1}`, the eight curves are the complete nef walls, and exact
cone enumeration gives five finite and four ideal vertices.  Hence it is a
finite-volume noncompact full-nef Coxeter polytope carried by a genuinely
one-dimensional surface family.  This is a project-new inverse-moduli
certificate; no worldwide novelty claim for the abstract `H4` polytope is
made.  See
[`looijenga_h4_inverse_carrier.md`](looijenga_h4_inverse_carrier.md).

The decoration is no longer an input in the bounded class
`C_triangle(4,8)`.  Starting only from primitive full-span wall vectors in
`I_(1,4)`, with at most eight square-minus-one/two walls, a finite compiler
recovers the canonical vector by adjunction and enumerates every
anticanonical triangle, blowdown model, period/root branch, and complete-nef
inequality.  Carrier nonemptiness is equivalent to one of exactly four
weighted Gram types `T6a,T6b,T7,T8`.  They give four finite-volume
noncompact `H4` polytopes, three with project-new full-nef carriers.  Exact
lineage identifies those three abstract chambers with known Tumarkin/Im Hof
types, while a new congruence theorem realizes all four groups as mod-two
preimages of indices `10,20,30,60` in `W(I_(1,4))`.  This is
an exact variable-class carrier theorem, but not a universality theorem for
all algebraic carrier species; a rejected input may still be an abstract
Coxeter polytope.  See
[`triangular_looijenga_h4_compiler.md`](triangular_looijenga_h4_compiler.md)
and
[`triangular_looijenga_h4_worldwide_lineage.md`](triangular_looijenga_h4_worldwide_lineage.md).

## The mixed-wall middle association: one-boundary Coble recognition

There is now a second exact middle theorem, with a deliberately narrower
scope.  Start only from an oriented integral rank-11 Coxeter chamber in
`I_(1,10)`.  The Coble compiler:

1. enumerates characteristic norm-minus-one facets `b`;
2. requires `b^perp` to be even unimodular of signature `(1,9)`, hence `E10`;
3. partitions the remaining walls into twelve norm-minus-two roots
   orthogonal to `b` and norm-minus-one classes `E` with `b.E=1`;
4. requires the primitive isotropic classes `f=E+b` to be exactly all ideal
   extreme rays of the complete root cone in `b^perp`; and
5. compares that complete root chamber with Kondō Type I or Type II.

Within the classified species of complex finite-automorphism one-boundary
Coble surfaces, these conditions are necessary and sufficient for the marked
full-nef carrier moduli to be nonempty.  This is `P -> surface`: the boundary,
Coble--Mukai lattice, roots, and cusp lifts are reconstructed from the
polytope rather than imported from a chosen surface.

Applied to the pre-existing `D16-D5-H10-30` wall certificate, the compiler
finds a unique characteristic facet, recovers `E10`, twelve roots, and all
seventeen cusp lifts, and identifies the Type-II moduli family.  This is a new
strong carrier realization of a known polytope, not a new abstract chamber.
Failure of the compiler obstructs only this Coble species.

An independent no-import implementation scans every available local
odd-unimodular `H10` chamber: six chamber types, 24 integral presentations,
and 154 characteristic candidates.  Only `P10/12` and `D16-D5-H10-30` pass
the raw incidence and cusp-regeneration gates; Kondō's classification excludes
the former only as a one-boundary finite-automorphism Coble carrier.  The
latter is the unique positive in this frozen denominator, with unique `E10`
saturation up to its root-diagram automorphism group.  See
[`polytope_first_coble_h10_hostile_audit.md`](polytope_first_coble_h10_hostile_audit.md).

Separately, in the opposite `surface -> P` direction, Kondō Types I and II
produce complete full-Neron--Severi finite-volume noncompact `H10` chambers
with `22` and `30` facets.  That positive construction must not be used as the
proof of the preceding inverse statement.  The exact split theorem and replay
are in
[`coble_one_boundary_full_ns_h10_theorem.md`](coble_one_boundary_full_ns_h10_theorem.md).

The surface-first calibration now also reaches the two-boundary Type-I Coble
surface.  Its complete full-`NS` nef cone is a 26-facet finite-volume
noncompact chamber in `H11`, with exact rays `508+42+0`, 3,641 good Vinberg
edges, and 106 good affine-neighbourhood audits.  Its arithmetic denominator
is `2 Pic~=U+E8(-2)+<-2>^2`.  This is a new project-level complete
realization, not a proof that a raw `H11/26` polytope reconstructs a Coble
surface and not evidence for broad `P -> surface` universality.  See the
[two-boundary Coble theorem](coble_type_i_two_boundary_full_ns_h11.md).

## Surface-only necessary-and-sufficient existence criterion

There is an even stronger formulation in the complete-root class which does
not assume that a finite-volume polytope has already been supplied.  Start
only with an even lattice `N` of signature `(1,n)`, `2<=n<=19`, and define
`M_N^K3` to be its exact-Picard marked K3 moduli problem.  Then

```text
there exists a finite-volume full-nef K3 Coxeter polytope with lattice N
  <=>
N embeds primitively in Lambda_K3
and Cox(X) is finitely generated for one (equivalently every)
exact-Picard carrier X in M_N^K3.
```

Equivalently, after the embedding gate,

```text
finite-volume complete (-2)-root chamber
  <=> finite Aut(X)
  <=> rational-polyhedral Eff(X)/Nef(X)
  <=> finitely generated Cox(X).
```

The left side is not placed in the definition of the surface locus.  Period
surjectivity first constructs `X` from `N`; the algebra of `X` then decides
whether its full root chamber is a finite-volume polytope.  When the criterion
passes, the finite irredundant rational-curve walls of `Nef(X)` extract the
polytope.  When it fails, no finite-volume complete-root chamber exists for
that lattice, although mixed-root or non-nef Coxeter polytopes may still do
so.

This is precisely a surface-classification-to-polytope necessary-and-
sufficient criterion.  In ranks 11--20 its classification has already been
carried out: it yields exactly nineteen finite-volume noncompact chamber types
in ranks 11--19 and none in rank 20.  Those nineteen polytopes were known
before; the new project content is the explicit association criterion,
classification transfer, and literal all-row verification.  The next goal is
to obtain an analogous criterion for a broader carrier category where the
output includes genuinely new polytope types.

## What the associated surface tells us about P

Once Level III holds, classical surface properties become equivalent to
polytope properties:

```text
finite hyperbolic volume
  <=> finite index of the (-2)-reflection group
  <=> finite Aut(X)
  <=> rational-polyhedral Eff(X)/Nef(X)
  <=> finitely generated Cox(X).
```

Moreover,

```text
P compact
  <=> N_a has no rational isotropic ray
  <=> X has no genus-one fibration.
```

Meyer's theorem therefore rules out compact full-nef K3 associations in
hyperbolic dimension at least four.  These are necessary-and-sufficient
criteria **inside the full-nef K3 association class**.  They become no-go
theorems for a broader abstract Coxeter class only after a universality
theorem proves that every member of that class must reach Level III.

The compactness obstruction is not K3-specific.  For the entire direct
full-nef surface class, integral facet reflections and finite volume make the
reflection group finite-index in `O^+(NS)`; compactness is rational
anisotropy, so Meyer again gives `n<=3`.  The bound is sharp at the K3 carrier
`S_(4,28)`.  Rational surfaces are stronger: a pulled-back ruling fibre gives
a nef isotropic class in every Picard rank at least two, so their full nef
projectivizations are never compact.  See
[`direct_full_nef_surface_compactness_obstruction.md`](direct_full_nef_surface_compactness_obstruction.md).

The general convex--Hodge theorem gives the forward volume implication in an
even broader form: on any smooth projective surface of Picard rank at least
three, a complete rational-polyhedral nef cone has finite hyperbolic volume
and is compact exactly when it has no isotropic extremal ray.  This is an
auxiliary consequence after strong association, not the main direction of
the present theorem.

## Totally-real K3 and hyperkähler eigencarriers

When `d>1`, the visible `H^n` is normally an eigenspace inside the full trace
lattice rather than the entire Picard space.  For a deformation type with
Beauville lattice `Lambda` and second Betti number `b_2`, the exact weak
association remains

```text
trace lattice N_a embeds primitively in Lambda
  <=> nonempty marked full-Picard eigencarrier moduli.
```

The period dimension is

```text
b_2 - d(n+1) - 2,
```

and necessarily

```text
n <= floor((b_2-2)/d) - 1.
```

Strong association additionally requires equality with the complete visible
nef/MBM root arrangement and an ample eigenframe.  The finite visible
wall-norm shells are exactly enumerable, so this is a discrete promotion
problem after the embedding gate.  In the admissible totally-real case the
field form is anisotropic: a complete invariant reflection chamber has
finite volume exactly when it is compact and its arithmetic chamber
stabilizer is finite.  Hence a finite-volume noncompact arithmetic polytope
cannot be the complete chamber of this particular eigencarrier mechanism.

The existing rank ladder is:

```text
golden H8  -> K3 trace rank 18 -> compact strong example;
silver H9  -> K3 trace rank 20 -> weak carrier; complete intrinsic
             field-norm-two chamber infinite-volume, but the full restricted
             K3 (-2)-wall arrangement remains open;
golden H10 -> OG10 rank 22     -> weak carrier, infinite complete chamber.
```

The qualifier on the silver row is essential.  Saturation adds a K3
`(-2)`-divisor of field norm `3-sqrt(2)` whose reflection is not integral over
`Z[sqrt(2)]`.  It is therefore absent from the intrinsic field-norm-two
reflection arrangement to which the stabilizer theorem applies, but may cut
the geometric eigen-nef section.  The tested Alexeev local-star restrictions
are non-Coxeter; neither those tests nor the intrinsic volume theorem close
the complete restricted K3 `(-2)` arrangement.

There is also a polytope-first ordinary OG10 calibration at the maximal
rank.  The raw 210-wall Gram of the classical Borcherds `D4` `H21` chamber
reconstructs `L=U+D20(-1)`.  Exhaustive discriminant gluing gives one
primitive embedding orbit in the OG10 lattice, so the exact-Picard period
locus is nonempty and all 210 facets are OG10 wall divisors.  Full-cone
promotion nevertheless fails: the additional primitive
`v=(-1,1,0^20)` has type `(-4,div1)` and its hyperplane contains a point
strictly interior to all 210 input halfspaces.  Thus the proposed full-nef
moduli locus is empty.  This is an exact one-datum promotion criterion, not a
polytope no-go or the required uniform inverse theorem.  See the
[OG10 carrier obstruction](borcherds_d4_h21_og10_carrier_obstruction.md).

## Uniform polytope-first automorphic arrangement carrier

There is now a broader AG association that does not begin by selecting a
surface deformation type.  Let `P` be compact arithmetic and let its even
primitive totally-real wall lattice `L_P` have signature `(1,n)` at one real
place and be negative definite at the others.  The coordinate-free
determinant line has quadratic form `2*(-1)^n*det(b_L)`, hence the unique sign
pattern needed to form a `(2,n)` orthogonal Shimura datum.  Its
reflection-stable even overlattices are the finite list of stable isotropic
subgroups of the discriminant module.

For every output, a canonical real slice pulls the special divisor of a
facet root back to exactly that facet mirror.  Negative-definite divisor
intersections give elliptic candidate spans; oriented halfspace feasibility
recovers the actual faces.  The full arithmetic DM stack retains still more:
the local inertia at a two-mirror crossing contains `I_2(m_ij)`, so its order
`2m_ij` recovers each finite Coxeter label.  Parallel/ultraparallel data are
retained by the marked real-slice Gram.  Compactness makes the extra mirrors in
each fixed Heegner shell finitely enumerable.  This is an unconditional,
uniform `P ->` projective AG-arrangement theorem.

On a neat cover, a general complete-intersection surface intersects the
special divisors in curves determined by `P`.  This records incidence and
orbifold crossing data, not an isometric Lorentzian Gram or automatic full
nef cone.  For a stronger Whittaker/Borcherds carrier, the prescribed
principal part exists if and only if its coefficients annihilate every
vector-valued Hilbert cusp form under Bruinier's finite pairing.  This iff is
carrier-specific: automorphic correction is not automatic for reflective
lattices.  The compact Bugaenko calibration makes this failure concrete:
its eleven facets give one standard Heegner column, the column also has 309
lower-face contacts and no interior cut, and the exact weight-`(5,5)`
Poincare estimate proves that the pure norm-two principal part is obstructed.
The complete globally minimal compensator shell, of norm `1/80`, is a
nonempty divisor wholly disjoint from the real slice; nevertheless its
target-plus-shell `2x2` Poincare matrix is invertible.  Thus the complete
first support is obstructed too.  The determinant safe-span theorem now
proves that some larger finite slice-disjoint correction exists; with the
stated Hilbert real/integral q-structure it is nonnegative/effective.  What
remains open is its least support and every intrinsic height, level, weight,
or degree bound.

Consequently this construction has not yet produced a compact dimension or
facet bound.  A dimension theorem requires either a uniform quantitative
bound on the automatic correction on a non-singleton compact class or a jointly exhaustive
non-Borcherds carrier branch, together with a uniform high-dimensional
invertibility estimate for every allowed cusp-pairing matrix.  A facet theorem further
requires a sheet/degree inequality converting Heegner or orbifold Chern data
into the number of walls of one real chamber.  See the
[`determinant-line transfer theorem`](polytope_to_automorphic_divisor_transfer_theorem.md),
the [`Bugaenko Whittaker gate`](bugaenko_h8_whittaker_gate.md), and the
[`first compensator theorem`](bugaenko_h8_whittaker_first_compensator.md).

## Exact decision interface

For a fixed exact arithmetic realization of `P`:

1. normalize the oriented wall roots and verify field/integral reflections;
2. form every admitted even trace/wall lattice and reflection-stable
   saturation;
3. compute signature, Smith data, and finite quadratic forms;
4. decide primitive K3/Beauville embeddings by complementary discriminant
   forms;
5. construct the marked period components;
6. test primitive `(-2)` facets and exhaust hidden roots for direct K3, or
   enumerate complete visible root/MBM shells for an eigencarrier, retaining
   an explicit stopping certificate for either enumeration;
7. extract surface properties—Cox ring, automorphism group, isotropic nef
   rays, fibrations—and translate them through the equivalences above; and
8. retain a literal Gram/cone audit to verify labelled recovery of `P`.

For a bare diagram or bounded arithmetic class, prepend step zero: prove a
finite exhaustive enumeration of all positive integral normalizations and
glue packages.  Without step zero and a terminating completeness certificate
at step six, failure is only package-relative rather than a decision for the
abstract polytope.

## The actual open middle problem

The current endpoints and first middle bridge are now clear:

```text
universal rational blow-up association
    = broad but too weak to constrain P;

complete-root K3 association
    = exact and classification-sensitive but too narrow.

fixed-decoration Looijenga association
    = exact character-theoretic middle theorem, with bounded decoration and
      complete-exceptional-wall enumeration still required.

one-boundary Coble association
    = exact mixed-norm recognition theorem inside a classified carrier
      species; exhaustive multi-species universality still required.
```

The main algebraic research objective is to construct a carrier category
between them for which:

1. every arithmetic polytope in a nontrivial surface-free class has an
   associated surface produced from its Coxeter data;
2. effectivity, irreducibility, and complete-cone equality have independently
   decidable algebro-geometric criteria;
3. the surface moduli are softer and classifiable; and
4. their classification yields a necessary-and-sufficient existence
   criterion for the original polytope class.

The first bounded Looijenga theorem, its worked example, and the exhaustive
variable triangular carrier class `C_triangle(4,8)` are now proved.  The
next targets are higher-rank bounded Looijenga classes, a jointly universal
finite carrier stack, full-Neron--Severi Coble surfaces, and bounded
Mori-dream rational surfaces.
The next broad theorem must be a **polytope-to-surface universality theorem**
for one such class, not merely another surface-first example.

## Frozen sources

- [`rational_surface_eigencarrier_universality.md`](rational_surface_eigencarrier_universality.md)
  proves Level I constructively.
- [`k3_weak_algebraic_carrier_theorem_audit.md`](k3_weak_algebraic_carrier_theorem_audit.md)
  proves Level II and the exact local embedding decision.
- [`k3_strong_inverse_carrier_transfer_theorem.md`](k3_strong_inverse_carrier_transfer_theorem.md)
  proves Level III and the surface-property dictionary.
- [`totally_real_eigencarrier_criterion.md`](totally_real_eigencarrier_criterion.md)
  proves the totally-real/hyperkähler extension.
- [`projective_surface_nef_hyperbolic_volume_theorem.md`](projective_surface_nef_hyperbolic_volume_theorem.md)
  supplies only the auxiliary strong-carrier volume implication.
- [`looijenga_h4_inverse_carrier.md`](looijenga_h4_inverse_carrier.md)
  proves the fixed-decoration Looijenga character criterion and its complete
  full-nef arithmetic `H4` calibration.
- [`triangular_looijenga_h4_compiler.md`](triangular_looijenga_h4_compiler.md)
  removes the supplied decoration, exhausts `C_triangle(4,8)`, and classifies
  its four positive full-nef carrier Gram types.
- [`triangular_looijenga_h4_worldwide_lineage.md`](triangular_looijenga_h4_worldwide_lineage.md)
  identifies the known abstract ancestors and proves the four-group mod-two
  congruence-preimage hierarchy.
- [`polytope_to_automorphic_divisor_transfer_theorem.md`](polytope_to_automorphic_divisor_transfer_theorem.md)
  proves the canonical determinant-line Shimura-arrangement carrier and
  isolates the exact Whittaker gate and the three missing bound lemmas.
- [`coble_one_boundary_full_ns_h10_theorem.md`](coble_one_boundary_full_ns_h10_theorem.md)
  proves the one-boundary Coble recognition criterion, its polytope-first
  reconstruction of the known `D16-D5-H10-30` chamber, and the separately
  booked surface-first Type-I/II full-nef chambers.
- [`polytope_first_coble_h10_hostile_audit.md`](polytope_first_coble_h10_hostile_audit.md)
  independently replays the raw `D16-D5` reconstruction, proves integral
  saturation uniqueness, and exhausts the six-type local `H10` denominator.
