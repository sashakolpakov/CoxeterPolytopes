# Cross-coface cap capacities and affine corner payments

This successor preserves the
[previous checkpoint](joint_port_facets_and_cusp_fibers_20261017.md).
It adds a geometric arithmetic constraint that excludes the preceding
restricted witness, an all-rank quadrilateral result for affine-A corners,
and an exact obstruction to one specific construction patch. The complete
internally checked unrestricted dimension-163 proof is unchanged; its
[status](dimension163_proof_status_20261017.md) is not merely numerical
evidence. The compact arithmetic ceiling21, all-simple-cusp ceiling61 and
compact-facet ceiling127 also remain unchanged. No new example is claimed.

## The arithmetic witness fails a cross-coface wall constraint

For a fixed maximal simplex S, several isolated A1 components of its star
give orthogonal immediate cofaces. A selected prism cap at one such root q
has normalized wall normal

    b = v - sqrt(tau) q/2,       (v,v) = 2 - tau/2.

Here complete contacts determine v in the nonsingular Lorentzian core.
For 4 < tau < 6, a second identical profile at a different isolated root
q' would have wall product (b,b')=(v,v) strictly between -1 and0. That
product is forbidden for distinct norm-2 Coxeter walls. At the same q the
normal is unique. Thus the whole selected profile can occur at most once
across all these cofaces of S.

The [new capacities](fresh_sqrt5_cross_coface_cap_capacity_20261017.md)
restrict to trivial core symmetry and raw stars with precisely one
internal A1 extra after base deletion. The exact census supplies241
columns,70 with positive load at the old masses, of which53 are violated.
All omitted cap types receive zero positive coefficient. No full cap
classification is asserted.

The simplest violating S-profile has nine isolated A1 cofaces and
tau=(7+sqrt(5))/2. Its old selected prism mass is7551/10^10 while its
simplex mass is839/10^10: the new capacity is violated by a factor of nine.
The [independent witness checker](review_sqrt5_cross_coface_full_base_witness_20261017.md)
reconstructs its core inverse, projection norm, automorphisms, actual
provenance and base map. It distinguishes the prism's other base type.
The former12.34494473291 ratio therefore does not certify a dual witness
for the enlarged model.

The [independent whole-family checker](fresh_sqrt5_cross_coface_cap_capacity_20261017_review.md)
reconstructs every one of the241 mappings and projection norms, including
all13 relevant core automorphism groups and their spherical deletions.
It confirms all70 positive loads and53 violations without importing the
producer or invoking an optimizer.

The [exact repair](fresh_sqrt5_cross_coface_cap_repair_20261017.md)
raises53 already positive simplex masses while keeping all355 selected
quad masses fixed. It rebuilds every inherited column, revalidates all593
local packing covers and adds240 explicit singleton terms. It also
recomputes every affected three-, four- and five-port charge, checking
the whole saved family in each case. The repaired ratio is
**10.247084111225892**, below12; a simpler uniform repair gives10.23143381779.
These are feasible dual witnesses, not optima. A low repaired witness
does not establish a feasible primal or a dimension theorem. The
[separate repair review](fresh_sqrt5_cross_coface_cap_repair_20261017_review.md)
checks the monotone transfers, changed local covers and common joint
provenance without claiming a second complete arithmetic replay.

One genuinely enlarged [finite-model solve](fresh_sqrt5_cross_coface_cap_full_model_20261017.md)
retains all82,806 face obligations,6,558 vertex budgets and428,310 old
variables. It adds the241 cap transfers and217 saved positive cuts:
31 local-packing cuts,24 triple cuts,47 four-port cuts and115 five-port cuts.
The resulting exact optimum bracket is

\[
\frac{171523868647821690}{13513502427307591}
\le T_*\le
\frac{3173194699319}{250000000000},
\]

or **12.69277669282928 <= T_* <= 12.692778797276**. Both endpoints
remain above12. The former finite-model bracket was near12.79397; the
new constraints improve it but a different obstruction survives.

The one bounded optimizer invocation takes69.153 seconds. Its rational
primal and dual recovery checks all retained matrix obligations; fresh saved
replay reconstructs the matrix and verifies both certificates. A
[separate source review](fresh_sqrt5_cross_coface_cap_full_model_20261017_review.md)
checks preservation, positive rational scaling, cap provenance and
monotone signed-column recovery. The new dual has3,447 positive face
masses and770 positive vertex masses. No complete local-cover certificate
or complete joint scan from the old masses is silently reused at this
new optimum: this solve includes exactly the stated finite cut family.
The geometric face/forest/field coverage is still partial, so even a
target12 primal here would not by itself be a dimension theorem.

## Why the preceding individual-coface search is exhausted

The [labelled scan](fresh_sqrt5_labelled_coface_base_packing_20261017.md)
retains simultaneous coface/subset orbits and the full spectator-copy
multiplicities. It checks164,219 floor-capacity columns with no violation.
Its complete core symmetry census has603 automorphisms on287 cores.

The [abstract packing-cone certificate](fresh_simplex_base_packing_cone_20261017.md)
then gives every nontrivial facet for pairwise-disjoint bases on two
through five vertices:1,1,8 and294 facets respectively. The five-vertex
case reduces to294 half-integral vertices of the Petersen fractional
stable-set polytope. Every facet's affine rank and all allowed packings
are checked exactly. The
[independent implementation](review_simplex_base_packing_cone_20261017.md)
reconstructs packings from set partitions and obtains the same vertices
by active-matrix rank, without using the producer's odd-cycle predicate.

Applying these complete cones to all34,520 exported coface slices checks
[510,656 inequalities](fresh_sqrt5_complete_coface_base_cones_20261017.md)
without a violation. There are870 positive slices and65,069 positive-slice
inequality checks. The separate
[review](fresh_sqrt5_complete_coface_base_cones_20261017_review.md)
states exactly which geometry, normalization and computations it covers.

In fact, the [singleton-cover certificate](fresh_sqrt5_singleton_coface_covers_20261017.md)
shows something stronger: every exported slice has total demand no larger
than its capacity. Taking each demanded base by itself, plus an empty
configuration, realizes its marginal vector. Additional within-slice
incompatibilities retaining each individual option cannot exclude this
particular witness. These are separate averaged slices, not a simultaneous
ambient realization. Coupling different cofaces is what makes the new
cap constraint effective.

## A complete all-rank affine-A quadrilateral case

The [mixed weight interface](fresh_simple_cusp_vinberg_weight_frontier_20261017.md)
extends the finite-vertex rule to large classical affine vertex diagrams.
Its total cost is at most3n per vertex; triangles are paid. Quadrilateral
payment is the outstanding step, so this interface alone is not a new
dimension theorem.

For an actual quadrilateral with an affine-A cycle corner, the
[new proof](fresh_simple_cusp_affine_a_quadrilateral_20261017.md)
settles that step in every n>=14. If the marked cycle roots are at
distance at most3, this corner is paid. Otherwise both arc lengths P,Q
are at least4. Classification leaves eight possible replacement contact
profiles at each neighboring corner. Small affine subdiagrams, forbidden
by the inherited all-simple-cusp exclusion, give uniform height bounds.

Under the assumption that the opposite marked roots are orthogonal,
eliminating the white paths gives a four-by-four Schur determinant. Its
normalized value is jointly concave in the two reciprocal heights.
All256 rectangle corners are nonnegative at1/P=1/Q=1/4, and strict
monotonicity makes the determinant positive at every actual P+Q>=14.
This contradicts the star's necessary corank. The opposite marked roots
must therefore be adjacent, paying that corner.

The [independent review](review_simple_cusp_affine_a_quadrilateral_20261017.md)
checks the classification, small-affine embeddings, height bounds, the
arc-length-four overlap, and the all-rank argument. A separate Fraction-pair
implementation recomputes all256 rectangle values and four positive kernels.
The [long tree-corridor lemma](fresh_simple_cusp_affine_tree_corridor_20261017.md)
also settles affine B/C/D corners when both marked roots are at least
distance4 apart and the two outer B/D blocks each have rank at least9.
The full finite/affine classification leaves four replacement profiles.
The all-rank inverse port diagonal of either endblock is1, and the Schur
formula then has a uniform positive gap at least1/2. A separate
[mathematical/source review](fresh_simple_cusp_affine_tree_corridor_20261017_review.md)
checks the classification, short affine embeddings and all-rank algebra.

The decorated-end tree cases and the full mixed quadrilateral theorem
remain open. In particular the inherited global exclusion forbids affine
orders below n; it must not be replaced by the assertion that every affine
subdiagram is forbidden.

## Construction: an all-label obstruction for twenty rings

The [central-patch test](fresh_nonarithmetic_release47_central_patch_20261017.md)
finds598 numerically coherent five-cap collections, all with33 finite
vertices, seven exterior intersections and15 target-face open ends. It
does not improve the frozen exact64-wall H9 control.

The [exact successor](fresh_nonarithmetic_central_ring_contact_gate_20261017.md)
reconstructs20 literal rings from their norm/contact recipes and rules out
one central closing cap for every finite contact label. The seven allowed
A2 profiles and all label tuples in{2,3,4,5}^4 give35,840 inverse-Gram norm
equations, each separated strictly from zero by rational intervals.
A label at least6 must be isolated within every required finite corner.
This leaves64 real branches determined by the other three contacts; every
purported high-contact product is strictly above -sqrt(3).

Generation and fresh saved replay pass. The
[source and proof review](fresh_nonarithmetic_central_ring_contact_gate_20261017_review.md)
checks the complete high-label reduction and the distinction between exact
recipe objects and ambient contact acceptance. The result does not exclude
other rings or several closing facets. Increasing only the central label
cutoff cannot help these20 rings; the next construction must change their
geometry or use more than one closing wall.

## Reproduction and archive

The producer packages use generation followed by fresh saved-data replay.
Independent arithmetic implementations and mathematical/source reviews
state their separate scopes; an internal review is not external refereeing
or proof-assistant formalization. Frozen predecessor files are preserved.

The [successor auditor](../scripts/audit_cross_coface_caps_and_affine_corner_payments_20261017.py)
and [saved audit](cross_coface_caps_and_affine_corner_payments_20261017_audit.json)
bind the new packages, recursively check input manifests including prepared
model sources, preserve the preceding archive's frozen hashes, check report
links and record the current README and plan. This is archive integrity,
not a substitute for the geometric arguments or exact certificates.
Preserved predecessor files are carried into the new cumulative hash map,
including files no longer reached by a current package's dependency graph.

```sh
python3 scripts/audit_cross_coface_caps_and_affine_corner_payments_20261017.py --verify
```
