# Bounds near known examples and complete arithmetic distance domains

This checkpoint preserves the [preceding archive](continuous_families_facet_incidence_and_weight_barriers_2026_10_16.md)
and records completed arithmetic coverage, explicit failures of proposed
weight rules, and broader construction experiments. **No new finite-volume
example or dimension bound is established.** The existing compact arithmetic
bound-21 and unrestricted finite-volume bound-163 claims remain internally
checked proof candidates, with external validation and priority unresolved.
Noncompact arithmetic has the separate published bound n<=21, n!=20;
Borcherds attains21. Unrestricted compact retains the published bound29.

The research priority is now a ceiling within roughly two to five dimensions
of known examples: about23–26 for unrestricted finite volume and10–13 for
compact polytopes. These are objectives, not conjectured theorems. The
[method-barrier report](fresh_dimension_gap_method_barriers_20261016.md)
explains why the current scalar estimates cannot reach those targets merely
by optimizing their remaining coefficients.

## Compact arithmetic: new coverage finds a genuine missing face

The [global target-40 proposal](fresh_sqrt5_global_discount_20261016.md)
passes its restricted exact vertex and face tests, but does not cover all
necessary quadrilaterals. The [first-extension screen](fresh_sqrt5_global_b1_first_extensions_20261016.md)
checks4,565 rows against3,067,805 common spectator cases. Its underpaid rows
include incomplete cores, so a failed row alone is not a geometric counterexample.

The [finite first-extension census](fresh_sqrt5_global_b1_finite_first_20261016.md)
contains1,646 classes, including196 singular stopping configurations.
One gives an [exact B6 quadrilateral](fresh_sqrt5_global_b6_discount_witness_20261016.md):
four actual finite corners, both required embeddings and integral Cartan
coefficients. Its payment is1-delta, disproving full coverage by the proposed
positive uniform discount. The spectator extension has the intended local
rank21 data; no global compact realization is asserted. This failure does
not affect the older undiscounted bound-21 candidate.

The [initial dotted completion screen](fresh_sqrt5_global_b1_dotted_first_completions_20261016.md)
retains14 necessary arithmetic completions with no further touching walls.
That restricted screen is superseded for distance coverage by the complete
[minimal dotted parameter census](fresh_sqrt5_global_b1_dotted_parameters_20261016.md).
All82 structural classes are handled:81 with one dotted pair and one with
two. Exact physical and conjugate inequalities leave377 necessary seeds,
including42 singular seeds. The two-dot case has17 surviving ordered pairs;
it imposes no spurious requirement that square roots of dotted invariants
belong to the field.

A [separate arithmetic reviewer](fresh_sqrt5_global_b1_dotted_parameters_20261016_review.md)
reconstructs all82 classes and377 seeds using independent field arithmetic,
different gauges and larger integer search rectangles. It checks2,578
integral cycle products. Both generation and saved replay pass. This is
complete distance-domain coverage at the stated minimal cores, not complete
coverage of all further attachments or a valid new weighting rule.

The finite and dotted arbitrary-attachment queues remain active outside
this archive. Their pruning is conditional on adjacent discounts at most
one half, zero distance-two discount and nonnegative bonuses. Their graph
alphabet is specific to Q(sqrt5). The separate target-40 repair is also
unfinished and excluded.

## Why the bound strategy must change

The [small-component certificate](fresh_sqrt5_discount_structural_floor_20261016.md)
exhibits an allowed rank17 forest3D5+A2 with baseline25 against threshold24.
Discounts restricted to components of rank at least six cannot change it.
The rank13 forest D5+2D4 already needs a uniform adjacent discount at least
four fifths in the tested extension. A rank11 example also defeats an
adjacent-only approach even when those adjacent weights vanish. These are
obstructions to specified weight families, not dimension theorems.

For unrestricted finite volume, even removing the entire loss in the
current balanced narrow-pair estimate would only exclude dimensions146
and above, giving ceiling145. The older Prokhorov scalar lower estimate
remains negative through dimension66. Thus coefficient refinement within
these fixed estimates cannot reach the twenties. New information about
simultaneous face completion, actual framing walls, or a different selected
narrow face is needed. The exact [scalar checks](fresh_dimension_gap_method_barriers_20261016.md)
state the premises and distinguish these limitations from universal barriers.

## Constructions: broader contacts and coordinated finite faces

The exact H20 chamber remains at158 walls, at least184 finite vertices and
five proper cusps; a strict escape proves infinite volume. Completed
[near-closure recovery](fresh_h20_near_closed_cluster_recovery_20261016.md),
[one-anchor families](fresh_h20_one_anchor_dotted_families_20261016.md),
[old-edge closure](fresh_h20_open_edge_closure_20261016.md), and
[new-edge continuation](fresh_h20_seed_edge_families_20261016.md)
produce larger numerical banks and compatible clusters, but all retained
whole-face tests leave exterior boundary. Closing old open edges creates
new endpoint obligations.

The [spherical attachment census](fresh_h20_spherical_attachment_profiles_20261016.md)
gives57 fixed finite profiles and one unbounded I2(m) family. In particular,
an E8 attachment joins A4, A2 and A1 with residual norm1/30. A
[four-wall E8 recipe](fresh_h20_e8_triangle_cluster_20261016.md)
is numerically compatible with all158 old walls and creates a compact
triangular two-face, conditional on exact wall verification. Exact diagram
deletions predict four new ambient finite corners. Six rational free
pairings and prescribed quadratic branches specify a route to exact
reconstruction; that reconstruction is not claimed here. Its target H3
still has14 exterior rays.

The exact H9 chamber remains at55 walls and infinite volume. The
[edge-completion search](fresh_nonarithmetic_h3_edge_completion_20261016.md)
and [complete seven-profile lemma](fresh_nonarithmetic_h3_a2_profiles_20261016.md)
expand the available wall contacts. A
[recentered coupled-facet search](fresh_nonarithmetic_h3_edge749_20261016.md)
finds261 numerical two-cap proposals with a finite triangular new facet.
Each cuts only one initial open end; the whole H3 still has exterior rays.
These are bounded numerical discoveries, with the exact profile lemma
separately identified. They do not establish finite volume.

## Verification and preservation

Included producers and fresh saved replays have passed. The
[archive checker](../scripts/audit_near_example_targets_checkpoint_20261017.py)
and [saved audit](near_example_targets_and_complete_dotted_domains_2026_10_17_audit.json)
preserve prior frozen sections and bind recursive manifests, reports and
local links. This proves archive integrity, not the mathematical validity
of a dimension theorem. Active queues, the incomplete weight repair and
new20261017 escape-partition work are excluded.

```sh
python3 scripts/audit_near_example_targets_checkpoint_20261017.py --verify
```
