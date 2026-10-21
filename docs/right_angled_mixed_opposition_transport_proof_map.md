# Opposition transport manuscript: proof and reproduction map

18 October 2026. New, unfrozen manuscript source:
[right_angled_mixed_opposition_transport.tex](../manuscripts/right_angled_mixed_opposition_transport.tex).
This map records a translation of the reviewed packets below. Writing the
translation did not run a checker, optimizer, Lean process or TeX build.
The coordinating agent owns the subsequent review and build queue.

This is the preserved predecessor's scope. The subsequent
[positive-part theorem](right_angled_dimension11_proof_map_20261019.md)
now proves the unrestricted mixed ceiling **11**, with independent
written reviews, exact generation/replay and conditional count-kernel
verification. That new theorem is not silently part of this source.
When this predecessor was written, the ceiling remained12. The manuscript
proves dimension **at most 8** under either of two additional hypotheses:
all three-faces have facet-plus-cusp sum at most ten; or there is a finite
vertex, no finite–finite edge anywhere, and all three-faces have sum at
most eleven. The latter is a global edge hypothesis. A single corner or
opposition component with ideal incident edges does not suffice.

## Mathematical coverage

|Manuscript part|Claim and hypotheses|Authoritative packet and separate review|
|---|---|---|
|Actual walls and small faces|Acute projection, distinct-normal gap, wall/triangle equivalence, actual P3 opposite, attachment, complete sum≤10 combinatorics, R7 finite-edge length|[weighted triangle proof](../results/right_angled_mixed_weighted_triangles_research_20261018.md), [low-facet rigidity](../results/right_angled_mixed_low_facet_rigidity_research_20261018.md), [independent rigidity review](../results/right_angled_mixed_low_facet_rigidity_independent_review_20261018.md), [nu10 transitions](../results/right_angled_mixed_nu10_edge_transition_research_20261018.md)|
|First conditional dimension theorem|All three-faces sum≤10 implies n≤8; actual undirected finite-edge/P3 moves with weights 1,1/20; complete causal tail and nine radical base rows|[small-three-face proof](../results/right_angled_mixed_small_three_face_dimension_research_20261018.md), [root written review](../results/review_right_angled_small_three_face_dimension_20261018.md)|
|Carrier-excess fraction|Only n12: at least f/385 nonsmall finite vertices and E≥f/5390; 178 exact interval boxes, no sampling|[fraction proof](../results/right_angled_mixed_carrier_excess_fraction_research_20261018.md), [independent source/calculus review and fresh replay](../results/review_right_angled_carrier_excess_fraction_20261018.md)|
|Actual opposition support|Cross edges preserved, true-twin incoming capacity, invariant h²=n−9, K9-free exit in every nontrivial component for n9–12, dimension12 boundary and finite-edge budgets|[actual-support proof](../results/right_angled_mixed_actual_move_support_research_20261018.md), [root written review](../results/review_right_angled_actual_move_support_20261018.md)|
|Reciprocal-square boundary input|Only n12: row cap 516992/3125, using an invariant-null frame per component; average degree and boundary inequalities|[moment identities](../results/right_angled_mixed_opposition_moments_research_20261018.md), [power-two proof](../results/right_angled_mixed_opposition_power_two_research_20261018.md); inherited as reviewed input by the actual-support review|
|Actual symmetric quads|Actual cosh3 wall criterion; a Q22 plus a common bad neighbor forces the complete D3 carrier; invariant common-neighbor multiplicity|[symmetric-quad proof](../results/right_angled_mixed_symmetric_quad_charge_research_20261018.md), [root full written review](../results/review_right_angled_symmetric_quad_charge_20261018.md)|
|Second conditional dimension theorem|Global absence of finite–finite edges and sum≤11; every carrier P3 or D3; matching complement; actual quad moves counted once with weight1/2; causal margins for n9–12|Last section of [actual-support proof](../results/right_angled_mixed_actual_move_support_research_20261018.md) and its [review](../results/review_right_angled_actual_move_support_20261018.md)|
|Universal quad error and maximum height|Unrestricted mixed corners: sum M≤2E12+binom(n−1,2)eFF/5; companion charge; local norm restrictions at a height-maximizing finite vertex|[symmetric-quad proof](../results/right_angled_mixed_symmetric_quad_charge_research_20261018.md) and [review](../results/review_right_angled_symmetric_quad_charge_20261018.md)|

The TeX reproduces the wall and convex-hull arguments instead of treating
formal candidate vertices as actual ones. It also reproduces the
small-face and attachment proofs, rather than assuming a complete
three-polyhedron census. The R7 length proof uses the rectangular cusp
link and three actual opposite normal products.

The two uses of quad weights are different and explicit. The unrestricted
error field uses common-neighbor multiplicity d_ij/2, which is invariant
under reversal. The second conditional theorem counts each quad once
with weight 1/2. Several three-face carriers never become several
different opposite vertices.

## External and inherited assumptions

The external geometric inputs are the standard acute-wall property for
Coxeter polytopes, orthant finite links, rectangular three-dimensional
cusp links, simplicity at a finite right-angled vertex, convexity and
ordinary polytope face incidences. The sole cited primary background paper
is Potyagailo–Vinberg, *On right-angled reflection groups in hyperbolic
spaces*, Comment. Math. Helv. 80 (2005), 63–73,
[publisher PDF](https://ems.press/content/serial-article-files/43032).

The manuscript uses the current published ceiling12 as context, not as a
new conclusion. Its status and reference provenance are in
[current_bounds_and_verification.md](current_bounds_and_verification.md).
It makes no new priority claim, universal geometry formalization, or
external refereeing claim.

Two finite arithmetic interfaces remain explicit:

1. The nine rows in the manuscript's radical table are compared in
   Q(sqrt2), and the higher-k and higher-m tails are proved symbolically.
   The table contains all nine coefficient pairs, so these comparisons
   can also be checked directly by rational squaring.
2. The fraction lemma gives every rational enclosure formula, the
   square-root grid10^-7, lower bound for sqrt2, all sixteen norm
   intervals for k1–11, and the two degenerate endpoints. The exact
   checker verifies all178 inequalities. This is the finite interval
   arithmetic step in the proof, not a claim that all directions were
   numerically sampled. The pointwise comparison is made only where
   the inflated variance is nonnegative; enlargement of the domain
   comes afterward.

The moment cap in the TeX is 165.43744, the constant used in the reviewed
actual-support packet. The later
[direction-polynomial sharpening](../results/right_angled_mixed_opposition_quartic_direction_review_20261018.md)
is a separate written result and was not silently substituted into the
frozen boundary constants. No optimization formulation, clique-hull
coverage assertion or aggregate witness is needed for this manuscript.

## Exact replay entrypoints

These are **argv descriptions for the coordinated proof queue**, not
instructions to start concurrent processes. Before any replay or TeX
build, obey [AGENTS.md](../AGENTS.md): shared exclusive proof lock, one
process across agents, the40GiB disk/6GiB memory start gates and8GiB/3GiB
abort gates, bounded output, and the cleanup-agent trigger below20GiB.
The coordinator must apply the bounded process wrapper used in the
corresponding review. Do not use generation flags or replace frozen
certificates merely to replay them.

    python3 -B scripts/check_right_angled_mixed_nu10_edge_transition_research_20261018.py
    python3 -B scripts/check_right_angled_mixed_small_three_face_dimension_research_20261018.py
    python3 -B scripts/check_right_angled_mixed_carrier_excess_fraction_research_20261018.py
    python3 -B scripts/check_right_angled_mixed_opposition_moments_research_20261018.py
    python3 -B scripts/check_right_angled_mixed_opposition_power_two_research_20261018.py

Each entrypoint compares its reconstructed exact record with the matching
JSON under results/, without the generation flag. The mathematical claims
of the actual-support and symmetric-quad packets have written reviews
rather than new computational certificates. Their immutable source/review
sets are recorded in:

- [small-three-face freeze](../results/right_angled_mixed_small_three_face_dimension_research_20261018_freeze.json)
- [carrier-fraction freeze](../results/right_angled_mixed_carrier_excess_fraction_research_20261018_freeze.json)
- [carrier root-review freeze](../results/right_angled_carrier_excess_fraction_root_review_20261018_freeze.json)
- [actual-support root-review freeze](../results/right_angled_actual_move_support_root_review_20261018_freeze.json)
- [symmetric-quad root-review freeze](../results/right_angled_symmetric_quad_charge_root_review_20261018_freeze.json)

No Lean claim is attached to the transport geometry or to the two
conditional dimension-eight theorems. A later formalization must use
the serial signal-safe Lean batch runner; unrestricted lake build is
not a reproduction step.

## Illustration and translation checklist

The source contains five TikZ figures:

- P3's two actual finite vertices and three ideal vertices, distinguishing
  opposition from an edge;
- the R7 finite-edge star and the causal cancellation mechanism;
- exact transport of cross edges between corner graphs;
- a symmetric quadrilateral and its induced-path direction graph;
- the K4-joined-independent8 maximum-height defect example.

Each caption identifies whether its dots denote hyperbolic vertices or
incident directions, and whether it is an incidence drawing rather than
a metric projection. The examples of abstract graph data are explicitly
not asserted to be geometrically realizable.

The unresolved steps remain visible: the full hypothetical quad field
can be spacelike (K8 disjoint union K2 disjoint union K2); current error
charges do not force a contradiction in all dimension12 systems; a local
height maximum does not imply a bound on the unweighted global triangle
average; the proposed triangle-count56 subproblem is unresolved.

The TeX has a [complete root source review](../results/review_right_angled_mixed_opposition_transport_tex_root_20261018.md)
and a separate [complete independent translation review](../results/review_right_angled_mixed_opposition_transport_translation_20261018.md).
Both are written mathematical reviews, with no fresh arithmetic execution
attributed to them. The PDF build is queued because available RAM has been
below the mandatory6GiB start gate; no build success is claimed here.
