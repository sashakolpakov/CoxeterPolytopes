# Contact-cost dimension bound and bounded cusp links

**Checkpoint: 2026-10-15.** The nonexistence lane now has an internally reviewed research proof candidate **n≤169** for arbitrary finite-sided, finite-volume real hyperbolic Coxeter polytopes. The H9 construction gains a certified connection between finite corners. The H20 construction gains three bounded cusp links. Both chambers still have exact infinite-volume escapes; no new finite-volume example or dimension is established.

The [preceding checkpoint](avoid_node_bound_and_joint_caps_2026_10_15.md) and its sources remain frozen. The [archive checker](../scripts/audit_contact_cost_cusp_links_checkpoint_20261015.py) and [integrity record](contact_cost_bound_and_cusp_links_2026_10_15_audit.json) bind the new artifacts and preserve the predecessor. Archive integrity is distinct from mathematical verification.

## Dimension bounds and algorithm scope

The [169 proof binder](fresh_nikulin_unit_charge_20261015_dimension169.md) combines the new affine-factor inequality with the unchanged geometric and compact-facet arguments. The [primary estimate](fresh_nikulin_unit_charge_20261015.md) lowers the charge per unit Schur cost from two to four-thirds. The [independent review](fresh_nikulin_unit_charge_20261015_review.md) reconstructs all 40 exceptional cases, all admissible contact supports, complete exceptional omission obstructions and 24 symbolic classical index classes. Generation and saved replay pass.

With `N=n−1`, the two endpoint costs total at most four; the single factor containing the narrow-facet normal incurs an additional deficit 2/9. The resulting bound is

    −(x,y) ≥ M(N−43)/9,    M≥1.

For `n≥170`, the pairing is at least 14, contradicting Nikulin's strict narrow-facet bound. The alternative compact-facet case gives 143. The local coefficient four-thirds is optimal for this fixed-slope, preset-omission lemma; that does not establish optimality of the global method or dimension bound.

The [fresh literature and scope check](fresh_frontier_bound_scopes_20261015.md) keeps the hypotheses explicit:

| Hypotheses | Bound or status |
|---|---|
| Compact hyperbolic Coxeter polytope | Published bound 29; known compact examples reach dimension 8 |
| Arithmetic finite-volume Coxeter polytope | Published bound 29 |
| Noncompact arithmetic Coxeter polytope | Dimension at most 21, excluding dimension 20 |
| General finite-volume Coxeter polytope | Published Prokhorov bound 995; Nikulin announced below 300; repository candidate 169 is internally reviewed, without external validation or priority claim |
| Hyperbolic manifolds without the reflection/Coxeter condition | These reflection bounds do not apply |

Thus compact H9 and nonarithmetic finite-volume H20 remain distinct construction targets. Arithmetic lattice searches in the excluded noncompact dimensions are not pursued. Vinberg's wall enumeration and Allcock's published lattice edgewalking are also distinguished from the repository's geometric endpoint search: borrowing the edgewalking strategy does not supply a complete unrestricted nonarithmetic algorithm.

## Compact H9: an all-label family exclusion, with the incidence still open

The [compact continuation report](fresh_h9_tremaining_20261015.md) starts with three literal T-incidence predicate binders. Only T10 is searched in this phase. Its first necessary-label survivor occurs after 9,008 visits and has internal block labels `(4,4,>=7,4)`.

The [high-label theorem](fresh_h9_tremaining_20261015_case10_high_bound.json) and [independent review](fresh_h9_tremaining_20261015_case10_high_bound_review.json) exclude all values at least 12 in the unbounded category. The [remaining metric proof](fresh_h9_tremaining_20261015_case10_rank_gate.md) and [separate review](fresh_h9_tremaining_20261015_case10_rank_gate_review.md) cover labels 7–11. Four actual-face budgets reduce the necessary projection patterns to 100 cases; 500 exact Chebyshev identities leave four cases, each excluded by a two-negative-direction residual plane. Omitted meeting labels and dotted distances remain covered rather than assigned a numerical cutoff.

The reviewed gate is composed through literal incidence automorphisms into 153 predicates. Resuming the saved cursor takes 14 further visits and reaches a different necessary-label survivor, with labels `(4,4,>=7,5)`. The state is saved at **9,022 cumulative visits with ten pending frames**. T10 remains open. T19 and T29 have initial binders but are unsearched. The catalogue total stays **164 of 240 recipes excluded**, with 18 unresolved T representatives; this phase does not claim another whole-incidence exclusion.

## Geometric H9: an exact bridge between finite corners

The [H9 report](fresh_h9_joint_corner_20261015.md) starts from the independently certified 32-wall joint-corner chamber. Its nine outgoing directions are reconstructed exactly; bounded numerical scouts supply proposals but no complete cap classification.

The selected cap has norm two, label 3 at wall 30, label 5 at wall 8 and seven specified right contacts. It requires a genuine degree-eight extension of the existing quartic field. The [independent review](fresh_h9_joint_corner_20261015_bridge_review.md) uses a separate primitive-element representation and verifies all 528 pairs and all 33 actual facets.

All 90 old finite vertices survive. Exhausting the nine possible new finite-vertex supports gives precisely two new vertices, so the complete finite census is **92 vertices, 253 finite edges and 322 open finite-corner directions**. The number of components of the finite-edge graph drops from eight to seven. The new five-vertex component joins the formerly isolated joint corner to an old two-vertex component; the main 81-vertex component remains separate.

The strict escape survives, proving infinite volume. No complete new exterior-ray or ideal-vertex census is claimed. The 2,267-record scout and 4,566-record old-bank rescore remain explicitly bounded numerical diagnostics.

## H20: three actual bounded cusps and a uniform obstruction

The [H20 report](fresh_h20_new_projection_20261015.md) adds an exact cap with future timelike source projection of square `−15/2`, followed by three mutually orthogonal null-source caps. The [independent review](fresh_h20_new_projection_20261015_review.md) reconstructs all **47 facets and 1,081 pair relations** using separate rational multiquadratic arithmetic.

The first cap cuts the preceding escape and closes a finite edge. Its 20 possible finite-neighbor supports are exhausted; exactly one gives a finite vertex. Four stored finite vertices and two finite edges survive the whole extension. This is not a complete ambient H20 vertex census.

Each new null-source cap closes the full link at an explicit source cusp. Its 23 incident walls have four positive-kernel affine components of ranks `8,8,1,2`, totaling 19. The physical normal span has rank 20. The resulting three ideal vertices therefore have bounded full Euclidean links, not merely selected bounded sublinks.

A renewed strict spacelike escape nevertheless proves infinite volume. A new all-profile theorem also shows that every compatible cap whose source projection lies on any of the **31 literal source parallel-pair rays** retains that escape. It combines the first nonright Coxeter contact's lower scale bound with exact Cauchy–Schwarz margins. The independent reviewer checks all 31 rays and all 465 distinct-ray products. This exhausts those parallel-derived directions, not every source cusp or every possible projection.

## Remaining work

The dimension argument now warrants external mathematical scrutiny and broader literature comparison. A [complete balanced-height obstruction](fresh_nikulin_balanced_height_20261015.md), with [independent review](fresh_nikulin_balanced_height_20261015_review.md), rules out the proposed unit-charge shortcut even when both heights coincide. All 36 distinct E8 omission pairs have a legal common-row counterexample; these witnesses force a coefficient at least 8/7. No matching upper estimate is proved here. Shared height ratios or stronger global accounting may still improve the four-thirds bound.

The user's subsequent priority is sharper dimension bounds, especially for arithmetic groups of simplest type. The arithmetic target is to improve 29 in the compact case; the noncompact ceiling 21 is sharp. New lanes will test the high-dimensional field restrictions against minimal indefinite diagrams and global multi-face constraints. They are outside this frozen checkpoint. Allcock's rational-lattice algorithm does not evade the already sharp noncompact arithmetic theorem.

Continue T10 from the saved label-5 family in a fresh namespace, using the coupled residual-block structure rather than restarting its closed label-4 family. In geometric H9, the remaining escape and disconnected finite-edge graph call for different coupled caps or changed retained geometry. In H20, completing more of the 31 parallel-derived cusps alone cannot remove the proved escape; use other source projections or replace retained walls. The two proved nonreflective arithmetic H9 hosts remain closed.
