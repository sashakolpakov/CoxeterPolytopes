# Dimension29 for large cusp factors and forced quadrilateral incidence

This checkpoint completes an internally checked proof of dimension
**at most29 when every irreducible affine cusp factor has rank at least17**.
There is no arithmetic hypothesis or limit on the number of factors.
Equivalently, any finite-volume example in dimension at least30 must
have a cusp factor of rank at most16. The result is restricted: the
unrestricted163, compact-arithmetic21, all-simple29 and
minimum-cusp-factor-rank9 bound85 are unchanged. No new finite-volume
example is established.

The [complete proof](fresh_cusp_minimum_factor_rank17_dimension29_20261017.md)
combines all-rank geometric payments with an analytic cusp-curvature
comparison. Its [separate internal review](fresh_cusp_minimum_factor_rank17_dimension29_20261017_review.md)
checks the written composition and independently derives the scalar
identities using rational polynomial arithmetic. No finite rank cutoff
or inherited163 tail is needed for this new theorem. Internal proof
review and successful exact checks do not constitute external refereeing
or proof-assistant formalization.

## Actual adjacent faces force a new quadrilateral

The [raw23/raw25 proof](fresh_sqrt5_core23_25_forced_quad_20261017.md)
shows that two selected prism profiles sharing a maximal simplex force
a compact quadrilateral with a native5 corner of type A4+3H3. Exact
projection fixes the cap cross product; all four spherical corner
diagrams give actual feasible vertices by projection of a strict
interior point. The pointed quadrilateral determines at most one
originating maximal simplex, by disjointness of its vertex ports.
The [separate review](fresh_sqrt5_core23_25_forced_quad_20261017_review.md)
records the written geometry review separately from exact computation.

For a common simplex with k isolated A1 spectators and residual forest
T, the [full represented family](fresh_sqrt5_core23_25_spectator_incidence_20261017.md)
gives the necessary actual-incidence inequality

` k a + b - k s <= h v_Y,   Y=(k-2)A1+T+2H3,   h=#H3(Y). `

At most one A profile occurs, B profiles occupy distinct cofaces, and
mixed profiles cannot occupy the same coface. Their mixed-pair count is
at least kA+B-k. Each target vertex has h eligible native5 pairs.
The certificate checks65 literal profile maps in52 represented contexts.
Missing selected profiles receive coefficient zero; this is not a
completeness claim for all geometric faces.

The original saved witness violates the first constraint. The
[first monotone repair](fresh_sqrt5_core23_25_incidence_repair_20261017.md)
adds3363761/3000000000 to vertex mass5572. The52-context scan finds one
further violation and adds50071/5000000000 to vertex mass5413. Every face
mass remains unchanged. All earlier packing, cap, signed and joint
constraints therefore persist by their monotone vertex budgets; their
2,711,432 inequalities are not needlessly rescanned.

The final exact partial-model ratio is **12.11531601432504>12**, down from
12.12870037679824. It is a repaired obstruction, not an optimum. Full
face, field, forest and global incidence coverage remain necessary
before a compact arithmetic dimension conclusion.

## Complete fine payments and the restricted bound

The [large-factor quadrilateral proof](fresh_large_affine_factor_quad_payment_20261017.md)
extends the ideal-corner payment to a varying affine factor of order
at least18, assuming no connected affine subdiagram of order at most9.
It permits arbitrary spherical spectators. Its additional cases include
inner affine extensions, fixed E7/E8 roles, and six new E8 scalar rows;
the positive-semidefinite spectator correction preserves the strict
determinant contradictions. It does not treat every finite-corner face.

The [all-rank curvature comparison](fresh_cusp_link_rank_monotonicity_20261017.md)
proves that, when all cusp-factor ranks are at least17, increasing any
one rank raises normalized link curvature by at least7/76. For n>=30
the curvature satisfies kappa/V>=3s, with equality only for one factor
of rank29. This is analytic and has no rank-enumeration cutoff.
Together with the completed fine payments this implies n<=29.

The [private affine-cycle lemma](fresh_private_affine_cycle_exclusion_20261017.md)
and the new [private affine-tree proof](fresh_private_affine_tree_exclusion_20261017.md)
exclude affine private diagrams U,V in an unpaid finite quadrilateral
star when no affine subdiagram has order at most17. For a tree, exact
resistance inequalities reduce the shared contact to a terminal
triangle. A unique connector has inverse diagonal zero; its Schur
reduction forces a second private affine component. Joining the two
remaining finite B/D arms then produces an affine diagram inside an
actual finite corner, a contradiction. The proof is all-rank; the
saved finite probes corroborate its identities.

The [cycle carrier extension](fresh_affine_cycle_carrier_exclusion_20261017.md)
also excludes large affine A cycles inside H-U and H-V, without assuming
the cycle exhausts the private diagram. Two external contacts must be
at the two marked black roots, and the actual connecting white path
would extend a doubled terminal of a large finite B diagram impossibly.

The [tree carrier proof](fresh_affine_tree_carrier_exclusion_20261017.md)
closes the remaining obligation for affine B/C/D subdiagrams inside
H-U and H-V. Grounded resistance forces the contact distance below8;
908 exact role cases reduce to11 surviving rows. The actual connecting
white path excludes two nonprivate contact patterns. Its only remaining
terminal connector is unique and forces the whole carrier to be private,
where the preceding theorem applies. This is an all-rank reduction.

It matters that the proof uses the actual connecting path. The
[exact counterexamples](fresh_large_affine_white_contact_counterexamples_20261017.md)
disprove the stronger claim that two spherical marked deletions alone
force every exterior contact onto black roots. A
[bounded continuation scout](fresh_large_affine_white_continuation_scout_20261017.md)
helped identify the missing path condition; the final all-rank argument
is independent of that scout's sampled orders and labels.

Both complements are now affine-free, which supplies the preliminary
hypothesis of the audited finite-quadrilateral proof. Its later carriers
have at most17 nodes and are excluded by the cusp-rank hypothesis.
Triangles receive at least2 and quadrilaterals at least1, with finite
vertex cost at most3n and ideal-vertex cost at most3sV. The shape
inequality's strictly positive5(n-1) term contradicts those budgets
for every n>=30, including the equality case of a local coefficient.

## Construction control and verification

The exact construction control remains63 actual facets,199 ambient
finite vertices and10 target-face open ends, with certified infinite
volume. The [changed-hinge diagnostic](fresh_nonarithmetic_release30_66_hinge_20261017.md)
checks the exact deletion costs and a bounded5,886-branch search.
Its546 numerically separated candidates have at least13 ends; the
remaining56 recipes have literal small affine circuits giving permanent
cusp obstructions. Numerical boundary warnings remain explicit. No
changed chamber is promoted by that diagnostic.

The [alternative end-cap search](fresh_nonarithmetic_alternative66_endcaps_20261017.md)
tests five changed versions of wall66 and158 coherent low-label caps;
every single addition retains10 open ends. The
[wall59 release](fresh_nonarithmetic_alternative66_release59_20261017.md)
admits25 right-side branches. In the frozen baseline66 control, deleting
wall59 loses seven exact ambient finite vertices and several compact
faces. Among766 opposite-cap pairs,756 are dotted and10 fail numerical
Coxeter-contact gates. The other66 alternatives have numerical boundary
checks only. The
[changed-chain diagnostic](fresh_nonarithmetic_alternative66_chain_patch_20261017.md)
replaces a four-ray path and tests43 endpoint caps plus432 pairs,
without reducing the10-end count. Its high-label checks remain numerical
and conditional. These are bounded construction diagnostics with exact
deletion costs, not new finite-volume controls or global impossibility
theorems.

Every included computational package has successful generation and a
fresh saved replay. Written internal reviews state their scope and do
not claim external refereeing or proof-assistant verification. The
[archive checker](../scripts/audit_forced_quadrilaterals_and_affine_carriers_20261017.py)
and [saved audit](forced_quadrilaterals_and_affine_carriers_20261017_audit.json)
rehash the frozen [previous checkpoint](spectator_field_gates_and_boundary_replacements_20261017.md)
and its cumulative sources, verify new dependency edges, and independently
compose the two rational incidence repairs. Previously verified dependency
graphs are inherited across that frozen hash boundary; the checker does
not repeat every historical mathematical calculation.
