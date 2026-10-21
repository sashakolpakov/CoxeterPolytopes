# Changed walls and shared face incidence

This checkpoint continues the [preceding archive](split_bridges_simplex_incidence_and_discriminants_2026_10_17.md).
**No new finite-volume example or lower dimension ceiling is established.**
The unrestricted finite-volume bound163 remains a complete internally checked
computer-assisted proof, with external review and formalization outstanding;
see its [status record](dimension163_proof_status_20261017.md). The compact
arithmetic bound21 is unchanged. The research priority remains construction
in a new dimension and upper bounds much closer to existing examples.

## Exact progress on the dimension-nine boundary

The [changed-wall construction](fresh_nonarithmetic_s6_release34_exact_20261017.md)
is a different65-wall chamber from the preceding split-pair branch. Starting
from the63-wall chamber, it deletes old34 and52, retains61 walls and adds a
replacement52 and three new caps. The previous three-cap angle mismatch is
resolved by an algebraic deformation of52. A rational root interval, opposite
endpoint signs and a positive interval derivative specify the exact parameter.
All250 new contacts, a common strict center and all65 facet witnesses pass.

The patch contains two whole compact quadrilaterals and a whole compact
pentagon. Its complete ambient finite-vertex census is196:188 preserved old
vertices and eight new ones. The old strict escape still lies beyond the
ideal boundary inside every wall, proving infinite volume.

The [exact intrinsic boundary census](fresh_nonarithmetic_release34_s6_boundary_20261017.md)
checks the whole H4+A2-supported three-face, including the openings created by
the deletions. Relative to the63-wall parent, its finite vertices increase
from19 to23, compact edges from22 to28, and **open ends decrease from19 to17**.
The one-open/two-open edge counts change from13/3 to13/2. No ideal vertices
or singular triples occur in this census. This is a net boundary improvement,
but the three-face remains noncompact.

The intrinsic calculation uses the walls with positive projected square;
future causal projected walls are redundant on the hyperbolic three-space.
Its nine exterior projective rays in each chamber are not a complete ambient
cone count. In particular, the whole ambient cone of this alternative has not
been enumerated. The [numerical scout](fresh_nonarithmetic_s6_release34_patch_scout_20261017.md)
and [numerical boundary ledger](fresh_nonarithmetic_s6_release34_boundary_20261017.md)
are preserved separately from the exact proofs.

## New incidence constraints for compact arithmetic bounds

The [enlarged flexible-simplex model](fresh_sqrt5_rank13_full_flexible_simplex_lp_20261017.md)
contains10,185 distinct simplex stars and1,870 whole vertex forests. Its old
pair-weight variables have an exact obstruction above12.0238, against the
dimension13 target12. This is a limitation of that partial model.

A [maximal-simplex port-packing theorem](fresh_simplex_port_incidence_20261017.md)
adds genuine incidence information. Two distinct maximal simplices through a
vertex share at most one incident edge direction. Thus suitable collections
of simplex corner patterns cannot all occur there, even when every individual
root-pair capacity permits their fractional masses. Exact tests reject the
old dual distribution at22 vertex types and supply22 valid new charge columns.

Keeping every old row and complete vertex budget, the
[augmented model](fresh_sqrt5_rank13_simplex_packing_repair_20261017.md)
has exact feasible weights with maximum
239780393396/19999999999, approximately11.98901967, **below12**.
This establishes useful additional power of the incidence method; it does
not close the face-coverage requirement.

The [all-saved-simplex separation](fresh_sqrt5_rank13_packing_all_star_separation_20261017.md)
passes all10,185 flexible stars but finds1,417 unpaid stars in128 of158 fixed
discriminant classes. The rule is defined on all6,789 rank13 forests by zero
extension, so absent forests have a valid zero vertex budget; missing face
payments remain genuine failures. A local arithmetic tetrahedron already
receives zero against demand8 with all four vertex types represented.

The [canonical quadrilateral separation](fresh_sqrt5_rank13_packing_quad_separation_20261017.md)
deduplicates1,399 marked dotted cores to1,335 connected cores and119,698
spectator contexts to108,805 full stars. Of81,183 flexible full stars,77,076
are unpaid;4,107 pass. The27,622 fixed-class stars are counted but not evaluated.
There are10,396 failures using only previously represented vertex forests.
One compact local arithmetic quadrilateral with spectator A2+A9 has three
corners of type2A1+A2+A9 and one of type2A2+A9. Its exact payment is zero
against demand1. The core has odd quotient rank and is determinant-flexible.
No missing-forest explanation removes that witness.

A further [mixed face-incidence theorem](fresh_mixed_face_port_incidence_20261017.md)
now supplies a shared simplex/quadrilateral budget. Each actual two-face at a
simple vertex occupies one pair of incident edge directions. Distinct canonical
face aggregates occupy disjoint pair sets. For any three edge directions,
at most one aggregate can occupy two of their three pairs. Averaging this
indicator over vertex-diagram automorphisms gives a valid charge that can pay
both kinds of aggregate. Literal features are checked on all12,420 saved
simplex stars and1,399 marked dotted cores. Individual triangles and individual
quadrilaterals get zero from this charge; the zero-payment quadrilateral above
still requires an actual repair of the weights.

The [dimension14 domain](fresh_sqrt5_rank14_core_spectator_domain_20261017.md)
contains178,723 stars obtained from310 connected cores and every permitted
orthogonal finite spectator. One [bounded solve](fresh_sqrt5_rank14_core_spectator_lp_20261017.md)
timed out without an exact primal or dual certificate; its domain replay
passed. This is unresolved computation. Even a dimension14 exclusion alone
would not establish a universal ceiling13 without excluding larger dimensions
or proving an applicable descent theorem.

## Dimension twenty and validation

The preceding exact171-wall H20 alternative remains infinite-volume. Its
[numerical boundary ledger](fresh_h20_replacement_boundary_20261017.md)
has17 finite vertices and11 open ends. Subsequent bounded contact searches
find compatible pairs and a three-cap roof, but
[complete boundary accounting](fresh_h20_replacement_three_cap_roof_20261017.md)
shows that closing a new compact facet can increase the remaining openings.
Similarly, [deleting old150](fresh_h20_replacement150_boundary_20261017.md)
exposes a compact triangular cavity and loses compactness of another facet.
These numerical diagnostics guide replacement-wall searches; they are not
exact new examples. Their saved generations and replays are preserved.

The finite and dotted attachment queues completed generation. Their full
saved replays remain active at this checkpoint. Their pruning remains
conditional on the stated discount and bonus assumptions; completing a replay
does not imply unrestricted face-catalogue coverage.

The [archive checker](../scripts/audit_changed_walls_shared_incidence_checkpoint_20261017.py)
and [saved audit](changed_walls_and_shared_face_incidence_2026_10_17_audit.json)
bind an explicit list of completed artifacts and preserve the preceding
archive. Exact certificates, bounded numerical diagnostics, written internal
geometric reviews and unresolved computations retain their distinct scopes.
Archive integrity alone proves no mathematical claim. New construction
scouts and the combined simplex/quadrilateral solve remain outside this archive.
