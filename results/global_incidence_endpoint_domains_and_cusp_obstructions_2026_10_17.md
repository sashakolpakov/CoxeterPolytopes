# Global vertex incidence, complete endpoint domains and cusp obstructions

This checkpoint preserves the [preceding archive](near_example_targets_and_complete_dotted_domains_2026_10_17.md).
No new finite-volume polytope or stronger dimension bound is established here.
The objective remains a new construction dimension and ceilings near known
examples, approximately23–26 for unrestricted finite volume and10–13 for
compact polytopes.

The unrestricted finite-volume bound163 has a **complete, internally checked
computer-assisted proof**, with no currently identified missing mathematical
step. It has not been externally refereed or formally verified. The
[status clarification](dimension163_proof_status_20261017.md) explains the word
*candidate* in frozen historical reports. The compact arithmetic bound21 and
the separate published noncompact arithmetic bound21, excluding20, retain
their previously stated scopes. Unrestricted compact retains the published29.

## A global incidence constraint, beyond per-vertex weights

The new [D4 argument](fresh_sqrt5_d4_short_face_incidence_20261017.md), with a
[separate mathematical and arithmetic review](fresh_sqrt5_d4_short_face_incidence_review_20261017.md),
applies to compact arithmetic polytopes over Q(sqrt5). In dimension d=4k,
if V_D counts vertices with full diagram kD4 and V_other all other vertices,

    (4k−9) V_D < (12k+2) V_other,  k≥3.

Same-type edge replacements are forced. There are no all-special triangles,
and an all-special quadrilateral can use only center–leaf corners. The
leaf–leaf exclusion uses an integral norm congruence and the inertness of2.
A global count of short faces then gives the frequency inequality. This is
a scoped incidence theorem, not a universal corner discount or dimension bound.

The [mixed D4/D5 census](fresh_sqrt5_mixed_d4_d5_migration_20261017.md) advances
the three actual small-component budget obstructions. For D5+2D4, 2D5+D4 and
3D5+A2 it finds24,30,38 possible same-type endpoint attachments and11,13,16
possible two-special triangle corners per vertex. All attachments and finite
framing labels are covered. Migrating roots can have framing products outside
the ground field even when their squares belong to it; labels4,6,10 are
retained accordingly. Free-parameter quadrilaterals remain unresolved.

The [zero-discount certificate](fresh_sqrt5_global_zero_discount_20261017.md)
also closes a proposed repair route:1,709 exact dual terms force budget at least
43232664297439/1055499993667, greater than40. This excludes the specified
expanded bonus family with both discounts zero, not a geometric dimension.
The [original repair stop](fresh_sqrt5_global_b1_first_repair_stopped_20261017.md)
and [bounded replacement timeout](fresh_sqrt5_global_b1_first_repair_ipm_timeout_20261017.md)
are preserved as unresolved solver outcomes, not infeasibility proofs.

## Complete H9 endpoint families and whole-region obstructions

The55-wall H9 chamber remains infinite-volume. The [E55 gate](fresh_nonarithmetic_a2_escape_partition_20261017.md)
excludes cutting the designated escape by any wall that induces a facet of
the chosen H4+A2 face. It does not exclude ambient cutters missing that face.

The [complete endpoint domains](fresh_nonarithmetic_h3_endpoint_domains_20261017.md)
cover all17 open edge ends of that three-face:296 finite contact profiles and
14 full high-label families, with continuous normal formulas and all55 old
wall constraints. The H4 cusp obstruction makes compactness of this face
necessary in any finite-volume completion.

Intersecting these profiles gives [three all-label single-wall obstructions](fresh_nonarithmetic_h3_region_truncations_20261017.md).
No one wall can close the entire triangular patch[7,30,49] or either
quadrilateral patch[8,30,44,52] and[30,34,47,52]. The proof exhausts finite
profiles and continuous high-label tails with exact outward intervals.
It leaves coordinated multiple-wall completion open. The next decisive
test is feasibility of the unavoidable right-to-H4² endpoint family at[7,30]
against all retained Coxeter constraints, before another finite cap bank.

## Esselmann's dimension20 gap and faster exact arithmetic search

The [D19 countercusp certificate](fresh_esselmann_d19_countercusp_20261017.md)
explains a concrete cusp obstruction behind the exceptional gap. A D5
complement in the Niemeier lattice with root system E6^4 lies in the genus
of D19 but its reflection roots span only rank18. The comparison D4
complement has20 spanning reflection directions. This is an exact example
illustrating the genus obstruction, not a reconstruction of Esselmann's
entire classification and not an exclusion of nonarithmetic H20.

A [conjugate-radical lemma](fresh_conjugate_radical_lemma_20261017.md) removes
a redundant exponential affine-subdiagram scan from the dotted attachment
queue. The [exact migration comparison](fresh_sqrt5_global_b1_dotted_migration_20261017.md)
freshly checks all6,537 processed predecessor cores and14,848 accepted
attachments. The retained set is unchanged. The stopped predecessor is
immutable; the faster successor and the separate finite queue remain active
and incomplete outside the archive. Their weight-pruning assumptions remain
explicit and do not cover arbitrary larger discounts.

H20 E8 cusp-link completion and coordinated boundary searches also continue
outside this checkpoint. No numerical wall recipe is promoted to a complete
finite-volume chamber here.

## Verification and outstanding scope

Included exact producers and saved replays pass; independent reviews are
identified individually above. The [archive checker](../scripts/audit_global_incidence_endpoint_checkpoint_20261017.py)
and [saved audit](global_incidence_endpoint_domains_and_cusp_obstructions_2026_10_17_audit.json)
bind the frozen predecessors, recursive sources and local report links.
Archive integrity does not itself prove the mathematical claims.

```sh
python3 scripts/audit_global_incidence_endpoint_checkpoint_20261017.py --verify
```

Pending work includes mixed-type quadrilateral incidence, completion and
replay of both arithmetic attachment queues, and compatible H9/H20 wall
systems that actually close all relevant boundary. A new construction still
requires full ambient finite-volume verification.
