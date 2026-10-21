# Continuous families, facet incidence and weight barriers

This checkpoint preserves the [preceding compact-face checkpoint](compact_faces_and_complete_weight_obstructions_2026_10_16.md)
and adds closed construction diagnostics and exact obstructions to restricted
weight models. **No new finite-volume example, construction dimension or
dimension bound is established.** The published noncompact arithmetic bound
21 remains distinct from the internally checked compact bound-21 candidate.

The exact H20 chamber still has158 walls, at least184 finite vertices and
five proper cusps; its strict escape proves infinite volume. The exact H9
chamber still has55 walls and a strict escape. Numerical extensions below
have not been promoted to new exact chamber stages.

## H20: continuous families and geometric incidence

The [continuous-family package](fresh_h20_two_anchor_families_20261016.md)
retains2,730 admissible conic families and41,984 samples with two prescribed
old-side contacts. Whole parameter intervals are isolated numerically, without
an arbitrary boost cutoff. Mutual-angle equations produce compatible pairs;
a subsequent all-acute relaxation and simultaneous angle recovery produce
16 compatible numerical clusters absent from the sampled bank. None closes
the target H3. The acute relaxation stops at its explicit100-iteration guard;
it is not an exhaustive exclusion.

The [intrinsic-facet and open-edge report](fresh_h20_intrinsic_facet_incidence_20261016.md)
explains the common obstruction. In a compact target face, every new polygon
facet needs at least three finite neighbors. Combining all the above roots
with the earlier three-contact bank gives50,314 rows; degree removal leaves
1,227. The complete finite-graph cover search on61 saved boundary points is
infeasible in52 states, and an independent310,129-row binary program agrees.
All56 previously tested clusters contain a facet with insufficient degree.
These are bounded numerical exclusions, not continuous-domain theorems.

The seven actual old facets have13 finite edges, of which five need a second
finite endpoint. This provides a concrete frontier for a continuous search
that adds walls at unfinished edges, including walls with one or zero finite
old-side contacts. The more general two-parameter and edge-driven searches
are outside this checkpoint.

## H9: an exact cusp obstruction and a bounded combined-bank failure

The [two-contact package](fresh_nonarithmetic_h3_two_contact_20261016.md)
computes1,144 initial conic branches and removes25 null-line branches by an
exact, all-parameter obstruction:23 affine-triangle recipes remain orthogonal
to a retained angle pi/5, forcing a cusp incompatible with finite volume.
The [affine gate](fresh_nonarithmetic_h3_two_contact_20261016_affine_gate.md)
binds the old exact contacts and positive null kernels.

The remaining1,119 conics are explored through mutual-angle equations and
whole-face intersections. A combined4,542-candidate repair bank is
numerically infeasible on24 accumulated boundary samples after all10,312,611
pair tests. This does not exclude the continuous families. A new search
will couple actual edge-endpoint obligations instead of selecting isolated
caps solely by point coverage.

## Compact arithmetic dimension21: model obstructions, not dimension bounds

The [A4 actual-count model](fresh_sqrt5_h3_a4_actual_count_split_20261016.md)
and [A3 multiplicity barrier](fresh_sqrt5_h3_a3_multiplicity_barrier_20261016.md)
show why changing one H3 compensation orbit is insufficient. The
[all-six-orbit count split](fresh_sqrt5_h3_all_actual_count_split_20261016.md)
still has an exact restricted-model lower bound above37.98, exceeding the
marked-facet target35.

A short [H4 triangle-compensation barrier](fresh_sqrt5_h4_triangle_compensation_barrier_20261016.md)
gives the exact lower bound285/8=35.625. Separating compensation by
[E8 presence](fresh_sqrt5_h4_e8_context_split_20261016.md) breaks that particular
coupling but still yields an exact restricted-model lower bound above37.78.
These are failures of the specified weight families. They do not exclude
dimension21 or invalidate the earlier compact bound-21 candidate, which
uses a different, undiscounted rule in dimension22.

The new global target40 model is an active successor outside this checkpoint.
Any viable discounted rule requires complete coverage of newly affected
faces as well as exact vertex inequalities before a dimension claim.

## Verification and preservation

All included producers and saved-certificate replays have passed. The
individual reports distinguish exact certificates, source reviews, numerical
intervals, finite-bank exhaustion and local optimizer failures. The
[archive checker](../scripts/audit_continuous_families_facet_incidence_checkpoint_20261016.py)
and its [result](continuous_families_facet_incidence_and_weight_barriers_2026_10_16_audit.json)
bind recursive input hashes, reports and local links while preserving the
preceding archive's frozen sections. Archive integrity is not mathematical
validation of a dimension theorem.

```sh
python3 scripts/audit_continuous_families_facet_incidence_checkpoint_20261016.py --verify
```
