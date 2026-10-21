# Prism-base capacities and an exact boundary replacement

This checkpoint continues the
[previous frozen archive](changed_walls_and_shared_face_incidence_2026_10_17.md).
No new finite-volume example or lower dimension ceiling is established.
The unrestricted finite-volume bound163 remains a **complete, internally
checked computer-assisted proof**, with no known missing mathematical step;
external refereeing and formalization remain outstanding. Its
[status note](dimension163_proof_status_20261017.md) distinguishes those
validation questions from mathematical completeness. The compact arithmetic
bound21 is unchanged; the published noncompact rational case is separate.

## New global incidence theorems

The [disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md)
strengthens the preceding simplex packing rule. At any simple vertex,
distinct maximal simplex faces cannot share even one edge direction.
Two adjacent triangles force a tetrahedral three-face by three-connectivity;
triangular adjacency of directions is therefore transitive. Each resulting
clique closes to a whole simplex by vertex degrees and connectivity.

Prisms also force [simplex-base incidences](fresh_sqrt5_prism_simplex_base_incidence_20261017.md).
A maximal d-simplex has `binomial(d+1,k)(n−d)` possible incidences with
k-prisms along simplex bases. Squares use all four edges as bases.
Exact reconstruction of every81,183 flexible dotted-prism provenance gives
133,346 base-star types. Of these,2,981 are represented in the old simplex
model;99,681 missing types are maximal edges and30,684 are higher-dimensional
simplexes. Their unknown masses are never treated as zero.

The stronger [shared-capacity theorem](fresh_prism_simplex_shared_capacity_20261017.md)
couples all prism dimensions at once:

\[
 \sum_{k,Q} k\,m(Q,S,k)f_Q\le(d+1)(n-d)f_S.
\]

Inside each immediate coface of S, bases of distinct maximal prisms have
disjoint vertex sets. If two bases meet, their common external direction,
the triangular base faces and the quadrilateral side faces force their
union to be one prism. The proof reconstructs its facets explicitly.
This applies to actual maximal prisms and the current unique dotted-star
aggregates, not automatically to overlapping all-finite product aggregates.
The square convention and maximality step received internal reviews.

Both the compact-simplex port rule and the
[compact-prism shared capacity](fresh_finite_volume_prism_capacity_20261017.md)
extend to finite-volume Coxeter polytopes when all counted faces are compact.
Their finite vertices are simple; nonsimple ideal vertices elsewhere do not
invalidate the connectivity or coface arguments.

## Exact limits of the enlarged arithmetic model

The [combined simplex/quad model](fresh_sqrt5_rank13_simplex_quad_combined_20261017.md)
preserves every old simplex constraint and includes68,150 quad rows
representing81,183 full-star provenances. It has6,558 vertex forests and
78,335 face rows. Its first bounded solve timed out without a certificate.

Adding6,528 disjoint-port charges gives a
[solved model with exact primal and dual certificates](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.md),
placing its optimum near12.78544, above the dimension13 target12. This is
a barrier to that weighting family on the saved domain, not a dimension
theorem. All source provenances remain available after equivalent old
constraints are merged.

Further exact dual certificates show that the obstruction persists under
each of these enlargements:

| Added charge family | Certified lower bound |
| --- | ---: |
| All shared three-port indicators | 12.7853100516 |
| Represented-base capacities, separately for each prism dimension | 12.5623648223 |
| Represented-base shared capacities across all prism dimensions | 12.5350583889 |

The [three-port checker](fresh_mixed_three_port_dual_20261017.md) verifies
47,725 orbit columns and every positive merged provenance. Only three
small vertex-mass corrections are necessary. The
[separate-base transfer export](fresh_sqrt5_prism_base_transfers_20261017.md)
has3,840 signed columns with zero vertex cost. The
[shared-base export](fresh_sqrt5_prism_shared_transfers_20261017.md)
has2,845 such columns and checks both capacity families. Each constructs
an exact repaired dual without another optimization run.

The last two certificates include the original pair, clique and disjoint-port
charges; they do not include the three-port charges, whose loads change when
new simplex masses become positive. These are separate enlargements, not a
certificate for their union. Missing base types and full face coverage also
remain outside these claims. The natural next extension is to add the forced
missing base stars and validate all desired charge families together.

## Exact progress in dimension nine

The [64-wall replacement](fresh_nonarithmetic_release54_replace12_20261017.md)
deletes old54 and replaces old12 in the preceding65-wall chamber. It has
197 finite vertices, all64 actual facets and a new whole compact heptagon.
No previously compact face of dimension at least two is lost. Only one old
finite vertex is replaced, and its one compact incident edge gets a new
endpoint. The exact [complete intrinsic boundary census](fresh_nonarithmetic_release54_replace12_boundary_20261017.md)
reduces open ends17 to16 and increases compact edges28 to30.

A strict ambient escape survives and proves infinite volume. The complete
ambient finite census is separate from the intrinsic exterior rays; a full
ambient exterior cone has not been enumerated for this alternative.

The [continuous cavity obstruction](fresh_nonarithmetic_release34_cavity_gate_20261017.md)
also excludes a whole box of two-cap contact parameters, including specified
large-label ranges and their parallel limits. The
[bounded scouts](fresh_nonarithmetic_release34_cavity_20261017.md)
are preserved separately from the exact construction and obstruction.

In dimension twenty, [replacement150 searches](fresh_h20_replacement150_triangles_20261017.md)
and [continuations](fresh_h20_replacement150_continuation_20261017.md)
reduce an exterior-radius diagnostic while leaving openings. A
[parabolic completion diagnostic](fresh_h20_replacement150_parabolic_cusp_20261017.md)
adds a proper cusp conditionally on its numerical base, but increases the
complete target boundary's open ends11 to12. These are numerical diagnostics,
not promoted exact examples.

## Finite-volume transfer and completed validation

The [finite-volume analysis](fresh_finite_volume_simplex_truncation_20261017.md)
gives exact cusp-truncation face counts and a Coxeter square pyramid showing
why the compact triangle-port theorem does not extend to arbitrary triangles
with ideal vertices. Artificial cusp faces exceed the usual uniform
simple-polytope payment budget at every product-of-simplexes cusp type.
Khovanskii's original inequality instead gives a direct original-face budget
with explicit cusp contributions. The remaining inequality needed to bound
those contributions is stated; no improved unconditional ceiling follows yet.

Both previously pending full attachment replays have finished successfully:
the [finite closure](fresh_sqrt5_global_b1_finite_closure_20261016.md) checks
all36,879 states, and the [dotted closure](fresh_sqrt5_global_b1_dotted_replay_status_20261017.md)
checks all25,535 states. Their pruning remains conditional on the specified
discount and bonus assumptions. Full replay does not imply unrestricted
coverage of every possible face or weighting rule.

The [checkpoint checker](../scripts/audit_prism_base_capacities_checkpoint_20261017.py)
and [saved audit](prism_base_capacities_and_boundary_replacement_2026_10_17_audit.json)
bind the completed packages and preserve the frozen predecessor. Exact
certificate replays, written geometric proofs, internal reviews and numerical
diagnostics retain their separate scopes. Archive integrity alone proves no
mathematical theorem.
