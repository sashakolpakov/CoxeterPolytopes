# Dominant residual regions and nonright escape cutters — 2026-10-16

No new finite-volume Coxeter polytope or construction dimension is established.
The compact arithmetic bound21 and unrestricted finite-volume bound163 remain
internally checked proof candidates, pending external mathematical validation.
The published unrestricted compact bound29 is unchanged. Published noncompact
arithmetic bounds21, with20 excluded, are a separate input; the proposed
arithmetic improvement is on the compact side.

This checkpoint preserves the
[preceding archive](completed_cusps_finite_edges_and_weight_barriers_2026_10_16.md).
Its [auditor](../scripts/audit_dominant_regions_nonright_escape_cutters_checkpoint_20261016.py)
and [manifest](dominant_regions_and_nonright_escape_cutters_2026_10_16_audit.json)
check recursive source hashes, predecessor preservation and local links.
Archive integrity is distinct from mathematical verification.

## H20: a compact triangle and coordinated cuts of large exterior regions

The [153-wall construction](fresh_h20_pi8_triangle_face_20261016.md) and
[independent review](fresh_h20_pi8_triangle_face_20261016_review.md) give an
entire compact (2,4,6) triangle face, two new finite vertices and preservation
of all five proper cusps. Exact vertices and convexity prove the whole face
is retained; the ambient chamber still has a strict spacelike escape.

The [fixed-chart numerical comparison](fresh_h20_residual_monotone_scores_20261016.md)
showed that this local success scarcely reduced the largest residual exterior
region. That evidence redirected the search. The
[154th wall](fresh_h20_dominant142_cap_20261016.md), with its
[independent review](fresh_h20_dominant142_cap_20261016_review.md), imposes an
exact z142 cutoff near0.430157 and creates two finite vertices.

The [coordinated two-axis cluster](fresh_h20_dominant_axis_cluster_20261016.md)
and [independent reconstruction](fresh_h20_dominant_axis_cluster_20261016_review.md)
then give **156 actual facets and at least180 finite vertices**. Both new
walls meet the154th wall at pi/3 and are mutually dotted. Their complete new
finite-vertex census is four; the five proper cusps and compact triangle
persist. Independent checks cover all309 new pairings and all156 incidence
signs at each new finite vertex.

In the same four-dimensional residual section, the numerical exterior
Euclidean volume estimate falls from1.289858 to0.0114417, about99.11%, while
the largest squared vertex radius falls from16.3345753 to1.41724998. These
are numerical search scores with no rigorous error bound, not hyperbolic
volumes or a claim about the whole20-dimensional exterior region. The same
strict E6 survives all156 walls, proving infinite volume. Exterior vertex
counts can increase while the actual region shrinks.

The [adjacent H3 diagnostic](fresh_h20_triangle_adjacent_h3_20261016.md) is
numerical only. No additional exact closure follows from it. Further walls
releasing the right contact with144, and coordinated closures of the joint
residual region, remain active work outside this checkpoint.

## H9: an all-real obstruction and a complete new cone

The [right-support projection theorem](fresh_nonarithmetic_h3_escape_projection_20261016.md)
and [independent review](fresh_nonarithmetic_h3_escape_projection_20261016_review.md)
show that no coherent acute wall perpendicular to S6={0,1,2,3,5,6} can cut
the renewed53-wall escape. Its projection is an exact positive combination
of four existing finite vertices. The argument uses any common future
interior, so it does not assume retention of the original chosen center.
It excludes that support restriction, not arbitrary completions.

The [two-wall construction](fresh_nonarithmetic_nonright_s6_20261016.md)
releases that restriction and has
[exact local certificates](fresh_nonarithmetic_nonright_s6_20261016_caps.json)
and an [independent local review](fresh_nonarithmetic_nonright_s6_20261016_caps_review.md).
Both cut the old escape, are mutually dotted and give55 actual facets.
All18 possible new finite supports are excluded, so no finite vertex is
added. Native wall norms are retained; the pi/4 contact to old wall12 has
product-1/2 because that old root has norm1/4.

The [complete cone certificate](fresh_nonarithmetic_nonright_s6_20261016_cone.json)
and [independent cone review](fresh_nonarithmetic_nonright_s6_20261016_cone_review.md)
give **179 finite vertices, zero ideal rays and2,509 exterior rays**. All179
old finite vertices persist. The two walls jointly remove257 original
exterior rays, but exhaustive clipping creates further exterior directions.
The independent computation checks443 rank-eight crossings and147,840
final inequalities. A renewed exact strict escape still proves infinite
volume. A separate coordinated closure of the actual H3 face remains a
construction target; the face is not assumed to have Coxeter dihedral angles.

The bounded scout and pair-selection data are discovery evidence only.
Their exhaustiveness applies only to their explicitly enumerated domain.

## Compact arithmetic21 remains open

The [facet-A1 payment model](fresh_sqrt5_marked_a1_facet_mixed_cross_20261016.md)
adds30 coefficients and the33 previously found nonpath triangle extensions.
Its exact restricted-model lower bound is
519578441788333/14786000133074, approximately35.1398915942, above the required35.
The larger triangle constraint domain means this increase does not measure
the effect of the facet split alone. This is a model obstruction, not a
dimension exclusion or a converged optimum.

Read the [diagnostic correction](fresh_sqrt5_marked_a1_facet_diagnostic_correction_20261016.md)
with that frozen report: the dominant D5 mark is its long endpoint, and
the corresponding third corner is H3+D6. A different short-leaf attachment
has H3+E6. The original source and numerical bounds were already correct;
the correction concerns the report's prose.

A [new restricted arithmetic gate](fresh_sqrt5_long_white_two_angle_quad_gate_20261016.md)
excludes one two-angle quadrilateral pattern with an arbitrarily long white
A-path, while leaving both opposite products free. It prevents that specific
pattern from imposing a proposed discount cap. It supplies neither complete
short-face coverage nor a feasible dimension21 weight rule.

The H3-central-root payment model and its solver improvements remain active
outside this checkpoint. Compact arithmetic21 over Q(sqrt5) remains unresolved.
None of these payment-model bounds is a new polytope or a global nonexistence
theorem.

```sh
python3 scripts/audit_dominant_regions_nonright_escape_cutters_checkpoint_20261016.py --verify
```
