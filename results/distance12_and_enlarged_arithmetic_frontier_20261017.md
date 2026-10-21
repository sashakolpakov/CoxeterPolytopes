# Distance12 and the enlarged arithmetic model

The new unrestricted nonexistence argument gives dimension at most136
for finite-sided, finite-volume real hyperbolic Coxeter polytopes, with
no arithmeticity or cusp restrictions. Its compact-facet consequence
gives ambient dimension at most103. Primary generation and saved replay
pass, as do a separate Fraction implementation and mathematical/source
review. Compact arithmetic remains21, and the minimum-cusp-factor-rank9
class remains29. No new finite-volume example has been established.

## Unrestricted146 to136

The [new proof](fresh_selective_distance12_dimension136_20261017.md)
pays pairs at distance at most12, together with all pairs in connected
components of order at most17. It exhausts the minimal-core configurations
of a quadrilateral that lacks a short marked path.

The complete inverse tests leave18 possible ports and115 disjoint
marked joins. Two underlying join types contain a wholly white affine
block; the other four force two negative directions. The overlap census
has162 consistent identifications:120 fail an inverse test, six contain
a white affine block, and36 give two isolated-core types.

The shared-vertex core was already present in the146 argument, now with
27 marked contexts. Its actual radical relation forces any new supporting
white contact to create a short path. The nine new shared-edge contexts
have a nonsingular18-root Lorentzian Gram. Preferred corners restrict
extra contacts to two end sets. Contacts at both ends give a short path;
at one end, the inverse Gram has positive entries and diagonal at least4,
so a norm-two additional root would have a negative orthogonal residual.
Thus the marked component is isolated and has at most17 nodes.

The resulting costs are12n at a finite vertex and
sum min(12,r_i/2) per cusp-link vertex. Cap24 concentration reduces the
remaining dimensions137–146 to1175 profiles. Every exact curvature
margin is positive; the minimum is172096/1953125. The inherited146
proof excludes the tail. The profile(24,24,24,24,24,15) has negative
margin at ambient136, so these same scalar weights do not give135.

The [independent review](fresh_selective_distance12_dimension136_20261017_review.md)
records the separate implementation and geometric review. These checks
are internal to the agent workflow, not external refereeing or
proof-assistant formalization.

## Arithmetic: reoptimization changes the obstruction

The [enlarged finite model](fresh_actual_corner_enlarged_model_20261017.md)
preserves every earlier exported column and adds1454 cap columns,
52 signed forced-incidence columns,294 actual-corner clique constraints
and98 saved joint cuts. One bounded optimization of430666 variables
terminates successfully. Exact primal and dual certificates give

    12.52181135130475 <= optimum <= 12.521820585116.

The finite system still cannot attain the target12. The earlier8.80486
result repaired one fixed witness; a new optimum uses different masses.
This is not a new arithmetic dimension theorem.

The [support reconstruction](fresh_actual_corner_enlarged_model_20261017_support.md)
finds2638 positive simplex/edge rows and2294 quadrilateral rows. Of these,
397 quadrilateral rows are included in the common explicit-provenance
map used by the new incidence constraints. The other1897 remain outside
that map. Dropping them and recomputing all finite loads gives a
restricted witness of11.903446304081223. No arbitrary provenance is
assigned to a merged row.

At the full new masses, [complete local simplex-packing covers](fresh_actual_corner_enlarged_model_20261017_local_cones.md)
cover all815 relevant forests and5011 block orbits. Their exact repair
lowers the fixed witness to12.440403231430313. A subsequent
[joint-family scan](fresh_actual_corner_enlarged_model_20261017_joint_ports.md)
checks3519737 inequalities using all simplex rows and the397 actual
quadrilateral classes; the1897 unexpanded classes have explicit zero
coefficients in those new constraints. The repair gives12.43082950275823.
This covers the stated families, not all geometry of the merged rows.

The next geometric cohort is small and explicit:36 unselected positive
rows of raw source25/canonical core21 each have exactly one saved
provenance. Together with19 selected positive rows of this source, they
give55 classes and440 finite corners for an actual-corner compatibility
calculation. No row splitting is needed for these36. Their IDs, exact
masses and aggregate objective contributions are recorded in the support
report.

## Cusp incidence: two shortcuts fail on known geometry

The [exact carrier census](fresh_nonincident_carrier_incidence_20261017.md)
reconstructs all3456 two-faces of the known H17 pyramid. It has1048
quadrilaterals unpaid by the fine rule, including987 compact ones. All
are carried by the two small E8 factors of the same nonincident cusp.
Its link has648 edges and81 vertices. Hence an unpaid-face bound of one
per link edge fails, as does six per link vertex per small factor:
even987 exceeds6 times2 times81.

An explicit all-N H3 strip family also has arbitrarily many compact
quadrilateral faces with the same fixed-type nonincident cusp carrier.
Those faces are all paid8. This rules out bounds on all carrier
incidences from the cusp-link combinatorics alone, while preserving
the possibility of a bound that uses actual payment deficits and more
geometry. Triangle fine weights are not claimed by this census.

## Construction: a different four-cap disk

The [new construction scout](fresh_new_cavity_construction_20261017.md)
uses the changed pentagon[7,8,44,63,30], with four boundary arcs of
lengths1,2,1,1 and either diagonal of the interior four-cycle. It avoids
the already excluded three-cap topology. Across94 coherent low-label R
seeds,242 low norm equations and780 continuous high-label samples, no
coherent mandatory neighboring cap A survives.

The [exact successor gate](fresh_new_cavity_construction_gate_20261017.md)
excludes the complete mandatory A–R corner domain:96 low R recipes,
10 symbolic high-label families and285 A–R cases. Exact subdivision
gives1295 certified leaves. Some continuous parameter boxes are
excluded because they lie strictly between consecutive Coxeter angles;
the result covers every finite label, without claiming that every real
contact parameter is impossible. Both central diagonals of this disk
are excluded. Other four-cap configurations remain outside the claim.

The exact construction control remains63 walls,199 finite vertices and
10 target-face open ends, with certified infinite volume. A genuinely
different boundary incidence, retained-wall geometry or deformation is
still needed before a finite-volume example can be promoted.

## Archive and next steps

The [incremental auditor](../scripts/audit_distance12_and_enlarged_arithmetic_frontier_20261017.py)
and [saved audit](distance12_and_enlarged_arithmetic_frontier_20261017_audit.json)
bind these packages and rehash the frozen predecessor files. Small
independent scalar checks recompute all1175 profile margins, the exact
finite-model bracket and the monotone fixed-witness repairs. They also
match all1175 independent-review margins,75 ports and36 corner domains,
and recount the285 construction interval trees. Large
historical censuses and optimization runs are not repeated by the archive.

Next: classify the next long-path exceptions; expand the unique-provenance
raw25 cohort before testing actual-corner compatibility; and replace the
excluded construction boundary incidence. Cusp-defect estimates must
survive the explicit unpaid H17 counts and distinguish them from paid
carrier incidences in the H3 family.
