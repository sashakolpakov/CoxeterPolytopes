# Noncompact section exclusion and changed-chamber constructions

No new compact or finite-volume Coxeter polytope has been found. The completed exact result is that **an H20 Coxeter hyperplane section meeting the strict interior of the literal Borcherds H21 chamber cannot be noncompact**. Any remaining interior Coxeter section of that chamber must be compact and avoid every source cusp. This is a source-specific theorem, not a general dimension-20 exclusion. The current goal remains open.

Three different compact-H9 incidence constructions and a changed-chamber H20 mirror construction were also attempted. They produced no metric candidate ready for finite-volume certification. Their scopes and evidence differ from the complete section theorem and are kept separate below.

A further exact result excludes one full quartic cusp family under arbitrary projected rational ambient caps, for every positive real shift. A separate compact-section audit preserves the open cases and adds 162 exact failed constructions; it does not establish a complete compact-section exclusion.

## The complete cusp-section obstruction

The [coverage proof](fresh_h20_cusp_coverage_20261014.md) closes the remaining logical gap between the archived local cases and an arbitrary section through a source cusp. It selects actual irredundant facets and uses bounded fibers to prove that an active rank-n affine simplex contributes either n or n+1 facets for n>=2. Coincident restrictions and inequalities touching only lower-dimensional faces are included. The single rank-one exception is covered by the earlier absorption proof. The earlier mixed-factor theorem already uses this same full scope.

Thus the eighteen complete cusp orbits, with 53 affine factors, are exhausted by the existing mixed-factor, rank-one, rank-two, all-retained, and one-omitted cases. The last missing cases were the 91 rank-three direction records. Their [new global certificate](fresh_h20_rank3_global_sections_20261014.md) now excludes every full shift interval: 75 continuous actual obtuse-ridge proofs and sixteen continuous forbidden-angle proofs. At the lower endpoint, 64 use the same witness, 22 need a separate ridge, and five are old-facet supporting hyperplanes outside the strict-interior scope.

The [independent verifier](fresh_h20_rank3_global_sections_20261014_review.md) uses a different algebraic-field implementation and checks all 226 face-point occurrences, all 305 parameter-polynomial positivity certificates, and every endpoint case. It imports no producer arithmetic. The local direction census and earlier exclusions remain separately audited dependencies. Both exact saved-result replays pass.

In the Klein model, a finite-volume polytope's only boundary-sphere points are its finitely many ideal vertices. A noncompact section must therefore contain an old cusp. The cusp classification rules this out, yielding the stated compactness conclusion. It does not cover source-facet boundary sections, subdivisions, changed ambient chambers, or arbitrary H20 polytopes. In particular, the source's arithmeticity alone would not justify excluding compact sections.

The [fresh known-example calibration](fresh_cusp_coverage_calibration_20261014.json) passes again: Bugaenko's compact H8 chamber has eleven facets, 41 finite vertices, and no cusps, with exact Vinberg completion verified.

## A changed set of mirrors for H20

The [quartic mirror-orbit construction](fresh_h20_quartic_orbit_completion_20261014.md) keeps the valid cusp link of each quartic direction but releases the other original chamber inequalities. It draws candidate caps from 29,043 integral roots obtained from the original wall set and one ambient reflection. This changes the construction and falls outside the source-section exclusion.

At five shift values for each of 44 quartic records, the numerical scout finds 8,136 cap proposals compatible with the fixed cusp walls, plus 96 near-angle remainders. Compatibility alone is insufficient: a subsequent necessary endpoint test finds no cap closing any of the 42,420 tested cusp edges, even when the near-angle proposals are included. This is a finite numerical mirror-bank result, not an exact exclusion of all mirrors or all parameters.

The [independent bank and scope review](fresh_h20_quartic_orbit_completion_20261014_review.md) passes. It checks every integer reflection recipe and norm, the quartic degree restriction, and all archived omission tuples; it does not certify the floating-point endpoint signs.

The [exact integral cap gate](fresh_h20_quartic_orbit_cap_gate_20261014.md), with its [independent arithmetic review](fresh_h20_quartic_orbit_cap_gate_20261014_review.md), first excludes every integral ambient norm-two/norm-four completion of local record0 at `u=1/2`. Its complete three-wall endpoint census leaves three necessary cap families, but none can end two of the three triangle-edge types. Thus 168 of the 252 cusp edges remain unclosable. This covers the whole specified root class, not merely the numerical bank.

The subsequent [rational-cap theorem](fresh_h20_quartic_rational_cap_gate_20261014.md) removes the ambient norm bound and extends to **every positive real shift of that one spatial direction**. Exact field arithmetic and rational polynomial gcds enumerate the necessary endpoint patterns for arbitrary rational ambient cap vectors. Any endpoint in the first triangle-edge class forces the shift to be rational; irrational shifts therefore leave its 84 edges unclosable. For rational shifts, the other two classes remain blocked, leaving 168 unclosable edges. No finite-volume completion by projected rational ambient caps can retain this fixed cusp family. Nonrational ambient caps, changed cusp walls, and other spatial directions remain outside this result.

The [independent rational-cap review](fresh_h20_quartic_rational_cap_gate_20261014_review.md) reconstructs all 21 rationality gcds, twelve necessary records, and three rational-shift identities using different field arithmetic. Its field-degree test for the absence of sqrt5 is also independent. Generation and saved replay pass.

## The remaining compact source sections

The [compact coverage audit and construction report](fresh_h20_compact_section_coverage_20261014.md) confirms that at most three source-wall contacts change at a finite section vertex. That is a local reduction, not an exclusion of every compact section. The audit also identifies the need to handle section vertices coinciding with old finite vertices separately from transverse open-edge parameterizations.

A fresh construction joins two isolated short A1 factors into A2 on nine selected edge supports, with three contact ratios, three positions, and both residual signs. All 162 recipes fail at exact actual finite ridges: 142 obtuse angles, nineteen squared-cosine gaps between0 and1/4, and one gap between1/4 and1/2. The unsampled continuous families remain open.

The report also proves a compact-only obstruction for any actual G2 ridge formed from two old orthogonal short roots, for all contact ratios and section positions. Its boundary analysis forces a compact face with an invariant rational integral reflection lattice, contradicting rational isotropy in the relevant dimension. This covers that specific local G2 pattern; it is not a general compact-section theorem.

The [separate compact-section review](fresh_h20_compact_section_coverage_20261014_review.md) passes the geometric argument and independently checks ten construction witnesses covering all nine chosen source edges. Its separate arithmetic verifies 2,100 source-wall signs. The full 162-case exact replay remains the primary implementation; the independent arithmetic review is explicitly a sample.

## Compact H9 from nonpolygon bases

The [nonpolygon construction report](fresh_h9_nonpolygon_construction_20261014.md) tries three actual convex base constructions, all different from the earlier product and polygon multiwedge families:

| Construction | H9 facets | Tested scope |
|---|---:|---|
| Six wedges over a dodecahedron | 18 | All 82 admissible symmetry types, from all 12,376 positive compositions |
| Six edge truncations of a cube followed by six wedges | 18 | Thirty explicitly ordered truncation recipes |
| Vertex-truncated dodecahedral prism followed by five wedges | 20 | Three specified symmetry classes |

All 115 label joins are empty with no state truncation. Their exact exclusions use complete mandatory Lannér placements and the actual vertex restrictions; arbitrary unassigned labels and dotted distances are relaxed. A separate incidence checker examines all 5,949,320 nine-facet subsets and verifies the saved vertex lists, minimal nonfaces, ridge degrees, recovered base dimensions, and nonproduct status. It is not a second implementation of the label joins. All primary and incidence saved-result replays pass.

These attempts yield no full Gram matrix to test for signature or compactness. Other bases, truncation orders, and weight choices remain outside the stated scopes.

## Search consequences

Do not continue searching the original Borcherds chamber's cusp-containing interior sections. The remaining original-chamber route is compact and cusp-avoiding; a proposed exclusion there would require additional arguments. For the changed-mirror route, do not enlarge rational ambient cap banks on local record0's fixed cusp pencil: the all-real theorem covers every rational coefficient and norm. Other directions and nonrational ambient caps require actual cusp-edge closures. Compact-H9 work must change the stated failed incidence families again.

The [previous checkpoint](finite_angle_gate_and_algebraic_cusp_directions_2026_10_14.md) remains the source for the arbitrary-normal finite-angle restriction and the finite B2 closure calculation. The new cusp theorem does not turn its broad numerical edge sweeps into exact exclusions.

The [checkpoint provenance audit](noncompact_section_exclusion_and_changed_chamber_search_2026_10_14_audit.json) checks the current source manifests and every local link in these reports. Its [replay script](../scripts/audit_noncompact_section_checkpoint_20261014.py) checks provenance, not mathematical correctness; the exact replays and independent reviews described above supply the mathematical evidence.
