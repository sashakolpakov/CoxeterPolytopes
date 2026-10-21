# Right-angled internal proof archive

The published dimension ceilings remain4 for compact,6 for entirely ideal,
and12 for mixed finite/ideal polytopes. The following new internal results
have separate mathematical, script and formalization scopes.

* [Second distinct ideal RA4 volume](../results/ideal_ra4_second_volume_20261018.md):
  volume8π²/3, attained by the facet double of the24-cell.
* [Equality structure](../results/ideal_ra4_second_volume_equality_20261018.md):
  every second-volume polytope is that reflected double.
* [Octahedral decomposition](../results/ideal_ra4_octahedral_block_decomposition_research_20261018.md):
  if every facet is isometric to an octahedron or its triangular-face
  double, the polytope is a chain of24-cells with opposite successive cuts.
  The direct proof has separate decomposition, opposition and TeX reviews;
  its [nine-page illustrated manuscript/PDF](ideal_ra4_octahedral_decomposition_reproduce.md)
  has five figures and a fresh exact D4 model check. The facet hypothesis
  is essential; no universal24-cell decomposition is claimed.
* [Compact facet capacity](../results/compact_ra4_perpendicular_facet_capacity_research_20261018.md):
  f>=2V+5 for each facet with V vertices, giving V<=48 strictly below the
  regular120-cell volume. The [seven-page manuscript/PDF](compact_ra4_facet_capacity_reproduce.md)
  has the complete geometric and ear arguments, two figures and a root
  mathematical/source review. Neither compact volume conjecture is proved.
* [Five-face carriers](../results/right_angled_five_face_carrier_20261018.md):
  the universal ideal5 sumκ48 and mixed incidence constraints.
* [Positiveκ construction](../results/ideal_ra4_positive_kappa_doubling_20261018.md):
  an exact16-cell ideal4 polytope withκ6, and unboundedκ under doubling.
* [One-layer attachment](../results/ideal_ra5_positive_carrier_one_layer_20261018.md):
  a universal ideal5 attachment with closed original cusps and an open
  ideal escape; its volume is infinite.
* [Mixed five-face density](../results/right_angled_mixed_five_face_density_20261018.md):
  ideal5 needs56 facets; ideal five-face fractions are less than1/10 in
  dimensions11/12, with explicit finite/ideal vertex-ratio constraints.
* [Ideal6 propagation](../results/right_angled_ideal6_carrier_propagation_20261018.md):
  at least227 facets and3859 vertices, and ideal six-face fraction below3/40
  in dimensions11/12.
* [Sharp cusp transport](../results/ideal_ra6_sharp_reciprocal_20261018.md):
  the universal ideal6 vertex minimum is strengthened to4462, with
  V>=ceil(2673F/136); independent proof reviews and the serial three-checker
  TeX package pass. The former3859 bound is a preserved predecessor.

* [Coupled cusp compatibility](../results/ideal_ra5_known_carrier_cusp_compatibility_20261018.md):
  within the stated24-cell/double/16-cell facet spectrum, at least30
  sixteen-cell facets and66 doubles are necessary, with a forced eight-carrier
  cusp. This is a restricted-spectrum theorem.
* [Octahedral blocks](../results/ideal_ra5_octahedral_blocks_20261018.md):
  the24-vertex ridge adjacencies partition those sixteen-cell facets into
  blocks of six. A Lorentz and coloring obstruction strengthens the necessary
  counts to36 sixteen-cell facets and84 doubles. The same explicit
  three-type restriction applies; no dimension ceiling follows.
* [Mixed weighted-face theorem](../results/right_angled_mixed_weighted_triangles_research_20261018.md):
  exact vertex-ratio bands and a forced canonical ten-facet four-face in
  dimensions11/12. Corrected exact replay and independent root geometric
  review pass. This has no restricted facet-spectrum hypothesis and does
  not yet exclude either dimension.
* [Low-facet rigidity](../results/right_angled_mixed_low_facet_rigidity_research_20261018.md):
  equality at facet-plus-cusp sum15 in dimension4 uniquely gives P4;
  a16-facet right-angled five-polytope uniquely gives P5. Exact model
  certificates and an independent geometric/source review pass.
* [Actual move support](../results/right_angled_mixed_actual_move_support_research_20261018.md):
  every nontrivial P3-opposition component in dimensions9--12 reaches a
  K9-free corner. The global ideal-edge, three-face sum<=11 class has
  dimension<=8. Independent written root review passes; no new ceiling.
* [Symmetric quad charge](../results/right_angled_mixed_symmetric_quad_charge_research_20261018.md):
  unrestricted actual quad moves, invariant common-neighbor weights and
  a bound on unbacked induced paths in terms of large-three-face excess
  and finite-edge incidences. Independent written root review passes;
  full-field causality and dimension12 exclusion remain open.

The mixed transport results now have a
[complete illustrated TeX exposition](../manuscripts/right_angled_mixed_opposition_transport.tex)
and [proof/reproduction map](right_angled_mixed_opposition_transport_proof_map.md),
with five figures and a
[full root mathematical translation review](../results/review_right_angled_mixed_opposition_transport_tex_root_20261018.md).
Its PDF build is queued behind the mandatory available-memory start gate.
The two conditional ceiling-eight statements retain their global hypotheses;
the unconditional mixed ceiling remains12.

The second-value and equality results are complete internal volume/structure
proofs. The octahedral-chain theorem is a restricted structural classification,
and compact capacity is a universal necessary inequality. The other
results are constructions or necessary conditions with their stated scope;
they do not supply new dimension ceilings. Their standalone manuscript
sources and reproduction maps are indexed in [the manuscript directory](../manuscripts/README.md).

The archive exporter recursively verifies every bound source hash and copies
the complete reports into TeX. It preserves historical status language.
It does not convert that prose automatically into a Lean proof. The
[formal coverage map](../formal/README.md) records the exact proved modules
and remaining internal geometric obligations. Published inputs are listed
separately from the internal arguments in each result-specific map.
