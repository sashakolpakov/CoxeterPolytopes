# Additional triangle constraints leave an exact gap above 35

The exact restricted-model lower bound is **519578441788333/14786000133074 ≈ 35.1398915942 > 35**. Rational weak duality excludes a budget-35 assignment in this payment family. The last numerical finite-cut lower bound was 35.1406927242 after three solves and 18,204 vertex cuts; a fully converged optimum is not claimed.

This model separates H4 compensation when the facet root is an isolated A1 from compensation when that root belongs to another non-H4 component. The distinction is made separately at each corner, preserving the literal facet root through the whole face. Thirty new coefficients are added, for 31,025 variables in total. All inherited bonus coefficients remain free and nonnegative.

The model also adds the 33 nonpath extensions found in the [preceding bounded scout](fresh_sqrt5_marked_mixed_a4_cross_free_bonus_20261016.md). Each attaches one finite component to node 0 of core 51 by a 3-edge and has all three spherical corners. The other 17 accepted extensions are A-path rows already covered. This is a larger constraint domain, so the increase from the preceding lower bound does not measure the value of the facet split alone.

The complete marked finite-forest dynamic program retains exact counts of A1, I2(5), H3, H4, A4, D4, and D5. The face obligations comprise 196 selected triangle families in 541,151 marked contexts and all 1,979 inherited baseline-zero quadrilateral rows in 227,408 common contexts, giving 53,389 distinct profiles. Before optimization, 87 full forest profiles were compared with direct literal root-pair summation.

The exact dual has 300 vertex terms, 83 triangle terms, and 831 quadrilateral terms. Each dual vertex is reconstructed and recounted by literal pairs in the producer replay. The largest vertex masses are 4H4+D5 marked at a short D5 leaf (about 35.64%) and 4H4+A5 marked at an A5 endpoint (14.30%). The dominant triangle is core51_attach_D5_root2 with three external H4 components and facet root 5. Its third corner contains H3+E6, for which the inherited H3-cross catalogue supplies no E6 source-root payment. Related A-tail rows have H3+A6 corners with the same gap. Targeted H3-central-root cross payments are a concrete next family, not part of this certificate.

The [source](../scripts/fresh_sqrt5_marked_a1_facet_mixed_cross_20261016.py) and [certificate](fresh_sqrt5_marked_a1_facet_mixed_cross_20261016.json) are verified by generation and a fresh saved replay. Same-type and mixed cross-pair normalizations are inherited unchanged and are recounted in the certificate.

**Scope:** a restricted payment-model obstruction. The added triangles cover only a bounded attachment domain; complete coverage of all newly affected triangles and formerly paid quadrilaterals remains open. No dimension bound, nonexistence theorem, or new polytope follows.
