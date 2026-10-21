# Correction to the facet-A1 diagnostic prose

The dominant triangle in the [frozen facet-A1 report](fresh_sqrt5_marked_a1_facet_mixed_cross_20261016.md) is core51_attach_D5_root2. Its attached root is the **long D5 endpoint**, and its corners are H4+D5, A1+D8, and **H3+D6**. Root 3 is a short leaf and instead gives H4+D5, A1+E8, and H3+E6. The dominant marked D5 vertex also uses its long endpoint.

The earlier report misidentified these endpoint and third-corner types in prose. Its frozen source and exact certificate already have the correct literal graphs and labels. The exact lower bound 519578441788333/14786000133074 and all scope restrictions are unchanged. Both D6 and E6 missing source-root orbits are included in the subsequent H3-center model.

The standard-library [source](../scripts/fresh_sqrt5_marked_a1_facet_diagnostic_correction_20261016.py) reconstructs connected components and identifies D/E diagrams by branch-arm lengths, then checks the dominant dual terms. The [certificate](fresh_sqrt5_marked_a1_facet_diagnostic_correction_20261016.json) passed generation and a fresh saved replay. This is a narrow diagnostic correction, not a new dimension or coverage claim. All earlier frozen files remain unchanged.
