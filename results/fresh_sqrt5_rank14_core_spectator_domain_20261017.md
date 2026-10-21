# Direct rank-14 simplex domain from connected Lorentz cores

This exact count pilot constructs 178,723 distinct full simplex stars in dimension 14. It uses all finite orthogonal spectator forests for every connected Lorentz component inherited from the complete saved rank-13 selected-endpoint triangle catalogue. No optimization or dimension theorem is claimed.

The [source](../scripts/fresh_sqrt5_rank14_core_spectator_domain_20261017.py) and [result](fresh_sqrt5_rank14_core_spectator_domain_20261017.json) bind the frozen [full triangle separation](fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017.json). The result includes every core, every old-to-core graph map, every finite spectator context, and literal simplex corner occurrences.

| Quantity | Exact count |
| --- | ---: |
| Inherited full stars | 12,420 |
| Connected Lorentz cores | 310 |
| Distinct full rank-14 stars | 178,723 |
| Triangle / tetrahedron / four-simplex stars | 144,328 / 26,266 / 8,129 |
| Whole rank-14 corner forests reached | 11,747 |
| All finite rank-14 forests in the field catalogue | 11,751 |
| Unordered root pairs at each vertex | 91 |
| Dimension-14 payment threshold | 14 |

The four finite vertex types absent from these stars are `A14`, `B14`, `D14`, and `A1 + B13`. The inherited connected cores have orders 3 through 14. In particular, connected order-15 cores are absent, so this domain does not cover all rank-14 triangles or simplex stars.

## Reconstruction and absence of duplicates

Each inherited Gram star is nondegenerate Lorentzian and has exactly one nonspherical connected component. All its other components are finite. The source recognizes the finite components against the exact weighted finite-tree catalogue and extracts the unique remaining connected component. Color refinement is only a candidate filter; every merger requires a complete weighted graph isomorphism, checked on every unordered pair, including nonedges.

All one-root deletions of each new core are independently recognized against the finite catalogue. Their list agrees with the inherited star's spherical deletions. These deletions are exactly the vertices of the maximal projected simplex. Every corner retains both its omitted vertex root and its pair of released roots. For `m` spherical deletions, the source records exactly `3 * binomial(m, 3)` corner occurrences and demand `2 * binomial(m, 3)`.

For a core of order `q`, all multisets of finite components of total rank `15-q` are enumerated. Orthogonal addition preserves its Lorentz signature, conjugate positivity and compact projected simplex; admissible finite-component gauges preserve arithmetic integrality. A one-root deletion of the full star is spherical exactly when it is one of the core's listed deletions. Deleting a spectator root leaves the Lorentz component intact.

A full-star isomorphism must preserve the unique Lorentz component. The other connected components are identified exactly by their finite types. Thus a core isomorphism class and a finite spectator multiset determine one full-star isomorphism class. The 178,723 contexts contain no duplicate full stars under this equivalence.

## Intended payment test and limits

Any subsequent optimization must rebuild complete 91-pair budgets for the rank-14 vertex forests, including the added spectator roots. It cannot reuse a rank-13 budget or transfer the preceding numerical upper bound directly. This pilot contains no weights and runs no optimizer.

The threshold comes from the same simple-polytope inequality: in dimension 14, `f2/f0 > 7*6/2 = 21`, while each vertex has `binomial(14,2) = 91` two-face corners. Hence `2T+Q >= 5 f2 - 91 f0 > 14 f0`. This explains the target 14 without extrapolating the rank-13 target 12.

There is no discriminant partition. The ambient form has odd rank 15; multiplying the whole form by a totally positive scalar can change its determinant squareclass. The earlier even-rank determinant classification therefore does not transfer unchanged.

The finite spectator enumeration is complete for these inherited cores. Missing connected cores, quadrilateral stars, other polygonal faces, and simultaneous global realization remain outside the claim. Local admissible simplex stars are not compact ambient reflection polytopes.

Generation and a fresh saved reconstruction replay passed.
