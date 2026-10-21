# Joint quadrilateral patterns at three small D-component forests

This complete finite-pattern census treats initial rank-13 forests D5+2D4, 2D5+A3, and 3D4+A1. Every other corner may be any of the [215 high-baseline ADE types](fresh_sqrt5_rank13_high_ade_endpoints_20261017.md). Opposite wall products remain free parameters.

| Initial forest | Tested attachment/bond choices | Literal finite-corner patterns |
|---|---:|---:|
| D5+2D4 | 83,266 | 41,678 |
| 2D5+A3 | 125,990 | 56,017 |
| 3D4+A1 | 172,098 | 81,669 |

The 179,364 literal patterns give 11,331 colored records. For each pair of released roots, every pair of complete finite endpoint attachments and both finite new–new products, zero and minus one, is tested. The fourth corner must also belong to the selected set. Each literal map records the two endpoint attachments, all three destination types, the full colored isomorphism and whether it swaps the free opposite products.

Every record contains the exact rational four-root Schur matrix, determinant and factorization. The fixed 3-edge graph forces both free products into the ground field in 6,246 records, exactly one in 3,995, and neither individually in 1,090. Their squares always remain in the field; a cycle may also force the mutual product there.

There are 9,601 coupled factors, 5,051 rational linear edge factors and 1,777 nonlinear univariate factors, counted across records. The rational-edge exclusion alone leaves all 78 initial root pairs for each forest. The nonlinear univariate factors are retained in this census and require the later finite-anchor and positive-null identity gates; no numerical parameter cutoff is used.

The [source](../scripts/fresh_sqrt5_rank13_small_d_joint_quads_20261017.py) and [saved output](fresh_sqrt5_rank13_small_d_joint_quads_20261017.json) passed generation and a fresh complete reconstruction replay:

```sh
python scripts/fresh_sqrt5_rank13_small_d_joint_quads_20261017.py
python scripts/fresh_sqrt5_rank13_small_d_joint_quads_20261017.py --verify
```

The implementation caches white inverses by their full matrix entries and symbolic geometry by all 16 ordered Schur entries. Ground-field flags are recomputed for each row. Deterministic color refinement and the multiset of corner types only select isomorphism buckets: every nonliteral merge still supplies a fully checked bijection. A separate narrow source read confirmed these cache and isomorphism properties. It was not an independent replay of the census or of the inherited rational-edge theorem.

The [exact arithmetic parameter calculation](fresh_sqrt5_rank13_small_d_quad_parameters_20261017.md) now applies the complete finite-anchor domain to these records. This census alone neither realizes the retained patterns in a compact ambient polytope nor proves a dimension bound.
