# Complete cone after the split endpoint pair

The complete 65-wall cone has **196 finite vertices, no ideal extreme rays, and 2,852 exterior rays**. The frozen 63-wall strict escape remains strictly feasible, so this chamber still has infinite volume.

Inputs are the [exact two-cap construction](fresh_nonarithmetic_s6_split_pair_20261017.md) and the [complete 63-wall cone](fresh_nonarithmetic_escape62_20261017_cone.md). The [checker](../scripts/fresh_nonarithmetic_s6_split_pair_20261017_cone.py) saves its crossing certificates and complete final inequalities in [the JSON](fresh_nonarithmetic_s6_split_pair_20261017_cone.json).

| Added wall | Removed exterior rays | New rays | Finite vertices after clipping | Exterior rays after clipping |
|---|---:|---:|---:|---:|
| 63 | 30 | 60 | 194 | 2,824 |
| 64 | 30 | 60 | 196 | 2,852 |

All 193 previous finite vertices persist. The three new supports agree exactly with the separate finite-neighbor census in the cap certificate. The previously compact tetrahedron, quadrilateral, and triangular prism retain their complete four-, four-, and six-vertex cones.

The computation starts directly from the preceding computational coordinates. It checks the basis identifiers and Gram against the cap source before clipping; no physical-coordinate conversion is performed. Each new wall has a certified strict sign on every current ray. For an opposite-sign pair, the full inherited ray list determines its common-zero face. It is an edge exactly when that face contains those two rays. Every accepted crossing also has an explicit rank-eight minor. The positive crossing combination then satisfies nine independent active equations; direct reconstruction with a fixed nonzero coordinate gives the same oriented projective ray with tighter rational bounds.

The checker verifies all 198,120 final wall inequalities and every final norm sign. There are no unresolved zero cut signs or ideal norm cases. The retained strict spacelike escape and its future chart are checked against all 65 walls. The [separate boundary ledger](fresh_nonarithmetic_s6_split_pair_20261017_edge_ledger.md) distinguishes completed edges from replacement open boundaries.

This is a primary clipping certificate using frozen arithmetic and inherited cone completeness. A separate read-only mathematical source audit by the parent agent accepted the basis checks, adjacency criterion, rank certificates, reconstruction, and complete finite-support comparison; it did not independently repeat the arithmetic.

Validation: generation and fresh saved replay both passed. Syntax, direct source manifests, and report links passed. Source, saved data, and this report are final and frozen.
