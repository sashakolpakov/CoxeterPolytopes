# Roulleau rank-5--8 complete-nef target-header stock compiler

## Scope and source

This audit compiles the 17 entries in Xavier Roulleau's *An atlas of K3
surfaces with finite automorphism group* ([arXiv:2003.08985](https://arxiv.org/abs/2003.08985))
whose `(Picard rank, complete irreducible (-2)-curve count)` header belongs to

```text
(5,5), (6,6), (7,7), (8,8), (8,9), (8,10).
```

Those are precisely the lower-rank headers that were empty in the 30-package
ADEG promotion ledger.  This is a header-targeted denominator, not a claim to
compile every rank-5--8 K3 lattice in the atlas.

The primary arXiv source snapshot is pinned by SHA-256:

- source tar: `3e445ed01446b47ded6c007b395308de34e566b736f0ba982470bc497adeb1e9`;
- main TeX: `e6eeaaa4a3f375f209762718e8ef9a361b0a95356cc8be2e9f99a6713abf5df3`.

Roulleau's final table explicitly says that its numeric column is the number
of `(-2)`-curves.  Each selected subsection supplies the graph and an ample
class, and identifies a wall subset generating the Neron--Severi lattice.

## Exact stock ledger

| stock | lattice | rank/walls | determinant | Smith diagonal | cusps |
|---|---|---:|---:|---|---:|
| `R5-UA3-NEF-H4-5` | `U+A3(-1)` | `5/5` | 4 | `1,1,1,1,4` | 1 |
| `R5-P4D4-NEF-H4-5` | `[4]+D4(-1)` | `5/5` | 16 | `1,1,2,2,4` | 3 |
| `R6-UA4-NEF-H5-6` | `U+A4(-1)` | `6/6` | 5 | `1,1,1,1,1,5` | 1 |
| `R6-U2D4-NEF-H5-6` | `U(2)+D4(-1)` | `6/6` | 16 | `1,1,2,2,2,2` | 5 |
| `R6-UD4-NEF-H5-6` | `U+D4(-1)` | `6/6` | 4 | `1,1,1,1,2,2` | 1 |
| `R7-UA5-NEF-H6-7` | `U+A5(-1)` | `7/7` | 6 | `1,1,1,1,1,1,6` | 1 |
| `R7-UD5-NEF-H6-7` | `U+D5(-1)` | `7/7` | 4 | `1,1,1,1,1,1,4` | 1 |
| `R8-UD6-NEF-H7-8` | `U+D6(-1)` | `8/8` | 4 | `1,1,1,1,1,1,2,2` | 1 |
| `R8-UA6-NEF-H7-8` | `U+A6(-1)` | `8/8` | 7 | `1,1,1,1,1,1,1,7` | 1 |
| `R8-UE6-NEF-H7-8` | `U+E6(-1)` | `8/8` | 3 | `1,1,1,1,1,1,1,3` | 1 |
| `R8-UA3x2-NEF-H7-9` | `U+A3(-1)^2` | `8/9` | 16 | `1,1,1,1,1,1,4,4` | 2 |
| `R8-UA2A4-NEF-H7-9` | `U+A2(-1)+A4(-1)` | `8/9` | 15 | `1,1,1,1,1,1,1,15` | 1 |
| `R8-UA1A5-NEF-H7-9` | `U+A1(-1)+A5(-1)` | `8/9` | 12 | `1,1,1,1,1,1,2,6` | 2 |
| `R8-UA2D4-NEF-H7-9` | `U+A2(-1)+D4(-1)` | `8/9` | 12 | `1,1,1,1,1,1,2,6` | 1 |
| `R8-UA1D5-NEF-H7-9` | `U+A1(-1)+D5(-1)` | `8/9` | 8 | `1,1,1,1,1,1,2,4` | 1 |
| `R8-UD4A1x2-NEF-H7-10` | `U+D4(-1)+A1(-1)^2` | `8/10` | 16 | `1,1,1,1,2,2,2,2` | 4 |
| `R8-UA2x3-NEF-H7-10` | `U+A2(-1)^3` | `8/10` | 27 | `1,1,1,1,1,3,3,3` | 9 |

## Certificates

For every stock the machine archive verifies:

- a symmetric complete wall intersection Gram with diagonal `-2`;
- exact Lorentzian inertia `(1,r-1,m-r)` and wall rank `r`;
- an integral published wall basis spanning the named Neron--Severi lattice;
- exact determinant and Smith diagonal;
- the published ample square and every wall pairing, all strictly positive;
- Vinberg's finite-volume edge/cusp criterion by exhaustive exact subset
  enumeration.

All 17 are **finite-volume noncompact**, not merely Lorentzian candidates.
The exhaustive Vinberg output has no bad elliptic edge and no bad parabolic
zero-neighborhood in any row.  A separate run of the project's pinned C++
Vinberg verifier reproduces the edge, finite-vertex, cusp, and endpoint
histograms for all 17 inputs.

Completeness is algebro-geometric: the atlas gives the complete irreducible
`(-2)`-curve count.  On a K3 surface those curve inequalities cut out the
full nef chamber inside the positive cone.  This is the no-hidden-wall fact
needed by the downstream polytope-to-surface promotion theorem.

## Reproducibility

- certifier:
  [`certify_roulleau_rank5_8_matching_full_nef_stocks.py`](../scripts/certify_roulleau_rank5_8_matching_full_nef_stocks.py)
- archive:
  [`roulleau_rank5_8_matching_full_nef_stocks.json`](../data/enumerations/roulleau_rank5_8_matching_full_nef_stocks.json)

The frozen canonical payload SHA-256 is
`e05af6a969ad9b27fc3929c8fe499cf4975f1f066409bfbff46df58dc8d8a201`.
