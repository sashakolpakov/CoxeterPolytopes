# Complete subgroup closure for Roulleau's 17-wall `U+E8+A1^3` stock

Roulleau's complete 17-curve nef chamber on
`U+E8(-1)+A1(-1)^3` has exact weighted-diagram group `S3`.  Its subgroup
conjugacy classes have orders `1,2,3,6`; the three proper classes produce
exactly two distinct fixed cones.  Both are finite-volume noncompact
hyperbolic Coxeter polytopes, and both full-reflection indices are exactly one.

| subgroup classes | chamber | saturated fixed divisor lattice | Vinberg volume audit | full reflection lineage |
|---|---|---|---|---|
| `C2` | `H11`, 15 facets | even rank 12, determinant `-8`, Smith `1^10 2 4` | 54 finite vertices, 3 cusps, 341 elliptic edges, no bad edge; all 7 affine-neighbourhood audits pass | the unique proper integral overclosure has index 2, is odd of determinant `-2`, and is explicitly isometric to `<2>+<-1>^11`; intrinsic all-root Vinberg stops at height 36 with exactly these 15 mirrors, hence index 1 |
| `C3` and `S3` | the same `H10`, 13 facets | exactly `U+E8(-1)+<-6>`, even determinant 6, Smith `1^10 6`; no proper integral overclosure | 31 finite vertices, 2 cusps, 176 elliptic edges, no bad edge; all 4 affine-neighbourhood audits pass | intrinsic all-root Vinberg stops at height 8 with exactly these 13 mirrors, hence index 1 |

The source transcription is independently pinned by more than the wall graph:
`A1,...,A13` reconstruct a determinant-eight lattice of signature `(1,12)`,
all four dependent curve coordinates are integral, the published ample class
has square 294 and the stated 17 intersections, and the published nef class
has square 2 and exactly the stated zero set.

The lineage answer is negative for abstract novelty.  The `H11/15` chamber is
Scharlau--Kirschmer List11 entry 2 and the classical determinant-two `qf2`
CoxIter chamber.  The `H10/13` chamber is List10 entry 9; it is absent from
the two pinned CoxIter lanes, but the explicit lattice isometry and
wall-by-wall intrinsic Vinberg conjugacy prove that it is nevertheless the
known full reflection chamber.  Thus this stock contributes two positive
divisor-lattice folds and **zero new abstract hyperbolic polytopes**.  Whether
these exact divisor-lattice extractions were previously recorded remains a
bibliographic question.

The surface-action assertion is separate.  Roulleau records `Aut(X)=C2` and
tag `t` for the general K3, so its nontrivial automorphism acts trivially on
`NS(X)`.  Every nonidentity diagram action here acts nontrivially on `NS`.
Consequently none of the folds is induced by an automorphism of the general
surface; no claim about special-period lifts is made or needed for the
abstract chamber proof.

Replayable files:

- [`roulleau_ue8a1x3_all_subgroup_folds.json`](../data/enumerations/roulleau_ue8a1x3_all_subgroup_folds.json)
- [`certify_roulleau_ue8a1x3_all_subgroup_folds.py`](../scripts/certify_roulleau_ue8a1x3_all_subgroup_folds.py)
- [`qf2_h11_vinberg.py`](../src/coxeter_polytopes/search/qf2_h11_vinberg.py)
- [`u_e8_m6_vinberg.py`](../src/coxeter_polytopes/search/u_e8_m6_vinberg.py)
- [`test_roulleau_ue8a1x3_all_subgroup_folds.py`](../tests/test_roulleau_ue8a1x3_all_subgroup_folds.py)

