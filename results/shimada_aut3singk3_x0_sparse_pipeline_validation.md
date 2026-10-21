# Shimada `X_0`: sparse fixed-section pipeline validation

## Status and terminal `X_0` decision

The exact `X_0` restriction has now been run from the pinned complete source.
The unique source-group rank-12 involution gives 1,734 wall orbits and 314
independently strict facets in `H11`, but it is **not Coxeter**.  The first
forbidden strict pair is `(5,37)`; its primitive normals have squares
`(-4,-6)` and inner product `4`, hence normalized squared ratio

\[
\frac{4^2}{4\cdot6}=\frac23,
\]

which is not a crystallographic Coxeter value.  The two exact projection
witnesses separately prove that both inequalities are genuine facets.  One
forbidden pair is terminal, so no Vinberg volume inference is made and no
full pair scan is needed.

The guarded run completed in 2.33 seconds with maximum process RSS 45.16 MiB.
Its exact output is
[`shimada_aut3singk3_x0_symplectic_h11_sparse_screen.json`](../data/enumerations/shimada_aut3singk3_x0_symplectic_h11_sparse_screen.json).
This closes the complete genuine-`Aut(X_0,a_0)` high-subgroup denominator;
the possibly larger abstract oriented weighted-wall group remains outside
this claim.

## Sparse exact algorithm

`scripts/certify_shimada_aut3singk3_sparse_fixed_section.py` performs the
following operations using exact Python integers:

1. Pin the published `GramS`, ample, wall, and `AutXa` source hashes.
2. Parse the wall stock twice, clearing one global denominator, while loading
   only the selected 20-by-20 involution matrix from the 1,440-element list.
3. Verify `A^2=1`, `A G A^T=G` in Shimada's row-vector convention, `aA=a`,
   and exact permutation of every source wall.
4. Form all size-one or size-two wall orbits, take their oriented sums, divide
   by positive content, and deduplicate the resulting oriented primitive
   rays.
5. For every negative ray `t`, certify strict irredundancy with the integral
   positive multiple

   ```text
   (-t^2) a + (a,t) t
   ```

   of its orthogonal-projection witness.  Timelikeness, zero target slack,
   and strict positive slack on every other restricted inequality are all
   tested exactly.
6. Stream the complete facet-pair census.  Each rational squared ratio is
   classified as `m=2,3,4,6`, parallel, ultraparallel, or forbidden.  The
   positive Gram SHA-256 is generated in row order without retaining the
   matrix.  If a forbidden pair occurs, its two primitive roots and two exact
   strict witnesses are embedded in the result.

Thus no `3278 x 3278` Gram or quadratic pairing cache is ever materialized.
Vinberg's endpoint and parabolic-neighbourhood tests remain a separate gate,
authorized only if every strict facet pair is Coxeter.

## Frozen replay validation

Fresh processes independently reconstructed all three completed folds from
the pinned primary supplement.  They agreed with the frozen certificates on
the orbit count, oriented-ray count and signs, every projection witness, the
complete streamed positive-Gram digest, and every relation count.

| replay | wall orbits | unique rays | strict facets | pair result | positive-Gram SHA-256 | elapsed | max RSS |
|---|---:|---:|---:|---|---|---:|---:|
| `X_2` rational | 242 | 207 | 37 | 0 forbidden; Coxeter pair survivor | `48a86bd712574d5fdd057f09d5f56a5c8ec41e8128f8ad82040575064202f303` | 0.808 s | 40.781 MiB |
| `X_2` symplectic | 266 | 266 | 90 | 352 forbidden | `023f4e46a6c582afdac7590e2db2eb7b537784e1d007107dd99317ba0d7a9bb7` | 0.807 s | 40.734 MiB |
| `X_1` symplectic | 594 | 580 | 162 | 915 forbidden | `cbd8b22f90a2dfdfa82cc5968a3b7e1dccbae681dc0e47b9b370bdc1b1eb2fbe` | 0.634 s | 41.188 MiB |

The reported RSS is the maximum resident set of the entire fresh Python
process through certificate construction, not an incremental estimate.  The
eight-test regression suite also contains parser, orbit/deduplication,
strict-witness, forbidden-angle, source-isometry, ample-fixing, and resource
fixtures.

## Exact `X_0` resource model and guard

Before knowing the actual involution orbit count, the worst-case model is:

- 3,278 walls and at most 3,278 unique restricted rays;
- 5,371,003 unordered facet pairs;
- 10,745,284 dense source-Gram entries avoided;
- 196,680 retained coordinate slots across walls, rays, and Gram covectors;
- at most 214,905,680 integer multiply-adds for all strict-witness slacks and
  the same bound for the row-streamed positive Gram;
- 21.05 MiB conservative incremental allocation model, with a 64 MiB hard
  operational budget.

The allocation figure is a transparent model, not measured X0 RSS.  Because
the machine was critically compressed during preparation, the CLI refuses an
X0 run unless both conditions hold:

- the operator explicitly supplies `--authorize-x0-full`; and
- immediately free memory is at least 512 MiB by default.

The completed `X_0` peak confirms that the sparse resource model was in the
right regime.  The terminal run used `--stop-on-first-forbidden`: one
forbidden pair between separately strict facets is already a complete
non-Coxeter obstruction.  Consequently an ungapped full-Gram pass and a
Vinberg volume job are neither required nor logically authorized.

## Reproduction

Memory-light estimate and regression suite:

```bash
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_sparse_fixed_section.py --case x0 --estimate-only
PYTHONPATH=.:src pytest -q tests/test_shimada_aut3singk3_sparse_fixed_section.py
```

Measured completed-stock replays:

```bash
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_sparse_fixed_section.py --case x2-rational --measure-resources --output /private/tmp/x2-rational-sparse.json
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_sparse_fixed_section.py --case x2-symplectic --measure-resources --output /private/tmp/x2-symplectic-sparse.json
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_sparse_fixed_section.py --case x1-symplectic --measure-resources --output /private/tmp/x1-symplectic-sparse.json
```

Terminal X0 first-obstruction replay, when the memory guard is clear:

```bash
PYTHONPATH=.:src python3 scripts/certify_shimada_aut3singk3_sparse_fixed_section.py --case x0 --authorize-x0-full --stop-on-first-forbidden --measure-resources
```

Artifacts:

- `scripts/certify_shimada_aut3singk3_sparse_fixed_section.py`
- `tests/test_shimada_aut3singk3_sparse_fixed_section.py`
- `results/shimada_aut3singk3_x0_high_source_group_checkpoint.md`
- `data/enumerations/shimada_aut3singk3_x0_symplectic_h11_sparse_screen.json`
