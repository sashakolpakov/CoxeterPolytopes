# Shimada `X_0`: high source-group checkpoint

## Result

Shimada's theorem-complete `X_0` chamber is now pinned as the next unmined
algebraic wall stock.  It has 3,278 walls on a rank-20 Néron--Severi lattice,
and its genuine ample stabilizer `Aut(X_0,a_0)` has order 1,440.  An exact
first-stage census reduces every proper source-group fixed section above
`H9` to one conjugacy class of involutions, of fixed rank 12 and hence
hyperbolic dimension 11.

This is a queue-reduction theorem, not a new Coxeter polytope and not a
no-go theorem for `X_0`.

## Complete source stock

The primary source is Shimada, [*The automorphism groups of certain singular
K3 surfaces and an Enriques surface*](https://arxiv.org/abs/1412.6904), with
its published computational supplement.  The checkpoint verifies the pinned
hashes and then replays:

- the rank-20 Néron--Severi Gram matrix, of determinant `-36` and signature
  `(1,19)`;
- the ample class `a_0`, of square 20;
- all 3,278 walls, in thirteen complete orbits of sizes
  `60,40,180,10,144,144,240,360,180,240,240,720,720`; and
- all 1,440 integral Gram isometries in `Aut(X_0,a_0)`, including their
  ample-fixing condition.

## Exact high-rank reduction

The fixed-rank histogram of all 1,440 matrices is

| fixed rank | 2 | 4 | 6 | 8 | 10 | 12 | 20 |
|---:|---:|---:|---:|---:|---:|---:|---:|
| elements | 144 | 584 | 450 | 170 | 46 | 45 | 1 |

The 45 rank-12 elements are involutions and form one conjugacy class: the
centralizer of a representative has order 32, so the class has size
`1440/32=45`.  Exhaustive multiplication of high-element pairs produces 30
closed `V4` subgroups, all with common fixed rank 8.  Therefore any proper
subgroup of this genuine source group whose common fixed lattice has rank at
least 11 is conjugate to this unique `C2` class.

## Deliberately open gates

The representative `H11` fixed section has not yet been constructed.  In
particular, the following remain open:

- restriction and deduplication of all 3,278 walls;
- exact facet irredundancy;
- the complete Coxeter pair-angle test;
- Vinberg's finite-volume test, which is authorized only if the angle gate
  first leaves a Coxeter survivor; and
- enlargement from the genuine surface group to the full abstract oriented
  weighted-wall group.

This stop is intentional: materializing the larger restriction or dense
wall Gram while the Record-31 compact search is active would compete for the
same constrained memory.  The next decisive computation is the sparse,
orbit-aware restriction and exact angle gate for the unique `H11/C2` class.

That computation now has an execution-ready integer-streaming driver,
independently replayed against all frozen `X_1` and `X_2` folds.  It has not
been launched on `X_0`: a 512-MiB immediately-free-memory guard preserves the
Record-31 priority.  See the
[sparse pipeline validation](shimada_aut3singk3_x0_sparse_pipeline_validation.md).

## Reproduction

The full source replay is:

```bash
PYTHONPATH=.:src python3 scripts/inventory_shimada_aut3singk3_x0_high_source_group.py --check-archive
```

The regression test is static and memory-light:

```bash
PYTHONPATH=.:src pytest -q tests/test_shimada_aut3singk3_x0_high_source_group_checkpoint.py
```

Frozen artifacts:

- `data/enumerations/shimada_aut3singk3_x0_high_source_group_checkpoint.json`
- `scripts/inventory_shimada_aut3singk3_x0_high_source_group.py`
- `tests/test_shimada_aut3singk3_x0_high_source_group_checkpoint.py`
