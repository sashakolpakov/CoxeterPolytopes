# Bounded rank9 triangle-extension pilot

The first bounded connected-white search produced343 exact colored graph
classes from all76 Q(sqrt5) minimal Lanner triangle markings. There are207
expanded classes,83 hereditary-safe stops and53 open rank9 boundary
classes. This is bounded evidence only; the boundary was not treated as
closure.

The [producer](../scripts/fresh_sqrt5_triangle_extension_frontier_20261016.py),
[literal export](fresh_sqrt5_triangle_extension_frontier_20261016_rank9.json),
and [fresh temporary-database replay](fresh_sqrt5_triangle_extension_frontier_20261016_replay.json)
are preserved. The [replay source](../scripts/fresh_sqrt5_triangle_extension_frontier_20261016_replay.py)
rebuilds the entire bounded run in a new database and compares its full
export exactly. Generation and a fresh saved replay pass. This is a
reproducibility check using the same producer, not independent arithmetic.

The pilot was extended in the separate
[full dimension21 connected-domain package](fresh_sqrt5_triangle_extension_complete_20261016.md).
Every earlier file is unchanged. The pilot itself does not include common
spectator/facet-mark payment coverage or imply a new dimension bound.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_sqrt5_triangle_extension_frontier_20261016_replay.py --verify
```
