# Exact thirty-one-wall H9 chamber: execution and remaining escape

19 October2026. The thirty-wall certificate now has a separate replay
PASS. The prepared thirty-first Pell mirror then passed generation and
a separately launched replay. These are serial exact arithmetic checks;
the chamber still has infinite volume.

## Source-bound evidence

|Artifact|SHA256|
|---|---|
|`scripts/check_h19_thirty_one_wall_seed_20261019.py`|`e81cd38a68c5262acc6421c16f8d8114e2480f5f453edde503a32de30710c76e`|
|`results/h19_thirty_one_wall_seed_20261019.json`|`5a3b25aa139793f304d535f56d6033ad103d88cffe07c20278bc7b42c2712d57`|
|`results/h19_thirty_one_wall_seed_generation_20261019.json`|`9e7b2716835a7851889c20c2df83f37db13fa178081f7746d2ef957ffd7e8e83`|
|`results/h19_thirty_one_wall_seed_replay_20261019.json`|`e04ae84c292804f392e9fb465e6582d72552330de12da646339ea1a79657cba3`|
|`results/h19_thirty_wall_seed_replay_20261019.json`|`7c4ee1ead10f5b8f225f517c6662d3e7a154dee97cc13df527b5c9d02a7ad5fe`|

The checker verifies465 signed native pairs,31 actual facet projections,
reflection integrality, the common strict center, the actual cut of the
previous escape, and the new E7+2A1 finite endpoint. Seven off-cap finite
vertices are checked; with the two inherited41-vertex compact cap counts,
the lower bound is89. Those cap counts remain inherited inputs rather
than a fresh exhaustive enumeration. The bounded walk itself discovers
28 finite vertices before its next spacelike endpoint, with24 queued.

The new endpoint is reached from vertex3, freeing old3 and first hitting
q. Its square is

    6227/4840 - (691/1210)*sqrt(5) > 0.

Every one of its31 wall products is nonpositive. Together with the
strict timelike center this gives a strict feasible spacelike vector
by a sufficiently small perturbation, and therefore an open ideal
boundary patch. The finite-wall chamber has infinite volume. This is
not inferred merely from the walk's failure to close its queue.

Generation and replay took1.234/1.419 seconds, with process-tree peaks
57,409,536/57,704,448 bytes. Start resources were182.35GiB disk and
7.00/6.62GiB available RAM. The standard shared lock,59-second limit,
1GiB process-tree cap and thread limits were retained. The saved
certificate is53,403 bytes; no temporary build files were produced.

To replay in a fresh bounded batch, use an unused ledger path:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock --ledger results/h19_31wall_recheck.json -- scripts/check_h19_thirty_one_wall_seed_20261019.py --verify
```

This is progress in completing an arithmetic first-type H9 chamber,
not a new finite-volume example, a reflectivity theorem, a complete
root census or a dimension-bound improvement. The next task is to cut
the new actual escape while preserving all current walls.
