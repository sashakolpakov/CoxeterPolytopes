# Separate arithmetic audit of the 93- and 94-wall H20 stages

Both exact stages pass the Fraction-based checker, including their strict spacelike escapes. Neither stage has finite volume.

| Stage | Added walls | New pair checks | Preserved inherited cusp signs | New finite vertices checked at insertion |
|---|---:|---:|---:|---:|
| 88 to 93 | 5 | 450 | 535 | 3 |
| 93 to 94 | 1 | 93 | 107 | 1 |

The checker reconstructs the new normals with arithmetic separate from the producer, verifies their norms, every new wall pairing, the strict center, and the saved target cuts. Finite vertices are certified by exact positive-definite rank-20 Gram matrices and feasibility against all walls present when that cap is inserted. The table counts vertices at their individual insertion steps; it does not assert that every earlier vertex survives subsequent insertions within or after the five-cap package.

Preserving the 107 inherited cusp links does not prove that newly created cusps have bounded links. The two new parallel pairs require a separate audit and completion. This check makes no claim that all cusps of these stages have been closed.

The subsequent [96-wall completion certificate](../data/enumerations/fresh_h20_edge_family_completion_20261013_cusp_repairs.json) repairs those two links and separately checks that all four vertices counted above survive every final wall. The separate Fraction audit on this page ends at 94 walls; the two-cap completion has its own source review and replay.

The 93-wall escape has square
`193562714077866193/128000000000000000000 > 0`.
The 94-wall escape has square
`210646311442057377/128000000000000000000 > 0`.
Each is future oriented and strictly satisfies every wall inequality in its stage, proving infinite volume. Deliberate changes to a normal, a declared norm and a target-cut value are all rejected.

The inherited 88-wall construction is a frozen dependency with a recorded hash. This audit does not re-enumerate that source or certify the large numerical proposal banks.

The [93-wall certificate](fresh_h20_edge_family_completion_20261013_independent.json) and [94-wall certificate](fresh_h20_edge_family_completion_20261013_four_independent.json) contain the full pair checks. Reproduce with:

```sh
python3 scripts/fresh_h20_edge_family_completion_20261013_independent.py --verify
python3 scripts/fresh_h20_edge_family_completion_20261013_four_independent.py --verify
```

Implementation: [shared separate checker](../scripts/fresh_h20_edge_family_completion_20261013_independent.py), [94-wall entry point](../scripts/fresh_h20_edge_family_completion_20261013_four_independent.py).
