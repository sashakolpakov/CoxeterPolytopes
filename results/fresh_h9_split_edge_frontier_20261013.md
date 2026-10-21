# Existing-wall edge frontier of the split-seed chambers

The [diagnostic](fresh_h9_split_edge_frontier_20261013.json) explicitly enumerates the outgoing edges at the known finite vertices of the saved 72-wall chamber, and then repeats the walk after adding the exact `walk2747` cap. Normals are reconstructed from the exact archives before their coordinates are converted to the saved physical frame.

| Displayed chamber | Finite vertices reached | Directed edge rays | End at an existing finite vertex | Reach the light cone before any existing wall |
|---|---:|---:|---:|---:|
| 72 walls | 4 | 36 | 4 | 32 |
| 73 walls, with `walk2747` | 5 | 45 | 6 | 39 |

The four closed directions in the first row are the two orientations of the U8 and V8 edges. The other 32 directions keep seven seed walls and one endpoint cap. They reach the light cone while every other currently supplied wall remains satisfied. For each, the archive also records the first wall hit farther along the affine ray and a numerical spacelike point before that hit.

With `walk2747`, a breadth-first walk starting only from the original four vertices discovers the cap's new finite vertex. Its new finite edge accounts for both additional closed directions. Eight outgoing directions at the new vertex remain open; one previously open direction has closed, so the open frontier grows from 32 to 39. This increase counts newly exposed work, not a deterioration of the chamber.

For each simple finite vertex the nine independent active equations determine its affine point. Dropping each one in turn gives all nine outgoing edge directions. Along each ray, the code compares the first positive existing-wall intersection with the positive solution of the light-cone equation. Finite hits are added to the queue, and closed edges must be found in both orientations. This exhausts the finite-vertex components reached from the four starting vertices. Other disconnected finite-vertex components are not excluded.

The ray ordering and frontier counts are **numerical diagnostics**. They do not enumerate undiscovered mirrors, certify the nearest arithmetic root, or implement Allcock's complete edgewalking algorithm. No numerically near-ideal ambiguity occurred in this run. Exact replay within the same script separately verifies all four original finite corners in both displayed chambers and the new finite corner against every wall inequality, including negative point norm and positive definite incident Gram. The exact new point matches the vertex discovered by the numerical frontier. These exact checks reuse the construction lane's arithmetic and finite-vertex routine; they are not a second independent implementation.

Source geometry is [the 72-wall scout](../data/enumerations/fresh_h9_split_endpoint_selected_geometry_20261013_scout.json). The extra root and its original exact certificate are in [the no-removal cap archive](../data/enumerations/fresh_h9_split_endpoint_walk_exact_20261013.json). Source hashes, each directed ray, all numerical intersection parameters, and the exact finite-corner certificates are retained in the diagnostic.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_split_edge_frontier_20261013.py --verify
```

The [script](../scripts/fresh_h9_split_edge_frontier_20261013.py) performs no cap search. A constructive continuation should select an open edge and solve for a compatible cap meeting it before the light cone; a complete arithmetic walk would additionally require a residual-coset census proving that the selected mirror is first.
