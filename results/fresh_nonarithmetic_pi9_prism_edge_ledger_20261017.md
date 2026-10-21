# Complete finite-edge and S6 boundary inventory for62 walls

The [complete cone](fresh_nonarithmetic_pi9_prism_cone_20261017.md) gives431
compact ambient geodesic edges and875 geodesic edges with exactly one finite
endpoint. Their finite endpoint incidences total2*431+875=9*193. This inventory
contains every edge incident to a finite vertex; it does not count all ambient
edges whose two ends are open.

For S6=[0,1,2,3,5,6], all42 projective edges are classified:

| Kind | Count |
|---|---:|
| Compact geodesic edge |22|
| Geodesic edge with one open end |13|
| Geodesic edge with two open ends |3|
| Edge wholly outside hyperbolic space |4|

Thus the actual target three-face still has16 open geodesic edges. Its28
projective vertices comprise19 finite and nine exterior vertices. The exterior
adjacency has two components; the certificate retains the distinction between
geodesic edges and edges wholly outside hyperbolic space.

The [checker](../scripts/fresh_nonarithmetic_pi9_prism_edge_ledger_20261017.py)
and [certificate](fresh_nonarithmetic_pi9_prism_edge_ledger_20261017.json) index
every eight-node support incident to a finite ray or an S6 ray. Each has exactly
two endpoints in the complete cone. No simplicity assumption is made for
unrelated exterior ambient vertices.

For two exterior rays a,b, their positive cone meets hyperbolic space exactly
when their pairing is negative and their two-by-two Gram determinant is
negative. The checker supplies the positive combination
(b²-a.b)a+(a²-a.b)b and verifies its negative norm and future orientation.
The saved determinant identity is
(a²*b²-(a.b)²)*(a²+b²-2a.b). Other projective edges are explicitly classified
as lying outside hyperbolic space.

Generation and a fresh saved replay passed. This is an exact extraction from
the primary complete cone, not an independent reconstruction of that cone.
It establishes no new finite-volume example.

```sh
python3 scripts/fresh_nonarithmetic_pi9_prism_edge_ledger_20261017.py --verify
```
