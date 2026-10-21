# Exact remaining edge boundary of the60-wall chamber

The [complete cone](fresh_nonarithmetic_pi8_quad_cone_20261017.md) gives422 compact geodesic edges and848 open geodesic edges incident to its188 finite vertices. This counts every edge incident to a finite vertex, but does not claim to count ambient edges with two open ends.

The [ledger generator](../scripts/fresh_nonarithmetic_pi8_quad_edge_ledger_20261017.py) and [exact output](fresh_nonarithmetic_pi8_quad_edge_ledger_20261017.json) also classify every projective edge of the target three-face S6={0,1,2,3,5,6}:

| Edge type | Count |
|---|---:|
| Two finite endpoints |21|
| One finite endpoint and one open end |12|
| Two open ends, meeting hyperbolic space |3|
| Entirely outside hyperbolic space |3|

Thus the whole S6 face has39 projective edges and15 geodesic edges requiring further closure. The twelve one-open edges have side pairs

    [7,33], [7,49], [8,44], [8,52], [30,33], [30,44],
    [30,47], [30,59], [34,47], [34,52], [49,58], [58,59].

The three two-open edges have pairs[30,49],[30,52],[33,49]. Pairs[44,52],[47,52],[49,59] lie entirely outside hyperbolic space. In particular a projective edge involving newly added walls need not be a hyperbolic edge. The next coordinated construction targets the remaining end of[58,59], whose finite Coxeter component is I2(9), together with neighboring open edges[49,58] and[30,59].

Every finite ray and every S6 ray has exactly nine active normals. Because the inherited cone certificate proves these are extreme rays of a full10-dimensional pointed cone, those nine normals are independent. Deleting one gives a rank-eight edge support. The complete ray list then has exactly two endpoints with that support, as checked literally. Other ambient exterior rays can be nonsimple, so the script makes no simplicity assumption about them.

For two spacelike endpoint vectors a,b, write A=a^2,B=b^2,C=a.b. Their positive span contains a timelike vector precisely when C<0 and AB-C^2<0. The explicit positive combination

    v=(B-C)a+(A-C)b

has square(AB-C^2)(A+B-2C). This identity gives strict rational negative bounds even when direct coordinate evaluation suffers interval cancellation. Its future orientation is checked separately. If C>0 or AB-C^2>0, every positive combination is spacelike. No zero-sign case remains in the ledger.

The finite incidence identity is2*422+848=9*188. The S6 identity is2*39=3*26. Generation and fresh saved replay pass. This is a new exact incidence extraction from the primary complete cone, not an independent reconstruction of that cone or a finite-volume result.

```sh
python3 scripts/fresh_nonarithmetic_pi8_quad_edge_ledger_20261017.py --verify
```
