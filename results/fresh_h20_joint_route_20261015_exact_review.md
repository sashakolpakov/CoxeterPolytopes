# Independent review of the43-wall H20 cap and finite edge

**PASS.** The exact new cap, all43 actual facets, all903 wall pairs, the complete new-cap finite-vertex census, the displayed finite edge and the renewed strict escape agree with a separate arithmetic implementation. No mathematical or scope issue was found.

The [reviewer](../scripts/fresh_h20_joint_route_20261015_exact_review.py) and [certificate](fresh_h20_joint_route_20261015_exact_review.json) bind the frozen [producer](../scripts/fresh_h20_joint_route_20261015_exact.py), its [JSON](fresh_h20_joint_route_20261015_exact.json) and the literal source geometry through recursively checked manifests. The reviewer imports no producer, helper, project arithmetic, symbolic package or floating-point routine. It implements the basis1,sqrt2,sqrt3,sqrt6 using four rational coefficients, Galois-norm inverses and rational square-root isolating intervals. All658 distinct nonzero signs required only64-bit interval endpoints; every conclusion is strict or an exact coefficient equality.

The physical form has signature(20,1). The first40 normals are independently matched to the source embedding and the next two to the positive coordinate directions. The cap's19 source equations are independent and uniquely determine its saved projection Y. They give Y²=-3/2 and future orientation. The transverse coordinates are exactly(-sqrt3/2,-1), giving cap norm2. The new pair counts are18right,3pi/3,1pi/6,1parallel and19dotted; the complete43-wall counts are649right,30pi/3,6pi/4,1pi/6,32parallel and185dotted.

Every saved facet point is independently checked to be future timelike, to lie on its own wall, and to satisfy every other inequality strictly. This proves actual facets rather than merely an acute Gram matrix. The common center is also timelike and strictly inside all43 inequalities.

For the vertex census, the reviewer reconstructs the cap's22 finite neighbors and tests all1,540 choices of19 other walls. It applies full exact matrix LDL to each20-normal Gram, independently of the producer's forest-Schur algorithm. Exactly one support is positive definite. A finite vertex of a hyperbolic Coxeter chamber has a spherical local chamber and is simple; hence every finite vertex on the new facet is included in this test. The unique positive support is explicitly realized by a timelike point with all nonincident inequalities strict. This is a **complete census only of finite vertices on the new cap**, not of every vertex in the43-wall chamber.

Both displayed endpoints have positive-definite20-wall supports and exact future normalization C·v=-1. Their common19-wall support is positive definite, while every other inequality is strict along the interior of their geodesic segment. The old endpoint is strictly below the new cap and the new endpoint strictly below the released old wall41. These checks establish an actual finite edge.

The source projection Y is strictly future timelike, so its pairing with every future timelike or nonzero future null source ray is negative. Thus the new cap does not cut the retained full H18 face; checking a single interior point alone would not establish this stronger claim. The reviewer checks the timelike norm and common orientation needed for the causal-cone argument.

The old escape has positive pairing with the new cap and is genuinely cut. The replacement escape retains its first19 coordinates and uses transverse coordinates(-3/4,-1/1000). Its norm is strictly positive, all43 wall pairings are strictly negative, and its center pairing is negative. Together with the strict timelike center this supplies an open ideal-boundary patch, so the new chamber still has infinite volume.

The arithmetic values for the conditional non-quasi-arithmetic argument are independently reproduced: Cartan two-cycle12+8sqrt2 and retained rational two-wall determinant-32. The external Vinberg criterion is an inherited theorem dependency; this review does not establish that criterion or any host classification. In particular it does not call the current infinite-volume chamber a lattice or a finite-volume example.

Generation and final saved replay both **PASS**. Syntax, source bindings and report links pass. The review source, certificate and note are frozen:

```sh
python scripts/fresh_h20_joint_route_20261015_exact_review.py
python scripts/fresh_h20_joint_route_20261015_exact_review.py --verify
```
