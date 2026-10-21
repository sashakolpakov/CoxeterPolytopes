# Exact endpoint walking from the changed H8-based H9 seed

Four coherent new caps close **six previously open edge germs**. The resulting chamber has **19 actual facets, exactly 53 finite vertices and 175 finite edges**, compared with 44 vertices and 167 finite edges in the [frozen seed](fresh_h9_h8_disjoint_cap_20261014.md). Three further new vertices involve two of the caps. The literal seed escape survives all four additions, so the chamber still has infinite volume.

## Complete starting frontier

The [exact frontier certificate](fresh_h9_h8_disjoint_walk_20261014.json) reconstructs all 62 open germs incident to the seed's finite vertices: 38 leave the bottom compact-H8 facet, and 24 leave the three top vertices. It inherits the seed's complete finite-vertex and finite-edge census, then checks every remaining release directly.

For a finite vertex `P` with incident Gram `G_I`, releasing wall `j` gives

`D=−Σ_i (G_I⁻¹)_{ij}N_i`.

Thus `P·D=0`, `D²>0`, `D·N_j=−1`, and the other eight support pairings vanish. The forward line is `P+tD`, with null crossing `t²=−P²/D²`. Every other old-wall intersection lies strictly beyond that crossing. The certificate stores an exact spacelike continuation on each edge with precisely the eight support walls active. The opposite null limit is future and has those same eight positive-definite incident normals, so its link is unbounded. This is a complete frontier at the 44 finite vertices, not a census of edges having no finite endpoint.

The normalized support diagrams give **1,076 finite spherical extension profiles and 64 symbolic isolated `I2(m)` tails**. Finite profiles use contacts 3, 4 and 5; the tails retain every integer `m≥6`. Completeness uses the finite Coxeter diagram classification: the new node has at most three neighbors in distinct old components; degree three permits only label 3; a label above 5 can occur only in an isolated rank-two component. The exact Schur residual is checked in every accepted profile.

A separate [local affine certificate](fresh_h9_h8_disjoint_walk_20261014_affine.json) checks 6,056 contact assignments on the 22 supports with no old label 5. It finds **14 profiles** that complete the eight-wall support to an irreducible affine rank-eight diagram. Each has Schur residual zero and a strictly positive kernel. The other 40 supports cannot give a bounded Euclidean link with one new wall because they already contain a label-5 edge. Disconnected extensions leaving a spherical component are not bounded cusp links.

## Bounded metric scout

For a finite profile, let `p` be its projection onto the support span and let `ρ=2−p²>0`. In a future orthonormal frame `U,V` for the complementary Lorentz plane, the outgoing family is

`r(s)=p+√ρ(sinh(s)U+cosh(s)V)`.

The scout samples every numerically detected all-dotted parameter interval and the stated finite-anchor events. It checks that the candidate endpoint

`P−(r·P)(r−p)/ρ`

is future timelike, lies forward on the designated edge, and satisfies every retained old-wall inequality strictly off that support. Wall-pair compatibility alone is not treated as edge closure.

The [numerical manifest](fresh_h9_h8_disjoint_walk_20261014_scout.json) records 1,652 sampled profiles, 365,680 real anchor roots and 23,926 deduplicated accepted endpoint proposals. None cuts the literal seed escape. Tail samples are `6,7,8,10,12,16,24,32,64`; anchor labels are `2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64`. Remaining nonsupport pairs are tested numerically against labels 2 through 256 or a strict dotted interval. These bounds and tolerances do not prove an all-parameter or all-label exclusion. The full proposal bank is stored as [deterministic gzip](../data/enumerations/fresh_h9_h8_disjoint_walk_20261014_scout_bank.json.gz), with its hash in the manifest.

The [affine-event scout](fresh_h9_h8_disjoint_walk_20261014_affine_scout.json) tests 1,666 finite/parallel anchor events for the 14 exact local affine profiles. It accepts 218 numerical events, none cutting the escape. Continuous tails and unsampled anchor labels remain outside that bounded test. No affine proposal is promoted to an exact added wall here.

## The exact coherent quartet

The [exact producer](../scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py) reconstructs its four caps from literal contact recipes, without importing the numerical scout or bank. All four are

`r_k=(p_k,−√K)`, where `K=6+16√5/5`.

Their spatial parts lie in `Q(√5)` and satisfy `p_k²=−2−3√5/2`, so every cap has ambient norm two. The extension arithmetic is exact in `Q(√5,√K)`.

| New wall | Starting germ | Nonright support contacts | Additional anchor |
| --- | --- | --- | --- |
| 15 | 13 | wall 9 at `π/3` | wall 2 at `π/5` |
| 16 | 27 | walls 0 and 10 at `π/3` | wall 1 right |
| 17 | 9 | wall 10 at `π/5` | wall 3 at `π/5` |
| 18 | 19 | wall 9 at `π/5` | wall 2 at `π/3` |

Each row supplies nine independent equations for `p_k`. The remaining old-wall products are then calculated, not prescribed. Among the four caps, `(15,16)` is dotted with product `−2φ`; `(15,17)` and `(16,18)` have product `−φ` and angle `π/5`; the other three pairs are right. All 66 new pairs comprise 31 right, six `π/3`, eight `π/5`, and 21 dotted pairs.

The original strict center remains strict. Explicit orthogonal projections of it give actual facet points for all nineteen walls. The complete finite-neighbor graph has 799 nine-cliques; exact positivity tests reject 746 and retain exactly 53. Each retained support is independently realized by projecting the strict center, with direct future, norm, zero-set and remaining-wall sign checks. This proves the complete finite-vertex census, including preservation of all 44 old vertices. The complete finite-edge census follows by common eight-wall supports.

The six closed old germs are `9,13,15,19,27,34`. Every saved endpoint has positive exact forward parameter and is distinct from the old endpoint. Two further finite edges connect new vertices `18–21` and `37–49` in the certificate's numbering. The three other new vertices, numbered 15, 25 and 26, involve two new caps and have no finite edge to another finite vertex in this chamber. The current remaining finite-corner count is `53·9−2·175=127`; no complete classification of those other ends is asserted.

The [exact certificate](fresh_h9_h8_disjoint_walk_20261014_exact.json) stores all normals, the full Gram matrix, the 66 new pairs, all nineteen actual facet points, all finite vertices and edges, and the six exact forward parameters. In this file `[[a,b],[c,d]]` means `a+b√5+(c+d√5)√K`; its radical differs from the frozen seed's legacy encoding, whose extension coefficients are zero.

## Global limit and verification scope

The old strict spacelike escape remains strict at every new wall. Its unchanged square is `−6837/256+6543√5/512>0`, so the 19-wall chamber has infinite volume. This construction changes actual edge geometry, but does not establish a finite-volume completion. The numerical absence of an escape cutter is not a universal endpoint obstruction. Nonright contact profiles, untested anchor labels and other residual directions remain separate possibilities.

The exact frontier, local affine profiles and quartet each have generation and saved `--verify` checks. The numerical scout replays its deterministic compressed bank; the affine scout replays its bounded events. Source hashes bind the frozen seed and arithmetic helpers. No frozen predecessor was edited.

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014.py --verify
python scripts/fresh_h9_h8_disjoint_walk_20261014_affine.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_h8_disjoint_walk_20261014_scout.py --verify
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h9_h8_disjoint_walk_20261014_affine_scout.py --verify
python scripts/fresh_h9_h8_disjoint_walk_20261014_exact.py --verify
```
