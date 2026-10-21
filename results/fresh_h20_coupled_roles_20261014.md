# A cusp-edge cutter coupled to a compact S75 hexagon

The [exact construction](fresh_h20_coupled_roles_20261014_exact.json) adds seven compatible walls to the frozen [100-wall chamber](fresh_h20_polar_pair_20261014.md). One wall cuts its saved strict escape and closes a selected cusp edge at a finite endpoint. Six other walls close the entire S75 hyperbolic plane as a compact right-angled hexagon. All 107 listed walls are actual facets.

There are exactly seven finite vertices incident to the new walls: the cusp-edge endpoint and the six hexagon vertices. Every old finite vertex survives strictly, giving **at least 141 known finite vertices**. The three previously bounded cusp links also survive unchanged. The other retained S61 face remains unbounded, and a renewed strict rational escape proves that **the 107-wall chamber still has infinite volume**.

## The cap on a different cusp edge

Use germ 0 of the [complete repaired-cusp walk](fresh_h20_e7_repaired_cusp_walk_20261014.json). Its positive-definite nineteen-wall support is

```
S = [2,3,4,5,6,7,13,18,31,65,79,83,84,88,89,90,91,92,94].
```

This is different from the four supports in the earlier escape-dependent endpoint-family obstruction. Let ν be the repaired cusp, and let W be the orthogonal projection of

`H0 = N0 + ··· + N19 + N24`, with `H0² = −2`,

to the support's Lorentzian two-plane. Put `h = W·ν < 0` and

```
Z = −2hW + W²ν,
R = Z − ν,
P = Z + ν.
```

The checker proves `Z² = ν² = 0`, `Z·ν = −2h²`, `R² = 4h² > 0`, `P² = −4h² < 0`, and `R·P = 0`. Both null directions have the chosen future orientation. Wall R, numbered 100, is right to precisely the nineteen support walls and strictly dotted to all other 81 old walls. The point P is future timelike, zero on the support, and strict at every other old wall. It is therefore an actual finite vertex, joined to the old proper cusp by an actual finite-to-ideal edge. The old escape pairs positively with R and is cut.

The support contains a B2 pair with product `−√2`. For the projection inverse only, the checker multiplies root 92 by `√2`, producing a rational positive-definite Gram without changing the subspace. It also verifies orthogonality to every original, unscaled support root.

## A compact right-angled hexagon on S75

The retained positive-definite support is

```
S75 = [1,2,3,4,5,6,7,8,9,10,12,13,14,15,17,31,75,87].
```

Write π for orthogonal projection to its Lorentzian three-plane. Choose

`C = π(2N0 + N65 + 2N79)`, with `d = −C² > 0`.

Exact Gram–Schmidt operations on π(N0) and π(N11) give vectors U,V with Gram `diag(−d,s,t)`, where `s,t > 0`. Thus `C/√d, U/√s, V/√t` is an oriented orthonormal Lorentz frame. The six new normals, numbered 101–106, are

`√2 C/√d + b U/√s + c V/√t`,

with `(b,c) = (2,0), (1,√3), (−1,√3), (−2,0), (−1,−√3), (1,−√3)`.

Every normal has square 2. Adjacent normals pair to zero; nonadjacent products are −4 or −6. These are the sides of a regular right-angled hexagon of inradius `asinh(1)`. All six normals are right to S75 and strictly dotted to every one of the other 83 walls of the intermediate 101-wall chamber, including R.

The certificate gives all six vertex coordinates in this frame. Their time coefficient is 1, their spatial square is 2/3, and hence each has square −1/3. Every vertex is strict at all old nonsupport walls and the four nonincident hexagon sides. The corresponding Klein halfplanes bound the hexagon with these six vertices strictly inside the unit disk. This proves compactness of the actual S75 face, rather than merely giving a cyclic list of compatible walls.

The two roles are distinct. R is strictly dotted to walls 8 and 15 and does not define a boundary on S75; the six hexagon walls do define its boundary and are strictly dotted to wall 61. Thus the construction does not restrict all additions to either of the previously excluded perpendicular families. The certificate also checks that all six hexagon walls retain the old escape, which is cut specifically by R.

The three positive normalization square roots are kept as explicit vector recipes. R, P, C, U, and V are stored in the inherited 32-term coefficient field; pairings with the old 98-wall tower use its 128-term representation. Exact coefficient identities establish all equalities; outward rational intervals establish strict signs and dotted inequalities. All interval decisions succeed with 256-bit bounds. No floating-point acceptance or expansion of the much larger compositum is needed.

## Facets, finite vertices, and retained cusps

The 721 new pair checks split as follows:

| Relation | Count |
|---|---:|
| Right | 133 |
| Strictly dotted | 588 |

A small inward perturbation `H = C + 2^−30 T`, where T is the inherited strict timelike center, is checked directly to be timelike and strict at all 107 walls. The [acute-extension preservation lemma](acute_extension_preservation_20261014.md) then proves actuality of all facets and strict preservation of all old finite vertices.

The complete new-wall finite-neighbor lists make the new vertex census small. R has precisely its nineteen old support neighbors, giving exactly P. Each hexagon wall has precisely its eighteen support neighbors and its two adjacent hexagon walls. Those two neighbors are mutually dotted, so a finite twenty-wall support uses just one of them. This gives exactly the six stated hexagon vertices and no others on a new wall. Each of the seven new vertices has fewer than twenty old active walls, so it is distinct from every old finite vertex.

All seven additions are checked strictly negative at each of the three old bounded cusps. Their full affine links are unchanged. Besides the new finite-to-ideal edge, the hexagon has six actual finite edges.

## The remaining global obstruction

The same construction leaves the old S61 face region unchanged. Its inverse support Gram is entrywise nonnegative, with diagonal entry 2/3 at wall 61. All seven new normals are strictly dotted to 61. For a normal n of square q with nonpositive support pairings, this gives

`(πS61 n)² ≤ q − (2/3)(n·N61)² < −q/3 < 0`.

For the hexagon walls, all other S61 pairings are zero; for R they are nonpositive. Projection of the strict center fixes the future orientation of these timelike projected normals. Therefore all seven inequalities are strict on the entire S61 hyperbolic plane and cannot alter its existing unbounded region. The two roles—escape cutting and S75 closure—have been coupled, but global closure has not been reached.

The certificate stores a new literal rational vector Enew with

`Enew² = 38017800175279 / 500000000000 > 0`.

Every one of its 107 wall pairings is strictly negative, and its pairing with the strict future center is negative. All these signs are rechecked by the exact producer. Thus it gives an open escape through the ideal boundary and certifies infinite hyperbolic volume. Its numerical discovery is not a dependency of the exact proof.

## Reproduction and scope

```
PYTHONPATH=src:. python scripts/fresh_h20_coupled_roles_20261014_exact.py --verify
```

The [primary checker](../scripts/fresh_h20_coupled_roles_20261014_exact.py) uses the [explicit inherited-field helper](../scripts/fresh_h20_coupled_roles_20261014_field.py). Generation and saved `--verify` both pass. It verifies every new pair, the common center, the seven vertices, the three old cusp signs, the renewed escape, and the complete new finite-neighbor lists. The old 100-wall mutual pairs and old finite-vertex lower bound remain inherited, hash-bound results. No old frozen artifact is edited. The separate high-precision changed-core diagnostic with the shorter prefix is numerical only and is not a dependency of this exact construction.

The [independent arithmetic and geometry review](fresh_h20_coupled_roles_20261014_exact_review.md) also passed generation and saved replay. Its separate field arithmetic reconstructs both projections, all 721 new pairs, the seven vertex points, the center and cusp signs, the S61 redundancy bounds, and the renewed escape. It independently verifies the finite-neighbor census on all seven new walls. The old 100-wall mutual-pair and old-vertex certificates remain inherited dependencies.
