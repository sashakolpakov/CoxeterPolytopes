# Independent review of the quartet polar extensions

**PASS.** The [independent checker](../scripts/fresh_h9_quartet_polar_extension_20261014_review.py) and [certificate](fresh_h9_quartet_polar_extension_20261014_review.json) reproduce the complete 19-, 23- and 25-wall cones and every exterior polar test. They also verify the literal extensions' wall pairs, actual facets, finite vertices, finite edges and surviving strict spacelike escape. Both extended chambers have **infinite volume**.

| Retained walls | Finite vertices | Ideal vertices | Hyperideal rays | Compatible exterior polars |
|---:|---:|---:|---:|---:|
| 19 | 53 | 0 | 285 | 4: 37, 38, 95, 186 |
| 23 | 69 | 0 | 413 | 2: 82, 85 |
| 25 | 81 | 0 | 513 | 4: 62, 95, 96, 128 |

There are no opposite-sheet causal rays at any stage. These results agree with the [first polar census](fresh_h9_quartet_polar_20261014.json), [four-cap extension](fresh_h9_quartet_polar_extension_20261014.json) and [second extension](fresh_h9_quartet_polar_second_20261014.json).

## Independent arithmetic and cone completeness

The review imports no producer scalar class or chamber-enumeration engine. It constructs a separate SymPy algebraic field
\[
K=\mathbb Q\bigl(\sqrt5,\sqrt{6+(16/5)\sqrt5}\bigr)
\]
and uses exact field operations, with outward rational intervals at 192 bits to certify signs. Its only imported helpers are previously independent number-field linear algebra and interval primitives. All source manifests are checked recursively.

The cone enumerator is separately implemented. It chooses ten independent wall functionals and works in their evaluation coordinates, where the initial chamber cone is a negative orthant. It then clips this pointed cone by each remaining wall. A crossing pair of extreme rays forms an edge precisely when its common active-wall face contains no third extreme ray. This complete-face test differs from the producer's cached rank-of-common-walls test. Each retained or newly formed ray is normalized by a positive scalar. At the end, the review checks all physical inequalities and independently verifies rank nine for each active set.

These steps give an inductive exhaustive extreme-ray census. They are not merely a replay of supplied ray coordinates. The independently obtained active sets agree exactly with all three saved complete censuses. Every saved polar normal is then compared projectively with its independently generated ray and checked in physical coordinates.

## Complete polar angle test

The relative radicand has norm \(-76/5\), with one positive and one negative real conjugate. Thus \(K\) is quartic with signature \((2,1)\). It has only the quadratic subfield \(\mathbb Q(\sqrt5)\): a second quadratic subfield would make it biquadratic and Galois, incompatible with mixed real and complex embeddings.

For a finite Coxeter angle, \(\cos^2(\pi/m)\) generates a totally real subfield. Degree four cannot occur here, and the rational and quadratic possibilities leave exactly
\[
m\in\{2,3,4,5,6,10\}.
\]
Consequently the angle test is exhaustive, not a finite numerical cutoff. The review independently checks all 1,211 exterior candidates, their full saved pairings, forbidden-index lists, norms, independent supports, center signs and escape pairings.

## Extension geometry

All 300 wall pairs of the final 25-wall chamber are independently reconstructed. They comprise 164 right pairs, 14 at \(\pi/3\), 10 at \(\pi/5\), six at \(\pi/10\), and 106 strictly dotted pairs. The intermediate 23-wall Gram and its 82 new pairs are also checked, followed by all 47 pairs introduced in the second extension.

The inherited center is strictly timelike and strict for every wall. All 23 and 25 displayed facet points are independently checked to be the orthogonal projections of that center, timelike and strict for every other wall. Thus every retained wall is an actual facet.

The finite-vertex active sets are exhaustive because of the independent cone census. Every finite point has exactly nine incident walls; the review also verifies its positive-definite incident Gram and all physical pairings. All 53 earlier finite vertices survive the four-cap stage, which adds 16. All 69 then survive the second stage, which adds 12. The complete finite-edge and unfinished-direction lists are independently reconstructed from the actual finite vertices:

| Walls | Finite edges | Unfinished finite-corner directions |
|---:|---:|---:|
| 23 | 213 | 195 |
| 25 | 247 | 235 |

The same exact spacelike point satisfies every inequality strictly and has negative pairing with the common future center. This supplies a strict infinite-volume escape for both extensions.

At the 25-wall stopping point, the four further compatible polars have all six mutual relations strictly dotted. They each miss the same escape. This is verified frontier data; those four walls are not added by this review.

## Replay and scope

Independent generation and saved replay pass:

```sh
python scripts/fresh_h9_quartet_polar_extension_20261014_review.py --verify
```

The result is complete for exterior polar hyperplanes of the three literal saved cones. It does not exclude nonpolar caps, changed retained walls, or arbitrary later completions. It does not assert a new finite-volume polytope or a new nonarithmeticity proof; any inherited nonarithmeticity discussion belongs to the earlier source construction.
