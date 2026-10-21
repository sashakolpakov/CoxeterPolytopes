# Independent review of the binary-cap finite-vertex obstruction

The local theorem is correct: **no finite vertex of the frozen chamber is incident to any of its 524,288 binary caps**. This does not exclude ideal vertices or finite vertices produced after adding different walls. It is a local structural obstruction, not by itself a finite-volume impossibility theorem.

Reviewed inputs are the [local producer](../scripts/fresh_h20_weighted_binary_20261014_local.py), [local certificate](fresh_h20_weighted_binary_20261014_local.json), and the [orbit certificate](fresh_h20_weighted_binary_20261014.json). The separate [review checker](../scripts/fresh_h20_weighted_binary_20261014_local_review.py) and [review record](fresh_h20_weighted_binary_20261014_local_review.json) use integer arithmetic and no producer or project imports. Independent generation, independent saved replay, and the primary local saved replay all **PASS**. Source and certificate hashes are bound in the review record.

## Independent finite-neighbor graph

The checker reconstructs the literal radius-one centers using twice their normalized coordinates, which are integers. Four times each squared distance is an integer combination `a+b*phi`. For equal-radius caps, the finite-contact distances are exactly `2`, `3`, and `2+phi`, corresponding to labels two, three, and five. The checker exhausts every nonzero sign mask and verifies that every other distance is at least four. Thus the finite-neighbor list is independently complete, without trusting the producer's bounded mask-generation loop.

There are exactly 95 binary neighbors of the reference cap. Its finite seed neighbors are the 19 coordinate walls on its own sign side: three right, seven at label three, one at label four, and eight at label five. The opposite seed walls are dotted. The central sphere is dotted because its normalized cosine is `(9+8phi)/16>1`. Every vertex-centered endpoint sphere is dotted: the minimum squared distance is `17/4+2phi>6>(1+sqrt(2))²`.

The checker independently reconstructs all pairs among these 114 neighbors from coordinates and seed gaps. It recovers exactly the saved **2,644 edges** and their labels, entry for entry. This verifies the full relevant local graph, while the complete chamber's global Coxeter and actual-facet assertions remain dependencies of the separate orbit certificate.

## The spherical restriction needs no classification catalog

Suppose a positive-definite incident diagram had two label-five neighbors of the reference root. Normalize the three roots to square one and put `c=cos(pi/5)=phi/2`; write the pairing between the two neighbors as `-d`, where `d>=0`. Their principal determinant is

```text
1 - 2c² - d² - 2c²d <= 1 - 2c² = (1-phi)/2 < 0.
```

This contradicts positive definiteness. Hence at most one label-five neighbor can occur. The premise cited in the producer follows directly from this minor calculation, independently of the classification of finite Coxeter diagrams.

There are 42 non-label-five neighbors and 72 label-five neighbors. The 73 certificates handle either no label-five neighbor, or one chosen label-five neighbor together with its non-label-five common neighbors. The independent checker checks every coloring's exact node set, every graph edge, every integer color, and the resulting bound. It does not run or trust the greedy coloring heuristic. A deliberately invalid monochromatic certificate is rejected.

The bounds, including the reference root and the optional label-five root, are:

| Full clique bound | Cases |
| --- | ---: |
| 18 | 1 |
| 13 | 64 |
| 19 | 8 |

At a finite H20 vertex, incident normals span the 20-dimensional positive-definite tangent space. If the reference root is incident, it can be extended to a basis of 20 incident normals. These form a positive-definite finite-pair clique containing that root. Every such clique has at most 19 roots by the certificates, a contradiction. This argument does not need an unproved simplicity assumption about vertices.

## Transitivity is an actual symmetry of the entire chamber

The two full-coordinate triples are `A={0,1,2}` and `B={0,3,4}`. Their choice depends on sign parity in the ten unit and eight golden coordinates; the heavy coordinate is omitted from parity.

An even-parity sign flip preserves the triple type and maps the cap family to itself. For an odd-parity target sign mask, first apply the unit-coordinate permutation `(1 3)(2 4)`, then the target sign flip. The permutation exchanges `A` and `B`, preserves the metric weights and parity, and the odd sign flip exchanges the required type again. Thus this map preserves the family rule and sends the reference cap to the specified target cap.

Both operations also preserve all signed coordinate walls, the central sphere, and the complete signed endpoint orbit. Consequently they are symmetries of the full chamber, rather than merely coincidences of pair-count statistics. The checker verifies the permutation identities and the type/sign rule for every target mask. The reference-cap obstruction therefore applies to every binary cap.

No mathematical or source issue was found. The useful next construction must supply genuinely new finite neighbors if it is to create finite vertices on these caps; the existing million-facet orbit alone supplies none there.

```sh
python3 scripts/fresh_h20_weighted_binary_20261014_local_review.py --verify
```
