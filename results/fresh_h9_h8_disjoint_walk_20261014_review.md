# Independent check of the four-cap extension

The [independent checker](../scripts/fresh_h9_h8_disjoint_walk_20261014_review.py) reconstructs the [exact quartet](fresh_h9_h8_disjoint_walk_20261014_exact.json) in a separate SymPy number field. It does not import the producer or its extension-arithmetic class. The [saved review](fresh_h9_h8_disjoint_walk_20261014_review.json) records PASS.

All 171 wall pairs have the claimed exact relations. The 66 new pairs comprise 31 right, six pi/3, eight pi/5 and 21 dotted pairs. A strict timelike interior point and nineteen strict facet points establish that every wall is actual.

The checker filters all 92,378 nine-wall subsets independently of the producer's clique generator. Exactly 799 are finite-neighbor cliques; 746 fail positive definiteness. Physical nullspaces of the remaining supports produce exactly 53 future finite vertices, with the saved incidences and projective points. Their common eight-wall supports give all 175 finite edges. Six exact positive forward parameters independently verify the old-edge closures.

The original strict spacelike escape also passes every new inequality strictly. Thus this chamber has infinite volume. The review establishes the complete finite-vertex and finite-edge censuses, not an ideal/exterior vertex census or a finite-volume completion.

```sh
python scripts/fresh_h9_h8_disjoint_walk_20261014_review.py --verify
```
