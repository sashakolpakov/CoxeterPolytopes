# Independent check of the exceptional-cap height family

The [independent checker](../scripts/fresh_h9_h8_cap_height_20261014_review.py) uses a separate SymPy quadratic field to reconstruct the [maximal-height seed](fresh_h9_h8_cap_height_20261014.json). The [saved review](fresh_h9_h8_cap_height_20261014_review.json) records PASS.

All 105 wall pairs have the claimed relations: 79 right, twelve pi/3, five pi/5 and nine dotted. Fifteen strict facet points establish actual facets. Independently filtering all 5,005 nine-wall subsets, then solving physical nullspaces, gives exactly 62 finite vertices and 220 finite edges. This differs from the producer's source-ray intersection census.

The exact spacelike escape has square `4363/256−3057√5/512>0` and pairs strictly negatively with every wall. Each normalized source direction has minimum nonzero contact one, forcing its positive cap norm `q≤2`. After normalizing the common vertical form coefficient to `D=(7+√5)/4`, a cap has last coordinate `−η`, with `η=√((a+q)/D)≤1` and `a=(√5−1)/4`. Its escape pairing differs from the maximal-height pairing by `5D(η−1)≤0`. The same strict escape therefore survives every independent choice of the three heights and every subfamily of these three directions. Additional projected directions remain outside this theorem.

The checker also verifies the source contact alphabet and exact pair margins used by the all-label height classification. The conjugate form is positive definite and `Norm(D)=11/4`; these facts alone make no claim of reflectivity.

```sh
python scripts/fresh_h9_h8_cap_height_20261014_review.py --verify
```
