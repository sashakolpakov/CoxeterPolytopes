# A common annular escape for the current small-radius construction families

Every finite selection from the families below leaves the same strict spacelike escape. This includes arbitrary replacement and selection of the present quarter-face and central-contact candidates. Their probe-cover searches therefore cannot produce finite volume. The obstruction does not exclude larger endpoint caps or other placements.

The [exact checker](../scripts/fresh_h20_intrinsic_caps_20261014_annular_gate.py) and [certificate](fresh_h20_intrinsic_caps_20261014_annular_gate.json) use the [golden seed coordinates](fresh_h20_intrinsic_caps_20261014_golden_seed.json). The 252 product-alcove vertices have squared norms exactly

`302+4phi`, `326+4phi`, or `446+4phi`.

Thus every vertex z satisfies `17<|z|<22`. Since the Euclidean alcove C is their convex hull, every point of C has norm below 22. Let `x=3z_0/4` for the first vertex. Exact arithmetic gives `13<|x|<14`; all 22 seed inequalities are strict at x.

| Cap family | Rigorous separation from x |
| --- | --- |
| Center any old vertex, radius at most 1 | Distance greater than 3; ball power greater than 8 |
| Center in `(1/2)C`, radius at most 1 | Center norm below 11; distance greater than 2; power greater than 3 |
| Center in `(1/4)C`, radius at most 3/2 | Center norm below 11/2; distance greater than 15/2; power greater than 54 |
| Radius at most 2, meeting the central radius-2 sphere at a finite acute angle or tangency | Center norm at most 4; distance greater than 9; power greater than 77 |
| Center zero, radius at most 2 | Power greater than 165 |

For the fourth row, if the radius is R and the normalized central contact cosine is c with `0<=c<=1`, then the center distance d satisfies

`d²=4+R²+4Rc <= (2+R)²`.

Therefore `d<=4` when `R<=2`. This includes every label and parallel contact, and the entire continuously parametrized central-contact scout, whose radii are `2(1-gamma)<=2`. It is not a bound on central-contact caps of arbitrary radius. The quarter-face assertion covers every point of `(1/4)C`, not merely the sampled barycenters.

In the null-coordinate model, put

`X=(x,(x²-1)/2,1)`.

Then `X²=1`, and the positive ball-power bounds imply strict negative pairing with every sphere normal in these families. The vertical seed pairings are also strictly negative. The known strict center `(0,9/2,1)` remains timelike and strict for every displayed sphere: its height is 3, larger than every permitted radius. The exact pairing `X.center<0` selects the same time orientation. Thus this is a common strict infinite-volume certificate for every finite choice of the displayed caps.

The bounded replacement MILP had reached its 40-round limit without a conclusion; it is archived as incomplete. This independent geometric obstruction supersedes that search as a route to finite volume. It does not retroactively turn its combinatorial feasibility status into a solver proof.

The next construction must cross this annular gap with larger caps, different centers, or a changed seed. Merely adding more small layers or central contacts from the excluded radius range cannot suffice.

Generation and saved replay pass:

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_annular_gate.py --verify
```
