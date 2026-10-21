# Independent exact review of the two separate octic caps

The [independent checker](../scripts/fresh_h9_octic_escape_caps_20261015_review.py) and [certificate](fresh_h9_octic_escape_caps_20261015_review.json) replay the [two cap certificates](fresh_h9_octic_escape_caps_20261015.json) using a separate SymPy quartic field, rational outward sign bounds and an additional quadratic extension. No producer arithmetic is imported. Generation and saved replay pass; no issue was found.

Both normals are independently reconstructed as `p-sqrt(delta)v` from the frozen nine-contact data. Their squares are exactly two. All 58 new pairings are checked with the appropriate [complete finite-label alphabet](fresh_h9_octic_angle_20261015.md), using the actual old wall norms. The first cap has seven right, two label-3 and twenty dotted relations. The second has seven right, one label-3, one label-5 and twenty dotted relations. Each retains the strict timelike center and cuts the previous escape.

For each separate 30-wall chamber, all 30 facet points are reconstructed and checked against all 30 walls: 1,800 exact signs in total. The saved endpoint is future timelike and has exact zero set `[0,1,3,5,6,7,10,27,29]`, with every other wall strict. A reversed-order nine-root Schur calculation is positive definite. Thus this is an actual finite endpoint on the specified old eight-wall face, rather than merely an accepted angle profile. No complete 30-wall vertex or cone census is claimed.

The two new walls cannot be retained together. Their mutual product is independently reconstructed as

`-8-4sqrt(5) + (-1/2+3sqrt(5)/10)sqrt(6115+2735sqrt(5)) > 0`.

Finally, the common renewed escape is checked in physical coordinates: its square is strictly positive, its center pairing is negative, and all 29 old-wall pairings and both separate added-wall pairings are strictly negative. It therefore proves infinite volume for each of the two 30-wall chambers. The construction gives actual local endpoint progress, with this global limitation explicit.

```sh
python scripts/fresh_h9_octic_escape_caps_20261015_review.py
python scripts/fresh_h9_octic_escape_caps_20261015_review.py --verify
```
