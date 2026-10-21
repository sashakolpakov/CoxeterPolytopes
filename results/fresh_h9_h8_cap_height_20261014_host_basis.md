# The arithmetic host is an E8 lattice plus a binary form

The [free-basis certificate](fresh_h9_h8_cap_height_20261014_host_basis.json) identifies the full fifteen-generator reflection module, without enlarging or saturating it:

`L = E8 ⊥ [[−2φ,−1], [−1,2]]` over `O=Z[φ]`.

In the source coordinates, a basis consists of the first nine coordinate vectors and `w=((φ−1)/2)e₀−e₉`. The [checker](../scripts/fresh_h9_h8_cap_height_20261014_host_basis.py) stores explicit integral matrices for both inclusions: every old generator is an `O`-linear combination of these ten vectors, and each basis vector is an `O`-linear combination of the old generators. Thus this is a basis of the actual module, not merely its vector space.

The determinant is `−3−2√5=−φ(3+φ)`, with field norm `−11`. A vector with coordinates `(x,y,t)` has norm

`−2φ x²−2xt+2t²+E8[y]`.

This provides direct integral coordinates for root enumeration. It does not establish reflectivity or add any wall to the chamber.

```sh
python scripts/fresh_h9_h8_cap_height_20261014_host_basis.py --verify
```
