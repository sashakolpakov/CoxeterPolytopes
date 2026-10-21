# All-scale right-contact endpoint gap at maximal cap height

Every forward cap perpendicular to the eight source walls at any original H8 vertex misses the [maximal-height seed's escape](fresh_h9_h8_cap_height_20261014.md). This holds for every real scale and every Coxeter label. Such caps may still be useful as intermediate extensions.

Normalize a source-vertex ray `v` to minimum nonzero source contact one. Write `b=−v²`, `g=v·E_source<0`, and `D=(7+√5)/4`. A norm-two forward cap perpendicular to those eight walls has the form

`r=(s v,−√((2+b s²)/D))`.

Compatibility with the remaining source walls forces `s≥1`. The [exact checker](../scripts/fresh_h9_h8_maxheight_corner_gap_20261014.py) checks all 41 strict margins

`g²−25D(2+b)>0`.

It follows that `r·E=s g+5√(D(2+b s²))≤s[g+5√(D(2+b))]<0`. The [certificate](fresh_h9_h8_maxheight_corner_gap_20261014.json) independently reconstructs each scalar margin in Fraction arithmetic and a separate SymPy number field. Completeness of the 41 source vertices is inherited from the reviewed source census.

This theorem concerns only right contact with eight original source walls. Nonright profiles and later supports remain open. It does not exclude adding such caps on the way to a larger completion.

```sh
python scripts/fresh_h9_h8_maxheight_corner_gap_20261014.py --verify
```
