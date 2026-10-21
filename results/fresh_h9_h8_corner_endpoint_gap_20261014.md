# All source-corner right endpoint caps miss the same escape

For the changed fifteen-wall H9 chamber, **every forward norm-two cap right to the eight walls at any source vertex misses the saved strict escape**, for every real scale compatible with Coxeter angles. This covers all 41 source-corner directions, including the 38 unfinished vertical edges. It does not exclude nonright support profiles or extensions on later supports.

The [checker](../scripts/fresh_h9_h8_corner_endpoint_gap_20261014.py) and [certificate](fresh_h9_h8_corner_endpoint_gap_20261014.json) verify all 41 scalar margins in two separate arithmetic implementations: rational quadratic/quartic arithmetic and a SymPy algebraic field. They inherit the independently reviewed source-vertex census from the [changed chamber](fresh_h9_h8_disjoint_cap_20261014.md).

Normalize each source vertex ray `v` so that its smallest nonzero contact with a source wall is one. Set `b=-v²>0`, `d=sqrt(5)/4`, and write the saved escape as `E=(w,-10)`. A norm-two forward cap right to those eight walls has the form

`r=(s v, -sqrt((2+b s²)/d))`.

Its nonzero source contacts include `s`. Since the smallest nonzero norm-two Coxeter contact is one, `s>=1`. Put `g=(v,w)<0`. For every one of the 41 source rays the exact calculation gives

`g² - 100 d (2+b) > 0`.

Consequently, for every real `s>=1`,

`(r,E) = s g + 10 sqrt(d(2+b s²)) <= s [g + 10 sqrt(d(2+b))] < 0`.

Thus this entire family leaves the same escape inside. No finite label bound, cap-height sampling or old-cap compatibility assumption is needed. Such caps may still be useful intermediate steps that create new supports; the theorem excludes only their ability to cut this particular escape directly.

```sh
python3 scripts/fresh_h9_h8_corner_endpoint_gap_20261014.py --verify
```
