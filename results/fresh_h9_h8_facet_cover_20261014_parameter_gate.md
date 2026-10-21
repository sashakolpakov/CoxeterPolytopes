# Exact diameter bounds along the cap segment

For the fixed mixed twelve-wall H9 chamber, let `z(t)=(1−t)H+tE` use the frozen center and rational escape. Consider its all-dotted caps, with a perpendicular copy of the fixed eleven-facet compact Bugaenko H8 polytope as the new facet. Acuteness and finite volume impose the exact necessary restriction

`0.99674949 < t < 0.996752404`.

The [checker](../scripts/fresh_h9_h8_facet_cover_20261014_parameter_gate.py) and [certificate](fresh_h9_h8_facet_cover_20261014_parameter_gate.json) use exact quartic arithmetic and outward-rounded 128-bit dyadic intervals. They reconstruct all 41 source vertices from their literal supports and bound all 820 vertex-pair distances. The source diameter is approximately 4.24510024762. The [numerical scout](fresh_h9_h8_facet_cover_20261014_parameter_scout.json) motivated these endpoints; it is not used by the proof.

Write `q=z²`, `g_i=(z,n_i)`, and `γ_i=−g_i/√(q n_i²)>1`. For every wall the numerator `2g_i′q−g_iq′` is linear in t and strictly positive from the exact rational lower bracket of the null crossing through t=1. Thus every γ decreases wherever the cap is all-dotted. All ball radii `ρ_i=asinh(1/√(γ_i²−1))` increase. For distinct walls,

`cosh d(v_i,v_j)=(γ_iγ_j+c_ij)/(√(γ_i²−1)√(γ_j²−1))`,

where `c_ij=−(n_i,n_j)/√(n_i²n_j²)≥0`. Its partial derivative with respect to γ_i has the sign of `−γ_j−c_ijγ_i`, so every center distance increases too. This proves the needed monotonicity on the whole all-dotted interval.

At the lower endpoint, all 78 quantities `d(v_i,v_j)+ρ_i+ρ_j` are strictly smaller than the distance between source vertices 32 and 33. A finite-volume attachment requires the source inside the closed ball union, whose diameter is the maximum of those 78 quantities. The entire earlier interval is therefore excluded. At the upper endpoint, centers 3 and 7 are farther apart than every pair of source vertices. Acuteness requires both centers inside the source, excluding that endpoint and all later all-dotted caps. Diameter is maximized at vertices of a compact hyperbolic polytope by convexity of distance along geodesics.

Every comparison uses algebraic formulas for the cosh of the sum; no floating-point inverse hyperbolic function enters acceptance. Arithmetic libraries are shared with frozen prior certificates. These are necessary placement conditions, not a constructed placement or a sufficient coverage test.

```sh
python scripts/fresh_h9_h8_facet_cover_20261014_parameter_gate.py --verify
```
