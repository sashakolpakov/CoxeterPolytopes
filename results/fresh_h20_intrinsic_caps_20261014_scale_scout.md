# Changing the three factor scales

The best elementary covering ratio among central labels `(m_A,m_T,m_E)` in `{3,4,5,6}³` is `(3,6,3)`, giving inradii `(1,√3,1)` when the central sphere has radius two. However, **rescaling this seed leaves the same quarter-vertex gap in the central/vertex/half-vertex construction**. Changing the long A11 factor is a more substantial next direction.

The [small exact checker](../scripts/fresh_h20_intrinsic_caps_20261014_scale_scout.py) reconstructs the factor vertex geometry from the frozen intrinsic seed and checks all 64 choices in Q(phi). Its [certificate](fresh_h20_intrinsic_caps_20261014_scale_scout.json) gives these unit-inradius data:

| Factor | Vertex squared distances from incenter | Minimum squared separation |
|---|---|---:|
| A11 | 286 at every vertex | 264 |
| A2 | 4 at every vertex | 12 |
| E6 | 12, 36, or 156 | 48 |

Write the inradii as `a,t,e`. They lie between one and √3. The minimum product-vertex separation is always `√12 t`. The nonpositive-pairing condition gives the necessary equal-radius bound `R²≤6t²`; it does not guarantee that all other pairings are Coxeter. Separately, the smaller radius `R²=3t²` gives an explicit Coxeter layer: nearest triangle-neighbor caps are parallel, every other cap pair is dotted, and all three omitted seed contacts are dotted. This is a feasible layer, **not a maximum**. Central compatibility is strict. The 252 radius-one half-vertex caps remain compatible with that parallel layer: their retained seed angles are exactly the chosen central labels; nearest same-factor pairs are pi/3 when `t=1` and dotted otherwise. All cross-layer pairs are dotted. These assertions follow from the archived factor distances and the sphere-pair formula, rather than a repeated full cap census.

The ratio `(286a²+4t²+12e²)/(6t²)` measures the nearest product vertex's squared distance relative to the necessary acute radius bound. Increasing `a` or `e` worsens it; increasing `t` improves it. This gives `(3,6,3)` uniquely among the 64 choices. It is a geometric benchmark, not a claim of feasible maximal radii or optimal global volume.

Let `z` be the vertex using omissions `(4,23,22)` and let `x=z/4`. The checker proves that its nearest half-vertex center is `z/2` and its nearest vertex center is `z`. Put `Q=(286a²+4t²+12e²)/16`. Its ball powers are at least `Q−4` for the central sphere, `Q−1` for every radius-one half-vertex sphere, and `9Q−6t²` for every equal-radius vertex layer satisfying acuteness. All are strictly above `1/8` in every one of the 64 cases. Every vertical seed inequality is strict at `x`. Assigning formal height squared `−1/8` therefore gives a strict future spacelike escape of square `1/8` for each finite template, even with this overinclusive radius bound.

At `(3,6,3)`, the explicit parallel vertex layer has radius three. At the larger necessary acute radius bound `R²=18`, the three powers are still at least `123/8`, `147/8`, and `1251/8`. Thus even the radius relaxation leaves a substantial gap. This certificate does not constrain different centers, arbitrary varying-radius selections, or the intrinsic cusp seed as a whole. In particular it does not exclude the interior cap mechanisms being studied separately. The initial draft confused the parallel-layer radius with the acute upper bound; this corrected certificate keeps them separate.

```sh
python scripts/fresh_h20_intrinsic_caps_20261014_scale_scout.py --verify
```
