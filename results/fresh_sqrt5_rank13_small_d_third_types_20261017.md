# Non-ADE third vertices adjacent to the small D-component types

The frozen rank-13 triangle catalogue identifies a small, explicit set of non-ADE vertices that can compensate triangles with a D5+2D4, 2D5+A3, or 3D4+A1 endpoint. All have baseline cost above the dimension-13 threshold of 12, so they cannot simply be treated as low-cost destinations.

| Initial forest | Marked triangle patterns | Distinct third forests | Third-vertex baseline range |
|---|---:|---:|---:|
| D5+2D4 | 19 | 11 | 14–21 |
| 2D5+A3 | 20 | 12 | 13–21 |
| 3D4+A1 | 11 | 6 | 13–19 |

The largest third costs include D4+D5+H4 at 21, A4+D5+H4 at 20, and A1+2D4+H4 and 2D5+H3 at 19. Every third forest has exactly one exceptional component: H4, H3, I2(5), or I2(10). At the actual triangle framing pair an H4 component already pays two baseline units; the other exceptional components pay one, leaving one unit to be shared by the two ADE endpoints.

The [exact source](../scripts/fresh_sqrt5_rank13_small_d_third_types_20261017.py) and [saved extraction](fresh_sqrt5_rank13_small_d_third_types_20261017.json) enumerate every matching record in the [complete third-corner catalogue](fresh_sqrt5_rank13_triangle_third_costs_20261017.md), preserving literal marked-forest codes, both ADE endpoint types, the third forest, its whole-vertex and corner costs, and original witness indices. The other ADE endpoint may be any of the 215 selected types; it is not restricted to the three initial forests.

Generation and a fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_small_d_third_types_20261017.py
python scripts/fresh_sqrt5_rank13_small_d_third_types_20261017.py --verify
```

This is an exhaustive extraction of inherited local necessary patterns. It identifies the concrete non-ADE types that a future joint payment system must address. It neither asserts their compact ambient realizability nor supplies quadrilateral coverage or a new dimension bound.
