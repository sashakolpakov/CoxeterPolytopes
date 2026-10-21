# A radius-four replacement cuts the annular gap

The [exact replacement certificate](fresh_h20_intrinsic_caps_20261014_radius4_replacement.json) changes one facet of the [587-wall quarter-face chamber](fresh_h20_intrinsic_caps_20261014_quarter_face_exact.md). It removes the radius-one cap at the first old vertex and inserts a radius-four endpoint cap. The new chamber again has **587 actual facets**, preserves the other 251 certified finite vertices, and has a new finite endpoint replacing the first one. It cuts the [annular escape](fresh_h20_intrinsic_caps_20261014_annular_gate.md), but the earlier eighth-vertex escape survives, so its volume remains infinite.

This is endpoint-bank recipe `[0,1,4,"finite"]` from the [bounded numerical proposal](fresh_h20_intrinsic_caps_20261014_annulus_endpoint_scout.json). The first number is transport index0, whose omitted walls are `[4,23,22]` and whose endpoint type is2; it is not endpoint type0. The exact checker does not certify the rest of that numerical bank.

Let `z` be the first product-alcove vertex. On its nineteen incident unit seed roots, solve

```text
(y,v_12)=1/2,      (y,v_i)=0 for the other eighteen roots.
```

The solution is rational in the spatial basis and has `y²=5/6`. Define

```text
c=z-4y,
r=(c,(c²-16)/2,1),      r²=16.
```

The removed cap is prior wall index23. Its pairing with r is `11/6>0`, so this is a genuine replacement, not a compatible addition. Against all 586 retained walls the exact pair census is:

| Relation | Pairs |
|---|---:|
| Right | 18 |
| pi/3 | 1 |
| Parallel | 1 |
| Dotted | 566 |

The pi/3 contact is old source wall12. The parallel contact is old source wall22, global index15. A [separate exact diagnostic](fresh_h20_intrinsic_caps_20261014_parallel_cusp_gap.md) proves that its new cusp has quotient rank16 rather than19, hence three-dimensional link lineality. This new cusp is unbounded.

The [Fraction checker](../scripts/fresh_h20_intrinsic_caps_20261014_radius4_replacement.py) binds the frozen construction inputs and replays every new relation. The common strict center is unchanged. On the original cusp edge, the new endpoint is

```text
X=(z,(z²+8/3)/2,1),      X²=-8/3.
```

It lies on exactly nineteen seed facets and r; all other inequalities are strict. The full rank20 incident Gram has positive rational LDL pivots. This proves an actual finite endpoint and actual new facet. The remaining 251 old endpoint points are strict on r. The old central, half-vertex and quarter-face facet points are also strict on r, and every seed facet still occurs at a retained finite vertex. Therefore all other 586 facets survive.

The norm-one annular escape at spatial position `3z/4` is cut with exact pairing

```text
115/48-phi/8 > 0.
```

The old norm-`1/2` eighth-vertex escape remains strict on every wall, with new-cap pairing

```text
-15173/192-49phi/32 < 0.
```

This replacement therefore escapes the radius-one vertex-cap restriction in the annular-family obstruction, but does not establish global closure. The sphere still has only a dotted contact with the central cap; changing the remaining covering components is a separate obligation.

```sh
python3 scripts/fresh_h20_intrinsic_caps_20261014_radius4_replacement.py --verify
```

Generation and saved replay pass. The [root's independent replacement review](fresh_h20_intrinsic_caps_20261014_replacement_review.md) checks the new geometry in separate arithmetic. There is no claim of a complete vertex or cusp census, or of finite volume.
