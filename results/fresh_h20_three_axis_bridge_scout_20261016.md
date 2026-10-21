# Eight three-axis bridge patterns after the E6 cutter

This is a bounded numerical diagnostic on the exact
[157-wall chamber](fresh_h20_release144_cap_20261016.md), not a new
promoted Coxeter extension. Of eight prescribed patterns, six pass all157
numerical angle/coherence gates. The best two alternatives have finite
corner patterns and trim more of the dominant residual region than the
all-pi/3 proposal. They are mutually obtuse and cannot be used together.

The [source](../scripts/fresh_h20_three_axis_bridge_scout_20261016.py)
and [data](fresh_h20_three_axis_bridge_scout_20261016.json) pass generation
and a fresh saved replay. No infinite-label completeness, exact new-cap
verification or proper-cusp completion is claimed.

## Finite prescribed family

Use the orthonormal residual axes `e141,e142,e143,e144,e0`, with e0
future timelike. The existing dominant-axis caps154,153,155 have boosts
`s_i=(s_star,s153,s_star)` and `c_i=sqrt(1+s_i^2)`. Choose
`m_i` independently in `{2,3}`, and put `t_i=cos(pi/m_i)`.
For a proposed cap

```
r=sqrt2*(-a141*e141-a142*e142-a143*e143+s*e0),
a_i=(s*s_i-t_i)/c_i,
sum_i a_i^2-s^2=1,
```

the last equation is quadratic in s. The positive root gives norm2 and
the three prescribed products. The cap is right to S16 and144. Both
failed patterns, `(2,2,2)` and `(2,3,2)`, have a non-Coxeter contact
with142, numerically at labels6.66198 and4.63301 respectively.

Every accepted pattern has negative common-center product and is
compatible numerically with the exact E6-cutting wall156. All remain
right144, so they do not cut E6 by the separate all-real gate; they are
intended to supplement wall156, not replace it.

## Common-region comparison

All scores use the same fixed C chart, all retained inequalities and the
same2^20 scrambled Sobol sphere directions. The starting157 estimate is
`.0114396643157`, with maximum squared vertex radius1.417249978.

| Labels at154,153,155 | New interior corners | New near-sphere corners | Estimated exterior4-volume |
|---|---:|---:|---:|
| 2,2,3 | 2 | 0 | .0100402298 |
| 2,3,3 | 1 | 1 | .0102649927 |
| 3,2,2 | 2 | 0 | .0100399688 |
| 3,2,3 | 2 | 0 | .0104649961 |
| 3,3,2 | 1 | 1 | .0102648664 |
| 3,3,3 | 0 | 2 | .0105054786 |

All six reduce the maximum squared residual radius to1.338894933.
The table's corner classifications are numerical and use a1e-6 sphere
band. In particular, the near-sphere points of the all-pi/3 proposal
come from affine A2 triangles with additional spherical components;
they are not certified proper cusps. A full bounded Euclidean link
would still be required before such a claim.

For `(2,2,3)`, the two prospective finite corners have residual supports
`{144,153,154,new}` and `{144,153,155,new}`; their diagrams are
`A2+2A1` and `A3+A1`, respectively. The `(3,2,2)` case interchanges
the two dominant axes. These are preferable finite-corner targets if
an exact bridge is pursued.

## The two strongest alternatives cannot be combined

The two roots for `(2,2,3)` and `(3,2,2)` have the same s. Their axis
coefficients differ by swapping a decrement `delta=1/(2*c_star)` between
the141 and143 coordinates. Their norm equations therefore give the
exact identity

```
<r223,r322>/2 = 1-delta^2 = 1-1/(4*c_star^2)>0.
```

The normalized product is approximately.8810145932, so the pair is
obtuse. The source also tests every pair among the six accepted roots;
all fifteen are numerically obtuse. No compatible batch is supplied by
this eight-pattern family.

For comparison only, intersecting both strongest halfspaces gives a
geometric exterior-volume estimate `.00996493893`, with the same maximum
squared radius1.338894933. **That union is inadmissible as an acute
Coxeter extension.** The source retains this score explicitly as a
diagnostic, not as a proposed construction.

These are Euclidean volumes in a fixed residual chart, not hyperbolic
volumes. No rigorous integration error bound or exact complete residual
cone census is asserted. The exact157-wall escape E7 remains a separate
certificate of infinite volume for the currently promoted chamber.

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_three_axis_bridge_scout_20261016.py --verify
```
