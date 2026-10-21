# Actual arithmetic local triangles obstruct uncompensated discounts

The [exact witnesses](fresh_sqrt5_marked_discount_triangles_20261016.json) and [checker](../scripts/fresh_sqrt5_marked_discount_triangles_20261016.py) exhibit two dimension-21 local chambers in which the proposed long-component discounts fail. Both rules satisfy the exhaustive marked-vertex bound 35, yet on an actual compact triangle they pay only **35/19**, below the required **2**. Generation and fresh saved replay pass.

These are scoped counterexamples to the two uncompensated payment rules. Each chamber has an explicit strict spacelike escape and infinite volume. They are not compact or finite-covolume arithmetic examples, and do not rule out compensated discounts.

## Long distance-two payments

Start with minimal triangle marking 46 in the inherited triangle catalogue. Its five-root diagram has label-3 edges `0–4`, `4–1`, `1–3` and label-5 edge `3–2`. Its black roots are `0,1,2`; roots `3,4` are common white roots. Attach a simply-laced path of 17 new white roots to root 0. There are now 22 roots in a Lorentzian space of dimension 22. Mark the far endpoint of the new path as the distinguished facet root.

The three triangle vertices have finite types

```
H4+A17,       A19+I2(5),       A21.
```

Their original corner payments are `1,0,1`. The `A21` corner has its two black roots at distance two across two label-3 bonds, so the proposed long-component rule reduces that payment by `3/19`. The `H4` component remains short, and the middle corner's black roots lie in separate components. The total is therefore `2-3/19=35/19`.

This identifies a concrete compensation channel: the unchanged `H4` corner sits beside an external `A17` component containing the marked facet root. Alternatively, a bonus could be assigned to the particular disconnected `A19`–`I2(5)` corner, retaining the facet mark and root positions. These contexts differ from the critical `A1+D20` and `D21` vertex types and have slack under the discounted vertex budget.

## Long adjacent payments

Start with black roots `0,1,2` having labels `(2,3,10)` on pairs `(0,1),(0,2),(1,2)`. Attach a simply-laced path of 19 white roots to root 0 and again mark its far endpoint. To keep a literal Gram matrix over `Q(phi)`, use

```
n0^2=n2^2=2,       n1^2=2*(2+phi),
(n0,n1)=0,        (n0,n2)=-1,
(n1,n2)=-(2+phi).
```

Every tail root has norm 2, with adjacent pairing `-1`. The normalized squared product on the label-10 edge is exactly `2+phi=4cos²(pi/10)`.

The three finite vertex types are

```
I2(10)+A19,       A21,       A20+A1.
```

Their baseline corner payments are `1,1,0`. Only the adjacent label-3 pair in the long `A21` component is discounted. Again the total becomes `35/19`. The `I2(10)` corner with an external marked `A19` tail is a specific possible compensation channel.

## Exact geometry and arithmetic scope

For each 22-by-22 Gram matrix, the checker verifies physical signature `(21,1)`, positive-definite conjugate under `phi -> 1-phi`, and integral Cartan entries in `Z[phi]`. The diagram is a tree; all cyclic products reduce to integral edge two-cycles, and a displayed edge cycle generates `Q(phi)`. This certifies the stated local V-arithmetic data.

All three corner principal matrices are positive definite. The literal vectors `-G^-1 e_i` for the omitted black roots are future timelike, mutually have negative pairings, and satisfy every wall inequality with precisely the intended zeros. The common 19-root white submatrix is positive definite. Its three-by-three Schur complement has Lorentzian signature with positive-definite two-by-two corner submatrices, so the selected face is the actual compact triangle bounded by those three remaining walls.

A saved rational perturbation of the sum of the three finite corner vectors supplies a common strict timelike center. All 22 walls are actual facets, verified by their projected-center witnesses. A separate saved perturbation of the last tail root's dual vector supplies a strict spacelike escape with every wall pairing negative. Thus the local chambers have infinite volume.

These witnesses explain why checking only the 76 minimal cores would miss a real obstruction. Compensation must be tested on connected white extensions, with the literal facet root retained. Neither these two counterexamples nor the vertex-budget calculation establishes full triangle or quadrilateral coverage for any repaired rule.
