# One cap cannot close the remaining two-ray B component

In the [198-vertex replacement47 chamber](fresh_nonarithmetic_replace47_close12_20261017.md), the remaining B component has the two exterior rays `[7,30,47]` and `[30,47,52]`. Its four geodesic boundary edges are7–30,7–47,30–52 and47–52. **No single Coxeter wall can supply finite terminal vertices on all four**, regardless of its coefficient field or finite labels.

At a finite terminal vertex on7–30, the old support `S6+[7,30]` has diagramH4+H4. The finite Coxeter classification permits no nonorthogonal one-root extension of either H4. Thus the proposed wall is automatically right to S6,7 and30. This also explains why allowing a nonzero A2 profile cannot evade this particular single-cap argument.

Normalize its square to2. In the Lorentz4 complement of S6, let `e7=n7+(n5+2*n6)/3`, of square4/3. The complete frame `(e7,n30,n47,n52)` has an exactly invertible Gram. Write `a=-r.n47` and `b=-r.n52`. Finite endpoints on7–47 and30–52 require

```
0 <= a² < 5/2,
0 <= b² < 2.
```

These follow directly from the two2×2 inverses: with the contact to7 or30 zero, the respective Schur costs are `4a²/5` and `b²`.

The four frame products determine r uniquely. Its forced norm is the quadratic form obtained from the lower-right2×2 block of the inverse frame Gram. Both diagonal coefficients are strictly positive, so it is convex in either variable separately. Its maximum on the relaxed closed rectangle `[0,sqrt(5/2)]×[0,sqrt(2)]` lies at one of the four corners. The exact outward enclosures give upper bounds approximately

```
0, 0.88492480, 1.75122132, 0.14281495.
```

Every corner is strictly below `219/125 < 2`. Thus r cannot have square2. The checker additionally verifies the inverse Gram of7–30 has positive entries and both diagonals greater than2; even without invoking H4², this forces its two contacts to be zero once rightS6 is assumed.

Generation and fresh saved replay passed for the [source](../scripts/fresh_nonarithmetic_replace47_single_cap_gate_20261017.py) and [certificate](fresh_nonarithmetic_replace47_single_cap_gate_20261017.json). The proof uses the frozen exact geometry and shared outward rational interval engine. It is a primary certificate, not a separate independent arithmetic review.

This excludes one finite-terminal wall for the whole indicated B patch. It does not exclude multiple caps, changed retained walls, or finite-volume completion by another mechanism. No common-center rejection is used in the norm argument.
