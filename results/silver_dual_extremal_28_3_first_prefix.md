# Silver `DualExtremal(28,3)` first necessary prefix

## Result

For

```text
L = <-(2*sqrt(2)-1)> + DualExtremal(28,3)
```

in `H28`, the complete mirror-norm alphabet has eight classes. Exact controls
through `H=93 + 66*sqrt(2)` contain 51
controls on 35 heights. Every earlier control
is empty in the rigorous spatial Gram superset. The endpoint has the following
two norm components and six necessary Gram types.

| mirror norm | time control | `(u^2,v^2,u.v)` | max `d.y` | margin for `d/2` |
|---:|---:|---:|---:|---:|
| 3 | 9 + 6*sqrt(2) | (138, 72, 99) | 24 + 18*sqrt(2) | 3 + 3*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (300, 132, 198) | 45 + 30*sqrt(2) | 3/2 |
| 6 | 12 + 9*sqrt(2) | (288, 138, 198) | 44 + 31*sqrt(2) | 2 - 1/2*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (276, 144, 198) | 44 + 32*sqrt(2) | 2 - 1*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (264, 150, 198) | 40 + 31*sqrt(2) | 4 - 1/2*sqrt(2) |
| 6 | 12 + 9*sqrt(2) | (252, 156, 198) | 33 + 27*sqrt(2) | 15/2 + 3/2*sqrt(2) |

The lattice has determinant `3^26`, minimum 8 and 1,260 minimum vectors. Its
mod-three pairing-radical preimage has index 9, exact minimum 12 and 65,520
minimum vectors; these data are independently replayed with PARI.

For any minimum vector `d`, the section point z=(1,d/2) for any spatial minimum vector d with d^2=8 has

```text
z^2 = 3 - 2*sqrt(2) > 0.
```

Every displayed margin is strict even against the full necessary augmented
Gram superset. Therefore the complete prefix through the stated height has an
open hyperideal patch and infinite volume.

## Scope

This is a bounded negative screen. It neither realizes every endpoint Gram
type nor decides later Vinberg walls, termination, reflectivity, or the final
chamber. It is not a compact example and not a global nonreflectivity theorem.

Source: [Nebe--Sloane catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/DualExtremal_28_3.html).
