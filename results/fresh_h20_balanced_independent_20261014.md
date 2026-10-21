# Independent replay of the balanced cusp endpoint obstruction

The [separate checker](../scripts/fresh_h20_balanced_independent_20261014.py)
and its [saved result](fresh_h20_balanced_independent_20261014.json) pass.
This is an independent arithmetic implementation of the
[endpoint census](fresh_h20_balanced_endpoint_20261014.md) and
[central-cap obstruction](fresh_h20_balanced_endpoint_20261014_central_gap.md).
It imports no geometry or arithmetic functions from either producer.

The checker reconstructs the three affine cycle Gram matrices, inverts the
finite rank-19 block with FLINT rational matrices, and enumerates all 723
reduced contact vectors. It obtains the same 37 positive-definite and four
semidefinite endpoint profiles. The inverse Gram is entrywise nonnegative;
its diagonal lower bound excludes three nonright contacts and every retained
label at least five. Thus the finite enumeration covers all retained labels.
It also checks all 392 literal cyclic transports.

Both scale choices are recomputed in the field with basis
`1,sqrt2,sqrt5,sqrt10`. All 82 profilewise radius bounds, including every
exceptional omitted-wall right or pi/3 radius, agree with the producer.
Exact coefficient arithmetic uses FLINT rationals; every nonzero sign is
certified by a 192-bit Arb interval. Algebraic zero is handled separately.
The final Cauchy bounds reproduce the strict norm-1/4 escape for the central
cap together with any endpoint-cap selection.

This proves an obstruction for the stated cap family. It does not exclude
the cusp geometry, boundary-face caps, changed walls, or general H20 polytopes.

```sh
python scripts/fresh_h20_balanced_independent_20261014.py --verify
```
