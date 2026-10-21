# Silver odd-unimodular family in H27--H29 through H1000

## Outcome

Three genuinely new rootless target-rank hosts are now exact:

| spatial lattice | hyperbolic dimension | norm-three shell | first packet walls |
|:--|--:|--:|--:|
| `T27` | 27 | 1,640 | 119,574 |
| `U28` | 28 | 2,240 | 98,280 |
| `dim29odd` | 29 | 1,856 | 93,786 |

For each host

\[
 L_M=\langle-(2\sqrt2-1)\rangle\oplus M,
\]

the complete six-class local arithmetic has no height-zero wall. The first
packet at `H=5+4sqrt(2)` is exactly the full norm-four shell shown above.
No compact example results: for any fixed `d^2=3`,

\[
 z=(1,(43/55)d)
\]

is strictly hyperideal against every possible wall through `H<=1000`.

## Host-specific exact audits

The literal full Gram matrices were parsed from normalized snapshots of the
primary Nebe--Sloane catalogue pages ([T27](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/T27.html),
[U28](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/U28.html), and
[dim29odd](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim29odd.html)). For every row, PARI independently gives
determinant one and minimum three, reproduces the catalogue kissing number,
and counts the exact norm-four shell. Each displayed Gram contains an
explicit diagonal norm-three basis vector.

## Uniform H1000 theorem

The universal arithmetic census has **134 controls at
104 heights**. Exactly
30 controls retain a necessary Gram type,
for 256 types at
24 heights. Its tight row is

\[
 H=455 + 322*sqrt(2),\quad q=2,\quad
 (u^2,v^2,u.v)=(460,226,322),
\]

with maximum `37+26sqrt(2)` and strict margin
`-51/55 + 37/55*sqrt(2)`. The witness is spacelike by
`8572/3025 - 2*sqrt(2)>0`.

This finite theorem depends only on integrality, unimodularity, minimum three,
and a norm-three direction—not on rank or on an automorphism-orbit sample.
The necessary binary and augmented Gram gates retain supersets, so the
pairing bounds apply to every actual wall.

Thus all three complete bounded prefixes have open hyperideal patches and
infinite volume. Walls above `H=1000` remain open; this is not reflectivity or
Vinberg termination.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/certify_silver_odd_unimodular_h27_h29_h1000.py --verify
PYTHONPATH=src pytest -q tests/test_silver_odd_unimodular_h27_h29_h1000.py
```
