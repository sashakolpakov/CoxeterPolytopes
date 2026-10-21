# Silver `N(A2^12)`: accepted H24 cutter and surviving hyperideal point

## Exact outcome

For both

\[
 \langle-\sqrt2\rangle\oplus N(A_2^{12})
 \quad\text{and}\quad
 \langle-2\sqrt2\rangle\oplus N(A_2^{12}),
\]

the complete first-cap candidate census through

\[
 H=24+17\sqrt2
\]

has only two nonempty post-affine levels:

| height | norm | complete first-cap candidates |
|---|---:|---:|
| \(8+6\sqrt2\) | 2 | 3,960 |
| \(24+17\sqrt2\) | 2 | 4,733,256 |

The first row is the complete accepted second packet.  An exact exhaustive
pairing audit finds 7,839,114 negative and 66 zero unordered pairings among
its distinct walls, with no positive pairing.

At the endpoint, the vector with twelve `A2` Dynkin-weight blocks

\[
 (2+2\sqrt2,2+2\sqrt2),\quad(1,4+3\sqrt2),
 \quad(0,0)^{10}
\]

is a primitive reflection-integral norm-two wall.  Its time coordinate is
\(4+3\sqrt2\) for the first coefficient and \(3+2\sqrt2\) for the second.
It satisfies every affine cap and its maximum pairing against all 3,960
second-packet walls is exactly \(-1\).  All intervening controls are empty,
so this is a chronologically accepted third-packet wall.  Its evident orbit
has \(12\cdot11\cdot2=264\) walls; all distinct orbit pairs are strictly
acute, with maximum \(-4-2\sqrt2\).

## Complete-superset hyperideal certificate

It is unnecessary to enumerate which of the 4,733,256 endpoint candidates
survive all 3,960 earlier inequalities.  A strict point survives the entire
candidate superset and therefore every possible accepted subset.

For `a=sqrt(2)`, place the rational Dynkin weight \(c=1/4\) in both
coordinates of all twelve blocks.  Its Lorentz norm is

\[
 24c^2-\sqrt2=\frac32-\sqrt2>0.
\]

The exact largest total Dynkin-weight sums are \(9+9\sqrt2\) at the second
packet and \(21+15\sqrt2\) at the endpoint.  Hence the maximum pairings are

\[
 -\frac74+\frac14\sqrt2<0,
 \qquad
 -\frac34-\frac14\sqrt2<0.
\]

For `a=2sqrt(2)`, use \(c=3/8\).  Its norm is

\[
 \frac{27}{8}-2\sqrt2>0,
\]

and its maximum pairings at the two levels are

\[
 -\frac58-\frac58\sqrt2<0,
 \qquad
 -\frac18-\frac38\sqrt2<0.
\]

Both points are also strict against the affine and height-zero walls.  Thus
both complete bounded prefixes have an open hyperideal patch and infinite
volume.  Later walls may cut these points: this is not a global
nonreflectivity theorem or a compact example.

## Reproduction

```sh
PYTHONPATH=src python3 scripts/audit_silver_a2_ternary_golay_h24_cutter.py
PYTHONPATH=src pytest -q tests/test_silver_a2_ternary_golay.py
```

The exact archive is
`data/enumerations/silver_a2_ternary_golay_h24_cutter.json`.
