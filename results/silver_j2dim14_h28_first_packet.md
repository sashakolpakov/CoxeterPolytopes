# Silver J2dim14 host in H28: exact first packet

## Outcome

The new rank-28 host

\[
L=\langle-(1+\sqrt2)\rangle\oplus J2dim14
\]

has no height-zero walls. Its exact first nonempty Vinberg packet is at
`H=14 + 10*sqrt(2)` and contains **374,976**
primitive walls. A strict hyperideal point survives the entire packet, so
this prefix has infinite volume and is not a compact example. This is a
bounded-prefix statement; later walls may cut the point.

## Exact lattice

The literal Nebe--Sloane catalogue Gram has determinant `6103515625`
and Smith invariants `1^14 5^14`. PARI recovers minimum eight and exactly
23,310 minimal vectors. The full preimage of the
mod-five radical has index `5^14`; its Gram divided by five is integral and
isometric to the original lattice, independently confirming the 5-modular
local gate.

PARI gives cumulative signed shell counts `23,310` through norm eight and
`398,286` through norm ten. Thus the norm-ten shell has exactly
**374,976** vectors.

## Eight coefficient lanes through H=300

All local alphabets and all height-zero exclusions are complete. Nonempty
counts here are necessary Gram packets; only the unit first packet is also
proved exactly realized.

| coefficient | norm classes | controls | necessary nonempty | first height |
|:---|---:|---:|---:|:---|
| `1 + 1*sqrt(2)` | 4 | 41 | 4 | `14 + 10*sqrt(2)` |
| `1*sqrt(2)` | 4 | 33 | 5 | `96 + 68*sqrt(2)` |
| `2*sqrt(2)` | 4 | 23 | 5 | `96 + 68*sqrt(2)` |
| `-1 + 2*sqrt(2)` | 8 | 55 | 1 | `124 + 88*sqrt(2)` |
| `5 + 5*sqrt(2)` | 4 | 41 | 4 | `14 + 10*sqrt(2)` |
| `5*sqrt(2)` | 4 | 33 | 5 | `96 + 68*sqrt(2)` |
| `10*sqrt(2)` | 4 | 23 | 5 | `96 + 68*sqrt(2)` |
| `-5 + 10*sqrt(2)` | 8 | 55 | 1 | `124 + 88*sqrt(2)` |

The four coefficients divisible by five reproduce their unscaled normalized
chronologies in the mirror-norm class carrying the inert-five factor.

## Exact first packet

There are 12 controls at
8 heights through the endpoint.
All 10 earlier controls have empty necessary
minimum-eight Gram censuses. The sole surviving row is

\[
q=2,\quad x=2+2\sqrt2,\quad
y=(1+\sqrt2)d,\qquad d^2=10.
\]

The binary Gram `(10,10,10)` forces the scalar form, so the packet is the
complete norm-ten shell and has 374,976 walls.

## Strict infinite-volume certificate

Let `h=G^(-1)e_13`, using zero-based index 13. Then `h^2=8/5` and
`h.d=d_13`. Direct enumeration of the complete norm-ten shell gives the
coordinate-absolute maxima

`5,5,3,14,13,7,6,4,9,10,10,7,4,2,7,5,4,2,4,6,4,4,2,2,2,2,2,2`,

hence `|h.d|<=2`. Put `z=(1,(3/2)h)`. Its Lorentz excess is

\[
\frac94\frac85-(1+\sqrt2)=\frac{13}5-\sqrt2>0.
\]

After dividing a packet root by `1+sqrt(2)`, the slab bound is
`2+2sqrt(2)`. The witness pairing is at most three, leaving exact margin

\[
-1+2\sqrt2>0.
\]

Thus an open hyperideal patch survives every first-packet wall.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_j2dim14_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_j2dim14_first_packet.py
```
