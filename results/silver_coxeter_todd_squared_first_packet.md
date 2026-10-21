# Silver Coxeter--Todd-squared host: exact first packet

## Outcome

The genuinely new rootless host

\[
  L=\langle-(12\sqrt2-13)\rangle\oplus K_{12}\oplus K_{12}
  \quad\text{in }\mathbf H^{24}
\]

has a complete first nonempty Vinberg packet at
`H=67 + 48*sqrt(2)`.  It contains **613,872**
primitive walls.  An exact strict hyperideal point survives the whole packet,
so this prefix has infinite volume and is not a compact example.  Later walls
may cut the point; this is not a nonreflectivity theorem.

## Why this lane is new

The spatial lattice is the rootless even 3-modular lattice `K12^2`, not a
Niemeier, Leech, odd-Leech-section, or `Lambda22` host.  The coefficient

\[
 a=12\sqrt2-13=(\sqrt2-1)(3+\sqrt2)(5-2\sqrt2)
\]

has norm `-119=-7*17`.  It is unit-square-equivalent to the archived Leech
direction `9+10sqrt(2)`, but the host is genuinely new: `K12^2` introduces
a rank-12 ternary pairing radical and a different complete first packet.

The authoritative Coxeter--Todd Gram matrix is from the
[Nebe--Sloane lattice catalogue](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/K12.html).  Exact integral checks give
`det(K12)=729`, Smith form `1^6 3^6`, positive leading minors, minimum four,
and 756 minimal vectors.  PARI independently recovers the shell counts
`756, 4032, 20412` at norms `4,6,8`.  The mod-three pairing radical has
minimum 12 and 756 minimal vectors.

## Complete local and chronological search

The primitive pairing-ideal argument leaves exactly **16** norm classes:

`(1 or 2) * (1 or 3) * (1 or (3+sqrt(2))) * (1 or (5-2sqrt(2)))`.

All 16 height-zero classes are empty.  Through the first packet there are
exactly **50** controls at
**35** heights.  All
**48** earlier controls are
empty after the exact odd-prime divisibility, evenness, minimum, Cauchy, and
short-combination gates.

## First packet

The simultaneous components are:

| root norm | control | spatial shell | primitive walls |
|---:|:---|---:|---:|
| `1` | `3 + 2*sqrt(2)` | `K12^2(4)` | 1,512 |
| `2` | `4 + 3*sqrt(2)` | `K12^2(8)` | 612,360 |

For norm one, the only Gram triple `(36,16,24)` forces
`u=3d,v=2d,d^2=4`.  For norm two, `(64,36,48)` is exactly a
`sqrt(2)`-multiple of a norm-one root, while `(68,34,48)` would give the
forbidden vector `(2u-3v)^2=2`.  The remaining triple `(72,32,48)` forces
`u=3f,v=2f,f^2=8` and is primitive.  Thus

\[
  2\cdot756+\bigl(2\cdot20412+756^2\bigr)=613872.
\]

Minimum four and integrality give strict same-height acuteness bounds
`q1/q1=-33 - 24*sqrt(2)`,
`q1/q2=-11 - 7*sqrt(2)`,
and `q2/q2=-32 - 24*sqrt(2)`.

## Exact infinite-volume certificate

Let `h=G_K12^-1 e1`.  Then `h` lies in `K12^*`, `h^2=4/3`, and `h.d` is
integral for every `d` in `K12`.  Put

\[
 z=(1,(9/5)h,0).
\]

Its Lorentz norm is

\[
 z^2=433/25 - 12*sqrt(2)>0.
\]

Cauchy plus integrality gives `|h.d|<=2` on the norm-four shell and
`|h.f|<=3` on every first-factor projection of a norm-eight vector.  Hence
the exact maximum pairings with the two packet components are respectively
`9/5 - 14/5*sqrt(2)` and
`-19/5 + 9/5*sqrt(2)`, both strictly negative.

The basepoint is also strict.  The segment from it to the spacelike `z`
therefore crosses the light cone away from every facet, producing an open
ideal patch.  The complete first-packet chamber has infinite volume.

## Reproduction

```bash
PYTHONPATH=src python3 scripts/audit_silver_coxeter_todd_squared_first_packet.py --verify
PYTHONPATH=src pytest -q tests/test_silver_coxeter_todd_squared_first_packet.py
```
