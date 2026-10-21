# Silver Leech H24 exact packet screen

Let `pi=sqrt(2)`, `epsilon=1+pi`, and let `Lambda` be the even unimodular
Leech lattice of minimum four.  This note records exact Vinberg-prefix
certificates for

\[
\langle-a\rangle\oplus\Lambda,
\qquad a\in\{\epsilon,\pi,2\pi,2\pi-1\}.
\]

It does **not** claim reflectivity, termination, compactness, or novelty.

Reproduce with

```console
python3 scripts/audit_silver_leech_h24.py
PYTHONPATH=src pytest -q tests/test_silver_leech.py
```

## Mirror norms and admissibility

All four coefficients are positive at `pi>0` and negative at the conjugate
embedding, so the form is hyperbolic at exactly one real place.  Local
reflection-integrality at the ramified prime `(pi)`, together with spatial
unimodularity and primitivity, gives the complete norm classes

| coefficient | canonical mirror norms |
|:--|:--|
| `epsilon` | `1, 2` |
| `pi` | `2, 2+pi` |
| `2pi` | `2, 4+2pi` |
| `2pi-1` | `1, 2, 3+pi, 6+2pi` |

For `a=2pi` and norm `4+2pi`, there is an additional indispensable condition
`y in pi Lambda`.  Applying it removes the tempting but nonreflective
`epsilon`-times-minimal-shell candidate at height `1+pi`.

## Unit coefficient: complete first packet

Write `y=u+pi v`.  If `y^2=A+B pi`, then

\[
u^2+2v^2=A,\qquad u\mathbin\cdot v=B/2.
\]

Evenness, minimum four, Cauchy--Schwarz, and the absence of norm-two Leech
vectors exclude every control before `h=14+10pi`.  At that height the only
surviving Gram triple is `(u^2,v^2,u.v)=(10,10,10)`, hence

\[
q=2,\quad x=2\epsilon,\quad y=\epsilon b,\quad b^2=10.
\]

The packet has `4,629,381,120` walls.  It need not be materialised.  For
distinct norm-ten vectors `b,c`, Leech minimum applied to `b-c` gives
`b.c<=8`, and therefore

\[
(r_b,r_c)=\epsilon^2(b\mathbin\cdot c-4\epsilon)<0.
\]

An explicit 24-tuple of norm-ten vectors in the integral ALCO MOG basis has
coordinate determinant two, so the affine section

\[
P=\{z: |z\mathbin\cdot b|\le2\epsilon\text{ for every }b^2=10\}
\]

is bounded.  Boundedness is **not** finite hyperbolic volume.  For any minimal
`u0^2=4`, minimum four gives `|u0.b|<=5`.  Thus the strict interior point
`z=4u0/5` lies in `P`, while

\[
z^2=64/25>\epsilon.
\]

It supplies an open hyperideal patch and decisively disproves finite volume
at the first packet.

## Unit coefficient through the complete H82 prefix

The next nonempty necessary shell is at `h=28+20pi`, with

\[
q=2,\quad x=4+2\pi,\quad
y=w+\epsilon v,\quad (w^2,v^2,w.v)=(4,14,6).
\]

The first boundary witness is cut, but the uniform replacement

\[
z=c u_0,\qquad c=(4+10\pi)/23
\]

is feasible for every such root and remains hyperideal:

\[
z^2-\epsilon=(335-209\pi)/529>0.
\]

For every arithmetically possible later root, an exact elementary oracle
chooses `w=u-kv` of minimal square and maximises `u0.y` over integral
pairings subject to the exact `3x3` Gram determinant and Leech minimum
constraints.  Through height `82+58pi`, the only further nonempty necessary
level is:

| height | control | necessary Gram types | rigorous max `u0.y` |
|:--|:--|--:|:--|
| `28+20pi` | `4+2pi` | 1 | `10+6pi` |
| `41+29pi` | `4+3pi` | 3 | `11+8pi` |
| `56+40pi` | `4+4pi` | 3 | `14+10pi` |
| `63+45pi` | `6+3pi` | 3 | `15+11pi` |
| `82+58pi` | `6+4pi` | 6 | `17+12pi` |

In every row `c max(u0.y) <= epsilon x`.  Empty necessary shells require no
wall-orbit or acuteness work.  Consequently the same explicit hyperideal ray
survives every possible wall, hence every accepted Vinberg wall, through the
complete `H=82` prefix.

## Nonunit coefficients

For both `a=pi` and `a=2pi`, the first reflective packet occurs at
`h=8+6pi`, has norm two, and is `epsilon` times the complete norm-six Leech
shell (`16,773,120` walls).

- For `a=pi`, its slabs are `|z.b|<=2pi`; `z=(2pi/3)u0` is hyperideal.
- For `a=2pi`, its slabs are `|z.b|<=4`; `z=(6/5)u0` is hyperideal.

Their next nonempty necessary packet is common at `h=24+17pi` and has

\[
y=w+\epsilon v,\qquad(w^2,v^2,w.v)=(4,12,5).
\]

Exact minimal-vector bounds `|u0.w|<=4`, `|u0.v|<=6` give explicit strict
hyperideal witnesses feasible through that entire packet for both seeds.
Thus none of the three silver Leech coefficients is a finite-volume survivor
at the frozen exhaustive prefix.

## Primitive norm-seven coefficient: `a=2pi-1`

The admissible coefficient `rho=2pi-1` has norm `-7`.  Local reflection
integrality gives exactly the norm ideals `O,(2),(rho),(2rho)`, represented
by the totally positive elements

\[
1,\quad 2,\quad 3+\pi,\quad 6+2\pi.
\]

For either class divisible by `(rho)`, the spatial component must lie in
`rho Lambda`; the exact control census enforces this before its Leech Gram
test.  All seven lower controls are empty.  The first packet is

\[
h=5+4\pi,\quad q=2,\quad x=2+\pi,\quad
y=\epsilon u,\quad u^2=4,
\]

so it consists of all `196,560` minimal Leech vectors.  Dividing each root by
`epsilon`, its affine section is

\[
P=\{z:|z\mathbin\cdot u|\le 4-\pi\ \text{for every }u^2=4\}.
\]

This is exactly a rescaling of the standard contact polar, whose reference
bound is two.  Dutour Sikiri\'c--Sch\"urmann--Vallentin classify its complete
vertex set into 232 `Co0` orbits.  The unique orbit of maximum squared norm
has reference square `8/3`, stabilizer `Co3`, and 552 incident minimal-vector
facets.  After rescaling, its square is

\[
\frac{8}{3}\left(\frac{4-\pi}{2}\right)^2
=12-\frac{16}{3}\pi.
\]

Table 2a gives the exceptional representative explicitly.  In integer MOG
coordinates scaled by `sqrt(8)`, its raw row is

```text
(2,2,2,2,-2,-6,2,-2,-2,6,2,2,-2,2,-2,2,2,2,2,2,-2,2,2,6)/3.
```

The raw coordinate square is `192`, so division by both `sqrt(8)` and three
gives `192/(8*9)=8/3`.  Equivalently, the exceptional direction is `2d/3`
for the norm-six Leech vector obtained by dividing the raw row by two.

Its Lorentzian norm is therefore

\[
13-\frac{22}{3}\pi>0,
\]

the last inequality following exactly from `39^2 > 2*22^2`.  Thus the full
first packet has a hyperideal vertex orbit and an open hyperideal patch: it is
not finite volume.  This conclusion uses the complete polar classification,
not boundedness of the affine section or a sampled direction.

The next nonempty necessary packet is at

\[
h=31+22\pi,\quad q=2,\quad x=4+3\pi,
\]

and its surviving Leech Gram types for `y=u+pi v` are
`(u^2,v^2,u.v)=(36,14,22),(32,16,22)`.  The first is genuinely realised:
the frozen scaled-`sqrt(8)` coordinates in the artifact have
`w^2=6,v^2=14,w.v=8` and `u=w+v`.  They also show that the simple first-packet
norm-six witness `z=3w/5` can be cut.  This is not a rescue of finite volume.
For any fixed norm-six Leech direction `d`, the replacement

\[
z=\frac{37}{67}d
\]

is strictly inside every first-packet slab and remains hyperideal, because

\[
(4-\pi)-\frac{111}{67}=\frac{157}{67}-\pi>0,
\qquad
z^2-(2\pi-1)=\frac{12703}{4489}-2\pi>0.
\]

An exact `3x3` positive-semidefinite Gram oracle, augmented by the Leech
minimum on sums and differences, gives the following rigorous pairing maxima
through the complete `H=82+58pi` prefix:

| height | norm | control | necessary Gram types | max `d.y` |
|:--|:--|:--|--:|:--|
| `31+22pi` | `2` | `4+3pi` | 2 | `14+8pi` |
| `56+40pi` | `6+2pi` | `12+8pi` | 1 | `21+15pi` |
| `62+44pi` | `2` | `6+4pi` | 4 | `17+13pi` |
| `82+58pi` | `6+2pi` | `14+10pi` | 2 | `21+15pi` |

For norm `1,2` the exact wall margin is
`(2pi-1)x-(37/67) max(d.y)`.  In a `(rho)`-divisible norm class, however,
the census uses `y=rho z`; cancelling `rho` in the wall inequality makes its
threshold `x`, not `(2pi-1)x`.  The four margins above are respectively
`(18+39pi)/67`, `(27-19pi)/67`, `(41+55pi)/67`, and
`(161+115pi)/67`, all strictly positive.  All other control levels in the
interval have empty necessary Gram censuses.  Hence a strict hyperideal open
patch persists through every possible (and therefore every Vinberg-accepted)
wall in this exhaustive prefix.

Primary source: M. Dutour Sikiri\'c, A. Sch\"urmann, F. Vallentin,
[*The contact polytope of the Leech lattice (complete version)*](https://arxiv.org/abs/0906.1427),
Discrete & Computational Geometry 44 (2010), 904--911.
