# Cubic-seven rootless-Leech `H24`: two packets and an exact `H=1000` obstruction

## Outcome

Let

\[
K=\mathbf Q(\theta),\qquad
\theta=2\cos(2\pi/7),\qquad
\theta^3+\theta^2-2\theta-1=0.
\]

For the actual Leech-spatial host

\[
L=\langle-2\theta\rangle\oplus\Lambda_{24}
\]

the former audit stopped at `H=10`, before the first wall.  The exact
continuation now reaches two completely realised and accepted simultaneous
Vinberg packets:

| packet | exact height | exact description | walls |
|---:|:--|:--|---:|
| 1 | `5+14theta+6theta^2` | the complete minimal Leech shell | 196,560 |
| 2 | `16+45theta+20theta^2` | ordered minimal pairs of inner product `-1` | 9,258,762,240 |

This is a genuine packet chronology, not a necessary-Gram census.  All
first/second and same-height pairings are nonpositive, so both entire packets
are inserted without choosing orbit representatives.

No finite-volume or compact example results.  If `d` is the explicit
norm-six Leech vector recorded below, then

\[
z=(1,2d/3),\qquad z^2=\frac83-2\theta>0,
\]

is strict against **every arithmetically possible wall through `H=1000`**.
The complete box has 18 controls.  Nine are empty after the exact
minimum-four Gram test; the other nine contain 47 necessary Gram types.  The
smallest exact wall margin is

\[
-\frac43+\frac23\theta+\frac43\theta^2>0.
\]

Thus an open hyperideal patch survives the complete `H=1000` prefix.  Later
walls can still cut it, so this is not a nonreflectivity or nonexistence
theorem.

Reproduce the frozen record with

```console
PYTHONPATH=src python3 scripts/audit_cubic7_leech_h24_h1000.py --verify
PYTHONPATH=src pytest -q tests/test_cubic7_leech_h24_continuation.py
```

## Four small coefficient directions

Put `delta=2-theta`, the totally positive generator of the ramified prime
over seven.  The four small coefficient directions supported at `(2)` and
`(delta)` are now treated under one exact local convention:

| timelike coefficient `a` | ideal support | complete canonical mirror norms |
|:--|:--|:--|
| `theta` | unit | `1, 2` |
| `2theta` | `(2)` | `2` |
| `theta*delta` | `(delta)` | `1, 2, delta, 2delta` |
| `2theta*delta` | `(2)(delta)` | `2, 2delta` |

The apparent unit coefficient `theta-1` is not another case:

```text
(theta^2-2)^2 theta = theta-1.
```

For a primitive root `r=(x,y)` of norm `q`, reflection integrality is

\[
(q)\mid2I(r),\qquad I(r)=(ax,(y,\Lambda_{24})).
\]

Even unimodularity excludes primes away from `2a`.  At the odd prime
`(delta)`, a positive norm valuation is exactly one and forces
`y in delta Lambda_24`.  At inert `(2)`, an odd coefficient permits
valuation zero or one, whereas a coefficient with `v_2(a)=1` forces
valuation one.  Class number one and `U^+=U^2` remove the remaining
totally-positive unit-square ambiguity.  These local alternatives give the
displayed lists, including the spatial gate

\[
y=\delta y_0,\qquad y_0^2=y^2/\delta^2
\]

for norms `delta` and `2delta`.

Every one of the four hosts is rootless at height zero.  The trace lattice
of `O_K tensor Lambda_24` has minimum 12.  The ungated targets `1`, `2`, and
`delta` have traces `3`, `6`, and `7`, while the only potentially larger
target `2delta` becomes nonintegral after division by `delta^2`.  The full
minimum-four coefficient-Gram enumerator consequently returns no type in
any height-zero norm class.

## First accepted packet for `a=2theta`

Write

\[
s=2\theta+\theta^2,qquad
x_1=1+2\theta+\theta^2.
\]

Every earlier control is empty.  At

\[
H_1=5+14\theta+6\theta^2
\]

the coefficient Gram equations have exactly the rank-one solution

```text
(u^2,v^2,w^2,u.v,u.w,v.w)=(0,16,4,0,0,8).
```

It forces `u=0`, `v=2M`, `w=M` with `M^2=4`.  Conversely every minimal
Leech vector gives

\[
r_M=(x_1,sM),\qquad
-2\theta x_1^2+4s^2=2.
\]

Both `s` and `x_1` are units, so these are primitive reflective roots.
There are exactly 196,560.  For distinct minimal vectors `M,N`, Leech
minimum gives `M.N<=2`; the exact maximum root pairing is

\[
-4-14\theta-6\theta^2<0.
\]

Hence the complete equal-height packet is accepted.

## The 9.25-billion-wall second packet

The next nonempty control is

\[
H_2=16+45\theta+20\theta^2,qquad
x_2=1+4\theta+2\theta^2.
\]

Its unique necessary coefficient Gram type is

```text
(u^2,v^2,w^2,u.v,u.w,v.w)=(4,40,10,10,5,20).
```

It forces `v=2w`.  Put `D=-u` and `E=D+w`.  Then

\[
D^2=E^2=4,qquad D\mathbin\cdot E=-1,
\]

and every root is uniquely

\[
r_{D,E}=\bigl(x_2,s(E-D)-D\bigr).
\]

Conversely every ordered minimal pair with inner product `-1` realises the
Gram type.  The complete minimal-shell distribution about one vector is

```text
inner product:  -4    -2     -1      0      1     2    4
count:           1   4600  47104  93150  47104  4600   1
```

The minimal shell is one `Co0` orbit, so the packet count is exactly

\[
196560\cdot47104=9,258,762,240.
\]

This very large packet is nevertheless accepted exactly:

- its maximum pairing against packet one is
  `-s^2=-3-7theta-3theta^2<0`;
- two ordered pairs are controlled by four cross inner products in
  `[-4,4]`.  Exhausting all `9^4` integral possibilities, all principal
  minors, and the exact minimum-four quotient leaves 494 joint Gram types;
  after excluding the identical root, their maximum root pairing is `-1`.

Thus neither orbit sampling nor a same-height choice is hidden in the packet
claim.

## Complete later-prefix obstruction

For a spatial field vector

\[
y=u+\theta v+\theta^2w,
\]

the exact coefficient equations for `y^2=A+Btheta+Ctheta^2` are

\[
A=U+2R-W,
\quad B=2P+4R-W,
\quad C=2Q+V-2R+3W,
\]

where `(U,V,W,P,Q,R)` is the Gram data of `(u,v,w)`.  The inverse trace
Gram bounds `U<=Tr(y^2)` and `V,W<=2Tr(y^2)/7`.  The enumerator retains
rank-one and rank-two cases and quotients every integral kernel before
testing minimum at least four.

Take the explicit scaled-`sqrt(8)` vector

```text
d=(-4,-2,2,0,0,0,2,2,2,0,2,0,2,0,0,0,0,0,0,0,0,0,0,2).
```

It is the sum of a type-`(4^2)` minimal vector and an even-sign Golay-octad
minimal vector with pairing `-1`, hence `d^2=6`.

For every necessary Gram type, append the three integral pairings
`(d.u,d.v,d.w)`.  Cauchy bounds them finitely.  All principal minors of the
resulting four-by-four Gram are tested, followed by the exact minimum of its
integral kernel quotient.  Maximising `d.y` over this complete necessary
superset gives:

| exact height | Gram types | max `d.y` | margin for `z=(1,2d/3)` |
|:--|---:|:--|:--|
| `5+14theta+6theta^2` | 1 | `6theta+3theta^2` | `2+2theta` |
| `16+45theta+20theta^2` | 1 | `3+12theta+6theta^2` | `2+2theta` |
| `25+70theta+31theta^2` | 1 | `3+18theta+9theta^2` | `2` |
| `33+94theta+42theta^2` | 2 | `5+16theta+8theta^2` | `(8+10theta+2theta^2)/3` |
| `52+146theta+65theta^2` | 4 | `8+25theta+11theta^2` | `(2+4theta+2theta^2)/3` |
| `64+180theta+80theta^2` | 1 | `6+24theta+12theta^2` | `4+4theta` |
| `81+227theta+101theta^2` | 8 | `12+32theta+14theta^2` | `(2theta+2theta^2)/3` |
| `126+353theta+157theta^2` | 13 | `16+40theta+17theta^2` | `(-2+4theta+2theta^2)/3` |
| `149+418theta+186theta^2` | 16 | `17+44theta+19theta^2` | `(-4+2theta+4theta^2)/3` |

The other nine controls are empty.  Every displayed margin is positive at
the distinguished embedding.  This checks all possible roots, so it is
strictly stronger than checking only the two accepted packets or a proposed
set of Conway orbits.

## Secondary small-coefficient screens

The same exact gate/Gram/fixed-direction calculation gives:

| coefficient | cutoff | controls | nonempty necessary levels / types | strict witness |
|:--|---:|---:|:--|:--|
| `theta` | 300 | 37 | `1 / 1` | `(1,(7/15)d)` |
| `theta*delta` | 300 | 51 | `0 / 0` | no wall exists in the prefix |
| `2theta*delta` | 300 | 13 | `2 / 2` | `(1,(3/5)d)` |

For `a=theta`, the witness excess and minimum margin are

\[
\frac{98}{75}-\theta>0,
\qquad
\frac65+\frac23\theta-\frac23\theta^2>0.
\]

For `a=2theta delta`, they are

\[
\frac{54}{25}-4\theta+2\theta^2>0,
\qquad
1+\frac45\theta-\frac35\theta^2>0.
\]

The zero-wall `theta*delta` prefix is sparse, not a chamber closure.

## Scope and files

This audit proves exact local mirror norms and ramified spatial gates for the
four displayed coefficients; exact height-zero emptiness; the complete two
accepted packets for `a=2theta`; all necessary Leech coefficient Gram types
through `H=1000`; and a strict hyperideal obstruction against every one of
those possible walls.  It does not enumerate walls above the cutoff, classify
all coefficient ideals, or prove global nonreflectivity.

- Backend: `src/coxeter_polytopes/search/cubic7_leech_continuation.py`
- Generator: `scripts/audit_cubic7_leech_h24_h1000.py`
- Tests: `tests/test_cubic7_leech_h24_continuation.py`
- Archive: `data/enumerations/cubic7_leech_h24_h1000_continuation.json`
