# Silver odd-Leech H24: exact H45 prefix for norm-23 coefficient

## Outcome

For

`L = <-(-3+4 sqrt(2))> + O_24`

over `Z[sqrt(2)]`, the first nonempty positive Vinberg packet occurs at

`H = 7/2 + 3 sqrt(2)`.

It consists of exactly **4096 primitive norm-two roots**
`r_b=(1+sqrt(2))(1,b)`, one for each norm-three vector `b` of the odd Leech
lattice.  All distinct roots in the packet have negative pairing, so the
whole packet is accepted.  The exact affine point

`(1,d)`, with `d=(4,4,0^22)/sqrt(8)` and `d^2=4`,

does **not** itself carry the H45 claim.  The exact replacement
`(1,z)` with `z=(5/6)d` satisfies every possible wall through
`H=45/2+16sqrt(2)` strictly and is
hyperideal.  Therefore this complete necessary-wall prefix has **infinite
hyperbolic volume**.  This is not a
global nonreflectivity result: later mirrors may cut the patch.

## Exact odd-Leech model

Use the [Nebe--Sloane catalogue basis](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/oddLeech.html) on the `sqrt(8)`
raw-coordinate scale.  With the catalogue extended Golay code `C` and
`Z_+^24={z: sum(z) even}`, the lattice is

`((2C+4Z_+^24) union (1+2C+4Z_+^24))/sqrt(8)`.

The encoded 24-row basis has integral Gram determinant **1**;
all 24 leading principal minors are positive.  Its Golay generators have
dimension 12 and weight enumerator
`1 + 759 x^8 + 2576 x^12 + 759 x^16 + x^24`.

The construction and basis are exactly equal, rather than merely similar
models.  Raw `4Z_+^24` has covolume `2^49`; adjoining the 12-dimensional code
reduces this to `2^37`, and the shifted coset halves it to `2^36`.  Dividing
24 coordinates by `sqrt(8)` gives covolume one.  Every catalogue basis row
lies in that construction, and its Gram determinant is also one, so the
containment and equal covolumes force equality.

The count 4096 is proved, not guessed from a sign-pattern shortcut.  In the
shifted coset all raw coordinates are odd, so norm three forces all 24 to be
`+/-1`.  Their residues modulo four are exactly one word of `C`, giving
precisely `|C|=4096` vectors.  The unshifted coset has minimum at least four:
a nonzero codeword has weight at least eight, while codeword zero uses a
nonzero even-sum integer vector.  Independently, PARI
`qfminim(G,3,,flag=2)` returned **4096** vectors from the determinant-one Gram
matrix (a `24 x 2048` stored half-shell).

An arbitrary collection of 4096 signs, a wrongly ordered Golay code, or an LP
that has not first proved compatibility with this basis is not a valid odd-
Leech shell model.  The certificate here avoids that failure mode entirely.

## Complete local mirror-norm candidates

Put pi=sqrt(2), a=-3+4*pi, and p=(a), of ideal norm 23. For a primitive r=(x,y), reflection integrality is (q(r))|2I(r), with I(r)=(a*x,(y,O24)). Away from (pi) and p, any prime dividing q would divide x and, by unimodularity, y, contradicting primitivity. At (pi), a is a unit and primitivity makes I(r) a unit, so v_pi(q)<=v_pi(2)=2; oddness supplies no even-lattice lower bound. Thus the dyadic exponents are among 0,1,2. At p, a positive q-valuation and q|2(y,O24) force y in p*O24. Then x is a p-unit, I(r) has valuation one, and q=-a*x^2 mod p^2 has valuation exactly one. Hence the p-exponent is 0 or 1, and the exponent-one classes have the spatial gate y in a*O24. Modulo totally-positive unit squares the six generators are 1,2+pi,2,5+pi,12+7*pi,10+2*pi.

At height zero the p23-supported cases contradict primitivity.  The remaining
canonical norms have traces `2,4,4`, whereas any nonzero
`y=u+sqrt(2)v` has `Tr(y^2)=2(u^2+2v^2)>=6`.  Hence there are no height-zero
walls.

## Exact control census through the first packet

The search uses the complete Minkowski box for `H<=8` and then applies the
exact endpoint `H<=7/2+3sqrt(2)`, crystallographic integrality, conjugate
positivity, and the p23 spatial gate.

| # | height | q | x | p23 gate | target after gate | integral? | necessary triples |
|---:|---|---|---|:---:|---|:---:|---:|
| 1 | `0 + (1/2) sqrt(2)` | `12 + (7) sqrt(2)` | `1 + (1) sqrt(2)` | yes | `61/23 + (43/23) sqrt(2)` | no | 0 |
| 2 | `1/2 + (1/2) sqrt(2)` | `10 + (2) sqrt(2)` | `1 + (1) sqrt(2)` | yes | `47/23 + (32/23) sqrt(2)` | no | 0 |
| 3 | `1 + (1) sqrt(2)` | `5 + (1) sqrt(2)` | `1 + (1) sqrt(2)` | yes | `36/23 + (25/23) sqrt(2)` | no | 0 |
| 4 | `1 + (1) sqrt(2)` | `10 + (2) sqrt(2)` | `2 + (1) sqrt(2)` | yes | `72/23 + (50/23) sqrt(2)` | no | 0 |
| 5 | `2 + (3/2) sqrt(2)` | `12 + (7) sqrt(2)` | `3 + (2) sqrt(2)` | yes | `183/23 + (129/23) sqrt(2)` | no | 0 |
| 6 | `2 + (2) sqrt(2)` | `5 + (1) sqrt(2)` | `2 + (1) sqrt(2)` | yes | `61/23 + (43/23) sqrt(2)` | no | 0 |
| 7 | `2 + (2) sqrt(2)` | `10 + (2) sqrt(2)` | `2 + (2) sqrt(2)` | yes | `122/23 + (86/23) sqrt(2)` | no | 0 |
| 8 | `7/2 + (5/2) sqrt(2)` | `10 + (2) sqrt(2)` | `3 + (2) sqrt(2)` | yes | `169/23 + (118/23) sqrt(2)` | no | 0 |
| 9 | `7/2 + (3) sqrt(2)` | `2` | `1 + (1) sqrt(2)` | no | `9 + (6) sqrt(2)` | yes | 1 |

The eight earlier controls are all p23-gated, and division by `a^2` makes
their reduced spatial targets nonintegral.  The endpoint has

`q=2`, `x=1+sqrt(2)`, `y^2=9+6sqrt(2)`.

Writing `y=u+sqrt(2)v` gives the unique necessary Gram triple
`(u^2,v^2,u.v)=(3,3,3)`.  Positive definiteness then gives `u=v=b`, with
`b^2=3`; conversely every one of the 4096 minimal vectors realizes a root.
Distinct Golay sign vectors have pairings `1,0,-1,-3`, producing root
pairings `-4-4sqrt(2), -7-6sqrt(2), -10-8sqrt(2), -16-12sqrt(2)`, all
strictly negative.

## Strict hyperideal witness

For every minimal `b`, the explicit norm-four vector `d` has
`b.d in {-1,0,1}` (counts `1024,2048,1024`).  The first walls in the affine
section are `b.z<=a`.  At `z=d`, the minimum wall margin is

`a-1 = -4+4sqrt(2) > 0`,

while the Lorentzian excess is

`d^2-a = 7-4sqrt(2) > 0`.

Thus a neighborhood of the segment's sphere-crossing point remains inside
the first-prefix chamber, proving an open hyperideal patch and infinite
volume for this certified prefix.

## Continuation through `H=45/2+16sqrt(2)`

The complete normalized arithmetic census through this endpoint contains
**22 controls**.  Apart from the first packet,
every possible level below the endpoint is empty after the p23 gate,
integrality, and all universal consequences of odd-Leech minimum three.  In
particular, the apparent intermediate target `(9,4,6)` is rank one and would
force a forbidden norm-one primitive direction.

At the endpoint,

`q=2`, `x=3+2sqrt(2)`, `y^2=47+32sqrt(2)`,

and exactly three necessary binary Gram types remain:

- `(27,10,16)`, equivalently `(b^2,w^2,b.w)=(3,10,-4)` for `b=u-2v`
  with `y=b+(2+sqrt(2))w`;
- `(25,11,16)`, equivalently `(d'^2,v^2,d'.v)=(4,11,5)` with
  `y=d'+(1+sqrt(2))v`;
- `(23,12,16)`, equivalently `(b^2,w^2,b.w)=(3,12,4)` with
  `y=b+(1+sqrt(2))w`.

There is no need to expand their enormous ordinary shells.  Extend each
binary Gram matrix by the fixed `d^2=4` row and column, require its exact
`3 x 3` Gram matrix to be positive semidefinite, and enumerate the finitely
many **joint** integral pairings.  Optimizing the full pair
`(first-vector.d, second-vector.d)`—not combining incompatible coordinatewise
bounds—gives for the unscaled `d`, respectively,

`y.d <= 10+6sqrt(2), 9+6sqrt(2), 9+6sqrt(2)`.

Consequently the unscaled `d` is explicitly **withdrawn as an H45 witness**:
all three necessary-superset maxima exceed the wall bound.  Scale along
the same ray and use the replacement `z=(5/6)d`.  Its maxima are

`y.z <= 25/3+5sqrt(2), 15/2+5sqrt(2), 15/2+5sqrt(2)`.

The wall bound is `a*x=7+6sqrt(2)`, leaving exact margins

`-4/3+sqrt(2), -1/2+sqrt(2), -1/2+sqrt(2)`.

All are strictly positive.  The tightest clearance in the entire prefix is
therefore `-4/3+sqrt(2)>0`.  The replacement remains hyperideal because

`z^2-a = 52/9-4sqrt(2)>0`.

Because the Gram scan deliberately optimizes over a
superset of realized odd-Leech pairs, the conclusion applies a fortiori to
every actual Vinberg wall through H45.

## Scope

Status: **H45_NECESSARY_PREFIX_INFINITE_VOLUME; EVENTUAL_CHAMBER_UNRESOLVED**.  The mirror-norm gates, height-zero exclusion,
control census, first shell, and the necessary-superset hyperideal witness are
exact and complete only through the stated H45 endpoint.  No claim is made
that later Vinberg walls preserve the patch or that the eventual reflection
chamber has infinite volume.
