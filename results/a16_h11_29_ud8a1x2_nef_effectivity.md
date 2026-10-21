# Nef framing and exact effectivity for the complex `A16 H11/29` carrier

## Positive theorem

The equation-level complex-K3 carrier of the finite-volume noncompact
`A16 H11/29` Coxeter chamber can be placed explicitly in Roulleau's complete
closed nef chamber for

```text
NS(X) = U + D8(-1) + A1(-1)^2.
```

One simultaneous marking transports all 29 oriented mirrors.  All 17/17
square-minus-two orientations are effective and have exact decompositions in
the complete stock of nineteen irreducible `(-2)`-curves.  The
irreducibility decision is exact: facet `23` is exactly `A16`, while the
other sixteen are reducible.  Thus one imported facet is a genuine wall of
the K3 nef cone.  The twelve square-minus-four, divisibility-two vectors are
kept separately, with no effectivity claim.

There is also an explicit common ample/Coxeter interior.  If `H` is the
nef-reduced carrier interior, then

```text
B1 = H + D136
   = (180,272,427,213,370,314,286,259,67,40,126,89)
```

is primitive and ample, has square `44258`, and is strictly positive on all
29 transported mirrors.  Its minimum imported-wall pairing is `30`, attained
at facets `5` and `11`.

This is a stronger algebraic realization of the known `H11/29` polytope.  It
is not a new abstract polytope, and the 29-wall chamber is not the full nef
chamber.

## Complete curve stock and the correct marking sign

The geometric input is the `U+D8+A1^2` subsection of Xavier Roulleau's
*An atlas of K3 surfaces with finite automorphism group*,
[EPIGA 6 (2022), article 19](https://doi.org/10.46298/epiga.2022.6286),
[arXiv:2003.08985](https://arxiv.org/abs/2003.08985).  It gives the complete
stock `A1,...,A19`; the first twelve curves form the displayed
Néron--Severi basis.  In that basis the published ample class is

```text
D136 = (2,5,9,4,10,12,15,19,9,9,6,0),   D136^2 = 136.
```

Its degrees on `A1,...,A19` are

```text
(1,1,1,1,1,1,1,1,1,1,7,2,2,12,18,18,12,18,18).
```

The old carrier uses a standard `U+D8+A1^2` basis, not Roulleau's curve
basis.  The certificate includes an explicit determinant-minus-one integral
change of basis and checks

```text
T^t Gram(Roulleau) T = Gram(U+D8+A1^2).
```

After this bridge, the archived square-`37418` strict carrier interior has
`D136`-degree `581903`.  Its negative has degree `-581903`.  Therefore the
correct positive-cone sign is `+`; no informal sign choice remains.

## The 87-reflection Weyl reduction

For a published curve `Ai` with `Ai^2=-2`, use

```text
s_i(v) = v + (v.Ai) Ai.
```

At each stage the deterministic rule fires the first curve having negative
pairing.  It produces the following 87-reflection word, acting from left to
right:

```text
1 4 3 5 3 6 5 3 7 8 7 11 8 7 6 5 3 4 10
8 7 6 5 11 8 7 12 1 13 1 12 15 18 13 1 2 18
13 1 2 3 5 6 7 8 11 12 1 2 3 5 6 7 8 13 1
2 3 5 6 7 12 1 2 3 5 6 13 1 2 3 5 12 1 2
3 4 13 1 2 3 5 6 12 1 2 13
```

Every firing has a strictly negative wall pairing and strictly lowers the
positive `D136`-height.  The product `W` is an integral determinant-minus-one
isometry of the Picard lattice.  Composed with the determinant-minus-one
standard-to-Roulleau bridge, it gives a determinant-one marking.

The endpoint is

```text
H = (178,267,418,209,360,302,271,240,58,31,120,89),
H^2 = 37418,                         H.D136 = 3352.
```

Its intersections with the complete curve stock are

```text
(0,62,0,0,0,27,0,0,124,178,0,0,178,240,294,62,418,116,240).
```

They are all nonnegative, so `H` is in the closed nef chamber; its zero walls
are `A1,A3,A4,A5,A7,A8,A11,A12`.  Simultaneously, its pairing with every
transported Coxeter normal is strictly positive, with minimum `27`.

## A primitive simultaneous ample interior

Consider

```text
B_k = kH + D136,                    k >= 0.
```

Because `D136` is ample and `H` is nef, every `B_k` is ample.  At `k=0`, the
pairings with imported facets `0,7,8` are respectively `0,0,-6`, so `D136`
alone is not strict for the imported chamber.  Already at `k=1`, however,

```text
B1 = H + D136
   = (180,272,427,213,370,314,286,259,67,40,126,89),
B1^2 = 44258,                       B1.D136 = 3488.
```

The coordinates are primitive.  Its degrees on all nineteen complete nef
walls are

```text
(1,63,1,1,1,28,1,1,125,179,7,2,180,252,312,80,430,134,258),
```

and its pairings with imported facets `0,...,28` are

```text
(178,182,64,31,65,30,132,54,118,370,65,30,134,199,190,
 69,67,195,136,39,184,202,197,80,390,68,203,376,380).
```

All entries in both lists are positive.  Hence `B1` is strictly inside both
the ample chamber and the transported 29-wall Coxeter chamber.  It is the
first nonnegative integral member of this family with the latter property.

## Exact effectivity and irreducibility decisions

Every row below is checked in Roulleau's `A1,...,A12` basis.  Each displayed
class has square `-2`, the coefficients are nonnegative integers, and its
support is connected.  The carrier pairing is the preserved strict pairing
with `H`.

| facet | exact effective divisor | `D136` degree | carrier pairing | decision |
|---:|---|---:|---:|---|
| `2` | `A1 + A2` | `2` | `62` | reducible |
| `3` | `A3 + A4 + A5 + A6` | `4` | `27` | reducible |
| `4` | `A2 + A3 + A5` | `3` | `62` | reducible |
| `5` | `A6 + A7 + A8` | `3` | `27` | reducible |
| `10` | `A2 + A3 + A4` | `3` | `62` | reducible |
| `11` | `A5 + A6 + A7` | `3` | `27` | reducible |
| `13` | `A1 + 2A2 + 2A3 + A4 + 2A5 + 2A6 + A7 + A8 + A11 + A12` | `21` | `178` | reducible |
| `14` | `A1 + 2A2 + 3A3 + A4 + 2A5 + 2A6 + A7` | `12` | `178` | reducible |
| `15` | `A1 + A2 + A3 + A4 + A5 + A12` | `7` | `62` | reducible |
| `16` | `A1 + A2 + A3 + A12` | `5` | `62` | reducible |
| `17` | `2A1 + 2A2 + 3A3 + A4 + 3A5 + 2A6 + A7 + A8 + A12` | `17` | `178` | reducible |
| `19` | `A3 + A5 + A6 + A7 + A8 + A11` | `12` | `27` | reducible |
| `21` | `A1 + 2A2 + 3A3 + 2A4 + 3A5 + 2A6 + 2A7 + 2A8 + A11` | `24` | `178` | reducible |
| `22` | `A1 + 2A2 + 4A3 + 2A4 + 3A5 + 2A6 + 2A7 + A8 + A12` | `19` | `178` | reducible |
| `23` | `A16` | `18` | `62` | **irreducible** |
| `25` | `A1 + A2 + 2A3 + A4 + A5` | `6` | `62` | reducible |
| `26` | `2A1 + 2A2 + 3A3 + 2A4 + 2A5 + 2A6 + 2A7 + A8 + A11 + A12` | `25` | `178` | reducible |

The replay finds each decomposition by negative-wall descent.  If a current
root `q` has `q.Ai=-m<0`, it replaces `q` by `q-mAi=s_i(q)` and records the
summand `mAi`.  The positive `D136`-degree strictly decreases at every step,
and the finite trace ends at one of the nineteen published curves.

Roulleau's stock is complete.  Consequently, equality with a stock vector is
an exact irreducibility test here: facet `23` is `A16`, while every other row
has an explicit nontrivial effective decomposition.  Component
multiplicities range from `1` to `18`, and the largest coefficient is `4`.

## The twelve square-minus-four vectors

These vectors all have square `-4`, divisibility `2`, and integral
Gram-preserving reflections.  They do not fall under the square-minus-two
effectivity argument.

| facet | transported coordinates | `D136` degree | carrier pairing |
|---:|---|---:|---:|
| `0` | `(-2,-2,-2,-1,-1,0,1,2,1,1,1,-2)` | `0` | `178` |
| `1` | `(-2,-2,-2,-1,-1,0,1,2,1,1,1,0)` | `4` | `178` |
| `6` | `(0,0,0,0,0,0,0,0,1,0,1,0)` | `8` | `124` |
| `7` | `(0,0,0,0,0,0,0,0,-1,1,0,0)` | `0` | `54` |
| `8` | `(0,0,0,0,0,0,0,0,1,0,-1,0)` | `-6` | `124` |
| `9` | `(0,0,0,-1,1,2,1,2,1,1,1,0)` | `14` | `356` |
| `12` | `(0,0,0,0,0,0,0,2,1,0,1,0)` | `10` | `124` |
| `18` | `(0,0,0,0,0,0,2,2,1,0,1,0)` | `12` | `124` |
| `20` | `(0,-2,-2,-1,-1,0,1,2,1,1,1,0)` | `6` | `178` |
| `24` | `(0,0,0,1,1,2,3,4,1,1,3,0)` | `34` | `356` |
| `27` | `(0,0,2,1,1,2,3,2,1,1,1,0)` | `20` | `356` |
| `28` | `(0,0,2,1,3,2,3,4,1,1,1,0)` | `24` | `356` |

The negative or zero `D136` degrees at facets `0,7,8` are a useful warning:
these are algebraic reflection vectors and strict imported halfspaces, not
curve-effectivity statements.

## Boundary of the theorem

This closes the earlier generic instruction to “Weyl re-mark and apply
Riemann--Roch.”  The basis bridge, sign, 87-reflection word, final nef vector,
all 29 transported reflection matrices, all 17 effective decompositions, the
one-versus-sixteen irreducibility split, and the common ample interior are
now exact and replayable.

It does not overturn the prior obstruction.  Twelve primitive imported
facets have square `-4`, and sixteen square-`-2` orientations are reducible
effective divisors rather than irreducible wall curves.  Therefore this is a
nef-framed mixed divisor-reflection carrier, not the full nef chamber and not
an automorphism-fixed chamber.

## Reproduction

- [certificate generator](../scripts/certify_a16_h11_29_ud8a1x2_nef_effectivity.py)
- [exact JSON certificate](../data/enumerations/a16_h11_29_ud8a1x2_nef_effectivity.json)
- [focused regression test](../tests/test_a16_h11_29_ud8a1x2_nef_effectivity.py)
- [underlying complex-K3 carrier](alexeev_h11_29_ud8_divisor_overlattice_carrier.md)
- [explicit BH `0.4.4.2` family](ud8a1x2_order4_weierstrass_family.md)

Replay from the repository root:

```bash
python3 -m scripts.certify_a16_h11_29_ud8a1x2_nef_effectivity --check-archive
python3 -m pytest -q tests/test_a16_h11_29_ud8a1x2_nef_effectivity.py
```
