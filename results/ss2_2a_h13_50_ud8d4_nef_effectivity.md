# Exact nef-chamber effectivity for the `SS2-2A H13/50` K3 carrier

## Positive theorem

The complex-K3 realization of the finite-volume noncompact `H13/50`
Coxeter chamber can be placed explicitly in the published closed nef chamber
of

```text
NS(X) = U + D8(-1) + D4(-1).
```

Under that single marking, all 50 mirrors remain strict at the transported
carrier interior.  More strongly, all 40/40 square-minus-two orientations
are effective, with exact nonnegative decompositions into Roulleau's complete
stock of 20 irreducible `(-2)`-curves.  Exactly 0/40 is itself one of those
irreducible curves: all forty are reducible.  The ten square-minus-four,
divisibility-two reflection vectors are listed separately, with no
effectivity claim.

There is also an explicit primitive class strictly inside both relevant
chambers.  If `H=v_nef`, then

```text
B2 = 2H + D154
   = (53,107,162,250,169,339,260,182,137,93,-18,14,38,16)
```

is ample, has square `12970`, and has positive pairing with all 50 imported
Coxeter normals, with minimum `4`.

The curve and ample-class input is the `U+D8+D4` subsection of Xavier
Roulleau's *An atlas of K3 surfaces with finite automorphism group*,
[EPIGA 6 (2022), article 19](https://doi.org/10.46298/epiga.2022.6286),
[arXiv:2003.08985](https://arxiv.org/abs/2003.08985).  Its curves
`A1,...,A20` form the complete irreducible `(-2)`-curve configuration;
`A1,...,A14` are the displayed Néron--Severi basis.  In that basis the
published ample class is

```text
D154 = (1,3,6,10,7,15,14,14,15,17,8,8,0,4),   D154^2 = 154.
```

Its degrees on `A1,...,A20` are

```text
(1,1,1,1,1,1,1,1,1,1,1,1,17,9,17,17,17,1,9,1).
```

## Signed marking and explicit Weyl transport

There is one essential component correction.  The old isometry `X` sends the
certified square-`2012` carrier interior to a vector of `D154`-degree
`-10008`, hence to the positive-cone component opposite the published ample
class.  Reflections cannot change that component.  We first choose the
equally valid signed marking `-X`; this negates both the interior and all 50
oriented normals, preserving every chamber pairing.  Its `D154`-degree is
then `10008`.

For `Ai^2=-2`, use

```text
s_i(v) = v + (v.Ai) Ai.
```

The following 60-letter Weyl word acts from left to right:

```text
1 2 1 3 2 4 3 2 1 5 6 4 3 2 7 6 4 3 2 1
8 7 6 4 3 2 1 9 8 7 6 4 3 2 11 13 10 12 14
10 9 12 13 10 14 16 1 20 1 2 3 4 6 7 8 16 1 2 3 20
```

At every letter the current pairing with the fired curve is strictly
negative, and the positive integer height against `D154` strictly decreases.
The exact replay multiplies the 60 integral reflection matrices to an
explicit determinant-one isometry `W`, checks `W^T Gram W = Gram`, and also
archives the combined isometry `W(-X)` from the doubled-dual carrier lattice.

The endpoint is

```text
v_nef = (26,52,78,120,81,162,123,84,61,38,-13,3,19,6),
v_nef^2 = 2012,                 v_nef.D154 = 1192.
```

Its pairings with the complete curve stock are

```text
(0,0,16,0,0,0,0,16,0,0,64,32,0,26,38,0,32,64,38,26).
```

They are all nonnegative, so this is in the closed nef chamber.  It lies on
the ten walls `A1,A2,A4,A5,A6,A7,A9,A10,A13,A16`; it is not claimed to be
ample.  Its pairing with every transported Coxeter normal is nevertheless
strictly positive, with minimum `6`.

## A simultaneous ample and Coxeter-strict interior

The nef endpoint can be moved off its ten genuine curve walls without
leaving the imported Coxeter chamber.  Consider the integral family

```text
B_k = kH + D154,                    k >= 0.
```

Every `B_k` is ample: `D154` is ample, `H` is nef, and the archived degrees
against the complete 20-curve stock are strictly positive.  The first two
translates give an exact boundary calculation:

```text
B_1 = H + D154,      B_1^2 = 4550,
```

is positive on 49 imported walls but fails facet `22` by exactly `-2`.  By
contrast,

```text
B_2 = 2H + D154
    = (53,107,162,250,169,339,260,182,137,93,-18,14,38,16),
B_2^2 = 12970,                 B_2.D154 = 2538.
```

Its degrees on `A1,...,A20` are

```text
(1,1,33,1,1,1,1,33,1,1,129,65,17,61,93,17,81,129,85,53),
```

so it is strictly inside the complete ample chamber.  Its pairings with the
transported facets `0,...,49` are

```text
(36,36,34,74,38,35,78,38,36,38,
 52,39,80,35,44,54,72,36,56,70,
 58,28,4,82,41,38,143,161,168,181,
 182,144,163,165,149,169,162,167,171,150,
 166,187,151,183,165,183,151,168,188,167).
```

All 50 are positive; the minimum is `4`, at facet `22`.  The class is
primitive.  Moreover, `D154` itself fails imported facets `8,14,22`, and
`B_1` fails only facet `22`, while `B_2` succeeds.  Thus `k=2` is the first
nonnegative integral member of this family that is simultaneously ample and
strictly inside the 50-wall Coxeter chamber.

## The 40 exact effective decompositions

The facet number is the zero-based index in the certified `H13/50` carrier.
Every displayed equality is checked in the `A1,...,A14` Néron--Severi
coordinate basis.  Coefficients are nonnegative integers; every support is
connected; and every class has square `-2` and divisibility `1`.

| facet | exact effective divisor | `D154` degree | carrier pairing |
|---:|---|---:|---:|
| `0` | `A2 + A3 + A4 + A6` | `4` | `16` |
| `1` | `A4 + A6 + A7 + A8` | `4` | `16` |
| `2` | `A3 + A4` | `2` | `16` |
| `4` | `A3 + A4 + A5 + 2A6 + A7` | `6` | `16` |
| `5` | `A8 + A9 + A10` | `3` | `16` |
| `7` | `A1 + A2 + A3 + A4 + A6 + A7` | `6` | `16` |
| `9` | `A1 + A2 + A3 + A4 + A5 + A6` | `6` | `16` |
| `10` | `A1 + A2 + A3 + A16` | `20` | `16` |
| `11` | `A4 + A5 + 2A6 + A7 + A8 + A9` | `7` | `16` |
| `13` | `A7 + A8 + A9` | `3` | `16` |
| `15` | `A6 + A7 + A8 + A9 + A10 + A13` | `22` | `16` |
| `17` | `A5 + A6 + A7 + A8` | `4` | `16` |
| `18` | `A4 + A5 + A6 + A7 + A8 + A9 + A10 + A13` | `24` | `16` |
| `20` | `A1 + A2 + A3 + 2A4 + A5 + 2A6 + A7 + A16` | `26` | `16` |
| `24` | `A4 + A5 + 2A6 + 2A7 + A8 + A9 + A10` | `9` | `16` |
| `25` | `A2 + A3 + A4 + A5 + A6 + A7` | `6` | `16` |
| `26` | `A1 + A2 + 2A3 + 2A4 + A5 + 3A6 + 2A7 + 2A8 + A9` | `15` | `64` |
| `27` | `A1 + 2A2 + 2A3 + 2A4 + A5 + 2A6 + 2A7 + 2A8 + A9 + A10 + A13` | `33` | `64` |
| `28` | `2A1 + 2A2 + 2A3 + 3A4 + 2A5 + 5A6 + 3A7 + 2A8 + A9 + A10 + A16` | `40` | `64` |
| `29` | `A1 + A2 + 2A3 + 3A4 + 2A5 + 4A6 + 2A7 + 2A8 + A9 + A10 + A13 + A16` | `53` | `64` |
| `30` | `2A1 + 2A2 + 2A3 + 3A4 + A5 + 3A6 + 2A7 + 2A8 + 2A9 + A10 + A13 + A16` | `54` | `64` |
| `31` | `A2 + 2A3 + 2A4 + A5 + 3A6 + 3A7 + 2A8 + A9 + A10` | `16` | `64` |
| `32` | `A1 + 2A2 + 2A3 + 3A4 + 2A5 + 3A6 + 2A7 + 2A8 + A9 + A16` | `35` | `64` |
| `33` | `A1 + A2 + 2A3 + 3A4 + A5 + 4A6 + 3A7 + 2A8 + 2A9 + A10 + A16` | `37` | `64` |
| `34` | `A1 + 2A2 + 2A3 + 2A4 + 2A5 + 4A6 + 3A7 + 2A8 + 2A9 + A10` | `21` | `64` |
| `35` | `A1 + A2 + 2A3 + 4A4 + 2A5 + 5A6 + 3A7 + 2A8 + 2A9 + 2A10 + A13` | `41` | `64` |
| `36` | `A2 + 2A3 + 2A4 + 2A5 + 3A6 + 2A7 + 2A8 + 2A9 + A10 + A13` | `34` | `64` |
| `37` | `A1 + 2A2 + 2A3 + 3A4 + 2A5 + 5A6 + 4A7 + 2A8 + A9 + A16` | `39` | `64` |
| `38` | `A1 + 2A2 + 2A3 + 4A4 + 3A5 + 6A6 + 4A7 + 2A8 + A9 + A10 + A13` | `43` | `64` |
| `39` | `A2 + 2A3 + 4A4 + 3A5 + 5A6 + 3A7 + 2A8 + A9 + A10` | `22` | `64` |
| `40` | `2A1 + 2A2 + 2A3 + 3A4 + 2A5 + 3A6 + 3A7 + 2A8 + A9 + A10 + A16` | `38` | `64` |
| `41` | `A1 + 2A2 + 2A3 + 3A4 + 3A5 + 5A6 + 3A7 + 2A8 + 2A9 + 2A10 + A13 + A16` | `59` | `64` |
| `42` | `A1 + 2A2 + 2A3 + 4A4 + 2A5 + 4A6 + 3A7 + 2A8 + 2A9 + A10` | `23` | `64` |
| `43` | `A1 + 2A2 + 2A3 + 3A4 + A5 + 3A6 + 3A7 + 2A8 + 2A9 + 2A10 + A13 + A16` | `55` | `64` |
| `44` | `A1 + A2 + 2A3 + 2A4 + 2A5 + 3A6 + 3A7 + 2A8 + 2A9 + 2A10 + A13` | `37` | `64` |
| `45` | `A1 + A2 + 2A3 + 3A4 + 2A5 + 4A6 + 4A7 + 2A8 + A9 + A10 + A13 + A16` | `55` | `64` |
| `46` | `A1 + A2 + 2A3 + 4A4 + 3A5 + 5A6 + 4A7 + 2A8 + A9` | `23` | `64` |
| `47` | `A2 + 2A3 + 4A4 + 2A5 + 5A6 + 4A7 + 2A8 + 2A9 + A10 + A13` | `40` | `64` |
| `48` | `2A1 + 2A2 + 2A3 + 3A4 + 3A5 + 5A6 + 4A7 + 2A8 + 2A9 + A10 + A13 + A16` | `60` | `64` |
| `49` | `A1 + A2 + 2A3 + 3A4 + 3A5 + 4A6 + 3A7 + 2A8 + 2A9 + A10 + A16` | `39` | `64` |

The replay obtains these decompositions by deterministic negative-wall
descent.  If a current root `q` has `q.Ai=-m<0`, it replaces `q` by
`q-mAi=s_i(q)` and records `mAi`.  Each step strictly lowers positive
`D154`-degree and the stored finite trace ends at one of `A1,...,A20`.
The exact coordinate sum is then checked afresh.  Total component
multiplicities range from `2` to `28`, and the largest coefficient is `6`.

None of the forty transported vectors equals any member of the complete
published 20-curve stock.  By adjunction, an irreducible curve in a
square-minus-two class on a K3 surface is a rational `(-2)`-curve.  Thus the
complete stock excludes any missing irreducible representative: the
irreducible count is exactly zero, not merely “not yet identified.”

## The ten square-minus-four divisor-reflection vectors

These ten vectors all have square `-4`, divisibility `2`, and integral
reflections.  Coordinates are again in the `A1,...,A14` basis.  They do not
fall under the square-minus-two effectivity argument.

| facet | transported coordinates | `D154` degree | carrier pairing |
|---:|---|---:|---:|
| `3` | `(0,0,-2,-2,-1,-2,-1,0,1,2,0,1,1,1)` | `22` | `26` |
| `6` | `(0,0,0,0,0,0,0,0,0,0,0,0,1,1)` | `26` | `26` |
| `8` | `(-2,-4,-6,-8,-5,-10,-7,-4,-1,2,2,1,1,1)` | `-16` | `26` |
| `12` | `(0,0,0,0,0,0,0,0,0,2,0,0,1,1)` | `28` | `26` |
| `14` | `(0,0,0,0,0,0,0,0,0,0,0,0,-1,1)` | `-8` | `26` |
| `16` | `(0,-2,-2,-2,-1,-2,-1,0,1,2,0,1,1,1)` | `20` | `26` |
| `19` | `(-2,-2,-2,-2,-1,-2,-1,0,1,2,0,1,1,1)` | `18` | `26` |
| `21` | `(2,4,6,8,5,10,7,4,1,-2,-1,-2,-1,-1)` | `16` | `6` |
| `22` | `(0,0,0,0,0,0,0,0,0,0,0,1,0,-1)` | `-8` | `6` |
| `23` | `(0,0,0,0,0,0,0,0,2,2,0,0,1,1)` | `30` | `26` |

The negative `D154` degrees in this table are another reason not to blur
these vectors into the effective `(-2)` theorem.  The statement here is
exactly that they remain algebraic lattice vectors defining integral
reflections and strict carrier halfspaces; there is no effectivity claim.

## What is new, and the boundary

This removes the previous generic “Weyl-transport and use Riemann--Roch”
handoff.  The marking, Weyl word, final nef vector, all fifty transported
normals, and all forty effective decompositions are now simultaneous and
exactly replayable.  The primitive ample class `2H+D154` additionally gives
an honest common interior point for the algebraic ample chamber and the
50-wall Coxeter chamber.  It is a stronger positive divisor realization of
a known finite-volume Coxeter polytope, not a new abstract polytope and not
a claim that the 50-wall chamber equals the full nef cone.

## Reproduction

- [certificate generator](../scripts/certify_ss2_2a_h13_50_ud8d4_nef_effectivity.py)
- [exact JSON certificate](../data/enumerations/ss2_2a_h13_50_ud8d4_nef_effectivity.json)
- [focused regression test](../tests/test_ss2_2a_h13_50_ud8d4_nef_effectivity.py)
- [underlying complex-K3 carrier](ss2_2a_h13_50_ud8d4_divisor_carrier.md)

Replay from the repository root:

```bash
python3 -m scripts.certify_ss2_2a_h13_50_ud8d4_nef_effectivity --check-archive
python3 -m pytest -q tests/test_ss2_2a_h13_50_ud8d4_nef_effectivity.py
```
