# Nef effectivity and an ample interior for the `A18-C06 H13/31` carrier

## Positive theorem

The equation-level `U+E8(-1)+D4(-1)` realization of the finite-volume,
noncompact `H13/31` Coxeter chamber can be placed explicitly against the
actual nef and ample cones of Roulleau's K3 surface.

After correcting the positive-cone component and applying one exact
235-reflection Weyl transport, the square-`385` carrier interior lies in
Roulleau's complete closed nef chamber.  All 31 mirrors move simultaneously
and remain strict there.  Among the 26 square-minus-two mirrors:

- exactly five are literal members of the published irreducible curve stock;
- the other 21 are reducible, with exact nonnegative decompositions into that
  stock; and
- the five square-minus-four, divisibility-two vectors remain separate, with
  no effectivity claim.

There is also an explicit primitive ample polarization inside every imported
halfspace.  If `H` is twice the transported rational carrier interior and
`D506` is Roulleau's published ample class, then

```text
B_7 = D506 + 7H
    = (578,1171,1765,882,1492,1220,963,707,522,338,169,1,-75,-77),
B_7^2=90666.
```

Its minimum degree on the complete 15-curve stock is `1`, and its minimum
pairing with the 31 imported Coxeter normals is `4`.  Moreover, `k=7` is
minimal in the chosen family `D506+kH`, for nonnegative integral `k`.

## The complete geometric chamber and the fiber-basis bridge

The geometric input is the rank-14 `U+E8+D4` subsection of Xavier
Roulleau's *An atlas of K3 surfaces with finite automorphism group*,
[EPIGA 6 (2022), article 19](https://doi.org/10.46298/epiga.2022.6286),
[arXiv:2003.08985](https://arxiv.org/abs/2003.08985).  Roulleau's direct
enumeration and completeness test gives all 15 irreducible `(-2)`-curves
`A1,...,A15`; `A1,...,A14` are a basis of `NS`.

The geometric elliptic basis is recovered directly from this graph.  Put

```text
F = A11 + 2A12 + A13 + A14 + A15,       O = A10.
```

Then `F` is the `I0*` fiber class, `F^2=0`, and `F.O=1`.  The columns

```text
F, F+O;
A1,A2,A3,A5,A6,A7,A8,A4;
A13,A12,A14,A15
```

have determinant `-1` in Roulleau's basis, and their Gram matrix is exactly

```text
U + E8(-1) + D4(-1).
```

This is the explicit bridge from the equation family's fiber basis to the
complete curve chamber.  Roulleau's ample class is

```text
D506 = (46,93,141,70,120,100,81,63,46,30,15,1,-12,0),
D506^2 = 506.
```

Its degrees on `A1,...,A15` are

```text
(1,1,1,1,1,1,1,1,1,1,1,1,25,1,1).
```

## Positive-cone correction and exact nef reduction

The old equation-to-curve marking sends the square-`385` carrier interior to
a vector of `D506`-degree `-15889`.  It therefore lies in the positive-cone
component opposite the ample cone.  The certificate negates both the marking
and all oriented normals.  This preserves every carrier pairing and changes
the degree to `15889`.

At each subsequent step, the deterministic reduction fires the curve `Ai`
with most negative current pairing, breaking ties by the smallest label.  For
`Ai^2=-2`, the reflection is

```text
s_i(v) = v + (v.Ai) Ai.
```

Every one of the 235 firings has negative wall pairing and strictly lowers the
positive integral `D506`-height.  The exact word, acting from left to right,
is:

```text
11 12 13 6 5 3 2 1 7 8 6 7 5 6 3 4 5 3 2 1 3 5 6 4 3
2 5 3 4 7 6 5 3 2 1 8 7 6 5 3 2 4 3 5 6 7 8 9 14 12
5 6 7 15 3 4 5 6 3 5 2 3 4 1 2 3 5 6 7 8 1 11 2 3 4
5 3 2 1 6 5 3 4 2 3 5 6 7 8 9 3 5 6 2 3 4 5 3 2 1
2 3 5 6 4 3 2 5 3 4 7 6 5 3 2 1 5 6 7 8 3 5 6 7 9
8 12 14 2 3 5 6 4 3 2 5 3 4 7 5 6 3 5 2 3 1 2 1 5 6
7 8 9 10 5 3 4 6 5 3 7 6 5 8 7 6 2 3 4 5 3 2 1 11 8
9 7 8 12 10 9 5 11 10 7 6 3 2 4 3 5 3 4 2 3 6 7 8 9 5
6 7 8 5 6 2 3 4 5 3 2 7 1 5 6 3 5 3 4 3 2 9 1 3 6
7 5 6 5 3 10 11 10 15 12
```

The product `W` is an integral Gram isometry of determinant `-1`; after the
signed fiber bridge, the combined isometry has determinant `+1`.  Its endpoint
is

```text
v_nef = (38,77,116,58,98,80,63,46,34,22,11,0,-9/2,-11/2),
v_nef^2 = 385,                  v_nef.D506 = 525.
```

The pairings with the complete curve stock are

```text
(1,0,1,0,0,1,0,5,0,1,0,1,9,11,0).
```

Thus `v_nef` is in the closed nef chamber.  It lies on six genuine curve
walls, but its pairing with every transported Coxeter normal is positive,
with minimum `1`.

Twice this rational point is the primitive integral big and nef class

```text
H = 2v_nef
  = (76,154,232,116,196,160,126,92,68,44,22,0,-9,-11),
H^2 = 1540,                     H.D506 = 1050.
```

Its degrees on the 15 curves are

```text
(2,0,2,0,0,2,0,10,0,2,0,2,18,22,0),
```

and its minimum imported-wall pairing is `2`.

## Exact effectivity and irreducibility of all 26 roots

Each oriented square-minus-two vector has positive `D506`-degree.  The replay
descends it through negative curve walls until it reaches one of
`A1,...,A15`.  If `q.Ai=-m<0`, then

```text
q = s_i(q) + m Ai,
```

and the `D506`-degree strictly drops.  The stored finite trace and the final
coordinate equality certify every decomposition below.

Because Roulleau's stock is complete, equality with one of its members decides
irreducibility.  The five irreducible imported mirrors are exactly

```text
facet 0 = A1,   facet 3 = A12,   facet 8 = A3,
facet 10 = A10, facet 14 = A6.
```

The remaining 21 are reducible.  In the table, the carrier pairing is against
`v_nef`.

| facet | exact effective divisor | classification | `D506` degree | carrier pairing |
|---:|---|---|---:|---:|
| `0` | A1 | irreducible | `1` | `1` |
| `1` | 2A1 + 4A2 + 6A3 + 3A4 + 5A5 + 4A6 + 3A7 + 2A8 | reducible | `29` | `22` |
| `2` | A9 + A10 + A11 | reducible | `3` | `1` |
| `3` | A12 | irreducible | `1` | `1` |
| `7` | A5 + A6 + A7 | reducible | `3` | `1` |
| `8` | A3 | irreducible | `1` | `1` |
| `9` | 2A1 + 4A2 + 6A3 + 3A4 + 5A5 + 4A6 + 3A7 + 2A8 + 2A9 | reducible | `31` | `22` |
| `10` | A10 | irreducible | `1` | `1` |
| `11` | A11 + A12 + A15 | reducible | `3` | `1` |
| `14` | A6 | irreducible | `1` | `1` |
| `15` | A3 + A4 + A5 | reducible | `3` | `1` |
| `16` | A2 + A3 + A5 | reducible | `3` | `1` |
| `17` | A2 + A3 + A4 | reducible | `3` | `1` |
| `18` | 2A1 + 4A2 + 5A3 + 2A4 + 4A5 + 4A6 + 4A7 + 3A8 + 3A9 + 2A10 + A11 + A12 + A15 | reducible | `36` | `29` |
| `19` | 2A1 + 4A2 + 5A3 + 3A4 + 5A5 + 4A6 + 3A7 + 3A8 + 2A9 + 2A10 + 2A11 + A12 | reducible | `36` | `29` |
| `20` | A2 + 2A3 + A4 + 2A5 + 2A6 + 2A7 + 2A8 + 2A9 + A10 | reducible | `15` | `15` |
| `21` | A2 + 2A3 + A4 + 2A5 + 2A6 + 2A7 + 2A8 + A9 + A10 + A11 | reducible | `15` | `15` |
| `22` | 2A1 + 4A2 + 5A3 + 3A4 + 5A5 + 4A6 + 3A7 + 3A8 + 3A9 + 2A10 + A11 + A12 + A15 | reducible | `37` | `29` |
| `23` | 2A1 + 4A2 + 7A3 + 4A4 + 6A5 + 6A6 + 6A7 + 5A8 + 3A9 + 2A10 + A11 + A12 + A15 | reducible | `48` | `43` |
| `24` | 2A1 + 4A2 + 7A3 + 3A4 + 7A5 + 6A6 + 5A7 + 5A8 + 4A9 + 2A10 + 2A11 + A12 | reducible | `48` | `43` |
| `25` | 2A1 + 4A2 + 5A3 + 2A4 + 4A5 + 4A6 + 4A7 + 3A8 + 2A9 + 2A10 + 2A11 + A12 | reducible | `35` | `29` |
| `26` | A1 + 2A2 + 2A3 + A4 + 2A5 + 2A6 + 2A7 + 2A8 + A9 | reducible | `15` | `15` |
| `27` | 2A1 + 4A2 + 6A3 + 2A4 + 6A5 + 6A6 + 6A7 + 6A8 + 5A9 + 4A10 + 3A11 + 2A12 + A15 | reducible | `53` | `50` |
| `28` | 2A1 + 4A2 + 7A3 + 4A4 + 6A5 + 6A6 + 6A7 + 5A8 + 4A9 + 2A10 + 2A11 + A12 | reducible | `49` | `43` |
| `29` | 2A1 + 4A2 + 7A3 + 3A4 + 7A5 + 6A6 + 5A7 + 5A8 + 3A9 + 2A10 + A11 + A12 + A15 | reducible | `47` | `43` |
| `30` | 2A1 + 4A2 + 6A3 + 4A4 + 6A5 + 6A6 + 6A7 + 6A8 + 5A9 + 4A10 + 3A11 + 2A12 + A15 | reducible | `55` | `50` |

All supports are connected.  The reducible total component multiplicities
range from `3` to `55`, and the largest component coefficient is `7`.

## The five square-minus-four reflections

These are integral divisor-reflection vectors of square `-4` and divisibility
two.  They are not curve roots and do not enter the effectivity argument.

| facet | transported coordinates | `D506` degree | carrier pairing |
|---:|---|---:|---:|
| `4` | `(2,4,6,3,5,4,3,2,1,0,-1,-2,0,-1)` | `26` | `9` |
| `5` | `(0,0,0,0,0,0,0,0,0,0,0,0,-1,1)` | `-24` | `2` |
| `6` | `(2,4,6,3,5,4,3,4,3,2,1,0,0,-1)` | `36` | `23` |
| `12` | `(-2,-4,-6,-3,-5,-4,-3,-2,-1,0,1,2,2,1)` | `24` | `9` |
| `13` | `(2,4,6,3,5,4,5,4,3,2,1,0,0,-1)` | `38` | `23` |

## The first simultaneous ample and Coxeter-strict translate

Consider

```text
B_k = D506 + kH,                    k >= 0.
```

Every `B_k` is ample because `D506` is ample and `H` is nef.  The exact
boundary is controlled by imported facet `5`:

```text
D506.facet5 = -24,             H.facet5 = 4,
B_k.facet5 = -24 + 4k.
```

For each `k=0,...,6`, facet `5` is the unique nonpositive imported facet; its
pairings are respectively

```text
-24, -20, -16, -12, -8, -4, 0.
```

Thus `B_6` reaches that hyperplane but is not strict.  The primitive class
`B_7` displayed in the theorem has complete-curve degrees

```text
(15,1,15,1,1,15,1,71,1,15,1,15,151,155,1)
```

and imported-wall pairings

```text
(15,337,17,15,152,4,358,17,15,339,15,17,150,360,15,17,
 17,17,442,442,225,225,443,650,650,441,225,753,651,649,755).
```

Consequently `B_7` is ample and strictly inside the `H13/31` chamber; the
minimum imported pairing is `4`, at facet `5`.

## Finite volume and exact scope

The geometric strengthening leaves the literal Vinberg certificate intact:

| invariant | value |
|---|---:|
| dimension | `13` |
| facets | `31` |
| finite vertices | `3,433` |
| cusps | `238` |
| elliptic rank-12 edges | `28,766` |
| connected affine-component audits | `529` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

This is a positive strengthening of an already known abstract Coxeter
polytope: it supplies an equation-level K3 marking whose imported chamber has
nonempty intersection with the ample cone, and it identifies the complete
effectivity and irreducibility pattern of its square-minus-two mirrors.

It is not the full nef cone.  The full K3 nef chamber has 15 irreducible curve
walls, while the imported chamber has 31 mixed divisor mirrors.  No surface
automorphism realization is claimed, and no worldwide priority claim is made.

## Reproduction

- [certificate generator](../scripts/certify_a18_c06_ue8d4_nef_effectivity.py)
- [exact JSON certificate](../data/enumerations/a18_c06_ue8d4_nef_effectivity.json)
- [regression test](../tests/test_a18_c06_ue8d4_nef_effectivity.py)
- [upstream equation-level carrier](a18_c06_ue8d4_divisor_overlattice_carrier.md)
- [explicit Weierstrass family](ue8d4_order6_weierstrass_family.md)

Replay from the repository root:

```bash
PYTHONPATH=. python3 scripts/certify_a18_c06_ue8d4_nef_effectivity.py --check-archive
PYTHONPATH=. pytest -q tests/test_a18_c06_ue8d4_nef_effectivity.py
```
