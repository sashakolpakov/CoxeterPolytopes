# The quartic E8+A1 host is nonreflective

The full root reflection chamber of

`L = <-3(T−1)> + E8 + <1>`, where `T=2cos(pi/8)` and `O=Z[T]`,

has **infinite volume**. This closes the host left open in the [previous two-edge report](fresh_h9_e8_first_20261013.md). The proof constructs two actual lattice automorphisms preserving one complete reflection chamber; their product has infinite order. It does not identify this lattice with D9 or its maximal overclosure, and makes no exclusion claim about arbitrary real H9 Coxeter polytopes.

The [global certificate](fresh_h9_e8_next_20261013_global.json) and [separate Fraction reconstruction](fresh_h9_e8_next_20261013_global_independent.json) both pass. The arithmetic search added three complete first edges, rather than increasing a height bound. Together with the two previous edges, there are five complete first steps, 71 integral-residual profiles, 46 interval representatives, six certified finite corners and fourteen coherent facets. Eight additional spherical profiles in the new scalar direction are rejected by the necessary integrality of their residual target before the norm-orbit census.

## Three new complete steps

| Starting corner and released direction | Integral-residual profiles | Interval points | Primitive first mirrors / new simple walls | New full corner |
| --- | ---: | ---: | ---: | --- |
| E7+A1+A1, release scalar A1 | 11 | 2 | 1 / 1 | E7+A1+A1, 65 mirrors |
| E7+A1+A1, release the E7 end node leaving D6+A1+A1 | 24 | 17 | 1 / 1 | D6+3A1, 33 mirrors |
| Initial E8+A1, release the E8 end node leaving D7+A1 | 14 | 18 | 2 / 1 | E8+A1, 121 mirrors |

The new simple norms are respectively `3(2+T)`, `2`, and `2`. On the D7 edge, two primitive reflective roots meet at the first parameter. An exact cone test identifies the unique simple one: both roots are nonnegative combinations of the retained eight walls and this root. The other primitive root is therefore a mirror of the same finite vertex stabilizer, not a second new facet.

The step certificates are [scalar](fresh_h9_e8_next_20261013.json), [D6](fresh_h9_e8_next_20261013_d6.json), and [D7](fresh_h9_e8_next_20261013_d7.json). Each has a separate Fraction census: [scalar replay](fresh_h9_e8_next_20261013_independent.json), [D6 replay](fresh_h9_e8_next_20261013_d6_independent.json), and [D7 replay](fresh_h9_e8_next_20261013_d7_independent.json). Every accepted root satisfies the actual E8 lattice condition, primitivity, and reflection integrality. All 91 pairs of the final fourteen walls pass the exact Coxeter-angle gate.

The fourteen-wall prefix still has the previous strict spacelike escape, of norm

`(1929321 − 2133504 T + 589824 T²)/634880 > 0`.

The global conclusion below is stronger than that finite-prefix observation.

## Completeness and the rational coordinate change

Write `beta=3(T−1)`. The complete primitive candidate norm alphabet remains the ten classes dividing `6`, as established in the previous lattice triage. All finite labels possible in this field are included through the exact `4 cos²(pi/m)` values for `m=2,3,4,6,8,16`.

At the E7 corner, the previously certified relative norm-one unit `(Z,X)` satisfies `Z²−beta X²=1`. Let `omega=(-1,1,0,…,0)` in the E8 spatial coordinates, with `omega²=2`. The rational Lorentz boost

`time' = Z*time − (X/sqrt(2))*(v,omega)`

`v' = v + ((Z−1)*(v,omega)/2 − beta*X*time/sqrt(2))*omega`

normalizes this corner. It is **not** an E8 lattice automorphism: its matrix in the actual lattice basis has three nonintegral entries. Its inverse nevertheless sends every lattice basis vector to integral time and scalar coordinates. The D6 direction also has integral pairing with `eta=(-1,-1,0,…,0)`. These necessary integrality facts justify the norm equations without admitting any false lattice symmetry. Every candidate is pulled back and checked in the original lattice.

For each finite spherical pairing profile, let `r0` be its fixed projection to the eight retained walls. The complete residual equations are:

- Scalar edge: `z²−beta*x²=q−r0²`.
- D6 edge: `z²−beta*x²=2(q−r0²)`.
- Initial D7 edge: `z²−beta*x²=4(q−r0²)`.

The D7 fundamental weight has square `4`; its time variable is twice the original time coordinate. The residual targets must be positive at all four embeddings. The first two new edges reuse the already certified absolute norm lists. The D7 edge has additional absolute targets, computed by the same certified PARI interface in the [new GP file](../scripts/fresh_h9_e8_next_20261013_d7.gp) and [raw transcript](fresh_h9_e8_next_20261013_d7_raw.txt).

A finite candidate provides an upper endpoint. All relative norm-one unit orbits are then enumerated through that endpoint, including the last preceding and first following representatives. The bounded unit-power proposal is only a way to obtain the finite upper bound; it is not the completeness argument. There is no old-wall compatibility pruning before firstness. The complete finite stabilizer at the start, the finite-contact spherical profiles, the full interval census, and the simultaneous-root cone test determine the next full-chamber sector.

## The two chamber symmetries

Use the final wall order: the original E8 walls `0,…,7`, original scalar wall `8`, earlier first mirrors `9,10`, and the new scalar, D6 and D7 walls `11,12,13`.

The D6+3A1 corner has incident frame

`[2,3,4,5,6,7,8,10,12]`.

The isometry `H` fixes its timelike point and simultaneously swaps the two D6 spinor leaves and the two norm-two A1 factors. In this frame its permutation is

`[0,1,2,3,5,4,6,8,7]`.

The isometry `J` transports the initial E8+A1 frame `[0,…,8]` to the new E8+A1 frame `[1,2,3,4,5,6,7,8,13]`. Its target permutation is

`[8,0,1,2,3,4,6,5,7]`.

Both maps send every incident simple normal to the indicated normal with positive scale `1`; they preserve the form, the actual E8 lattice, and its inverse image. Both satisfy `H²=J²=I`. The positive scale on the transported timelike point is also checked exactly. Thus the argument uses oriented chamber sectors, not just an unoriented diagram match or an ambient rational isometry.

The H corner is joined to the initial corner by the complete initial E7 step followed by the complete D6 step. The J corner is joined directly by the complete initial D7 step. The global check replays all three path orientations, the retained walls, strict departure/arrival signs, and containment in all fourteen half-spaces. It closes the full finite mirror systems at the relevant corners, obtaining `121`, `33`, and `121` mirrors. These are complete incident sectors of one full `W(L)` chamber.

A lattice automorphism preserves the full reflective root arrangement. Mapping one complete oriented corner sector to another in the same chamber therefore preserves that chamber. Consequently both H and J stabilize it. Their product has trace

`−9866664830356008865430581914690`

`−5339800527638310253817310787632 T`

`+16843450439595559276468453333752 T²`

`+9115609690914863918898611662416 T³ > 10`.

The sign is certified exactly. A finite-order real ten-dimensional matrix has trace at most ten, so HJ has infinite order. The admissible integral orthogonal group is discrete. An infinite chamber stabilizer yields infinitely many disjoint translates of a sufficiently small interior ball, forcing infinite chamber volume. Hence this full lattice reflection group is nonreflective; it cannot yield the sought finite-volume H9 example.

## Reproduction and audit scope

The separate Fraction implementations regenerate the new profile and interval censuses, lattice tests, finite closures, pair gates and geometry. The global Fraction checker reconstructs both matrices directly from the source/target frames, rather than trusting producer matrices. The certified PARI absolute-norm and full-unit completeness interface is a shared, explicit dependency; this is not a second implementation of PARI. A [separate source and mathematical review](fresh_h9_e8_next_20261013_global_review.md) also passes both global replay commands and confirms the complete-sector argument.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_d6.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_d6_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_d7.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_d7_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_global.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_next_20261013_global_independent.py --verify
```

This concludes the arithmetic search on this host. Possible indefinite identifications with other lattices remain unproved and are unnecessary for the exclusion.
