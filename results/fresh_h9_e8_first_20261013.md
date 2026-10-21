# Two complete first edges in the quartic E8+A1 host

**Superseded status:** the [subsequent complete-edge and chamber-symmetry proof](fresh_h9_e8_next_20261013.md) establishes that this full E8+A1 host is nonreflective. The two-edge checkpoint below records its earlier open status; its exact artifacts remain unchanged.

The lattice `L = <-3(T−1)> + E8 + <1>`, with `T^4−4T²+2=0`, now has **two certified first edge steps, three actual finite corners and eleven compatible facets**. Separate Fraction implementations replay both complete norm-orbit censuses. The eleven-wall prefix still has an exact strict spacelike escape. The full lattice's reflectivity and possible indefinite conjugacy with earlier hosts remain open; this is not a new finite-volume polytope or a claim of a new full reflection group.

## Structural comparison

The [preceding exact triage](fresh_h9_tripled_d9_20261013_e8_triage.md) establishes the maximal integral lattice, safe primitive norm pool `q | 6` with ten representatives, and full zero-height E8+A1 chamber. Its spatial root system has 121 mirrors. Its displayed reflection set differs from D9's in both directions, but that alone does not settle indefinite conjugacy.

The new [exact comparison](fresh_h9_e8_first_20261013_comparison.json), replayed by its [checker](../scripts/fresh_h9_e8_first_20261013_comparison.py), gives a stronger caution against transferring the closed D9 theorem. The first new binary reflection preserves both L and the standard maximal lattice `M0=<-3(T−1)>+O^9`, but fails the D9 module test on three basis vectors. After a coordinate permutation preserving D9 and M0, it strictly cuts the previously certified first D9 endpoint. Thus M0 has extra mirrors inside that old full D9 chamber. Saturation does not identify the full root arrangements or transfer the nonreflectivity proof.

L and M0 remain distinct maximal integral neighbors in the displayed coordinates, with common ` <-3(T−1)>+D8+<1>` and index ideal (2) in each. Whether another indefinite isometry identifies them has not been decided. Neither has the reflectivity of M0. This round searches a still-open full arrangement rather than assuming that every maximal overclosure of the closed D9 lattice is already excluded.

## The E8-fixed edge

Retain all eight E8 simple walls and release the scalar A1 wall. There is an explicit finite candidate endpoint, so any earlier first contact would also be strictly inside hyperbolic space. A finite rank-nine spherical diagram extending the retained induced E8 diagram can only be E8+A1. Consequently every possible first facet is orthogonal to E8. This finite-contact reduction received a separate mathematical review; it is not an assertion about an unbounded edge or an ideal endpoint.

Write `beta=3(T−1)`. The complete equations are

`z² − beta*x² = q`, for all ten primitive candidate norms `q | 6`.

The [first-edge certificate](fresh_h9_e8_first_20261013.json) contains five corrected relative-norm orbits. Its complete interval has two algebraic points: one is imprimitive and the other is the unique primitive reflective first mirror. It is

```text
r = (x, 0,0,0,0,0,0,0,0, z),
x = −85−48T+143T²+78T³,
z = −135−72T+228T²+123T³,
r² = 3(2+T).
```

The endpoint has scalar coordinate

`beta*x/z = (−3−45T+11T²+12T³)/17`

and is timelike. Both endpoints have complete E8+A1 stabilizers with 121 mirrors. Their distinguished primitive A1 norms have different square classes: the ratio `3(2+T)` has rational norm 162, which is not a square. A K-isometry therefore cannot identify these two pointed full chamber sectors directly.

The [producer](../scripts/fresh_h9_e8_first_20261013.py) and [separate Fraction replay](../scripts/fresh_h9_e8_first_20261013_independent.py) both pass. The [independent result](fresh_h9_e8_first_20261013_independent.json) checks every corrected orbit and interval endpoint, root primitivity, reflection integrality, all embeddings, both finite stabilizers, 45 pair gates and the exact escape.

## The initial E7+A1 edge

Release the E8 simple root at index 1 in the archived basis. Its fundamental weight is `omega=(-1,1,0,...,0)` in the E8 spatial coordinates, with `omega²=2`. This keeps an E7+A1 edge and uses the same relative field, because `sqrt(2)` already belongs to K.

There are exactly **12** integral spherical profiles in the complete field-compatible angle alphabet. For a profile projection r0, the residual equation is

`z²−beta*x² = 2(q−r0²)`,

with normal `r=(x/sqrt(2), r0−z*omega/2)`. All necessary norms, nonsquare residuals, conjugates and relative-unit orbits are included. The complete first interval has **seven** points and one primitive reflective survivor.

Let

```text
Z = −29689−16068T+50682T²+27429T³,
X = −18616−10074T+31780T²+17199T³.
```

The new first mirror is `(sqrt(2)*X, −Z*omega, 0)`, of norm two, using `Z²−beta*X²=1`. Its endpoint has full type **E7+A1+A1**, with **65** mirrors. It is compatible with the preceding binary cap. This gives eleven actual facets and three finite vertices in a single full-chamber prefix.

The [E7 certificate](fresh_h9_e8_first_20261013_e7.json), [producer](../scripts/fresh_h9_e8_first_20261013_e7.py), [separate Fraction checker](../scripts/fresh_h9_e8_first_20261013_e7_independent.py) and [independent result](fresh_h9_e8_first_20261013_e7_independent.json) all pass. The independent code regenerates all profiles in reversed order and replays the complete intervals, 65-mirror closure, all 55 combined pair gates, strict interior and strict escape.

## Completeness dependency and live frontier

Both searches use certified PARI absolute norm representatives and full unit groups in `K(sqrt(beta))`, whose signature is (2,3). The relative norm map on the four fundamental units has exponent vectors `(0,2,0)`, `(-2,0,0)`, `(0,0,2)`, `(0,0,0)` in the base unit group. Hence the norm-one kernel is generated by the fourth unit and sign. The raw certified transcripts and exact unit corrections are archived. A 100-power scout only supplies a finite upper endpoint; every orbit representative between the origin and that endpoint is then exhausted. No old-wall compatibility filter is used before firstness. The separate implementations share this explicit certified PARI completeness dependency rather than pretending to replace the number-field algorithms.

The [frozen frontier](fresh_h9_e8_first_20261013_frontier.json) checks all 55 nine-wall subsets and finds exactly the three stated finite vertices. Its [checker](../scripts/fresh_h9_e8_first_20261013_frontier.py) records all **23** remaining outgoing directions and tests every other existing wall for a finite endpoint. Their full-root first mirrors remain uncomputed.

The common strict escape has square

`(1929321−2133504T+589824T²)/634880 > 0`.

It satisfies all eleven wall inequalities strictly and lies on the future side of the common timelike interior. Thus the current prefix has infinite volume. No further root search is hidden in the frontier census, and no conclusion about the full host follows from this finite-stage escape.

```sh
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013_e7.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013_e7_independent.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013_comparison.py --verify
PYTHONPATH=src:. python scripts/fresh_h9_e8_first_20261013_frontier.py --verify
```
