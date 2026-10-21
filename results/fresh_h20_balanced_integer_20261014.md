# Balanced H20: integer-gap caps and an exact family obstruction

An exact 810-facet Coxeter chamber extends the earlier 418-facet chamber by
392 jointly compatible radius-two caps. Each new cap cuts one of the 392
previously escaping probes `2z/5`. **The chamber still has infinite volume.**
A stronger exact enumeration proves that the entire radius-two integer-gap
family cannot complete the retained baseline to finite volume, regardless
of how many such caps are chosen or whether their mutual compatibility is
relaxed.

The baseline consists of the 22 affine `A6 × A6 × A7` seed walls, the
central radius-two cap, and 392 vertex caps of radius `sqrt(28)`. The 810-wall
construction also retains the three-wall motif from the
[418-wall chamber](fresh_h20_balanced_20261014_bridge.md).

## Construction and exact verification

The seed inradii are one. For a radius-two cap, prescribe its 22 nonnegative
integer seed gaps, whose sums in the three factors must be `7,7,8`.
Gap zero gives a right angle; gap one gives `pi/3`; gap two gives parallel
walls; larger gaps give disjoint walls. The factor composition counts are
`1716,1716,6435`. The center is uniquely reconstructed from these gaps using
the inverse finite seed Gram matrix.

The [numerical scout](../scripts/fresh_h20_balanced_integer_caps_20261014.py)
found 555 caps cutting `2z_000/5`, all compatible with the baseline. Exactly
two also retain the old three-wall motif at this particular probe. Cyclic
seed permutations transport these recipes to all 392 vertices. The
[joint selection](../scripts/fresh_h20_balanced_integer_joint_20261014.py)
found 392 mutually compatible caps covering all 392 probes while retaining
the motif. Its [saved proposal](fresh_h20_balanced_integer_joint_20261014_motif.json)
is numerical; it is not itself a proof of a chamber.

The separate [exact producer](../scripts/fresh_h20_balanced_integer_exact_20261014.py)
reconstructs all 810 rational Lorentzian normals and verifies:

| check | exact result |
|---|---:|
| distinct actual facets | 810 |
| wall pairs | 327,645 |
| right pairs | 10,405 |
| `pi/3` pairs | 6,318 |
| parallel pairs | 410 |
| disjoint pairs | 310,512 |
| retained finite endpoints | 392 |
| removed `2z/5` probes | 392 |

Every facet has a strict timelike relative-interior witness. Every old
finite endpoint retains exactly 20 incident walls. The full root list,
selected gap recipes, source hashes and strict escape are in the
[exact certificate](fresh_h20_balanced_integer_exact_20261014.json).
Generation and saved replay passed:

```sh
python3 scripts/fresh_h20_balanced_integer_exact_20261014.py --verify
```

The [independent chamber review](fresh_h20_balanced_integer_exact_20261014_review.md)
reconstructs all normals and all pair relations with standard-library
integer and rational arithmetic, without importing producer code or FLINT.
It also diagnoses an actual affine `A2 + 5 A1` cusp link and explains why
this deficient crystallographic link admits a local repair. Such a repair
does not provide global finite volume.

The spacelike point with spatial coordinate `3z_000/5`, Lorentz square one,
and final coordinate one is strict on all 810 inequalities. Together with
the strict timelike interior point it yields an open patch at infinity,
proving infinite volume. The construction has no claim of a complete new
vertex or cusp census.

## The whole integer-gap family is insufficient

The [integer-only enumeration](../scripts/fresh_h20_balanced_integer_obstruction_20261014.py)
independently rebuilds the inverse seed matrices from the type-A formula.
It does not import the scout's numerical geometry. All spatial squared
distances lie in `(1/28)Z`; distances to `pz/5` lie in `(1/700)Z`.
The factor distances are nonnegative, so discarding a factor whose distance
already exceeds four cannot lose any radius-two cap meeting the probe.

The enumeration considers every composition whose cap contains or is
tangent to either probe. Exact results are:

| probe | strictly containing caps | compatible with baseline | tangent caps |
|---|---:|---:|---:|
| `2z_000/5` | 555 | 555 | 0 |
| `3z_000/5` | 555 | 0 | 0 |

The first row exactly reproduces the scout's 555 gap recipes. For every
rejected cap in the second row, the certificate records an incompatible
retained endpoint wall. The endpoint gate covers all finite Coxeter labels:
its cosine has the form `r/sqrt(7)` with rational `r`; its square is rational.
For an angle `pi/m` with integer `m≥2`, rational cosine square permits only
`0,1/4,1/2,3/4`. Only zero can be a rational square divided by seven.
Thus right or disjoint contact is necessary; an untested larger finite label
cannot evade the gate.

Every baseline-compatible cap in this family therefore has power at least
`1/700` at `x=3z_000/5`. In coordinates with quadratic form `x²−2DU`, the point

```text
E = (3z_000/5, 98783/2800, 1),       E² = 1/1400
```

is strict on **every** such cap. Its central-cap power is `1664/25`, every
old endpoint-cap power is at least `84/25`, and its minimum seed gap is
`2/5`. Hence the entire family, even without mutual-compatibility constraints,
leaves a strict spacelike escape. Cyclic symmetry gives the same conclusion
at each old vertex direction.

The [whole-family certificate](fresh_h20_balanced_integer_obstruction_20261014.json)
passes generation and saved replay:

```sh
python3 scripts/fresh_h20_balanced_integer_obstruction_20261014.py --verify
```

This rules out completion using only radius-two integer-gap caps while
retaining the baseline. Other radii, noninteger gaps, changed endpoint walls,
and changed seed chambers are outside the result. Adding more caps from the
same dictionary cannot solve the finite-volume problem.
