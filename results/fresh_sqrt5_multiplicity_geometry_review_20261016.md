# Independent H3 multiplicity and local-geometry review

**PASS.** The [primary multiplicity obstruction](fresh_sqrt5_h3_multiplicity_barrier_20261016.md)
has the correct normalization and proves the exact restricted-model bound
`T >= 4997/139 > 35`. The face behind the obstruction is an actual compact
quadrilateral in a local admissible arithmetic chamber. This does not assert
that the local chamber has a compact ambient completion.

The [independent checker](../scripts/fresh_sqrt5_multiplicity_geometry_review_20261016.py)
and [certificate](fresh_sqrt5_multiplicity_geometry_review_20261016.json) passed
generation and fresh saved replay. They import no payment-model or barrier
producer. Exact field and Gaussian operations are shared only with the
standalone Fraction `Q(sqrt(5))` helper in the earlier
[4202 geometry witness](fresh_sqrt5_adjacent_discount_witness_20261016.md).

## Literal face and arithmetic validity

Class7661 has ten roots, white roots4 through9, and nonright bonds

```
0--6:5, 1--9:5,
2--3:3, 2--5:3, 3--7:3, 3--8:3,
4--5:3, 4--6:3, 8--9:3.
```

All roots have norm2, with products `-1` on3-bonds and `-phi` on5-bonds.
The exact full Gram has physical inertia `(8,1,1)` and conjugate inertia
`(9,0,1)`. Deleting root1 gives a nonsingular Lorentz basis. A strict future
center has square `381-421 sqrt(5)/2 < 0`.

The six white roots are positive definite. All four adjacent black-pair
supports are positive definite and have actual future timelike vertices,
strict at the other two black inequalities. The checker also tests all six
black pairs and both orientations of their one-dimensional kernels. Exactly
the four adjacent pairs are feasible. Because the projected side normals
span the Lorentz three-space, these rays exhaust its pointed side cone.
The whole actual white face is therefore compact.

Adjoining the orthogonal spectator forest `A1+4H3` realizes the exact common
context in dimension21. The resulting 23-root Gram has physical inertia
`(21,1,1)` and conjugate inertia `(22,0,1)`; all529 Cartan entries are integral
over `Z[phi]`. A strict timelike center is checked explicitly. Thus neither
the field nor local arithmetic signature excludes this face/context.

## Count and normalization audit

The four corner forests before adjoining spectators are:

| Black pair | Corner forest | Relevant role |
|---|---|---|
| 0,1 | A1+H3+H4 | H3 and H4 are both excluded destination hubs |
| 0,3 | A4+H4 | No H3 black root |
| 2,1 | A1+A4+H3 | H3 black1 is the5-edge leaf; A4 black2 is an endpoint |
| 2,3 | E8 | One adjacent3 baseline payment, eligible for the discount |

At corner `(2,1)`, the common external four H3 copies give **five actual H3
components**; the external A1 and core A1 give two actual A1 components.
The designated white facet root4 is internal in A4, so it is neither black
root and does not remove the black-pair payment. The old H3 cross rule pays
each H3-root/source-root pair `beta/(3c)`, hence this corner pays `beta/15`.
Its H3 black root is a leaf, so the newer central-root compensation is absent.

A direct source audit of `compile_corner`, `add_bonus`, and the E6 exemption
confirms that the other bonus channels do not apply: the first two corners
have no eligible allocated cross type, while the last has a directly
adjacent E8 pair and no allocated long-pair bonus. The E6 exemption does not
change this E8 edge. The frozen full sparse profile is therefore exactly

```
1 - delta_adjacent3 + beta/15.
```

The reviewer checks this full profile against the primary JSON, but does not
reimplement the entire roughly32,000-variable payment compiler. The graph,
root roles, actual multiplicities, scalar counts, and geometry are separately
reconstructed.

For the marked-A1 forest `2A1+H3+4A4`, one isolated A1 root is omitted from
the marked angle count. It removes no intra-component payment or H3/A4
pair. The remaining baseline is `3+4*5=23`. There are eight A4 endpoints and
three H3 roots, each contributing `beta/3` per endpoint. Hence this vertex
has the necessary lower estimate `T >= 23+8 beta`.

For marked `A1+D20`, omitting the A1 leaves the entire D20 component. An
independent breadth-first distance count gives19 adjacent pairs and19 pairs
at distance2. With the assumed distance2 discount zero, its estimate is
`T >= 38-19 delta`. All omitted bonus terms in both vertex estimates are
nonnegative. The same beta is selected because the current coefficient
depends on A1 count2, I2(5) absence, and H3 presence, **not actual H3 count**.

## Exact consequence and permissible repair

The face requirement gives `beta >= 15 delta`. Therefore

```
T >= 23+120 delta,
T >= 38-19 delta.
```

Multipliers `19/139` and `120/139`, respectively, cancel delta and give

```
T >= (19*23+120*38)/139 = 4997/139 > 35.
```

This proves a model obstruction once this literal local face domain and the
two forest bounds are imposed. It is not a new dimension exclusion or a
geometric nonexistence theorem. Giving the H3-to-A4-end coefficient separate
values at actual counts1 and5 is mathematically legitimate: the rule can
depend on the complete finite vertex forest, while each face corner uses
its own actual count with common external spectators. It breaks this
specific dual coupling. It does not prove that the successor optimization
is feasible or covers all affected faces.

```sh
python3 scripts/fresh_sqrt5_multiplicity_geometry_review_20261016.py --verify
```
