# Exact five-vertex colored-tree denominator for JK-H12

Status: **PASS/FREEZE; exact Burnside calculation** (2026-08-08).

The frozen 384-element seed action on the 112 pairwise-disjoint doubling
walls gives the following complete search denominator before ordinary
Coxeter-matrix reconstruction:

| unlabelled tree | raw proper colorings modulo tree automorphisms | seed-and-tree orbits |
|---|---:|---:|
| `P5` | 76,587,336 | 201,747 |
| `K1,4` | 6,210,820 | 16,872 |
| degree-3 fork | 75,897,360 | 199,362 |
| **total** | **158,695,516** | **417,981** |

These are colored-tree orbits, not yet Coxeter-polytope isomorphism classes.
Different rows could in principle reconstruct isomorphic ordinary Coxeter
matrices.  The separate shape certifiers perform that stronger comparison.

## Proper-coloring conventions

For `P5`, an edge word `(a,b,c,d)` satisfies

```text
a != b,  b != c,  c != d
```

and is identified with `(d,c,b,a)`.  There is no reversal-fixed proper word,
so the raw count is

```text
112 * 111^3 / 2 = 76,587,336.
```

For `K1,4`, all four incident colors are distinct and unordered:

```text
C(112,4) = 6,210,820.
```

For the degree-3 fork, let `a` be the central edge on the long arm, `d` its
extension, and `{b,c}` the unordered short-arm colors.  Then `a!=d`, the
three central colors `a,b,c` are distinct, and `d` may equal `b` or `c`.
Thus

```text
112 * 111 * C(111,2) = 75,897,360.
```

## Burnside formulas

For a seed-action element `g`, write

```text
F  = number of fixed colors,
T  = number of colors lying in 2-cycles,
F2 = F + T = number of colors fixed by g^2.
```

The `P5` identity-on-tree contribution is `F(F-1)^3`.  Combining `g` with
path reversal contributes `T(F2-1)`: the first and second colors must be
fixed by `g^2`, while the second lies in a genuine 2-cycle.  Therefore

```text
sum_g [F(F-1)^3 + T(F2-1)] = 154,941,696,
154,941,696 / (2*384) = 201,747.
```

For `K1,4`, a fixed four-subset is the coefficient of `x^4` in

```text
product over cycles C of g: (1 + x^|C|).
```

The sum of this coefficient over the group is `6,478,848`, hence

```text
6,478,848 / 384 = 16,872.
```

For the fork, the identity-on-tree contribution on ordered short arms is
`F(F-1)^2(F-2)`.  Combining `g` with the short-arm swap contributes
`F(F-1)T`.  Hence

```text
sum_g [F(F-1)^2(F-2) + F(F-1)T] = 153,110,016,
153,110,016 / (2*384) = 199,362.
```

The machine archive includes the full ten-row seed-action cycle-type table
and all fixed-coloring contributions.

## Scope

This closes the exact combinatorial denominator for all five-vertex tree
shapes.  It does not infer finite volume, root integrality, or pairwise
nonisomorphism of reconstructed Coxeter matrices.  Those gates belong to the
shape census.  In particular, the denominator explains why `K1,4` is the
tractable first complete positive layer: the other two shapes each have
about 200,000 seed orbits.

## Reproduction

```bash
python3 scripts/certify_jacobian_kummer_h12_five_vertex_tree_denominators.py --check-archive
PYTHONPATH=src:. pytest -q tests/test_jacobian_kummer_h12_five_vertex_tree_denominators.py
```

Machine archive:
[`jacobian_kummer_h12_five_vertex_tree_denominators.json`](../data/enumerations/jacobian_kummer_h12_five_vertex_tree_denominators.json).

Certifier:
[`certify_jacobian_kummer_h12_five_vertex_tree_denominators.py`](../scripts/certify_jacobian_kummer_h12_five_vertex_tree_denominators.py).
