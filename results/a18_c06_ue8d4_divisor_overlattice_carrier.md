# Equation-level `U+E8+D4` carrier for the `H13/31` chamber

## Positive theorem

The finite-volume noncompact `A18_4_0_C06` Coxeter chamber with 31 facets in
`H13` has an equation-level K3 divisor carrier.  Its archived rank-14 wall
lattice `M` has signature `(1,13)`, determinant `-64`, and a unique
nontrivial even overclosure invariant under all 31 facet reflections.  This
overclosure is

```text
[N:M] = 4,        det(N) = -4,
N = U + E8(-1) + D4(-1).
```

The certificate supplies an explicit integral 14-by-14 inclusion matrix
`I`, of determinant `4`, satisfying

```text
I^T Gram(N) I = Gram(M).
```

Every one of the 31 facet reflections extends through `I` to an integral
unimodular reflection of `N`.

## The surface family

The ambient lattice is the exact Néron--Severi lattice of a very general
member of the already certified three-parameter Weierstrass family

```text
y^2 = x^3 + t^3(t-1)(t-a)(t-b)(t-c),
```

where `a,b,c` are nonzero, pairwise distinct, and different from one.  Its
reducible fibers are `I0* + II*`, giving `D4(-1)+E8(-1)` together with the
fiber and zero section `U`; the archived moduli argument proves that the very
general Picard rank is exactly 14.

This produces a full-rank mixed divisor-reflection subsystem on an explicit
algebraic K3 family.  It is a new equation-level realization in this project
of the already known `H13/31` chamber.

## Primitive mirrors

Primitive reduction in the ambient divisor lattice is essential.  In the
wall-lattice presentation the normal squares are

```text
14 of square -2,       5 of square -4,       12 of square -8.
```

The twelve square-`-8` presentation vectors acquire coordinate content two
in `N`.  After dividing every embedded vector by its content, the actual
primitive divisor mirrors are

```text
26 of square `-2` and divisibility 1,
 5 of square `-4` and divisibility 2,
 0 of square `-8`.
```

For each wall the certificate checks the square law, equality of the two
reflection formulas before and after primitive reduction, integrality,
involutivity, preservation of `NS`, and intertwining with the original
reflection on `M`.

After Weyl transport of the marking, a strict chamber-interior class is
ample.  Its oriented pairing with every primitive wall is positive, with
minimum one; K3 Riemann--Roch then selects the effective sign of all 26
square-`-2` divisor classes.  The five primitive square-`-4` mirrors are
algebraic divisor reflections, not asserted to be irreducible curves.

## Literal Vinberg finite-volume gate

The inherited exact audit is the full Vinberg test, not merely a list of
finite-looking subdiagrams:

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

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
affine component has the required full-rank zero neighborhood.  Therefore
the chamber is finite-volume and noncompact.

## Exact boundary of the result

This is not a new abstract polytope: the same chamber already had an exact
arithmetic certificate and a special-period order-four K3 realization.  The
new gain is a simple equation-level realization on a very general explicit
K3 family with the complete reflection bridge written down.

It is not the full nef cone.  The five primitive square-`-4` mirrors are not
K3 curve walls, and additional square-`-2` classes in the index-four ambient
Picard lattice may refine the 31-wall chamber.  Nor is it claimed that the
family's order-six automorphism preserves this chamber.

## Reproduction

- [certificate generator](../scripts/certify_a18_c06_ue8d4_divisor_overlattice_carrier.py)
- [exact JSON certificate](../data/enumerations/a18_c06_ue8d4_divisor_overlattice_carrier.json)
- [regression test](../tests/test_a18_c06_ue8d4_divisor_overlattice_carrier.py)
- [explicit Weierstrass family](ue8d4_order6_weierstrass_family.md)
- [original K3 lift](a18_c06_h13_31_k3_order4_lift.md)

Replay from the repository root:

```bash
python3 -m scripts.certify_a18_c06_ue8d4_divisor_overlattice_carrier --check-archive
python3 -m pytest -q tests/test_a18_c06_ue8d4_divisor_overlattice_carrier.py
```
