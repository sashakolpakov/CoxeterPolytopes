# An equation-level `U+D8+A1` K3 family

## Positive theorem

For parameters `r,a,b,c,d,e,f,g,h`, put

```text
P(t) = -3t^6 + 6t^5 + at^4 + bt^3 + ct^2 + dt - 3r^2,
Q(t) =  2t^9 - 6t^8 + (3-a)t^7 + (1+a-b)t^6
        + et^5 + ft^4 + gt^3 + ht^2 - rdt + 2r^3.
```

On an explicit nonempty Zariski-open subset of this nine-parameter space,
the minimal resolution of

```text
y^2 = x^3 + P(t)x + Q(t)
```

is a Jacobian elliptic K3 surface with fibers

```text
I4* + I2 + 12I1.
```

Its very-general Neron--Severi lattice and Mordell--Weil group are

```text
NS(X) = U + D8(-1) + A1(-1),       MW(X) = 0.
```

This gives a literal equation-level surface stock for a rank-11 divisor
lattice of hyperbolic dimension ten.

## Exact open-set witness

The discriminant has the forced form

```text
4P(t)^3 + 27Q(t)^2 = t^2 R(t),      degree R = 12.
```

At zero and infinity the first nonzero coefficients are respectively

```text
9r^2(12cr^2-d^2+12hr),
-9(a^2+6a+12b-12c-12e+9).
```

Their nonvanishing gives exact `I2` and `I4*` fibers.  Squarefreeness of
`R`, together with `gcd(R,P)=1`, gives twelve `I1` fibers.  The specialization

```text
(r,a,b,c,d,e,f,g,h) = (1,0,0,0,1,0,0,0,0)
```

has nonzero exact discriminant of `R` and nonzero resultant `Res(R,P)`, so
the required open set is genuinely nonempty rather than formal.

## Why the very-general lattice is exact

The two reducible fibers and the zero section produce
`U+D8(-1)+A1(-1)`, hence Picard rank at least 11.  All nine normalized
parameters are rationally recoverable from the coefficients.  The unique
`I4*` and `I2` fibers fix infinity and zero; the leading coefficients and
the `6t^5` normalization leave only a finite coordinate stabilizer.  Thus
the family has a nine-dimensional moduli image.  The Picard-rank-at-least-12
Noether--Lefschetz loci have dimension at most eight, so the very-general
rank is exactly 11.

Shioda--Tate then gives Mordell--Weil rank zero.  For a nonzero section the
largest possible local height corrections are `2` at the `D8` fiber and
`1/2` at the `I2` fiber.  The height is therefore at least

```text
4 - 2 - 1/2 = 3/2,
```

which excludes torsion.  Consequently the displayed trivial lattice is
saturated and is the full Neron--Severi lattice.

## Published geometric marking

[Roulleau's atlas](https://doi.org/10.46298/epiga.2022.6286), subsection
`U+D8+A1`, publishes the same lattice with an affine-`D8` fiber, an `I2`
fiber, a section, and a curve basis.  The certificate supplies an explicit
integral determinant-minus-one change of basis from

```text
[F, F+O, eight D8 components, the nonidentity I2 component]
```

to Roulleau's curve basis and verifies the Gram identity exactly.  The
lattice and fibration frame are classical.  The normalized short-Weierstrass
chart above is a project-level equation package; no worldwide priority is
claimed for it.

## Boundary of this standalone result

This file certifies the algebraic K3 family and its full divisor lattice.  It
does **not yet** attach the separate `A18-C04 H10/28` mixed Coxeter chamber.
That composition is a distinct theorem so that the classical surface family
and the new chamber attachment cannot be conflated.  No surface-automorphism
claim is made.

## Reproduction

- [certificate generator](../scripts/certify_ud8a1_weierstrass_family.py)
- [exact JSON certificate](../data/enumerations/ud8a1_weierstrass_family.json)
- [focused regression test](../tests/test_ud8a1_weierstrass_family.py)

```bash
PYTHONPATH=src:. python3 -m scripts.certify_ud8a1_weierstrass_family --check-archive
PYTHONPATH=src:. python3 -m pytest -q tests/test_ud8a1_weierstrass_family.py
```
