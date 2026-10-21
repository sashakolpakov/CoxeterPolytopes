# Equation-level `A18-C04 H10/28` K3 carrier

## TOOT — positive theorem

The known finite-volume noncompact 28-facet `A18-C04` Coxeter polytope in
`H10` now has a literal equation-level K3 realization.

For

```text
P(t) = -3t^6 + 6t^5 + at^4 + bt^3 + ct^2 + dt - 3r^2,
Q(t) =  2t^9 - 6t^8 + (3-a)t^7 + (1+a-b)t^6
        + et^5 + ft^4 + gt^3 + ht^2 - rdt + 2r^3,
```

the very-general minimal resolution of

```text
y^2 = x^3 + P(t)x + Q(t)
```

on the certified nonempty nine-parameter open set has

```text
fibers:  I4* + I2 + 12I1,
NS(X):   U + D8(-1) + A1(-1),
MW(X):   0.
```

After choosing the marked geometric basis

```text
[F, F+O, eight nonidentity I4* components, the nonidentity I2 component],
```

the certificate gives 28 explicit primitive divisor mirrors and all 28
integral reflection matrices.  The mirrors generate the full Neron--Severi
lattice with index one.  Their types are

| square and divisibility | count |
|---|---:|
| `(-2,1)` | 16 |
| `(-2,2)` | 8 |
| `(-4,2)` | 4 |

Every mirror coordinate has absolute value at most `6`.

## The algebro-geometric anchor

The initial scaled-dual isometry had unwieldy coefficients.  A deterministic
reduction through 226 reflections in Roulleau's 14 published `(-2)` curve
walls produces the small geometric coordinates above.  More importantly, it
sends the sum of the 28 oriented mirrors to

```text
H = 24F + 24(F+O) - 8D1 - 16D2 - 24D3 - 32D4
    - 37D5 - 42D6 - 26D7 - 21D8 - 12C.
```

This class has square `508`, is big and nef, and remains strictly inside the
imported Coxeter chamber, with minimum mirror pairing `3`.

Let `A` be Roulleau's published square-`208` ample class, transported to the
same fiber basis.  Then

```text
A + H = 39F + 38(F+O) - 15D1 - 29D2 - 42D3 - 54D4
        - 62D5 - 69D6 - 40D7 - 35D8 - 18C
```

has square `1408`.  Its minimum pairing with the 14 complete nef walls is
`1`, and its minimum pairing with the 28 imported Coxeter walls is `7`.
Thus it is an explicit ample polarization lying strictly inside the imported
polytope: the imported Coxeter chamber has nonempty intersection with the
K3 ample cone.

This is the most geometric part of the result.  The polytope is not merely
matched to the name of an abstract lattice; one of its interior points is an
honest ample divisor on the equation-level surface.

## Reusable ample-intersection criterion

This calculation isolates a useful gate strictly weaker than “the imported
chamber is the whole nef cone.”  Let `C` be an oriented rational Coxeter cone
in `NS(X)_R`, cut out by `(x,r_i)>=0`.  Then its relative interior meets the
ample cone if and only if its relative interior meets the nef cone in the
positive cone.  The nontrivial implication is constructive: if `H` is nef,
has positive square, and satisfies `(H,r_i)>0` for every imported wall, then
for any ample `A` the class `mH+A` is ample and remains strict for every
imported wall once `m` is large enough.

When a complete rational-polyhedral nef chamber is known, this becomes a
finite exact feasibility test:

```text
H^2 > 0,
(H, every complete-nef wall) >= 0,
(H, every imported Coxeter wall) > 0.
```

Thus a failed full-nef promotion need not end the example track.  One can
instead search for a big-nef interior witness and promote it to an ample
interior witness.  In the present case no scaling is needed: `H+A` already
works.

## Literal Vinberg gate

The composition preserves the oriented 28-by-28 Gram matrix exactly.  The
pinned Vinberg audit therefore transports without a finite-volume guess:

| invariant | value |
|---|---:|
| hyperbolic dimension | `10` |
| facets | `28` |
| finite vertices | `1,087` |
| cusps | `74` |
| elliptic rank-9 edges | `6,956` |
| connected affine components | `204` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Hence every elliptic edge has exactly two finite-or-ideal endpoints and every
affine component has the required full parabolic neighborhood.  The chamber
is finite-volume and noncompact.

## Exact lineage and boundary

The ingredients are kept separate.

- Classical surface input: Roulleau published the
  `U+D8(-1)+A1(-1)` K3 lattice, fiber configuration, curve basis, ample
  class, and complete nef chamber.

- Previously certified polytope input: the abstract `A18-C04 H10/28`
  chamber and its finite-volume audit were already known in this project.

- New positive gain: the separately replayable equation-to-28-wall
  attachment, the small divisor coordinates, and the ample polarization
  inside the imported chamber.

This mixed 28-wall arrangement is **not the complete nef chamber**.  Four
mirrors have square `-4`; the 24 square-minus-two classes are not asserted to
be irreducible curves.  No automorphism of the surface is claimed to realize
the reflections or preserve the chamber.  This is not a new abstract
polytope, and no worldwide-priority claim is made.

Primary surface lineage: [Xavier Roulleau, *An atlas of K3 surfaces with
finite automorphism group*](https://doi.org/10.46298/epiga.2022.6286),
subsection `U+D8+A1`.

## Reproduction

- [composition generator](../scripts/certify_a18_c04_ud8a1_weierstrass_carrier.py)
- [exact composition JSON](../data/enumerations/a18_c04_ud8a1_weierstrass_carrier.json)
- [focused regression test](../tests/test_a18_c04_ud8a1_weierstrass_carrier.py)
- [standalone equation theorem](ud8a1_weierstrass_family.md)
- [standalone scaled-dual carrier](a18_c04_ud8a1_scaled_dual_carrier.md)

```bash
PYTHONPATH=src:. python3 -m \
  scripts.certify_a18_c04_ud8a1_weierstrass_carrier --check-archive
PYTHONPATH=src:. python3 -m pytest -q \
  tests/test_a18_c04_ud8a1_weierstrass_carrier.py
```
