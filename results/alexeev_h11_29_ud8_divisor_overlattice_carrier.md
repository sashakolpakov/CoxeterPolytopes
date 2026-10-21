# TOOT: an explicit K3 divisor-overlattice carrier for `H11/29`

## Positive theorem

The finite-volume noncompact 29-wall Coxeter chamber in `H11`, previously
obtained as the transposition section of Alexeev's dual 118-wall chambers,
has a genuine equation-level K3 divisor carrier.

Let `M` be its archived rank-12 even reflection lattice.  Then

```text
signature(M) = (1,11),       det(M) = -256.
```

Among every integral overlattice of `M` invariant under all 29 facet
reflections, the possible indices are `1,2,4,8`; the even ones have indices
`1,4`.  Consequently there is exactly one nontrivial reflection-stable even
overclosure, and it satisfies

```text
[N:M] = 4,                   det(N) = -16,
N = U + D8(-1) + A1(-1)^2.
```

This is an explicit isometry, not a match of determinant or discriminant
histograms.  The certificate gives a 12-by-12 integral matrix `I`, with
columns equal to the chosen basis vectors of `M` in the displayed basis of
`N`, such that

```text
det(I) = -4,                 I^T Gram(N) I = Gram(M).
```

Every one of the 29 facet reflections of `M` extends to an integral
unimodular reflection of `N` and intertwines through `I`.

## The actual surface family

The ambient lattice `N` is the exact very-general Néron--Severi lattice of
the already certified four-parameter Weierstrass family

```text
y^2 = x^3 + P(t^2)x + t Q(t^2),

P(s) = -3s^3 + alpha s^2 + beta s + (-3r^2+3-alpha-beta),
Q(s) = 2s^4 - alpha s^3 + gamma s^2 + q1 s + q0,
q1   = 9r-2alpha r-beta r-r^3-8+3alpha-2gamma,
q0   = 2r^3-2+alpha-gamma-q1.
```

On the recorded nonempty parameter open its fibers are
`I4* + 2I2 + 10I1`.  For a very general member, the Mordell--Weil rank and
torsion are zero and

```text
NS(X) = U + D8(-1) + A1(-1)^2.
```

The order-four automorphism

```text
sigma(x,y,t) = (-x, i y, -t)
```

is Brandhorst--Hofmann row `0.4.4.2`.  Its Picard involution is the identity
on `U+D8(-1)` and exchanges the two `A1` roots.  Under the explicit embedding
above it preserves `M` integrally.  It does **not** preserve the selected
29-wall Coxeter chamber: only 14 of the 29 oriented simple walls return to
that oriented simple-wall list.  Thus the theorem is a stable divisor-
overlattice realization, not an automorphism-fixed-section claim.

## Literal hyperbolic-volume gate

The inherited Vinberg audit is complete:

| invariant | value |
|---|---:|
| dimension | `11` |
| facets | `29` |
| finite vertices | `3,876` |
| cusps | `253` |
| elliptic rank-10 edges | `28,798` |
| connected affine-component audits | `668` |
| bad edge endpoints | `0` |
| bad parabolic neighborhoods | `0` |

Every elliptic edge has exactly two finite-or-ideal endpoints, and every
affine component has a full rank-10 zero-neighborhood.  Hence the chamber is
finite-volume and noncompact by Vinberg's literal criterion.

## Primitive ambient mirrors and algebro-geometric scope

There are two normalizations which must not be conflated.  In the archived
`M` presentation the 29 primitive normals have squares

```text
11 of square -2,       12 of square -4,       6 of square -8.
```

After applying `I`, 23 normal vectors remain primitive in `NS(X)` and six
have coordinate content two.  Dividing each embedded normal by its ambient
content gives the actual primitive `NS(X)` mirror split

```text
17 of square -2,       12 of square -4 and divisibility 2,
0 of square -8.
```

Among the 17 square-minus-two classes, 16 have divisibility one and one has
divisibility two.  Scaling a normal does not change its hyperplane or its
reflection.  The certificate nevertheless checks this rather than assuming
it: for every facet, the presentation normal is exactly its content times
the primitive ambient normal, both reflection formulas agree, and that
common reflection extends the archived action on `M`.

Choose a generic point in the strict Coxeter chamber away from all ambient
`(-2)` hyperplanes, and change the marking by the K3 Weyl group so that this
point is ample.  Riemann--Roch then chooses the effective signs of all 17
primitive square-minus-two divisor classes.  The 12 primitive square-minus-
four mirrors are algebraic divisor reflections; they are not asserted to be
irreducible curve walls.

Extra roots of the index-four ambient Néron--Severi lattice can refine the
29-wall chamber.  Therefore this is a **full-rank mixed divisor-subsystem
chamber**, not the full nef cone of the surface.  It adds a new geometric
realization in this project; it does not add a new abstract Coxeter
polytope, arithmetic class, or compact example.

## Bounded primary-source lineage audit

The geometric bridge appears unrecorded in the three directly relevant
primary sources audited here:

- Alexeev's [2-elementary reflection-group paper](https://arxiv.org/abs/2209.09110),
  Section 3, records the 118-root `(16,6,1)` chamber and its `S6` symmetry,
  but not the transposition section or its even index-four overclosure.
- Brandhorst--Hofmann's [classification paper](https://arxiv.org/abs/2112.07715)
  and [database version 2.1](https://zenodo.org/records/10410278) record row
  `0.4.4.2`, but not this 29-mirror chamber, its root list, or the displayed
  Weierstrass equation.
- Roulleau's [EPIGA K3 atlas](https://doi.org/10.46298/epiga.2022.6286),
  Section 12.2, treats `U+D8(-1)+A1(-1)^2` through its 19 curve walls, not
  through this mixed 29-mirror subsystem.

Thus the precise lineage statement is: this equation-to-overlattice-to-
chamber realization **appears unrecorded in these audited primary sources**.
This is bounded evidence, not a worldwide novelty theorem.

## Reproduction

- [certificate generator](../scripts/certify_alexeev_h11_29_ud8_divisor_overlattice_carrier.py)
- [exact JSON certificate](../data/enumerations/alexeev_h11_29_ud8_divisor_overlattice_carrier.json)
- [regression test](../tests/test_alexeev_h11_29_ud8_divisor_overlattice_carrier.py)
- [original `H11/29` theorem](alexeev_118_dual_h11_29_positive_example.md)
- [explicit `U+D8+A1^2` family](ud8a1x2_order4_weierstrass_family.md)

Replay from the repository root:

```bash
python3 -m scripts.certify_alexeev_h11_29_ud8_divisor_overlattice_carrier \
  --check-archive
pytest -q tests/test_alexeev_h11_29_ud8_divisor_overlattice_carrier.py
```
