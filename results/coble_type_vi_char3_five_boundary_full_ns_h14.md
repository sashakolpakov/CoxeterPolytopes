# TOOT: the five-boundary Type-VI Coble surface has the complete `H14/60` nef chamber

## Outcome

Kondō's Type-VI Coble surface in characteristic three, with five boundary
components, has a complete sixty-wall nef cone in its full Néron--Severi
space.  Projectivizing the positive cone gives a finite-volume noncompact
Coxeter polytope in `H14`:

| datum | exact value |
|---|---:|
| boundary `(-4)`-curves | 5 |
| `(-2)`-curves | 10 |
| `(-1)`-curves | 45 |
| total facets | 60 |
| finite vertices | 27,845 |
| ideal vertices | 1,027 |
| elliptic edges | 244,330 |
| hyperideal extreme rays | 0 |

The arithmetic mirror chamber is exactly Alexeev's published
`(r,a,delta)=(15,7,1)` chamber, equivalently the already certified central
section of Shimada's `D5` chamber.  Thus this is a **new project-level full-nef
Coble realization of a known abstract polytope**, not a new abstract
`H14/60` polytope or a new arithmetic reflection group.

The exact replay is
[`certify_coble_type_vi_char3_five_boundary_full_ns_h14.py`](../scripts/certify_coble_type_vi_char3_five_boundary_full_ns_h14.py),
and the frozen certificate is
[`coble_type_vi_char3_five_boundary_full_ns_h14.json`](../data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json).

## Do not confuse the two Type-VI characteristics

Kondō's classification has two different Type-VI Coble rows:

- characteristic `3`, with **five** boundary components;
- characteristic `5`, with **one** boundary component.

The theorem here concerns the first row.  The five fixed points of the
Cremona involution in characteristic three become five fixed curves on the
resolution and hence five boundary components on the quotient.  The earlier
Type-VI Enriques-degeneration discussion that produces the one-boundary
Coble surface belongs to the characteristic-five row and cannot be used as
the source of this chamber.

The source tarball and main TeX file are pinned in the certificate by SHA-256.
The used input is Kondō's classification table, Example `TypeVI3`, Figure
VI, and his quoted Coble--Mukai nef, Vinberg, primitive-isotropic/fibration,
and Halphen results.

## From twenty Coble--Mukai roots to the full divisor lattice

The twenty effective irreducible roots in Kondō's Figure VI admit a
uniform five-set description.  Index them by the ten two-subsets of
`{1,2,3,4,5}`:

- the ten `e` roots have Petersen adjacency, so two pairs meet when they are
  disjoint;
- the ten `f` roots have triangular-graph adjacency, so two pairs meet when
  they share one point;
- `e_ij . f_kl = 2` precisely when `{i,j}={k,l}`.

Their `20 x 20` Gram matrix has rank ten.  An explicit ten-root minor has
determinant `-1` and inertia `(1,9,0)`, so these roots recover the complete
Coble--Mukai lattice

```text
CM(S) = II_(1,9) = U + E8(-1).
```

Exact double description of their Coble--Mukai chamber gives 657 rays:
600 finite and 57 primitive isotropic.  There are no hyperideal rays.

Write `h_i=B_i/2`.  The five `h_i` are mutually orthogonal with square `-1`
and are orthogonal to `CM(S)`.  The ten visible boundary-pair exceptional
curves are

```text
E_ij = (f_ij - h_i - h_j)/2.
```

Adjoining these classes to the twenty-root data determines the integral
glue.  The reconstructed Néron--Severi lattice is odd unimodular of signature
`(1,14)`, hence

```text
NS(S) = I_(1,14).
```

This also follows geometrically from `B=-2K`, the five disjoint
`B_i^2=-4`, and rationality: `K^2=-5` and `rho(S)=15`.

## Why the list of sixty negative curves is complete

Kondō proves that the displayed twenty classes are all effective
irreducible Coble--Mukai roots.  They consist of ten actual `(-2)`-curves and
ten roots arising from `(-1)`-curves meeting two distinct boundaries once.
The remaining `(-1)`-curves can be recovered without a bounded curve search.

For every `(-1)`-curve `E`, adjunction and `B=-2K` give

```text
B . E = 2.
```

If this incidence is split as `1+1` between two boundaries, then
`2E+h_i+h_j` is an effective irreducible root and therefore is already one
of Kondō's ten boundary-pair roots.  If both incidences lie on `B_i`, then

```text
f = E + h_i
```

is a nef isotropic Coble--Mukai class.  The primitive-isotropic/fibration
correspondence and the Halphen blow-up description identify the primitive
ray of `f` with one of the 57 cusps.  Conversely, an integral class
`E=f-h_i` has `E^2=K.E=-1`; the corresponding boundary fibre, equivalently
Riemann--Roch together with the fibre description, makes it an exceptional
curve.  The exact lattice glue performs the only remaining test.

Of the 57 primitive cusps:

```text
22 lift at no boundary;
 7 lift only at B1;
 7 lift only at B2;
 7 lift only at B3;
 7 lift only at B4;
 7 lift only at B5.
```

Thus there are exactly 35 one-boundary exceptional curves, in addition to
the ten boundary-pair exceptional curves.  The resulting negative-curve
stock is exactly

```text
5*(-4) + 10*(-2) + 45*(-1) = 60 curves.
```

The exact Gram check shows nonnegative intersection between every two
distinct listed curves.  Together with the preceding exhaustion this also
rules out a hidden decomposition of a lifted square-minus-one class: any
negative irreducible component would have to be one of the same listed
classes, while all distinct intersections are nonnegative.

Since every negative curve on a Coble surface is a boundary component, a
`(-2)`-curve, or a `(-1)`-curve, these sixty inequalities cut out the entire
nef cone.  Curves of nonnegative square lie in the future light cone and add
no further inequality there.  This is therefore a full-nef result, not a
weak divisor carrier.

## Literal finite-volume proof

All sixty wall pairings are acute and of Coxeter, parallel, or
ultraparallel type, and the wall vectors span all of `NS(S)`.  The literal
Vinberg audit checks the two parts that must not be conflated with merely
finding some cusps:

1. every elliptic rank-thirteen edge diagram has exactly two finite or ideal
   endpoint extensions;
2. every affine zero-neighborhood has full parabolic rank thirteen.

The audit enumerates 244,330 elliptic edges and 2,704 connected affine
components.  Its endpoint histogram is

```text
  9,660 edges: 0 finite + 2 ideal endpoints
 79,510 edges: 1 finite  + 1 ideal endpoint
155,160 edges: 2 finite  + 0 ideal endpoints
```

Both the bad-edge count and the bad-parabolic-neighborhood count are zero.
There are 27,845 finite vertices and 1,027 cusps, so the conclusion is
**finite-volume noncompact**.

As an independent cone-level check, exact rational cdd enumeration gives
28,872 extreme rays, split as

```text
27,845 finite + 1,027 ideal + 0 hyperideal.
```

All sixty inequalities are irredundant and every facet has rank fourteen.
The full ray list and inequality matrix are hash-pinned; the approximately
ten-minute enumeration is available as an optional slow replay.

## Exact arithmetic lineage

Adjoin the boundary halves and double the form on
`Pic~(S)=<Pic(S),B_1/2,...,B_5/2>`.  The resulting integral arithmetic lattice
has determinant `128`, Smith diagonal `1^8 2^7`, and exact identification

```text
2 Pic~(S) = U(2) + E8(-1) + A1(-1)^5,
(r,a,delta) = (15,7,1).
```

The sixty primitive mirror normals split as

```text
15 of square -2 and divisibility 1,
35 of square -2 and divisibility 2,
10 of square -4 and divisibility 2.
```

Every reflection is integral, and the sixty mirrors span the whole
arithmetic lattice.  Exact weighted complete-Gram isomorphism identifies all
3,600 entries with the frozen Shimada `D5 H14/60` certificate.  The existing
Shimada-to-Alexeev certificate then gives the transitive equality

```text
Type-VI Coble chamber = Shimada D5 H14/60 = Alexeev (15,7,1).
```

This equality is stronger than commensurability and closes abstract-polytope
lineage.  The positive gain is the complete Coble divisor/nef realization.

## What this contributes to the inverse track—and what it does not

The marked physical Gram contains exactly five square-minus-four boundary
walls.  Their common orthogonal is `II_(1,9)`, and restricting the other 55
inequalities to it recovers exactly

```text
20 negative primitive Type-VI root walls
+35 distinct isotropic redundant cusp restrictions.
```

That is valuable training data for a `POLYTOPE -> SURFACE` recognizer: the
surface leaves a rigid boundary-frame and cusp-lift fingerprint inside the
full chamber.  In this certificate, however, the physical divisor
normalization and boundary marking are inputs inherited from the surface.
It is **not** yet an unmarked theorem saying that an arbitrary arithmetic
polytope with this Coxeter diagram produces a Coble surface.  An inverse
theorem must reconstruct the frame from polytope data alone, prove the
integral de-scaling/glue conditions, and then prove algebraic existence.

## Replay

```bash
PYTHONPATH=src python3 -m scripts.certify_coble_type_vi_char3_five_boundary_full_ns_h14 --check-archive
python3 -m pytest -q tests/test_coble_type_vi_char3_five_boundary_full_ns_h14.py
```

The independent full cone enumeration is deliberately optional:

```bash
PYTHONPATH=src python3 -m scripts.certify_coble_type_vi_char3_five_boundary_full_ns_h14 \
  --check-archive --full-cone
```
