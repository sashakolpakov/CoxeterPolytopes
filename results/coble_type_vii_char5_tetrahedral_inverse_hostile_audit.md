# Hostile audit: the characteristic-five Type-VII tetrahedral inverse

Status date: 2026-08-09.

## Verdict

**Green.** An implementation that imports no production helper independently
reconstructs the full chain

```text
raw H10/72 P
  -> five intrinsic plane frames and a forced 4+6 design
  -> an explicit strong Coble surface over algebraic closure of F5
  -> the same marked 72-wall finite-volume chamber.
```

The audit recovers all five frames, tests all `5*24=120` colored markings,
matches every one of the twenty root vectors coefficient by coefficient,
transports all 72 walls integrally, and recomputes the twenty-root cone by a
fresh exact cddlib conversion. The same-`P` conclusion is not inferred from a
`12+4+4` degree profile.

This remains an **L3 one-datum inverse calibration**. It is a project-new
explicit plane realization and strict polytope-to-surface-to-same-polytope
proof for the known `H10/72` chamber. It is not a new abstract hyperbolic
polytope, a new dimension, an L4 theorem for a non-singleton class, or a
worldwide priority claim for the equations.

## The material warning survives

The boundary sextic is not ten-nodal. Its singularities are exactly

```text
four A1 ordinary nodes at the special centers;
six A2 ordinary cusps at the ordinary centers.
```

Accordingly its saturated Jacobian scheme has length

\[
                 4\cdot 1+6\cdot 2=16,
\]

whereas its radical is the reduced degree-ten center scheme. Any description
of this plane sextic as having ten nodes, or of its Jacobian as reduced of
degree ten, is false. One blowup at each center nevertheless makes its strict
transform smooth: one blowup resolves both an ordinary node and an ordinary
cusp as a strict plane curve.

## 1. Independent finite geometry

The audit does not read the production finite-geometry output. It begins with
the six points of the complete `F5` conic

\[
                         xz-y^2=0
\]

in the parametrization

```text
o0=(0:0:1)  o1=(1:1:1)  o2=(1:3:4)
o3=(1:2:4)  o4=(1:4:1)  o5=(1:0:0).
```

It enumerates all 15 perfect matchings of `K6`, forms their three chord
lines, and tests concurrency directly. Exactly ten matchings are concurrent,
at ten distinct poles. An exhaustive factorization of the fifteen edges into
five perfect matchings gives six one-factorizations, with concurrency census

```text
0 concurrent matchings: 1 factorization
4 concurrent matchings: 5 factorizations.
```

The lexicographically first usable factorization, selected without consulting
the target certificate, is

```text
(01)(23)(45), (02)(14)(35), (03)(15)(24), (05)(12)(34)
```

with poles

```text
s0=(1:1:4), s1=(1:3:1), s2=(1:2:2), s3=(1:0:2),
```

and residual matching `(04)(13)(25)`.

From these ten centers the independent linear algebra gives:

- twelve distinct chord lines with exact support `s_i+o_j+o_k`;
- four unique smooth conics, each six-point evaluation matrix having rank
  five and nonzero symmetric-matrix determinant;
- four unique quartics through all six ordinary points and double at three
  special points, each 15-row condition matrix having rank 14;
- no cubic through all ten centers, since the cubic evaluation rank is 10;
- a unique sextic double at all ten centers, since its 30 gradient rows have
  rank 27 in dimension 28.

Adding the ten sextic value rows leaves rank 27, as it must because degree six
is nonzero in characteristic five and Euler's relation applies. The resulting
point coordinates and all conic, quartic, and sextic coefficient vectors agree
exactly with the production certificate only after this independent
reconstruction is complete.

## 2. Exact curve schemes and local types

A separate Singular program factors the four quartics and the sextic over
`F5`, saturates their Jacobian ideals, computes radicals and minimal primes,
and compares with the explicit homogeneous point ideals.

For every quartic it returns

```text
one nonunit irreducible factor, exponent one
projective Jacobian degree       3
radical degree                   3
minimal primes                   3
Jacobian = radical = the three prescribed special points.
```

The independently translated local equations have nonzero quadratic
discriminant at all twelve prescribed quartic double points, with histogram
`2:6, 3:6`. Thus these are exactly ordinary nodes. Arithmetic irreducibility
also implies geometric integrality. If a quartic split geometrically,
Frobenius would act transitively on its components. Every rational node would
then lie on all components; its two branches permit only two conjugate
conics. Their three transverse node intersections use three of the four
Bézout intersections, forcing a fourth singular point, contrary to the exact
singular scheme.

For the boundary sextic Singular returns

```text
one nonunit irreducible factor, exponent one
projective Jacobian degree       16
radical degree                   10
minimal primes                   10
radical = the ten center ideals
Jacobian != radical.
```

Independent affine expansions give nonzero quadratic discriminants at the
four special points, histogram `2:2, 3:2`. At every ordinary point the
quadratic tangent cone is a square, has exactly one `F5` tangent direction,
and the cubic coefficient along that direction is nonzero, with histogram
`1:1, 2:3, 3:1, 4:1`. These are precisely four `A1` points and six `A2`
points. Their total delta invariant is ten, the arithmetic genus of a plane
sextic, so the resolved strict transform is smooth and rational.

Geometric integrality of the sextic follows without a hidden smoothness
assumption. If its single `F5`-irreducible factor split geometrically into a
Frobenius orbit of `r>=2` equal-degree components, every one of the ten
rational singular centers would lie on every component. Any two components
would share ten distinct points, but their common degree is `6/r<=3`, so
Bézout permits at most nine.

## 3. Root/boundary separation

The audit computes local affine expansions rather than importing tangent
records. At the 84 simple root-center incidences, the root tangent direction
does not vanish on the boundary quadratic tangent cone:

```text
line tangencies      1:18, 4:18
conic tangencies     1:12, 4:12
quartic tangencies   1:12, 4:12.
```

At the twelve quartic nodes, the resultant of the quartic and boundary
quadratic tangent cones is nonzero, with histogram `1:6, 4:6`. Hence no root
branch shares a boundary branch. The assigned local intersections exhaust
Bézout:

```text
line       3*2       =  6 = 1*6
conic      6*2       = 12 = 2*6
quartic    6*2+3*4   = 24 = 4*6.
```

There are no residual intersections. The blowup therefore has twenty smooth,
geometrically integral, rational `(-2)`-curves disjoint from the smooth
rational strict transform of the sextic. Since

\[
 |-K|=\varnothing,\qquad |-2K|=\{\widetilde F\},
 \qquad \widetilde F^2=-4,
\]

this is an explicitly constructed strong Coble surface.

## 4. Five raw frames and all 120 exact markings

Starting only with the raw Gram and 72 vectors, the lattice audit finds

```text
wall norms                         -1:52, -2:20
intrinsic boundary b               unique -1 wall orthogonal to all 20 roots
other norm-minus-one walls         51
maximal disjoint clique sizes      7:80, 9:25, 10:5.
```

The five ten-curve frames are, in zero-based wall indices,

```text
21 22 23 24 50 51 56 59 62 65
21 28 29 30 50 51 56 59 62 65
22 28 34 35 50 51 56 59 62 65
23 29 34 39 50 51 56 59 62 65
24 30 35 39 50 51 56 59 62 65.
```

For each frame the raw lattice forces

\[
                     H={b+E_1+\cdots+E_{10}\over3}.
\]

After a deterministic parity correction, every basis
`[H,E_1,...,E_10]` has determinant one and Gram
`diag(1,-1^10)`, and `b=3H-sum(E_i)`. The root coordinates force the same
`12 line + 4 conic + 4 quartic` multiplicity design, including the intrinsic
four-special/six-ordinary split.

For each frame a colored incidence isomorphism search has exactly 24
solutions. Under every one of the `5*24=120` solutions, every raw root vector
equals its geometric divisor class coefficient by coefficient. Every raw
wall also has integral coordinates in every determinant-one frame. The five
full-wall digests are

```text
b5db5263b925f79155f52e36bd20b201ae64de5ebadf4a6f0e2be208bf033c9b
21210102a7a9ecef5e7c12847a519059966d8440f2f545f53ab78e93e0f1b11b
7e9e885aec84541ad39389b9c32891f670242c3e6d5233995dfebcc03fad6cdb
ac284d6dc98f378d4311b1f9af041e6674deb4e5778646f5574d6b963676d06c
ae9ebc49a8df971d7942a40abbd1fe7c373c90fc5f204fdfc70cad3520c54758.
```

Their aggregate digest is

```text
74be27cfc3a7e6d21264b3e85d08b2d1c956e9cfcbaf12fdb1e9f212453d0c49.
```

## 5. Fresh root-cone and cusp computation

The twenty raw root inequalities are converted independently by exact
`cddexec_gmp`. Their rank-ten Coble--Mukai cone has

```text
450 finite rays + 51 ideal rays + 0 hyperideal rays = 501 rays.
```

The root span has rank ten in the determinant-minus-one unimodular CM
lattice, its pairings have census `0:135, 1:30, 2:25`, and all twenty root
reflections are integral involutive lattice isometries. Every root inequality
is irredundant: its active rays span a rank-nine facet. Thus the finite-volume
cone is genuinely the Coxeter chamber of the constructed root reflections,
not a redundant inequality presentation used to infer finite index.

The sum of the primitive rays has square `26,564,430`, minimum root slack
`1,882`, and minimum pairing with every causal ray `4,705`. Thus the rays lie
on a common closed positive-cone sheet and the chamber has finite volume.
The exact ray digests are

```text
all 501   80e1e81adc197bb92e818f52f3e695aa1d4d61fc2f4aecf6a550ea1b99fbb5ec
51 ideal  8621adee95e09a0213356e044f923a0fd902e22aa55ead6ffc518959fbfc73d7.
```

For each of the 51 other norm-minus-one raw walls `E`, the audit forms
`f=E+b`. Every `f` is primitive isotropic in `b-perp`, and exact set equality
holds:

```text
{all 51 cdd ideal rays} = {E+b for all 51 raw E-walls}.
```

This is stronger than a matching count and is the essential bridge for the
same-wall-set return.

## 6. Why the same-72-wall return is not circular

The explicit surface and its twenty effective `(-2)`-curves exist before any
Type-VII classification theorem is invoked. The subsequent handoff uses only
general Coble results from [Kondō's paper](https://arxiv.org/abs/2102.00591):

1. The fresh finite-volume root chamber makes the integral reflection group
   generated by the twenty actual roots finite-index in `O(CM)`.
2. The general finite-index theorem gives finite automorphism group, and the
   `Vinberg2` proposition then says these are all effective irreducible CM
   roots. This does not assume that the constructed surface was already the
   known Type-VII example.
3. Every primitive nef isotropic cusp gives a genus-one fibration. With one
   boundary, the general Halphen description forces an index-two `I1` or `II`
   blown fibre

   \[
                          B+2E=2f,
                 \qquad E=f-b.
   \]

   Conversely any `(-1)`-curve gives `f=E+b`. After division to the primitive
   fibre direction it is vertical; the single-boundary fibre has exactly one
   exceptional component, which forces multiplicity one and excludes both a
   nonprimitive multiple and an extra hidden `(-1)`-curve.
4. The general negative-curve trichotomy for a Coble surface leaves only the
   boundary, `(-2)`-curves, and `(-1)`-curves.

Therefore the complete negative wall stock of the constructed surface is

```text
b + 20 constructed roots + 51 cusp lifts = 72 walls.
```

Because each cusp lift is literally `E=f-b` and the set of `f`'s is literally
the set of raw `E+b` classes, these are the same 51 raw walls. Together with
the exact root markings and boundary, this returns the same marked wall set.

Only then is the already-independent full-`P` archive used for its complete
Vinberg census:

```text
dimension              10
facets                 72
finite vertices    18,761
ideal vertices      1,731
hyperideal rays          0
elliptic edges      143,565
bad edges                 0
bad cusp neighbourhoods   0.
```

Thus the returned chamber is finite-volume and noncompact. Surface existence
or the Type-VII classification is not used to manufacture the carrier; the
classification label is attached only after the exact return.

## Scope and lineage

The abstract chamber is the known project `COBLE-VII-p5-b1-H10-72` datum,
and Kondō's characteristic-five Type-VII surface is the known after-the-fact
classification target. This audit establishes a strict explicit plane
realization and an independently checked same-marked-`P` transition. It does
not establish worldwide novelty for the plane equations.

Most importantly, this is still not the requested general inverse theorem:
it treats one rigid chamber rather than an arbitrary member of a non-singleton
surface-free class and derives no AG classification or dimension bound for
that class. Its correct rating is L3, not L4.

## Reproducibility

- [hostile archive](../data/enumerations/coble_type_vii_char5_tetrahedral_inverse_hostile_audit.json)
- [independent hostile auditor](../scripts/audit_coble_type_vii_char5_tetrahedral_inverse_hostile.py)
- [hostile regression tests](../tests/test_coble_type_vii_char5_tetrahedral_inverse_hostile_audit.py)
- [production theorem](coble_type_vii_char5_tetrahedral_inverse.md)

The auditor imports no production helper. It calls Singular for exact
factorization and saturated singular schemes and `cddexec_gmp` for a fresh
rational H-to-V conversion. Six tests include byte-for-byte hostile archive
replay.
