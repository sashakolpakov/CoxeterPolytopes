# TOOT: a strict polytope-first Type-VII plane realization

Status date: 2026-08-09.

## Verdict

The raw `H10/72` chamber now has an explicit, exact

```text
P -> intrinsic finite design -> strong Coble surface -> the same marked P
```

return.  No known Type-VII surface is imported to obtain the carrier.  The
raw lattice itself produces five and only five plane frames; finite geometry
over `F5` realizes the forced `4+6` design; exact curve and singularity gates
give a strong Coble surface; and all five raw frames return coefficientwise to
the constructed divisor classes.

This is a positive L3 inverse result and a new explicit realization in the
project.  It is not a new abstract hyperbolic polytope, a new dimension, a
worldwide priority claim for the equations, or the requested class-wide L4
polytope-to-surface theorem.

## Hostile warning: the boundary has cusps

The plane sextic does **not** have ten nodes.  It has

```text
4 ordinary nodes at the four special centers;
6 ordinary A2 cusps at the six conic centers.
```

Consequently its saturated Jacobian scheme is nonreduced and has length

\[
                 4\cdot1+6\cdot2=16,
\]

although its radical is exactly the reduced ten-center scheme.  Any
certificate saying “ten nodes” or “reduced Jacobian of degree ten” for this
model is false.  One blowup resolves an A1 node and also resolves an A2 cusp,
so the strict transform is nevertheless smooth, as required.

## 1. What the raw polytope forces

Start only with the 72 primitive wall vectors in `I_(1,10)`.  Their norms are

```text
-1: 52          -2: 20.
```

The boundary half `b=-K` is intrinsic: it is the unique norm-minus-one wall
orthogonal to all twenty roots.  The zero-intersection graph on the other 51
norm-minus-one walls has maximal-clique census

```text
size 7: 80       size 9: 25       size 10: 5.
```

Thus the five plane blowdown frames are exactly

```text
21 22 23 24 50 51 56 59 62 65
21 28 29 30 50 51 56 59 62 65
22 28 34 35 50 51 56 59 62 65
23 29 34 39 50 51 56 59 62 65
24 30 35 39 50 51 56 59 62 65
```

in zero-based raw wall numbering.  For every frame `E_1,...,E_10`, the raw
lattice forces

\[
                        H={b+\sum E_i\over3}.
\]

After a deterministic parity swap in the frame order, every basis
`[H,E_1,...,E_10]` has determinant `1` and Gram
`diag(1,-1^10)`, and

\[
                         b=3H-\sum E_i.
\]

The twenty roots have the same full profile in every frame:

```text
12  H - E - E - E
 4  2H - E - E - E - E - E - E
 4  4H - six simple E - three double E.
```

The double coefficients intrinsically split the ten centers into four
special and six ordinary centers.  The twelve lines give four disjoint
perfect matchings of `K6`, one at each special center.  The unused fifth
matching controls the remaining three conics.  The four quartics omit the
four special points one at a time.

This extraction precedes and does not assume any algebraic surface.

## 2. The finite-geometric realization over F5

Put the six ordinary points on

\[
                            Q_0:xz-y^2=0
\]

as

\[
 p_t=(t^2:t:1),\quad t=0,1,2,3,4,
 \qquad p_\infty=(1:0:0).
\]

In the labels used by the certificate,

```text
o0=(0:0:1)  o1=(1:1:1)  o2=(1:3:4)
o3=(1:2:4)  o4=(1:4:1)  o5=(1:0:0).
```

Among the fifteen perfect matchings of the six points, exactly ten have their
three chord lines concurrent.  There are six one-factorizations of `K6`.
Five contain exactly four concurrent matchings and one contains none:

```text
number of concurrent matchings in a one-factorization:  0:1, 4:5.
```

Choose

```text
M0=(01)(23)(45),  pole s0=(1:1:4)
M1=(02)(14)(35),  pole s1=(1:3:1)
M2=(03)(15)(24),  pole s2=(1:2:2)
M3=(05)(12)(34),  pole s3=(1:0:2)
M4=(04)(13)(25),  residual and nonconcurrent.
```

The twelve root lines are the three chords in each of `M0,...,M3`.  Each has
exact center support consisting of its pole and its two conic endpoints.

The four root conics are

\[
\begin{array}{ll}
C_0=&xz-y^2,\\
C_1=&x^2+2xy+2xz-y^2,\\
C_2=&x^2+2xy-xz-2y^2+yz-z^2,\\
C_3=&xz+2y^2-2yz+2z^2.
\end{array}
\]

Here `C0` contains all six ordinary points; `C1,C2,C3` contain all four
special points and respectively the pairs of `M4`.  All four six-point
evaluation matrices have rank five, their center supports are exact, and
their symmetric matrices have nonzero determinant.  They are therefore
smooth and geometrically integral.

For each special point `s_i`, impose passage through all six ordinary points
and a double point at the other three special points.  The resulting fifteen
quartic conditions have rank fourteen and a one-dimensional kernel.  With
`R_i` denoting the quartic omitting `s_i`, the four equations are

\[
\begin{aligned}
R_0={}&x^3y+x^3z+x^2y^2+x^2yz-x^2z^2-2xy^3+2xy^2z+xyz^2\\
     &-xz^3-y^4-2y^3z-y^2z^2+yz^3,\\
R_1={}&x^3y+2x^3z+2x^2y^2-x^2yz-2x^2z^2+2xy^3-xy^2z+xyz^2\\
     &-2xz^3-2y^4-2y^3z-2y^2z^2-yz^3,\\
R_2={}&x^3y-x^3z-x^2y^2+2x^2yz+x^2z^2+xy^3-2xy^2z+xyz^2\\
     &+xz^3+y^4-2y^3z+y^2z^2+2yz^3,\\
R_3={}&x^3y+2x^2yz+x^2z^2+xy^3-2xy^2z+xyz^2+y^4-2y^3z+2yz^3.
\end{aligned}
\]

Singular proves each `R_i` irreducible over `F5`, with saturated Jacobian
equal to the three prescribed reduced special-point ideals.  Every one of
the twelve double points is an ordinary node.

These quartics are geometrically integral.  Indeed, arithmetic irreducibility
makes any geometric components a Frobenius orbit.  An ordinary node permits
at most two components, so the only possible split is into conjugate conics.
They meet transversely at the three rational nodes, contributing three to
Bézout's total four; the forced fourth intersection would be another
singular point, contradicting the exact three-point singular scheme.  Their
strict transforms are smooth rational `(-2)`-curves.

## 3. The unique boundary sextic

The cubic evaluation matrix at all ten centers has rank ten, hence

\[
                              |-K|=\varnothing.
\]

The thirty first-derivative conditions on sextics have rank 27 in the
28-dimensional sextic vector space.  Their kernel is generated by

\[
\begin{aligned}
F={}&x^4y^2+x^4yz-x^4z^2-x^3y^3+2x^3yz^2+x^3z^3\\
 &+x^2y^4-x^2y^3z+x^2yz^3+x^2z^4-xy^5-2xy^4z\\
 &+2xy^3z^2+2xyz^4-y^6+2y^5z-y^4z^2-2y^3z^3+y^2z^4.
\end{aligned}
\]

Euler's identity in characteristic five makes the value equations
consequences of the gradient equations; adding them leaves rank 27.  Thus on
the blowup

\[
                         |{-2K}|=\{\widetilde F\}.
\]

Exact factorization says `F` is irreducible over `F5`.  Its saturated
Jacobian has length 16 and ten minimal primes; its radical equals the
intersection of the ten explicit point ideals.  The local tangent cones give

```text
special centers: discriminant histogram 2:2, 3:2  -> four A1 nodes;
ordinary centers: one tangent and nonzero cubic tangent term
                  with histogram 1:1, 2:3, 3:1, 4:1 -> six A2 cusps.
```

Geometric integrality follows by the same Frobenius--Bézout argument used in
the Type-VI inverse.  If the irreducible sextic split geometrically into
`r>=2` equal-degree components, every rational singular center would lie on
all components.  Two components would share ten points, but their degree is
`6/r<=3`, so Bézout permits at most nine.

A plane sextic has arithmetic genus ten.  Both A1 and A2 have delta invariant
one, so the ten singularities have total delta ten.  One blowup at every
center makes the strict transform smooth and rational.  Its class is

\[
             6H-2\sum E_i=-2K,\qquad \widetilde F^2=-4.
\]

The blowup is therefore a strong Coble surface constructed directly from the
raw design.

## 4. Boundary/root transversality

The audit separates smooth-root incidences from quartic-node incidences.

At the 84 centers where a root is smooth, its tangent direction evaluates
nontrivially on the boundary tangent cone:

```text
lines       value 1:18, value 4:18
conics      value 1:12, value 4:12
quartics    value 1:12, value 4:12.
```

At the twelve special centers where a quartic is nodal, the resultant of the
quartic and boundary binary tangent cones is nonzero:

```text
resultant 1:6, resultant 4:6.
```

Thus no root branch shares a boundary branch.  Local intersections attain
their multiplicity products, and Bézout is exhausted:

```text
line       3*2 =  6 = 1*6
conic      6*2 = 12 = 2*6
quartic    6*2 + 3*4 = 24 = 4*6.
```

There are no residual intersections away from the centers.  All twenty root
strict transforms are smooth rational `(-2)`-curves disjoint from the
boundary.

## 5. Exact marking for every raw frame

For each of the five raw frames, the colored incidence graph records point,
line, conic and quartic vertices and distinguishes simple from double center
incidences.  Each frame has exactly 24 isomorphisms to the explicit finite
geometry—the stabilizer expected for a tetrahedral frame.  The audit tests all

\[
                             5\cdot24=120
\]

markings.  Under every marking, every one of the twenty raw root vectors is
equal coefficient-for-coefficient to its geometric divisor class.  This is
strictly stronger than agreeing on the `12+4+4` degree profile.

The same determinant-one bases transport every one of the 72 raw wall
vectors integrally.  The per-frame full-wall coordinate digests are

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

## 6. Root-cone completeness and the same P

To avoid smuggling surface existence through the old Type-VII label, the
twenty raw root inequalities are independently sent to exact rational
cddlib.  Their rank-ten Coble--Mukai cone has

```text
450 finite rays + 51 ideal rays + 0 hyperideal rays = 501 rays.
```

The sum of its primitive rays is a strict timelike interior vector: its square
is `26,564,430`, its minimum root slack is `1,882`, and its minimum pairing
with any ray is `4,705`.  Thus the twenty-root chamber is finite volume on one
Lorentz time sheet.

Set `f=E+b` for each of the 51 remaining norm-minus-one raw walls.  Every
such `f` is primitive, isotropic, and lies in `b`-perpendicular.  Exact set
comparison—not merely a count—gives

```text
{51 cdd ideal rays} = {f=E+b for the 51 raw walls}.
```

The complete primitive-ray and cusp digests are

```text
all 501  80e1e81adc197bb92e818f52f3e695aa1d4d61fc2f4aecf6a550ea1b99fbb5ec
51 cusps 8621adee95e09a0213356e044f923a0fd902e22aa55ead6ffc518959fbfc73d7.
```

Now apply only general Coble geometry, not a Type-VII existence theorem.  The
twenty constructed effective roots generate a finite-index reflection
subgroup.  The general Coble finite-index and root-exhaustion propositions
therefore exclude hidden effective irreducible roots.  For a one-boundary
strong Coble surface, every primitive nef cusp gives the `(-1)`-curve

\[
                         E=f-b=f+K,
\]

and every `(-1)`-curve returns `f=E+b`.  Together with the general negative
curve trichotomy this yields the complete wall stock

```text
1 boundary + 20 (-2)-curves + 51 (-1)-curves = 72 walls.
```

Under the exact integral marking this is the frozen raw wall set.  Its
independent P-side Vinberg certificate supplies

```text
dimension                 10
facets                    72
finite vertices           18,761
ideal vertices             1,731
hyperideal extreme rays        0
elliptic edges            143,565
volume status             finite, noncompact
Vinberg failures                0.
```

Hence the construction returns the same marked `H10/72` polytope.

## 7. Imported versus new

Imported:

- the frozen raw 72-wall lattice and its already-independent full-cone and
  literal Vinberg census;
- general theorems about Coble reflection groups, root exhaustion, negative
  curves, and the one-boundary cusp lift.

Constructed or recomputed here:

- all five intrinsic contraction frames and determinant-one plane markings;
- the `F5` conic, the five usable factorization frames, and all ten centers;
- the twelve lines, four smooth conics, four geometrically integral nodal
  quartics, and the geometrically integral nodal/cuspidal sextic;
- `|-K|=empty`, uniqueness of `|-2K|`, every local singularity, every tangent
  gate, and the strong-Coble conclusion;
- all 120 exact root-vector markings, all five integral 72-wall transports,
  the direct `450+51+0` root cone, and exact identification of all 51 cusps.

Kondō's known characteristic-five Type-VII surface can be attached as an
after-the-fact classification label.  His printed construction uses a
`P4` Cremona quotient and does not print this plane certificate.  This local
source comparison does not establish worldwide priority for the equations.

The result is therefore:

> a project-new strict P-first explicit realization of a known positive
> `H10/72` datum, at inverse level L3.

It is not L4 because it treats one rigid datum rather than a variable class of
arithmetic Coxeter polytopes with necessary-and-sufficient surface-carrier
criteria.  It nevertheless supplies a second complete P-first calibration,
now with nodal quartics and a mixed node/cusp boundary, against which such a
class-wide theorem can be designed.

## Reproducibility

- [exact certificate](../data/enumerations/coble_type_vii_char5_tetrahedral_inverse.json)
- [certifier](../scripts/certify_coble_type_vii_char5_tetrahedral_inverse.py)
- [regression tests](../tests/test_coble_type_vii_char5_tetrahedral_inverse.py)
- [raw finite-design archive](../data/enumerations/coble_remaining_seven_contraction_designs.json)
- [raw `H10/72` archive](../data/enumerations/coble_type_vii_char5_one_boundary_full_ns_h10.json)

The certifier calls Singular for exact factor and singular schemes and
`cddexec_gmp` for an independent exact root-cone conversion.  Six regression
tests include byte-for-byte archive replay.
