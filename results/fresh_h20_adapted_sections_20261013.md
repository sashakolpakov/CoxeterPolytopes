# Simultaneous cusp-adapted repairs after ambient doubling

**No Coxeter H20 section survives in the 32 continuous families tested here.**
This result covers every real parameter in four specified cusp pencils and
every subset of three specified, mutually orthogonal additional seams. It
does not exclude other cusp directions or other chamber assemblies.

The constructive attempt first repairs a persistent ridge to pi/3 at four
new exact quadratic section normals. Each section has finite volume, a
complete rank-19 Coxeter cusp, and a separate actual obtuse ridge. The
simultaneous local analysis then excludes the entire declared parameter
families, rather than only those four samples. Ambient unions contain up
to 64 Borcherds chambers, with up to 5,577 true ambient facets. These are
ambient-facet counts, not a census of section facets.

[Exact construction and data](../scripts/fresh_h20_adapted_sections_20261013.py)
are preserved in the
[certificate](../data/enumerations/fresh_h20_adapted_sections_20261013.json).
The [independent Fraction replay](../scripts/fresh_h20_adapted_sections_verify_20261013.py)
and its [result](../data/enumerations/fresh_h20_adapted_sections_verify_20261013.json)
check the ambient unions, cusp relations, continuous ridge proofs, unique
angle parameters, and actual residual ridges. The checker uses only the
standard library and earlier independent verification helpers; it imports
no producer, NumPy, SymPy, or numerical search code.

## The simultaneous construction

Start with the previously certified eight-chamber union obtained by
successive ambient doublings across original walls 169, 194, and 71.
The first two mirrors are parallel and generate an infinite subgroup.
Thus these ambient unions are outside the earlier finite spherical-star
class.

Use the four changed cusp pencils from the
[preceding construction](fresh_h20_multidouble_sections_20261013.md), with
normal `e(t)=e0+t*u`, where `e0^2=4`, `u^2=(e0,u)=0`, and `u` is a cusp
vector. All four satisfy `(e,r1)=2`. The five short neighbors of `r1` are
original walls `22,23,24,49,50`. They have norm two, pair `-1` with `r1`,
and are mutually orthogonal.

Two of those neighbors are cusp-incident. The three nonincident neighbors
provide a simultaneous repair: double across any subset, giving eight
choices for each pencil. The available seams are:

| Cusp bridge defining the pencil | Additional seams | Parameter interval admitting an interior |
|---|---|---|
| 1/22 | 23, 49, 50 | `(-7/15,8/15)` |
| 1/24 | 23, 49, 50 | `(-7/15,8/15)` |
| 1/49 | 22, 23, 24 | `(-4/9,5/9)` |
| 1/50 | 22, 23, 24 | `(-4/9,5/9)` |

Every selected seam is a current true facet and has nonpositive pairing
with every other current boundary normal. Exact reconstruction therefore
certifies each doubled union as convex: a reflected normal is a
nonnegative combination of its predecessor and the seam normal, and the
seam sign identifies the two halves. Each resulting boundary plane
contains a full facet of a preceding half. Exact duplicate removal gives
the listed true ambient facets. Every doubling doubles the number of
ambient chambers.

All seams are strictly nonincident to the corresponding cusp. The exact
set of incident boundary normals stays unchanged. Their projected Gram
has the same complete rank-19 Euclidean Coxeter link throughout each
open interval. The interval endpoints are forced by the positive affine
kernel relations: if a component's null multiple `mu(t)` is nonpositive,
a putative strict future timelike section point would pair negatively
with a positive sum of its wall normals, but nonnegatively with
`mu(t)*u`. Thus parameters outside the stated interval, including its
endpoints, cannot give a full-dimensional section of these unions.

## An actual ridge for every parameter

Put `L` for the lower endpoint and `x=(e(t),r23)`. In every pencil,

```
x=-2(t-L),             L<t<L+1,             -2<x<0.
```

The certificate contains two rational future timelike points `pL,pR` in
the original face `r1=r23=0`. They lie in the sections at parameters `L`
and `L+1`, respectively. Their other zero-wall sets have intersection
exactly `{1,23}`. With `A=-(pL,u)>0` and `B=-(pR,u)>0`, the vector

```
p(t)=(L+1-t)*B*pL + (t-L)*A*pR
```

is future timelike, lies in `e(t)^perp`, and is strict on every other
original wall for every interior parameter. The checker additionally
evaluates both endpoint points against every boundary normal in every
declared union. Their common zero set is exactly the displayed pair of
boundary normals, so this is a genuine H18 dihedral ridge of the H20
section. A positive-definite two-normal Gram also permits a small inward
perturbation, proving full dimension throughout the interval.

This verification matters because the ambient unions and the sections
are not assumed acute. A positive-definite raw pair alone would not prove
the existence of a ridge.

If wall 23 was **retained**, the ridge remains between `r1` and `r23`.
Its projected Gram, multiplied by four, is

```
[[4,       -4-2x],
 [-4-2x,   8-x^2]].
```

Its squared angle cosine is `(x+2)^2/(8-x^2)`, strictly between zero and
one half. Consequently the only possible Coxeter angle is pi/3. The
necessary equation is

```
5x^2+16x+8=0.
```

The independent Sturm count finds exactly one corresponding parameter in
the admissible interval.

If wall 23 was **doubled**, the two surviving normals are `r1` and
`s23(r1)=r1+r23`. Their scaled projected Gram is

```
[[4,     -2x],
 [-2x,   8-(2+x)^2]].
```

The off-diagonal entry is positive throughout `-2<x<0`. The determinant,
`16-16x-8x^2=16+8*(-x)*(2+x)`, is positive. The same actual ridge point
therefore proves an obtuse angle for **every** admissible parameter.
Doubling other selected neighbors does not remove this certified ridge.

## The only remaining angle candidates

For the assemblies retaining wall 23, the exact pi/3 equation leaves the
following four section normals, all in `Q(sqrt(6))`:

| Cusp bridge | Unique parameter | Separate actual obtuse ridge, in original wall labels |
|---|---|---|
| 1/22 | `(5-3*sqrt(6))/15` | 1/45 |
| 1/24 | `(5-3*sqrt(6))/15` | 1/11 |
| 1/49 | `(16-9*sqrt(6))/45` | 1/3 |
| 1/50 | `(16-9*sqrt(6))/45` | 1/5 |

Every row has an exact strict timelike section interior and an actual
pi/3 ridge between original walls 1 and 23. The saved residual ridge point
is future timelike, lies on exactly its two listed walls, and is strictly
inside all the others. The residual pair has positive projected pairing.
Its two walls are never seams in the declared assemblies. Direct exact
evaluation against every reconstructed ambient boundary confirms that
this actual obtuse ridge survives all eight additional-seam choices.

The sections have finite volume: a full-dimensional hyperplane section
of the finite-vertex Klein polytope has only finite or inherited ideal
vertices. Intersections created in open ambient edges lie strictly inside
the Klein ball. The resulting H20 ideal polytope therefore has finite
volume; its certified cusp makes it noncompact. This geometric inheritance
does not turn the forbidden angles into Coxeter angles.

## Replay and limits

```sh
python scripts/fresh_h20_adapted_sections_20261013.py
python scripts/fresh_h20_adapted_sections_verify_20261013.py
```

The producer reconstructs sixteen declared ambient states, yielding fourteen
distinct boundary lists because the two seam families overlap. This is not
a count up to isometry. The independent replay passed: it rebuilds every
state, checks its complete boundary digest, validates all 32 continuous
section families, and rejects a forged endpoint and a forged actual-ridge
point. The quadratic angle equation is counted independently using
Fraction Sturm chains, with seven small Sturm regressions also passing.
An additional agent reviewed the endpoint interpolation, actual-ridge gates,
and angle interval argument at source level and found no issue; that review
was not another producer run.

The obstruction is local and scoped: retaining wall 23 fixes one angle
parameter but leaves another actual obtuse ridge; doubling it creates an
actual obtuse ridge throughout the whole admissible interval. It excludes
the four stated pencils in the declared simultaneous-neighbor assemblies.
It does not exclude further seams that remove the residual ridges, changes
to the cusp direction, nonreflection gluings, or a Coxeter polytope in H20.
