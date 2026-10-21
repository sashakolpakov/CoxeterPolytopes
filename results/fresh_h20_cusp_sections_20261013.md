# H20 sections through the complete cusp census

The new construction preserves an entire Coxeter cusp, but does not produce
a Coxeter polytope. An exact all-parameter obstruction excludes **59 distinct
normal planes**, represented by **118 oriented cusp pencils**. These come
from equal-length orthogonal twin walls in 14 of the 18 recorded ideal-vertex
representatives of the Borcherds H21 chamber. Every interior section in these
families has finite volume and a complete Euclidean Coxeter cusp link; an
explicit actual ridge always has its angle strictly between pi/3 and pi/2.

This is a different construction from the previous adjacent-bridge cusp
cuts. It merges two nonadjacent cusp walls and removes a globally redundant
inequality. The complete ideal census supplies new directions, including
60 oriented pencils at the newly recorded `D4^5` cusp, orbit 58. The result
does not exclude general sections through a cusp, mixed-component directions,
or a Coxeter polytope in H20.

## The construction and its full parameter interval

Write Q with spacelike sign positive. Let delta be one of the recorded
future null cusp vectors. Choose incident short walls a,b satisfying

```text
r_a^2=r_b^2=2,   (r_a,r_b)=0,
(r_a,r_i)=(r_b,r_i) for every other wall incident to delta.
```

All such equal-length orthogonal twins in this literal 18-cusp census are
short. Set

```text
e0=(r_a-r_b)/2,       e(t)=e0+t*delta.
e(t)^2=1,            (e(t),delta)=0.
```

The old incident-wall pairings with e(t) are constant: 1 at a, -1 at b,
and zero elsewhere. Orthogonal projection to the section gives

```text
Pi_t(r_a)-Pi_t(r_b)=-2t*delta.
```

For `t>0` and a future timelike point, the a inequality is strictly stronger
than the b inequality. Thus b is globally redundant in the hyperbolic
section, not merely absent from the formal cusp diagram. Removing b leaves
the old incident Gram with the norm of a changed from two to one. Every
remaining pair is exactly right, label 3, label 4, label 6, or parallel.

Let the component containing a,b have positive kernel marks k_i, with
`k_a=k_b=m`. The frozen primitive cusp conventions give

```text
sum k_i*r_i = m*delta.
```

After projection and removal of b, double the mark at a. The positive
relation becomes

```text
2m*Pi_t(r_a) + sum_{i != a,b} k_i*Pi_t(r_i)
    = m*(1-2t)*delta.
```

The other affine components keep their positive multiples of delta.
All projected component Grams are connected, have nonpositive off-diagonal
entries, and have a positive kernel. Congruence by the kernel diagonal
turns each into a weighted graph Laplacian, proving positive semidefiniteness
and rank one less than its number of walls. Their ranks sum to 19. For
`0<t<1/2` the positive null multiples therefore give a bounded rank-19
Euclidean Coxeter link.

The interval is exact. At `t=1/2`, the normal is the positive combination

```text
e(1/2)=r_a+sum_{i != a,b} k_i*r_i/(2m).
```

It cannot vanish at a strict chamber interior point. For larger t, the
additional multiple of delta is another nonnegative wall combination.
Reversing a,b covers `t<0`. Hence a genuine interior section of any of
these normal planes has exactly `|t|<1/2`, up to overall scaling. The
remaining projective normal, delta itself, is null and does not define a
hyperbolic H20 section. At zero, the projected a,b normals coincide, so
one duplicate is removed.

Examples of the local change are `D4^5` to `B3 x D4^4` at cusp 58, and
the rank-20 fork component at cusp 14 to the affine `C19` chain. These
names describe the verified local Coxeter diagrams, not new global examples.

## An actual forbidden ridge for every parameter

For each oriented pencil, the checker finds a short wall j attached only
to one of the twins within their original affine component. In 112 cases
it is attached only to the omitted leaf b; in the other six, only to the
retained leaf a. It is nonincident to the cusp and satisfies
`(delta,r_j)=-1`.

Two exact ambient face points certify incidence throughout the interval:

* z0 lies on exactly `{a,b,j}`.
* zT lies on exactly the original component with b deleted, together with j.

Both face Gram matrices are positive definite. The checker projects the
strict ambient center onto each face, verifies negative norm and future
orientation, and checks all 210 inequalities, with precisely the stated
zero sets. Their zero sets intersect in `{a,j}`. Put

```text
A=-(delta,z0)>0,       B=-(delta,zT)>0,
x(t)=(1/2-t)*B*z0+t*A*zT.
```

For every `0<t<1/2`, x(t) is future timelike, belongs to `e(t)^perp`, and
vanishes on exactly a,j. Every other ambient inequality is strict. The
two projected normals have positive-definite Gram matrix, so this is a
genuine dihedral ridge. A small inward perturbation also supplies a strict
section interior. This avoids inferring a ridge from a raw projected pair
alone.

The retained normal a has projected square one. In the two attachment
cases, the ridge's squared cosine is respectively

```text
(1/2-t)^2 / [2-(1/2-t)^2],
(1/2-t)^2 / [2-(1/2+t)^2].
```

Both are positive. Subtracting four times the numerator from the denominator
gives, respectively,

```text
3/4+5t-5t^2,          3/4+3t-5t^2.
```

These concave quadratics are strictly positive at both ends of `[0,1/2]`,
and hence throughout. Thus `0<cos(theta)^2<1/4`: no integer Coxeter label
is possible. At `t=0`, z0 is still an actual ridge point after identifying
the duplicate a,b normals, and its squared cosine is `1/7`. The reversed
leaf order covers every negative parameter. The collapsed interval ends
do not give interior sections and introduce no missing Coxeter exception.

## Scope, finite volume and replay

The full-rank positive cusp relations prove that delta remains an actual
ideal vertex with bounded link. Finite volume of the section follows from
the known finite-volume H21 chamber: in the Klein model, a hyperplane
section of its finite/ideal vertex polytope has only finite vertices and
inherited ideal vertices. A finite triangulation therefore has finite
hyperbolic volume. No arithmetic-reflectivity assumption is used.

The census is complete for equal-length orthogonal twins among the
incident walls of the **18 saved literal cusp representatives**. Its 118
oriented rows give 59 distinct rational two-dimensional normal planes;
both orientations and every real parameter are covered. The same proof
transports under an ambient chamber symmetry. This is not a classification
of all Coxeter sections of the Euclidean links, all cusp-adapted normal
planes, or all H20 sections. A constructive continuation must change the
local direction, rather than only translating one of these leaf-merging
planes along its cusp.

* [Exact checker](../scripts/fresh_h20_cusp_sections_20261013.py) and
  [exact certificate](fresh_h20_cusp_sections_20261013.json).
* [Numerical proposal scout](../scripts/fresh_h20_cusp_sections_20261013_scout.py)
  and [sample data](fresh_h20_cusp_sections_20261013_scout.json), used only
  to explore the mechanism. Its missing ridge detections are not exclusions.
* [Complete ambient corner census](fresh_h21_source_corner_walk_20261013.json).
* [Independent source review and Fraction ridge replay](fresh_h20_cusp_sections_20261013_review.md),
  with [its checker](../scripts/fresh_h20_cusp_sections_20261013_review.py)
  and [PASS record](fresh_h20_cusp_sections_20261013_review.json).

```sh
python3 scripts/fresh_h20_cusp_sections_20261013.py --verify
```

The producer and saved replay both pass. The checker uses only standard
library integers and rational arithmetic, imports no numerical proposal
code, and reconstructs the twin census, local links, positive identities,
both actual face points, and angle inequalities from the frozen source.
The separate review reconstructs all 118 archived ridge identities, face
signs, local links and zero-parameter duplicates with its own Fraction
operations; it also reran the primary checker. Its detailed report states
the distinction between this independent certificate replay and the
source-level review of continuous geometric implications.
