# Removing invisible polar cuts from the H9 construction

**Outcome (2026-10-13): polar redundancy is real, but it does not rescue a
compact Coxeter candidate in the previously searched family.** An exact
geometric criterion identifies every retained polar facet for the simple,
nonideal vertices under consideration. Replaying the same numerical
templates with that reduction removes 1,036,376 proposed walls. All 67,344
noncusped skeletons still fail a pairing between a retained polar facet and
an original facet. This is a bounded numerical result, not a theorem that
compact H9 Coxeter polytopes do not exist.

This report resolves the polar-redundancy limitation in the
[previous Gram construction](fresh_bugaenko_gram_extension_20261013.md).
The original script and transcript are preserved; this experiment does
not add labels, templates, or enumeration height.

## Exact retention criterion

Let the original normals span Lorentz space of signature `(n,1)`, have
positive square, and satisfy `(r_i,r_j)<=0` for different indices. Assume
their Gram diagram is connected and that its negative eigenspace has
dimension one. Consider the polyhedron obtained by adding the polar
inequality of every hyperideal vertex of the original cone. For a simple
hyperideal vertex `p`, the proposed polar wall is a hyperbolic facet of
this intersection **if and only if the Gram diagram of the old walls
incident to `p` is connected**.

Here a hyperbolic facet must have a relative-interior point inside the
open Klein ball. A Euclidean face outside or on the sphere does not count.
The criterion is not a claim about nonsimple or ideal vertices.

### Centered Klein model

The smallest Gram eigenvalue has a strictly positive eigenvector `w`.
The vector `h=sum w_i*r_i` is timelike and satisfies all old inequalities
strictly. For every proposed polar normal `p`, the cone inequalities give
`(p,r_i)<=0`, with at least one strict inequality, so `(h,p)<0` as well.
Thus the old and truncated intersections have a common strict timelike
interior point.

Choose the centered Klein model with this point at the origin. Rescale
the old normals positively to `r_i=(1,u_i)`. The old Klein polytope is

\[
 P=\{x:u_i\cdot x\leq1\text{ for every }i\}.
\]

It is bounded and contains the origin in its interior: the positive
combination giving `h` yields a positive relation among the `u_i`, which
span the spatial space. Its Euclidean polar is
`P^\circ=conv{u_i}`. A vertex `p` of `P` contributes the polar inequality
`p.x<=1`. When `|p|<=1`, this inequality is automatically strict for
`|x|<1`. Consequently, inside the open ball, imposing all the hyperideal
polar cuts is exactly the intersection `P intersect P^circ`.

For a simple hyperideal vertex `p`, its incident set `I` contains `n`
old walls. Their Gram matrix is nonsingular of signature `(n-1,1)`, since
their normals span the orthogonal complement of the spacelike vector
`(1,p)`. The corresponding face of the Euclidean polar is the simplex

\[
 F_p=\operatorname{conv}\{u_i:i\in I\}.
\]

### Disconnected incident diagram: the cut is invisible

The incident Gram decomposes into one Lorentzian block and one or more
positive-definite blocks. If a point `x` belongs to `F_p` and satisfies
all the old inequalities, write
`x=sum a_i*u_i`, where `a_i>=0` and `sum a_i=1`. Its incident pairings are

\[
 \widetilde G_I a\leq0,
 \qquad (\widetilde G_I)_{ij}=u_i\cdot u_j-1.
\]

For each positive-definite block `E`, nonnegative coefficients imply
`a_E^t G_E a_E<=0`, hence `a_E=0`. Every feasible point on this polar wall
therefore lies in a proper face of `F_p`, of dimension at most `n-2`.
The wall has no relative-interior hyperbolic facet point.

### Connected incident diagram: the cut is visible

The negative Perron eigenvector of the incident Gram is strictly
positive. Its combination `q=sum b_i*r_i`, with `i in I`, is future
timelike, lies on the polar plane of `p`, and is strict for each incident
old inequality. Every other old pairing is nonpositive. For any other
old cone vertex `v`, at least one pairing `(r_i,v)`, `i in I`, is strictly
negative: otherwise `v` and `p` would be the same ray. Thus `(q,v)<0`, so
`q` is strict for every other polar cut.

To make any remaining old inequalities strict, use the projection

\[
 h_p=h-\frac{(h,p)}{(p,p)}p.
\]

This is timelike, lies on the same polar plane, and satisfies all old
inequalities strictly because `(h,p)<0` and `(p,r_j)<=0`. A sufficiently
small perturbation of `q` toward `h_p` keeps all other polar inequalities
strict. It gives the required relative-interior hyperbolic facet point.

### Simultaneous removal and the original walls

All distinct invisible cuts may be removed simultaneously inside the
open ball. If their removal admitted a new point, take a generic segment
from the common strict interior point to such a point. At its first
crossing of the removed hyperplanes, exactly one is crossed and all
other inequalities are strict. That would give a facet point for a wall
proved invisible, a contradiction. The hyperplanes are distinct because
the original cone rays are distinct.

Every old wall remains a facet. For a square-two old root `r_i`, the
projection `h_i=h-(h,r_i)*r_i/2` is timelike and lies on that wall. For
every other old or polar normal `s`,

\[
 (h_i,s)=(h,s)-\frac{(h,r_i)(r_i,s)}2\leq(h,s)<0.
\]

Thus neither an invisible polar cut nor a redundant old wall can now
cause a false rejection at the polar/old angle gate, subject to the
numerical scan correctly identifying its stated nondegenerate input.

## An exact example where redundancy occurs

Take chain labels `[4,3,3,3,3,3,3,3]`, attach nodes 9 and 10 to old
positions 0 and 5 with label 6, and set their dotted weight to
`t=sqrt(115)-6*sqrt(2)>2`. The nine-node chain has all leading principal
minors equal to 2. The remaining Schur complement is

\[
 \begin{pmatrix}-23/2&-\sqrt{115}\\-\sqrt{115}&-10\end{pmatrix},
\]

which has rank one and is negative semidefinite. The full Gram therefore
has inertia `(9,1,1)`. An exact kernel vector has signs
`[+,+,+,-,-,-,-,-,-,+,-]`. Omitting old walls 0 and 3 gives a feasible
simple cone vertex because those kernel signs are opposite.

Its incident diagram splits into the elliptic component `{1,2}` and
the Lorentzian component `{4,5,6,7,8,9,10}`. The latter consists of a
positive A5 chain with remaining Schur complement
`[[2,-t],[-t,-2]]`. Thus the incident Gram has inertia `(8,1)` and the
vertex is hyperideal. Its polar cut is invisible by the criterion above.
The [exact checker](../scripts/polar_reduction_exact_example_20261013.py)
verifies these algebraic facts. This is an example of redundancy, not
a compact Coxeter polytope.

An independent numerical SOCP checks every hyperideal polar cut of the
same template. For each cut it maximizes a common strict margin for the
other inequalities and for the unit ball while imposing equality on the
cut under test. All 24 checks agree with the graph criterion: 5 visible
facets and 19 invisible cuts. Visible margins range from `0.0167` to
`0.0631`; invisible optima are within `6e-11` of zero. This numerical
cross-check supplements the proof; it is not the proof.

## Reduced rescreen

The [new scanner](../scripts/polar_reduction_h9_20261013.py) imports the
previous template generator, replaces only its polar analysis, and
caches connectivity by the old diagram topology. It considers the same
`45*4^10 = 47,185,920` templates. Labels on the original chain and its
attachments remain in `{3,4,5,6}`; the subsequent polar-angle check has no
imposed maximum Coxeter label.

| stage | count |
|---|---:|
| rank-ten Lorentz skeletons | 113,099 |
| parked: zero kernel entry | 2,465 |
| parked: numerically ideal ray | 43,290 |
| noncusped skeletons receiving polar reduction | 67,344 |
| invisible polar cuts removed, summed over these skeletons | 1,036,376 |
| visible polar facets retained, summed over these skeletons | 268,980 |
| reduced skeletons rejected by retained polar/old angle | 67,344 |
| candidates requiring exact compactness verification | 0 |

Each of the 67,344 evaluated skeletons loses at least one invisible cut;
each still has a retained polar wall. None reaches the polar/polar gate.
No new compact example or new-dimensional example was found.

The retention theorem and exact redundancy example are exact results.
The broad rescreen still uses floating-point rank, signature, vertex,
and angle decisions with tolerance `1e-7`. Singular bases, degenerate
Schur complements, zero kernel entries, and cusp-bearing cases retain
their previous parked status. Other graph topologies, finite labels,
nonpolar cuts, and nonsimple constructions remain open. The scan is not
an exact no-go certificate even for its bounded template family.

Replay and validation:

```sh
python3 scripts/polar_reduction_exact_example_20261013.py
OPENBLAS_NUM_THREADS=1 python3 scripts/polar_reduction_socp_check_20261013.py
OPENBLAS_NUM_THREADS=1 python3 scripts/polar_reduction_h9_20261013.py
python3 -m py_compile scripts/polar_reduction_exact_example_20261013.py scripts/polar_reduction_socp_check_20261013.py scripts/polar_reduction_h9_20261013.py
```

Transcripts:

- [Exact redundancy example](../data/enumerations/polar_reduction_exact_example_20261013.json)
- [Independent SOCP check](../data/enumerations/polar_reduction_socp_check_20261013.json)
- [Reduced numerical rescreen](../data/enumerations/polar_reduction_h9_20261013.json)
