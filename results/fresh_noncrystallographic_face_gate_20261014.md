# A noncrystallographic ridge requires a large compact Coxeter face

An actual `pi/5` ridge in a finite-sided, finite-volume Coxeter H20
polytope would itself be a **compact Coxeter H18 polytope**. More generally,
an actual finite ridge label outside `{2,3,4,6}` in dimension `n` forces
a compact Coxeter face of dimension `n-2`. This is a consequence of
published face theorems and the cusp condition, not a claimed new
dimension theorem or a nonexistence proof for all H20 polytopes.

This changes the search priority. The weighted golden H20 seed does not
merely ask for an unknown noncompact H20 example: retaining its central
`pi/5` ridges also asks for compact H18 examples. The current H18-hybrid
lift, with only finite angles `pi/2`, `pi/3`, and `pi/4`, avoids this
particular requirement. Irrational divergent-wall distances can supply
nonarithmeticity without introducing these noncrystallographic ridges.

## Two separate ingredients

Borcherds' face theorem says that a spherical subdiagram with no
component of type `A_n` or `D5` defines a Coxeter face. Apply it to the
rank-two diagram `I2(m)`, where `m` is outside `{2,3,4,6}`. The resulting
face has dimension `n-2`. Its induced angles need not equal the original
angles; an `H4` extension of `I2(5)` can, for example, yield `pi/10`.
See [Allcock, Theorems 2.1 and 2.2](https://arxiv.org/pdf/0903.0138),
which give the geometric statement and proof, and
[Borcherds, Example 5.6](https://math.berkeley.edu/~reb/papers/coxeter/coxeter.pdf)
for the original underlying result.

The face cannot contain an ideal vertex. At such a vertex the two
reflections defining the ridge would act on a bounded Euclidean cusp
link. In its `d=n-1` dimensional linear representation their product
has trace `d-2+2cos(2pi/m)`. A cocompact Euclidean reflection group
preserves a full translation lattice, so that trace must be an integer.
For integer `m>=2`, integrality occurs exactly for `{2,3,4,6}`. Thus
every vertex of the ridge is finite. A face of a finite-sided
finite-volume polytope has only finite or ideal vertices, so this face
is compact. The [preceding cusp argument](fresh_h20_weighted_binary_20261014_dotted_completion_review.md)
gives the trace proof and its Euclidean-lattice reference in detail.

Combining this with Vinberg's compact dimension bound gives the general
necessary condition `n<=31` whenever such a ridge is retained. This
does not improve the unrestricted finite-volume dimension bound: a
polytope may have only crystallographic finite labels. More immediately,
an H20 `pi/5` route entails compact H18, and an H22 `pi/7` route entails
compact H20. Those compact dimensions remain unknown in the literature
checked; they are not proved impossible by this argument.

## The original weighted ridge has twenty divergent facets

Use physical spatial coordinates and unit vertical normals. Let `a` be
the central normal of square four, and `b` the positive golden seed of
square one. Their Gram matrix and determinant are

```text
S = [[4,-phi],[-phi,1]],       det(S)=k=3-phi.
```

For either side of a unit coordinate, restriction to `S-perp` gives
square `1-1/k>0`. For a heavy side the square is `1-2/k<0`; for another
golden side it is `1-phi^2/k<0`; the opposite target-golden side also
has negative square. The latter walls have no hyperbolic boundary on
the ridge. The twenty unit sides are actual facets there, since their
three-root extensions with S are positive definite and have strict
projected interior points.

The normalized negative pairings between their restricted normals are
`phi^2` for different coordinates and `3+2phi` for opposite sides.
Both exceed one: all twenty facets are mutually divergent. The
[exact checker](../scripts/fresh_noncrystallographic_face_gate_20261014.py)
and [certificate](fresh_noncrystallographic_face_gate_20261014.json)
verify these projections and the actual-ridge/absent-core-vertex signs
over `Q(phi)`. The cited general theorems are mathematical dependencies,
not conclusions of this finite arithmetic calculation.

## A genuine whole-family stopping rule

Suppose a proposed finite-volume completion retains this core and the
chosen ridge. Assume every additional wall that bounds the ridge is
orthogonal to both a and b, and that all its finite contacts with other
such walls and with the twenty old unit walls are right. These include
the small orthogonal-cap constructions and additions whose relevant
relations remain right, parallel, or dotted.

Restriction leaves each new normal unchanged. Its old right contacts
remain right. Its old parallel or dotted contacts with a unit seed
become divergent because that seed's projected norm decreases.
Consequently the induced H18 face would be right-angled. The first
argument makes it compact. Compact right-angled hyperbolic polytopes
have dimension at most four, a contradiction. The bound and its proof
appear in [Potyagailo–Vinberg, §2](https://ems.press/content/serial-article-files/43032).

This excludes arbitrary finite additions within the stated right-angled
face class, regardless of local vertex counts. It does not exclude
nonright contacts, `H3/H4` attachments, or replacement of the golden
ridge. A finite edge closure in this class can test a geometric
edgewalking implementation, but iterating it cannot solve the volume
problem. Further H20 work should prioritize a seed without this compact
H18 requirement, or explicitly undertake that harder compact-face problem.

```sh
python3 scripts/fresh_noncrystallographic_face_gate_20261014.py --verify
```
