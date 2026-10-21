# Dotted-only additions cannot complete the binary-cap chamber

**Conclusion.** Retain the frozen weighted binary chamber, and fix any one of its binary facets. No finite, full-dimensional acute-compatible Coxeter completion can have finite volume if every added wall is dotted to that facet. Thus a proposed family dotted to all binary caps cannot repair the chamber, regardless of its size, labels, or interactions among the added walls.

The hypotheses are retained facets and acute-compatible additions. This does not exclude deleting or replacing the chosen binary facet, changing its old neighbors, or adding walls which meet it. The conclusion is stronger than saying the current facet lacks finite vertices, but remains a scoped completion obstruction.

## A finite-volume Coxeter ridge at angle pi/5 has finite vertices

Let `P` be a finite-sided finite-volume hyperbolic Coxeter polytope of dimension at least three, and let `F` be an actual ridge with dihedral angle `pi/5`. In the Klein model the projective closure is a polytope whose vertices are finite or ideal. The nonempty face `F` therefore has a vertex of one of these two kinds.

An ideal vertex on `F` is impossible. A sufficiently small horospherical section at such a vertex is a bounded Euclidean Coxeter polytope. The two reflections in the facets defining `F` induce Euclidean reflections meeting at `pi/5`. On the `d=n-1` dimensional horosphere their product has linear trace

```text
(d-2)+2cos(2pi/5) = d-3+phi,
```

which is irrational. On the other hand, a cocompact discrete Euclidean isometry group has a full translation lattice by Bieberbach's theorem. The full translation subgroup is normal, so every linear part preserves this lattice and has integer trace. This is a contradiction. See [Gelander, *Lectures on Lattices*, Theorem 2.9, PDF p.5](https://www.math.utah.edu/pcmi12/lecture_notes/gelander.pdf). Equivalently, the classification of bounded Euclidean Coxeter chambers as products of affine simplices has no label five: [Felikson, §1.2 and Table 1.1, PDF pp.10–12](https://www.maths.dur.ac.uk/users/anna.felikson/talks/cl_cox.pdf).

Thus every vertex of `F` is finite, and at least one exists. No assertion that hyperbolic edges have finite length, or that a facet automatically inherits Coxeter angles, is used. The ridge itself is a face of the original Coxeter polytope; the horospherical angle is its original dihedral angle.

The same argument applies to every finite ridge label `m` outside `{2,3,4,6}`. The relevant trace is `d-2+2cos(2pi/m)`. If the last term were rational, its algebraic integrality would make it an integer in `[-2,2]`; for integer `m>=2`, this gives precisely `{2,3,4,6}`. Consequently a retained facet with no finite vertices and an actual ridge of any other finite label requires at least one added wall meeting that facet in every finite-volume completion. Parallel or divergent pairs are not finite-angle ridges and are not covered by this remark.

## Dotted walls leave the retained facet unchanged

Use outward normals, with chamber inequalities `(x,r)<=0`. Let `a` be the chosen retained facet normal, `q_a=a²>0`, and let `r` be any added normal dotted to it:

```text
(a,r)<-sqrt(q_a*r²).
```

Let `c` be a strict timelike interior point of the proposed completion, not necessarily the original chamber's center. Its projection

```text
c_F = c - (c,a)*a/q_a
```

is timelike. Acute compatibility and strictness imply that it lies in the relative interior of the retained facet: for any other normal `w`, the correction to `(c,w)` is nonpositive.

The restriction of `r` to `a-perp` is represented by

```text
r_F = r - (r,a)*a/q_a,
r_F² = r²-(r,a)²/q_a < 0.
```

Moreover `(r_F,c_F)<0`, so `r_F` has the same time orientation as `c_F`. It pairs strictly negatively with every future timelike point of `a-perp`. Therefore the added inequality is strictly satisfied on the entire hyperbolic supporting plane of the chosen facet. It creates no boundary there, removes none of the old facet, and cannot create a finite vertex on it.

This is a direct half-space argument, with no need to assume a general facet-volume theorem. Applying it to every added wall shows that the old facet, including its old ridges, remains unchanged.

## Application to the binary chamber

The [independent local certificate](fresh_h20_weighted_binary_20261014_local_review.md) proves that no finite vertex of the frozen chamber is incident to a binary cap. Its finite-neighbor inventory also includes eight retained golden seed walls at angle `pi/5` to that cap.

These are actual ridges, not merely finite-angle intersections of supporting planes. For any such pair, its two-root Gram matrix is positive definite with inverse having nonnegative entries. Projecting the common strict timelike center to their common orthogonal subspace produces a timelike point strict on every other wall: the projection adds a nonnegative combination of the two outward normals, whose pairings with all other walls are nonpositive. Hence the common codimension-two face is nonempty and actual.

Under dotted-only additions the chosen facet and these ridges remain unchanged. If the resulting chamber had finite volume, the first lemma would force a finite vertex on one of these ridges, hence on the chosen facet. The independent local obstruction says that facet has none, and dotted additions cannot create one. This contradiction proves the stated no-completion theorem.

## A simpler obstruction at the retained central facet

For the weighted 39-wall core, the same obstruction already applies to its central radius-two facet, independently of the binary finite-neighbor graph. Use normalized spatial coordinates with weights `1` on ten axes, `2` on one axis, and `phi^2` on eight axes; their sum is `N=20+8phi`. The vertical seed inequalities are `-1<=x_i<=1`. The central normal is `(0,-2,1)` in coordinates `(x,D,U)` with norm `x^2-2DU`.

On the central plane, the chart `U=1` gives `D=2`. A finite vertex on this 19-dimensional facet would require nineteen independent active vertical seed walls. At most one wall from each opposite pair can be active, so all nineteen coordinates would have to equal `+1` or `-1`. But such a point has Lorentz norm `N-4=16+8phi>0`, not a finite or ideal norm. Thus the central facet has no finite vertices.

Its positive golden seed ridges are actual and have angle `pi/5`. For a golden coordinate `g`, the point `(e_g,2,1)` is orthogonal to the central normal and the positive `g` seed normal, is strict on every other seed wall, and has norm `phi^2-4=phi-3<0`. The two normal norms are `4` and `phi^2`, their pairing is `-phi^2`, and their normalized cosine is `phi/2=cos(pi/5)`.

Consequently any finite acute-compatible completion retaining this core, with every additional wall strictly dotted to the central facet, still has infinite volume. The full binary orbit and the new outer-code, hole, and small-cap families all satisfy this dotted hypothesis, as their exact pair certificates verify. This simpler application does not require the binary 114-node local graph. A productive completion must include a wall meeting the central facet or replace that facet; adding only central-dotted layers cannot suffice.

This note is a mathematical proof review based on the frozen local graph/certificate and the whole-chamber strict-center/pair hypotheses. It adds no numerical search and does not alter the frozen source or certificate.
