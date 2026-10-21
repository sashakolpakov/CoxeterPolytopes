# Right-angled projection/height kernel: independent source review

20 October 2026. Reviewer: `/root/resource_cleanup`.

**PASS_SOURCE_ONLY; no mathematical or interface blocker found.** This is
an independent statement/proof-intent and installed-API read, not a Lean
elaboration, axiom audit or geometric instantiation. No Lean, TeX or
scientific checker was executed or imported. No primary source, catalogue,
ledger or existing verified module was edited. This review does not change
the existing four verified kernels/36 actual-print status.

## Binding

* `formal/CoxeterProofs/RightAngledProjectionHeightKernel.lean`:
  `2fac8028e432a2104a58e423c65fc38c9e214fed5f52210e408aabe3f15f0e90`,
  282 lines,12333 bytes.
* `docs/right_angled_projection_height_logic_20261020.md`:
  `736146d3c9354c4126c8eb0bc22b3f75d8aa024e2d555689b9ad545bdea61f2a`.
* Imported finite-coordinate kernel:
  `formal/CoxeterProofs/RightAngledCubicScalarKernel.lean`,
  `f55568035b8852a6284f1bd404c897758f30b301b1d85f38f09e4d9de365e6ee`.

The initial read used projection source `4c5715bd…`. I independently
decoded its preserved gzip and checked that its full suffix from
`namespace CoxeterProofs.RightAngledProjectionHeightKernel` is byte-identical
to the held source. The change is a status-header comment only; the map
updates the source hash/byte count. The actual map hash above supersedes
the earlier abbreviated map notices.

## Mathematical checks

The coordinate form is explicitly `-x.time*y.time + sum(i<d,x_i*y_i)`.
Only the first `d` spatial entries are used; the unbounded function tail
does not create a claimed finite-dimensional geometric realization.
Symmetry and bilinearity follow from the imported transparent finite sum.

For `v=q+t*n`, the three primitive equations `q.q=-1`, `n.n=1`,
`q.n=-t` give `v.n=0` and `v.v=q.v=-1-t^2`. Positive `lambda` with
`lambda^2*(1+t^2)=1` gives the derived projected unit norm and negative
pairing with `q`. No future-point conclusion is an input. For every
**other** wall with `q.m<0` and `n.m<=0`, adding `t*n`, then scaling
positively, preserves strict negativity. The selected wall instead has
derived equality zero; no assertion of strict feasibility on that wall
is needed or justified.

When `n.m=0`, the exact new height is `lambda*h`. Hence
`(1+t^2)*(h')^2=h^2`; `0<h<=t` gives the local squared recurrence and
positive new height. This uses only a maximum of old heights, not an
improved bound for the new one.

The finite-descent identity is correct:

```text
(1+a)*(1-(k+1)*b)
  = (1-k*a) + (k+1)*(a-(1+a)*b).
```

With the previous strict invariant, the right side is positive. Failure
of the next invariant makes the left side nonpositive. The one-step
statement remains true without assuming `b>=0`. In `descent_strict`, the
induction hypothesis retains its `j<=steps` antecedent, and each successor
uses only the recurrence at `j<steps`. The zero case is exactly `0<1`.
Consequently the arbitrary finite-chain statement neither presupposes
an infinite chain nor the desired terminal height.

For `radial_gap`, the two hypotheses give
`p*hi^2<ri^2` and `p*hj^2<rj^2`. Multiplication using the stated signs
gives `(p*hi*hj)^2<(ri*rj)^2`; positivity removes the squares. Scaling
`ri*rj*z<=hi*hj` by positive `p` then proves `p*z<1`.
This works for every `p>0`, not just2 or4. The all-wall product inequality
is a primitive unnormalized Gram bound, not the asserted strict radial
gap. Substituting the one-/three-step height bounds gives the intended
strict gaps without assuming them.

The hypotheses are nonvacuous. Over rational scalars, take spatial
dimension2, `q=(5/3;4/3,0)`, `n=(0;-1,0)`, `t=4/3`, `lambda=3/5` and
`m=(3/4;0,5/4)`. All unit/incidence, positivity, orthogonality and
`0<h<=t` hypotheses hold, with projected point `(1;0,0)` and old/new
heights `5/4,3/4`. Independently, a chain `a_j=1/(j+1)` satisfies the
cleared recurrence and strict invariant at every finite index.

## API and scope

I checked the finite-sum signatures in the imported source and the named
`OrderedRing` multiplication, square and Nat-cast lemmas in installed
Lean4.32.2 source. Their argument order and hypotheses match the uses.
The finite induction binder placement and cast rewrites are appropriate.
Actual tactic elaboration and runtime remain untested; this source read
cannot replace root's bounded compilation.

Exactly eight axiom requests name raw projection identities, normalized
projection, normalized strictness, orthogonal height identity, maximal
squared step, one-step descent, finite descent and radial gap. The raw
strictness helper is included through the normalized theorem's dependency.
No actual printer output is claimed by this review.

The map correctly retains existence of positive normalizers, actual unit
walls/interior points, maximum selection, right-angle all-other-wall signs,
intrinsic ridge-neighbor completeness, future-sheet interpretation and
face inheritance as geometric inputs. Actual radial vectors, their
indexing/injectivity and tensor realization are also outside this module.
The source does not yet construct the old `ProjectionGeometry` record or
formalize a hyperbolic polytope. None of these stated inputs assumes the
new strict height/radial bounds or a final dimension conclusion. The
remaining geometric boundary is explicit rather than hidden in a
conclusion-shaped field.
