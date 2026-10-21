# Root review: actual mixed replacements, including finite targets at cusps

20 October 2026. **PASS, full written review.** No scientific job, build,
enumeration or Lean execution. This is a necessary-condition tool, not a
dimension theorem.

Reviewed in full:

* `compact_farthest_vertex_mixed_replacement_inequalities_20261020.md`,
  SHA256 `3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84`;
* the separate outside-author review
  `review_farthest_mixed_replacement_resource_20261020.md`,
  SHA256 `f869aca8774e6584bc9693187ecb6c298d86d7337b6bbb9894cef0b8748c8d6e`.

Paths are under `results/`. I independently rederived the projection,
Lorentz normalization, cap-subset criterion, weighted inequality and the
two-cap diagnostic, rather than inferring validity from the review.

For an actual acute-wall Gram G which is positive definite, the
disjoint-positive-parts argument proves `G^-1 h>=0` when `h>=0`.
Thus `q0=p+M G^-1 h` satisfies every other actual wall inequality;
`q0^2=p.q0=-(1+D)` proves both timelikeness and future orientation.
Full selected rank n gives an actual finite vertex. This does not
establish realizability of abstract wall data.

Compactness bounds all actual projected points by the vertex maximum:
normalizing a convex combination of future unit vertices divides by a
number at least one. Writing the projected point as
`(p-v)/sqrt(1-R)` consequently gives
`x.v>=1-sqrt(1-R)`. The strict bound `x.v<sqrt(R)<1` justifies the
squaring step `R<=2s-s^2`, with the zero replacement treated separately.

For isolated actual caps the Gram is the stated diagonal congruence of
`K-r r^T`. Its positive definiteness is exactly `r^T K^-1 r<1`.
Adding the retained old walls gives a full-rank actual replacement;
solving its equations gives `v=sum (K^-1 r)_i u_i`. No positivity of
these spatial coefficients is necessary. Nonnegative combinations of
the replacement inequalities give the strict weighted certificate by
Cauchy--Schwarz and `|x|<1`.

The displayed B20 diagnostic has `0<R<1` and `R>0.9048`, whereas
`0<s<0.618` implies `2s-s^2<0.854076`. It therefore fails the
two-cap replacement requirement with the advertised margin. Its
separate native-alphabet failure is disclosed; this is not a realized
Coxeter example or a global exclusion.

## Finite-volume extension with its precise boundary

The outside review's useful extension is also correct. Let a
finite-sided, finite-volume acute polytope have a nonempty set of finite
vertices. Choose any interior future unit z and select p maximizing
`-z.p` over those finite vertices. The decomposition still has `|x|<1`.
Every **full-rank positive-definite actual-wall replacement** is itself
a finite vertex, so its comparison and all cap-subset/weighted
inequalities remain valid, irrespective of the number of cusps.

This proves no comparison with an ideal endpoint or an arbitrary
lower-rank projected point. It supplies no finite vertex in an entirely
ideal polytope, and no adequate replacement-coverage theorem. These
remain genuine obstacles to a universal finite-volume dimension bound.
The compact rank-at-most-four exclusion is not automatically extended
to cusps by this observation.

**Disposition:** retain within the existing farthest-vertex research
package; no additional primary manuscript or formalization until it
contributes to a material dimension result.
