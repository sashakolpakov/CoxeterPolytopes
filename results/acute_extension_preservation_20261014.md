# What an acute extension preserves

This elementary lemma applies to the 95-to-96-wall construction and to future additions. It does not establish finite volume.

Let the old facet normals be spacelike and pairwise nonpositive. Add spacelike normals that are nonpositive against every old normal and against one another. Suppose the enlarged system has a common strict future timelike point `C`.

**Every listed wall is an actual facet.** For a normal `n`, the point

`P = C − (C·n)n/n²`

is timelike, lies on its hyperplane, and is future because `C·P = P² < 0`. Its pairing with any other normal is strictly negative: the correction adds a positive multiple of `n`, whose pairing with the other normal is nonpositive. Thus the wall contains a relatively open facet.

**Every old finite vertex is strictly retained by every new wall.** Choose the independent incident normals at such a vertex `v`, with positive-definite Gram matrix `G`. Its inverse is entrywise nonnegative. One direct proof writes `G⁻¹ = a⁻¹ Σₖ(I−G/a)ᵏ` for sufficiently large positive `a`; the summands are entrywise nonnegative and the series converges.

Write `C = λv + Σcᵢnᵢ`. Then `λ > 0` by time orientation and `c = G⁻¹(C·nᵢ)` is entrywise negative. For a new normal `r`,

`λ(r·v) = r·C − Σcᵢ(r·nᵢ) < 0`,

because `r·C < 0`, `cᵢ < 0` and `r·nᵢ ≤ 0`. All its old nonincident inequalities remain unchanged. No old finite vertex is removed.

**Every old bounded full cusp is strictly retained.** At a proper cusp `ν`, the incident normals span `ν⊥` and their Euclidean affine components supply a relation `ν = Σaᵢnᵢ` with every `aᵢ > 0`, after positive rescaling. The positive signs follow by pairing each affine relation with a strict future center. Acuteness gives `r·ν ≤ 0`. Equality would force every incident pairing to vanish; then `r` would lie in `(ν⊥)⊥ = span(ν)`, contradicting that `r` is spacelike. Therefore `r·ν < 0`, so the cusp and its existing bounded link are preserved locally.

For the present construction, these conclusions require the exact all-wall acuteness, positive normal squares and common strict center. They cannot be inferred from a numerical angle screen alone. They also do not turn newly created open edges into bounded cusps or establish a complete vertex census.
