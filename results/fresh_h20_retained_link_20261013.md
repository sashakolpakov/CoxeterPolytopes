# Mixed-factor cusp sections of the Borcherds H21 chamber

No interior Coxeter H20 section through any of the 18 archived cusp representatives can have its section normal active in two or more irreducible affine factors of that cusp link. This includes arbitrary real mixing coefficients and arbitrary translation along the null cusp ray. It leaves single-factor directions open; it is not a nonexistence theorem for H20 polytopes or for sections avoiding these cusps.

The [exact certificate](fresh_h20_retained_link_20261013.json) and [standard-library checker](../scripts/fresh_h20_retained_link_20261013.py) verify the literal data used below: all 18 cusp rays, their 53 affine factors, 411 positive-definite one-wall deletions, the 247 connected deletions, and all 2,067 possible attachment nodes in those connected deletions. Of the attachment tests, 2,039 use an old angle π/3 and 28 use π/4. The only rank-one factor is the pair of walls 162 and 164 at cusp 0; its other factors have ranks 5, 7, and 7.

The argument uses the usual acute-angle property of Coxeter chambers: outward normals of every pair of actual facets have nonpositive pairing, including nonadjacent facets. A bounded Euclidean Coxeter chamber is a product of simplexes. These are classical facts about acute-angled polyhedra; see [Vinberg, *Hyperbolic reflection groups*, §2](https://www.mathnet.ru/eng/rm2140). Everything specific to this chamber is rechecked with rational arithmetic.

## Reduction using actual facets

Work in the positive-definite quotient of the cusp's null-orthogonal space. Its dimension is 20, and its chamber is the product of the archived affine simplexes. Normalize the tangent section normal `e` to square one. A factor is *active* when the projection of `e` into that factor is nonzero. Adding an arbitrary multiple of the null ray to the original hyperbolic section normal changes the affine placement of the section, but none of these projected Gram calculations.

Assume the section meets the strict interior and its bounded 19-dimensional Euclidean link is Coxeter. Take an irredundant set of actual facet inequalities, choosing only one representative when projected inequalities coincide. In an active factor of rank `n >= 2`, fix all other factor coordinates at a strict section point. The remaining fiber is a bounded, full `(n-1)`-dimensional slice. Consequently at least `n` of that factor's `n+1` original inequalities must remain: inequalities from other factors are constant on this fiber and cannot bound it. Thus at most one inequality from that factor is omitted. This argument also covers possible duplicate projected inequalities; the chosen irredundant representatives still have to bound the fiber.

An active rank-one factor might instead lose both bounds. There is only one such factor in the entire census, at cusp 0, and this exceptional case is already excluded by the frozen [A1-absorption theorem](fresh_h20_cusp_sections_20261013_absorption.md). For the rest of this proof every active factor retains at least all but one of its inequalities.

The positive affine relation in each factor forces both signs among the original pairings `(e,r_i)`. At least one retained pairing is nonzero, because the roots remaining after one deletion span the factor. If all inequalities of an active factor survived, both signs would remain there. A nonzero pairing in any other active factor would then have opposite sign to one of them. Since roots in different factors are orthogonal, their projected pairing would be

`(r_i,r_j) - (e,r_i)(e,r_j) = -(e,r_i)(e,r_j) > 0`,

contradicting the acute-angle property. Hence exactly one inequality is omitted from each active factor.

Let `N` be the sum of the active factor ranks. The active part of the bounded section then has exactly `N` actual facets in dimension `N-1`: it is a simplex. Its Gram matrix is irreducible, positive semidefinite, and has a strictly positive kernel. Inactive factors split off orthogonally.

A three-node negative Coxeter triangle has normalized off-diagonal magnitudes `a,b,c >= 1/2`. Its determinant is

`1 - a^2 - b^2 - c^2 - 2abc <= 0`,

with equality only at `a=b=c=1/2`. In a positive-semidefinite Gram matrix this exceptional triangle has a positive null vector and must be orthogonal to every other node. Therefore a larger irreducible Euclidean simplex cannot contain such a triangle. The three-facet exception cannot occur here: there is at most one rank-one factor, and its possible partners have ranks 5 and 7; all other factors have rank at least two.

All retained nonzero pairings in different factors must have the same sign, to avoid the positive projected pairing above. Three active factors would therefore give a negative triangle. Two nonzero retained pairings in one factor, together with one in another, would also give a negative triangle: the old same-factor pairing is already nonpositive. We conclude that there are exactly two active factors and exactly one retained root with nonzero pairing in each.

For either factor, delete its omitted root and call the resulting positive-definite Gram matrix `B`. This retained matrix must be connected. Otherwise a component not containing its sole changed root would stay isolated and positive definite, contradicting the positive kernel of the active simplex. Connectedness is thus a consequence, not an assumption about deleting a node of an affine diagram. The checker explicitly distinguishes 247 connected deletions from 164 disconnected ones.

## Two factors of rank at least two: all labels excluded

Let the sole attachment root have old square `L`, pairing `u` with `e`, and inverse diagonal `A=(B^-1)_{ii}`. The normal's squared contribution from this factor is

`q = A u^2`.

Set `D=A L`. Choose any unchanged neighboring root; one exists because the retained matrix is connected and has rank at least two. In this source its old squared cosine `c0` is either `1/4` or `1/2`. Projection only decreases the attachment root's squared norm, so its new squared cosine is

`c = c0/(1-u^2/L) > c0`,  and  `q = D(1-c0/c)`.

The two-node inverse bound gives `D >= 1/(1-c0)`. The checker independently verifies this bound for every one of the 2,067 attachment nodes, preferring a π/3 neighbor when one is available. The new angle must have a strictly larger integer label than the old one. Therefore:

- an old π/3 contact changing to any label `m >= 4` costs `q >= 2/3`;
- an old π/4 contact changing to `m >= 5` costs `q >= 2(sqrt(5)-2) > 1/3`;
- an old π/4 contact changing to `m >= 6` costs `q >= 2/3`.

Since the two factor contributions must sum to one, the only case not already excluded by these inequalities has both old π/4 contacts changing to π/5. In that case

`q_1 + q_2 = (D_1+D_2)(sqrt(5)-2)`.

Both `D_i` are positive rational numbers, so this irrational number cannot equal one. This is an all-label exclusion; no finite list of angle labels is imposed.

## The retained-one-A1 case

At cusp 0 the other active factor is D5 or A7, with only π/3 edges. Let `s=u^2/L`. Its changed neighboring angle implies `s >= 1/2`. The literal inverse census gives `D >= 2` for every connected D5 deletion, and `D >= 7/4` for every connected A7 deletion.

D5 would already use squared normal length `D s >= 1`, leaving no nonzero A1 component. For A7, the A1 contribution is `1-Ds > 0`. The squared cosine of the cross-factor contact is exactly

`(1-Ds)/(D(1-s))`.

It is positive and at most `1/7`, because `D >= 7/4` and `s >= 1/2`; the expression decreases in both variables on the admissible domain. It cannot be a nonright Coxeter angle, whose squared cosine is at least `1/4`. A right angle is impossible because both factors are active. This closes the case retaining exactly one A1 inequality. The case retaining both was already excluded by the sign argument, and the case retaining neither is covered by the cited absorption theorem.

## Verification and remaining scope

Run:

```sh
python scripts/fresh_h20_retained_link_20261013.py --verify
```

The checker reconstructs all wall pairings from the original integer Lorentz form, checks every cusp's full incident set, positive affine kernels and null-ray relations, and verifies every inverse bound used in the all-real argument. It binds the complete 18-cusp census and the previously replayed absorption certificate, including that certificate's source hashes. The geometric reduction from actual facets to the finite list of inverse bounds is the proof above; the checker is not a formal proof assistant.

Saved-result replay passed. A separate [SymPy reconstruction](../scripts/fresh_h20_retained_link_20261013_review.py) independently rechecked all 53 factors, 411 positive-definite deletions and 2,067 attachment bounds without importing this producer; its [certificate](fresh_h20_retained_link_20261013_review.json) passed. The reviewer also read the actual-facet reduction and all-label proof. This distinguishes the independent exact matrix census from the mathematical proof review.

A surviving constructive route must use a single irreducible affine factor at its chosen cusp, or choose a section avoiding these cusps. In particular, this result does not exclude arbitrary single-factor high-support normals, gluing, changed ambient chambers, or Coxeter polytopes unrelated to these sections. At cusp 0 the pure A1 plane is separately excluded by the earlier actual-ridge proof, leaving its D5 and A7 factors as the remaining local directions.
