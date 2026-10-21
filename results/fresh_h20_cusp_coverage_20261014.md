# Coverage audit for H20 sections through a source cusp

The archived local cases are exhaustive after the mixed-factor obstruction: **there is no additional single-factor facet-count case hidden by redundant or coincident inequalities.** An active affine factor of rank n>=2 leaves either n or n+1 distinct actual facets. The unique rank-one factor is covered separately by the A1 absorption proof.

The last **91 classified rank-three, one-omitted records** and all their interior lower endpoints are now [globally excluded](fresh_h20_rank3_global_sections_20261014.md), with a [separate arithmetic and proof review](fresh_h20_rank3_global_sections_20261014_review.md) passing. Combining the cases proves: **no Coxeter H20 hyperplane section meeting the strict interior of the literal Borcherds H21 source can contain a source cusp.** Every remaining interior Coxeter section would have to be compact and avoid all old cusps. Source-facet supporting hyperplanes are boundary sections and remain outside this particular theorem.

The [source-bound index checker](../scripts/fresh_h20_cusp_coverage_20261014.py) and [certificate](fresh_h20_cusp_coverage_20261014.json) verify the literal factor partition and coverage mappings. Generation and saved replay pass. The geometric completeness bridge is the proof below, not a finite collection of sampled slices. The checker does not duplicate all underlying exclusions or the complete corner walk.

## From a hyperbolic section to a bounded Euclidean section

Let P be the known finite-volume H21 chamber, H a hyperbolic hyperplane meeting its strict interior, and delta a source cusp contained in H. A sufficiently small cusp neighborhood has a bounded Euclidean link

```text
C = C1 × ... × Ck,
```

where each C_i is an affine Coxeter simplex of rank n_i and `sum(n_i)=20`. Intersecting H with a horosphere gives an affine hyperplane E in this product. Its intersection with C is bounded and full dimensional in E. Full dimensionality is not an extra assumption: the geodesic from an interior point of H intersect P to delta stays in the strict interior away from its ideal endpoint and eventually enters every sufficiently small cusp neighborhood.

Write the spatial normal as the sum of its mutually orthogonal factor components. A factor is active exactly when its component is nonzero. Adding a multiple of delta changes the position of E, but not the active factors or projected normal pairings. The normal has positive square, so at least one factor is active.

If P intersect H is Coxeter, its Euclidean link is a bounded Coxeter chamber. Such chambers are products of affine simplexes; the standard classification is given in [Felikson's lecture notes, section 1.2 and Table 1.1](https://www.maths.dur.ac.uk/users/anna.felikson/talks/cl_cox.pdf). Apply this statement only to the **irredundant actual facets** of the link, not to every input inequality.

## The fiber lemma handles redundancy and coincidence

Choose one original inequality representing each distinct actual facet of C intersect E and discard every redundant inequality. The resulting representation still defines the entire section. Each selected inequality is assigned to the original factor from which its chosen representative came.

Fix a point p in the strict interior of C intersect E. For an active rank-n factor, freeze all other factor coordinates at those of p. The remaining fiber is a bounded full-dimensional `(n-1)`-polytope in an affine hyperplane inside that factor. Every selected inequality from another factor is a **strict constant** on this fiber; it cannot bound any direction. Hence the selected inequalities assigned to this factor alone must bound it.

A bounded full-dimensional d-polytope with d>=1 requires at least d+1 inequalities with nonconstant normals. Indeed, at most d normals cannot positively span the d-dimensional translation space: if their rank is deficient there is lineality, and if exactly d are independent there is a direction strictly negative on all of them. Either gives an unbounded ray. Thus an active rank-n factor with n>=2 contributes at least n selected inequalities. Since it originally had n+1, it loses at most one.

This also addresses duplicate restrictions from different factors. If an inequality from factor i has the same restriction to E as one from factor j, it is constant on the factor-i fiber, because the factor-j expression is constant there. It cannot replace one of the n nonconstant inequalities required to bound that fiber. An arbitrary choice of irredundant representatives therefore cannot conceal a second omission in a factor of rank at least two.

For an inactive factor the fiber has dimension n, so all n+1 of its inequalities remain actual and distinct. It splits off orthogonally. Rank one is the only exception to the active-factor lower bound: its fiber is a point, and it can contribute zero inequalities. The source has only one rank-one factor, at cusp0 with walls162 and164.

## The mixed-factor theorem already uses this full scope

The [mixed-factor proof](fresh_h20_retained_link_20261013.md) does **not** assume that all old inequalities survive. Its first step is exactly the irredundant-representative fiber reduction above.

For every active factor, the positive affine relation forces both signs among its original normal pairings. If all n+1 inequalities of one factor survived, one of its pairings would have opposite sign to a nonzero retained pairing from another active factor. The original roots from different factors are orthogonal; their projected pairing would be positive. Both are actual facet normals, contradicting the acute-angle property of a Coxeter chamber. The retained roots of every other rank>=2 factor span that factor after at most one deletion, so at least one such pairing is nonzero. The same is true for a retained A1 bound.

Accordingly, after excluding absorption of both A1 bounds, every active factor in a mixed section contributes exactly n facets. If their total rank is N, the active section has dimension N-1 and N facets: it is a simplex, and its Coxeter diagram is irreducible affine. The existing negative-triangle argument then forces exactly two active factors with exactly one changed retained root in each, followed by the archived all-label inverse-norm obstruction.

There is no unhandled three-facet simplex exception in these data. The only rank-one factor has partners of ranks5,7,7; all other factors have rank at least two. Every potentially mixed irreducible affine simplex therefore has more than three facets, as required by that argument.

If both A1 inequalities can be omitted from an irredundant representation, removing them from the original description also leaves the section unchanged. Projection to the other factors is an affine isomorphism onto their **entire product**. This is the hypothesis of the [A1 absorption theorem](fresh_h20_cusp_sections_20261013_absorption.md), which includes arbitrary mixing with the other factors and closes that exception. It also covers the pure A1 normal plane. It is not limited to a preferred set of retained original normals.

## Exhaustive single-factor alternatives

After the mixed obstruction, inactive factors split off and exactly one affine rank-n simplex is sliced. For n>=2, the preceding lemma gives precisely the following alternatives:

| Actual local facets | Geometry | Archived coverage |
|---|---|---|
| n+1 | Product of two affine simplexes, if Coxeter | [Complete all-retained classification](fresh_h20_single_factor_classification_20261014.md) |
| n | A bounded affine Coxeter simplex | One-omitted classification, by rank |

No section with fewer than n distinct facets is bounded and full dimensional. A coincident pair of input inequalities counts once and belongs to the one-omitted case. An inequality touching only a lower-dimensional face is redundant and likewise counts as omitted. At rank2 the section is a segment, so there is no all-retained case; at rank1 the section is a point and the absorption theorem applies.

The complete eighteen-cusp census has 53 irreducible affine factors:

| Rank | Factors | Status |
|---|---:|---|
| 1 | 1 | Pure A1 and absorption globally excluded |
| 2 | 4 | [All real normal spans excluded](fresh_h20_single_factor_links_20261013.md), two C2 and two A2 |
| 3 | 3 | All-retained excluded; all [91 exact one-omitted records](fresh_h20_rank3_global_sections_20261014.md) globally excluded, including interior endpoints |
| At least4 | 45 | All-retained excluded; [387 omissions and their175 surviving local records globally excluded](fresh_h20_omitted_factor_classification_20261014.md) |

The index checker matches these literal factor sets against the complete cusp walk, checks all387 higher-rank omissions, and verifies all120 rank-three target systems: four omissions times ten labeled affine-triangle placements for each of three factors. The three ranks and omissions are not inferred from a matching count alone.

At a one-omitted lower endpoint, several minima `k_i/m_i` may coincide. This merely makes several coefficients in the redundancy identity vanish; it does not create a new spatial direction or a third facet-count case. For a global ridge proof one must still check the actual endpoint zero set and redundant-wall identity. If all retained coefficients vanish, the section is the old omitted facet hyperplane and misses the strict source interior. The remaining interior endpoint cases belong to the already stated one-omitted interval and must be tested, not dropped.

## Why this excludes noncompact interior sections

In the Klein model the closure of the finite-volume source is the convex hull of its finitely many finite and ideal vertices. Strict convexity of the boundary sphere implies that its only points on that sphere are its old ideal vertices. Intersecting with a hyperplane cannot create another boundary-sphere point. Therefore every noncompact interior section contains an old source cusp. Conversely, a section avoiding all old cusps is a closed subset of the source closure separated from the sphere and is compact in the hyperbolic metric.

The [complete corner walk](fresh_h21_source_corner_walk_20261013.md) supplies all eighteen cusp orbits under explicit wall-permuting isometries. Transport to one of the literal representatives preserves the source, active factors, angle conditions and redundancy. The classification is consequently not restricted to eighteen isolated rays; those representatives cover every source cusp.

The formerly conditional step is now supplied by the completed91-record global proof:75 continuous obtuse-ridge and16 forbidden-gap certificates,64 direct endpoint failures,22 endpoint follow-ups, and5 old-facet degeneracies outside the strict-interior scope. The coverage index binds that final certificate and its independent review, rather than relying on a reported success count.

This result does not exclude compact interior sections avoiding all cusps, source-facet boundary sections, different ambient chambers, gluing, subdivisions, or arbitrary H20 Coxeter polytopes. A compact section cannot be dismissed merely because the original source is arithmetic.

```sh
python scripts/fresh_h20_cusp_coverage_20261014.py --verify
```
