The stated A1-absorption obstruction is sound. An interior section through the specified cusp cannot be Coxeter if both A1 inequalities are redundant in its Euclidean link. This is a conditional obstruction at this cusp, not a theorem excluding arbitrary H20 sections or sections that retain an A1 link facet.

I reviewed the [proof](fresh_h20_cusp_sections_20261013_absorption.md) and [source](../scripts/fresh_h20_cusp_sections_20261013_absorption.py), and freshly replayed its exact certificate. A [separate Fraction checker](../scripts/fresh_h20_cusp_sections_20261013_absorption_review.py) and [review certificate](fresh_h20_cusp_sections_20261013_absorption_review.json) additionally reconstruct the signed supports, all D5 star norm forms, the 16 A7 width calculations, and the pure-A1 ridge interpolation. It imports no absorption-producer arithmetic. The source and data hashes are bound; the prose report is not hash-bound so its review link may be added without a provenance cycle.

The essential scope reduction is justified. If the normal has zero A1 component, the independent interval's two bounds remain necessary. Otherwise the section is a graph over the other three factors. Redundancy of both interval bounds means the graph is defined over the whole D5 × A7 × A7 product. Consequently every remaining facet survives, and every cross-factor pair is an actual ridge. Two active factors would each have both signs in their root pairings, by the positive affine relations; choosing opposite signs in different factors gives an actual obtuse ridge. Thus at most one remaining factor is active. This is not an assumption that an arbitrary section activates only one factor.

Within that factor, all old simplex facet pairs meet. Acuteness therefore requires every positive-support node to be adjacent to every negative-support node in the original diagram. The absence of a four-cycle in the D5 tree and A7 eight-cycle forces a star. The independent enumeration gives 32 signed supports for each A7 and 26 for D5. A changed node with an unchanged neighbor has squared cosine \(1/[2(2-v_i^2)]>1/4\). The next allowable finite Coxeter angle is at most π/4, so \(v_i^2\geq1\). The projected facet is nonzero because the active component has square less than one. This bound treats every integer label, without a finite label cutoff.

The D5 obstruction admits a short independent norm check. Center the star at wall 29 or 48, write its outer-leaf magnitudes as \(p,q\), and the other-center magnitude as \(b\). Reconstructing the quotient metric gives the following exact squares:

| Support | Active component square |
| --- | --- |
| One outer leaf and its center | \(p^2/2\) |
| Both centers | \(b^2\) |
| Two outer leaves and their center | \((p^2+q^2)/2\) |
| Outer leaf, its center, and the other center | \(p^2/2+b^2\) |
| All three neighbors of a center | \((p^2+q^2)/2+b^2\) |

In the first case the central magnitude is \(p/2\geq1\), giving square at least two. In the third it is \((p+q)/2\geq1\), giving square at least one. In all other cases the other-center magnitude has an unchanged neighbor and is at least one. Every case contradicts the strict norm budget below one. These formulas were rederived for all 14 nonempty star choices at the two trivalent nodes, rather than copied from the producer.

For A7, the adjacent opposite pair is forced to magnitudes one; the two-leaf star violates the norm budget. The surviving update has active square 7/8. On the horosphere \((\delta,x)=-1\), its functional width is 7/8. The A1 component of square 1/8 has width \(\sqrt2/4\): \(r_{162}-\delta\) has norm four and ranges from −1 to 1. The squared width ratio is therefore exactly 49/8. A translation cannot make both A1 bounds redundant. The factors' independent full ranges matter here and are supplied by the graph-over-the-full-product premise.

Finally, the pure-A1 plane covers every remaining normal, including its arbitrary affine shift by δ. Only \(-1/2<t<1/2\) can meet the strict interior. The two saved ambient face points have exactly the asserted three zero wall pairings and are strict at every other one of the 210 walls; the separate replay checks all these signs. Their positive interpolation is coefficientwise orthogonal to \(e(t)\), is future timelike, and has precisely walls 45 and 63 active. The projected pairing \(1/4-t^2>0\) and determinant \(3-4t^2>0\) therefore describe an actual obtuse ridge throughout the whole interval. The boundary and exterior parameters are excluded by the supporting wall identities, and the pure-null direction is not a hyperbolic section normal.

All checks passed:

```sh
python scripts/fresh_h20_cusp_sections_20261013_absorption.py --verify
python scripts/fresh_h20_cusp_sections_20261013_absorption_review.py
python scripts/fresh_h20_cusp_sections_20261013_absorption_review.py --verify
```
