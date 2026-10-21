# Independent review of all91 rank-three global exclusions

All91 local records are excluded throughout their full interior section intervals, including every admissible lower endpoint. The global geometry and algebra pass independent reconstruction. No mathematical issue was found in the actual-ridge proof, field-root selection, or degeneracy treatment.

The [separate checker](../scripts/fresh_h20_rank3_global_sections_20261014_review.py) imports no producer or project arithmetic. It uses SymPy's rational algebraic-field arithmetic and real-root isolation, whereas the producer uses its own standard-library Fraction polynomial/Sturm implementation. The [saved review certificate](fresh_h20_rank3_global_sections_20261014_review.json) binds the frozen [global archive](fresh_h20_rank3_global_sections_20261014.json), its producer, literal source data, and the independently audited local census. Generation and saved replay pass.

| Independently checked | Count |
|---|---:|
| Continuous families | 91 |
| Obtuse / forbidden-gap families | 75 / 16 |
| Ambient face-point occurrences with all210 signs | 226 |
| Positive polynomial certificates | 305 |
| Endpoints excluded by the same pair | 64 |
| Separate endpoint witnesses | 22: four obtuse, eighteen gaps |
| Old-facet endpoint degeneracies | 5 |
| Rational / quadratic / quartic fields | 33 / 14 / 44 |

## Independent arithmetic

Each archived minimal polynomial is checked to be irreducible, and its designated rational interval contains exactly one real root. All saved scalar coefficients are interpreted in that field with their stated ascending/descending conventions. Algebraic signs use rational Horner enclosures and exact root-interval refinement; no floating-point sign decision is accepted.

For every record the review reconstructs the normal from its three retained roots and kernel coefficients, verifies its positive square q, its orthogonality to the affine null shift D, and the positive constant L defining the full shift interval. It reconstructs the affine positive relation and the omitted-wall redundancy coefficients. The local completeness theorem remains the separate [91-record saturation audit](fresh_h20_rank3_omitted_sections_20261014_review.md), whose source hashes are checked here; this global audit does not independently re-enumerate the conics a third time.

The review then reconstructs the pairing and norm polynomials from the old Gram and normal. These are scaled by q, so the angle quotient is still `pairing²/(first_norm*second_norm)`. Each positivity certificate is rechecked by restoring endpoint-zero factors, transforming the reduced polynomial to every archived rational subinterval, and recomputing all Bernstein coefficients. Those coefficients are strictly positive at the designated algebraic embedding, and the subintervals cover [0,1] without gaps. Thus the strict conclusions hold for every real `0<u<1`, not merely sampled parameters.

Every gap is either `0<cos²<1/4` or `1/4<cos²<1/2`, between consecutive Coxeter values. Neither the uniform crystallographic-angle theorem nor an arithmeticity assumption is used.

## Actual-ridge proof

The two rational ambient face points for each interval have exactly their stated zero walls, with every other old inequality strict, and are future timelike. After normalization they satisfy `(D,p0)=(D,p1)=-1`, `(e0,p0)=0`, and `(e0,p1)=-L`. Their positive interpolation lies in the section with normal `e0-u*L*D` and has precisely the common zero set K.

The K roots have positive-definite old Gram. Their section projections are independent because the section normal pairs with p0 as `uL`, whereas every K root is orthogonal to p0. They lie in the positive-definite tangent space at the timelike interpolation point, so they define a local orthant. A small perturbation keeping the claimed pair zero and making other K inequalities strict proves an actual ridge. This uses full geometry, not merely a positive pair Gram.

At u=1, the omitted wall is the sole additional zero source inequality. Its exact ambient identity makes it a nonnegative combination of K inequalities on the section; all coefficients outside K vanish. It is therefore redundant. Multiple minimizing kernel/mark ratios simply yield extra zero coefficients and do not invalidate the argument. Independence is still certified by `(e1,p0)=L`.

Sixty-four endpoints inherit a strict forbidden angle from the interval's pair. For the remaining22 interior endpoints, the checker reconstructs a different pair of ambient face points, the exact endpoint pairing/norms, and a strictly positive pair determinant, with either an obtuse angle or one of the same forbidden gaps. All226 face-point occurrences are checked against all210 walls.

In five cases every retained redundancy coefficient vanishes. The endpoint normal is then a nonzero multiple of the old omitted facet root, and its hyperplane does not meet the strict source interior. These are correctly marked boundary degeneracies, not silently accepted interior sections.

## Scope and consequence

The [main report](fresh_h20_rank3_global_sections_20261014.md) correctly claims all91 archived omission/target branches, not91 projectively distinct planes. The local census and its independent saturation check establish completeness of those branches; this review supplies the independent global geometry and arithmetic.

Together with the [coverage bridge](fresh_h20_cusp_coverage_20261014.md) and the separate mixed-factor, rank1/rank2, all-retained and rank-at-least-four exclusions, the result closes every **interior H20 Coxeter section containing a source cusp**. In particular, a remaining interior Coxeter section would have to be compact and avoid all old cusps. Supporting hyperplanes of source facets are boundary sections and need their separate earlier facet-section results if one removes the word "interior" from that conclusion. No claim about arbitrary H20 polytopes, changed ambient chambers, gluing, or subdivisions follows.

```sh
python scripts/fresh_h20_rank3_global_sections_20261014_review.py --verify
```
