# The direct A2 coupling is also obstructed

Every coherent direct A2 coupling of the fixed complete40-facet H18 source has infinite volume. The proof covers all finite inducing-contact labels, including π/5, and all assignments of those profiles. It does not cover additional walls that induce no source facet, changed source metrics, or general H20 completions.

A focused repository search found no earlier direct A2 coupling theorem for this source. The [B2 result](fresh_h20_b2_h18_coupling_20261014.md) cannot simply be transferred: A2 permits different profiles and permits special endpoints on an individual source π/3 edge. The new obstruction instead uses the connected source diagram and an exact common escape.

The [producer](../scripts/fresh_h20_a2_h18_coupling_20261014.py) and [certificate](fresh_h20_a2_h18_coupling_20261014.json) pass generation and saved replay. They reconstruct all1600 source Gram entries and check the source diagram, exact source-point gaps and three scalar profile margins. They do not enumerate profile assignments or expand a large ambient chamber.

## Complete inducing profiles

Normalize the two A2 support roots to norm2, with Gram `[[2,−1],[−1,2]]`. Their fundamental weights satisfy

`ω0²=ω1²=2/3`, `(ω0,ω1)=1/3`.

Let n be a fixed source root of norm q and normalize an inducing ambient root to norm2. Its possibilities are exactly

`r=√(2/q)n`,

or

`r=√(h_m/q)n−a_m ω_c`, for `c∈{0,1}` and `m∈{3,4,5}`,

where `a_m=2cos(π/m)` and `h_m=2−2a_m²/3`. Thus there are seven profiles: ordinary, and two colors at each of π/3, π/4 and π/5.

Indeed two nonzero support contacts already cost at least2 in the positive A2 projection, leaving no positive source projection. A single contact has positive source projection exactly when `cos²(π/m)<3/4`, giving m=3,4,5. At m=6 the projection is null; larger finite labels and parallel/dotted support contacts cannot induce a facet of the fixed H18 face. In particular omitting π/5 would make this classification incomplete.

For profiles i,j with projected norms h_i,h_j, the exact pair rule is

`(r_i,r_j)=√(h_i h_j/(q_i q_j))G_ij+a_i a_j J_ci,cj`,

where `J_same=2/3`, `J_opposite=1/3`, and the last term is zero if either profile is ordinary. This formula applies to every source pair and supplies the complete metric compatibility rule without an angle cutoff.

On a source π/3 edge the only admissible pairs are:

| First profile | Second profile | Colors | Ambient relation |
|---|---|---|---|
| Ordinary | Ordinary | Irrelevant | π/3 |
| π/3 | π/3 | Same | Right |
| π/3 | π/4 | Opposite | Right |
| π/4 | π/3 | Opposite | Right |

An ordinary/special pair has product `−√(h_m/2)∈(−1,0)`, which is not any acute Coxeter product for two norm2 roots. The remaining special pairs are positive except for opposite-color π/3–π/3, whose product is `−1/3`, also forbidden. The certificate checks the complete color-collapsed transition table by exact squared comparisons.

## Why24 source facets must be ordinary

The source's27 π/3 edges form one connected component on24 vertices. Ordinary/nonordinary status must agree across every such edge by the preceding rule. This component contains the path `0–1–2`, while source roots0 and2 are orthogonal.

If one node of the component were special, every node would be special. But any two special profiles on orthogonal source facets have strictly positive ambient product: both A2 inverse-Gram entries are positive. Roots0 and2 would violate acuteness. Therefore all24 component vertices are ordinary. This uses the whole source diagram, not the false claim that a single π/3 edge forces ordinary endpoints.

The remaining16 source facets retain all seven profiles in a relaxed bank. Including both A2 support walls gives `24+16·7+2=138` normals. Some choices in this bank are mutually incompatible; keeping the entire bank only strengthens the escape argument.

## One exact point misses the entire bank

Use the literal source point x from the earlier B2 certificate, whose coordinates and all40 pairings are rechecked here. It satisfies

`x²=(3857−2800√2)/529<0`,

and `(x,n_i)≤−√q_i` for all16 potentially special source facets. For every other source facet `(x,n_i)≤0`. Let h be the sum of all40 source normals. Its norm is `−2952−1280√2`, it is future timelike, and all its source-wall pairings are strictly negative.

Take the exact rational constants

`s=3207/1000`, `ε=1/10^7`,

and the physical point

`E=(s x+εh,−ω0−ω1)`.

The transverse coordinates in the simple-root basis are `(−1,−1)`, its transverse square is2, and its pairing with either support root is−1. The certificate verifies

`s²h_m>a_m²` for m=3,4,5,

as well as `s²x²+2>0` and the positive square after the strict perturbation. In particular the difficult π/5 inequality is

`s²(2−2φ²/3)>φ²`.

For every special profile,

`(E,r)<−s√h_m+a_m<0`.

Every ordinary profile is also strict because of εh. Thus E is strict for all138 normals simultaneously. Its square is exactly

`509083214552184799/6612500000000000 − (281226472383113/5166015625000)√2 > 0`.

The separate point `C=(100h,−ω0−ω1)` is strict timelike for the same whole bank, and `(E,C)<0` verifies the future orientation. The segment from C to E crosses the future null cone inside every wall halfspace, giving an open ideal-boundary region. Every coherent direct assignment therefore has infinite volume. This conclusion does not need a finite-vertex or cusp census.

The complete finite-volume source is inherited from its frozen H18 construction. The arithmetic here reuses the standalone exact source-field helper from the B2 review; it is a new producer, not an independent replay of itself. This closes the direct A2 alternative at the feasibility stage. A constructive continuation would need additional noninducing walls or a changed source, rather than another assignment of these inducing profiles.

The [independent review](fresh_h20_a2_h18_coupling_20261014_review.md) expands all 138 normals in a separate multiquadratic field, checks their norms and all strict escape/center signs, reproduces the transition table and verifies the ambient inertia. Its different exact escape perturbation also works. Generation and saved replay pass.
