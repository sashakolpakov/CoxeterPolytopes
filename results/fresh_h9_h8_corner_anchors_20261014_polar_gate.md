# No low-polar anchor can complete the raw H8 attachment

For the frozen mixed-D4 twelve-wall chamber and the fixed compact 11-facet Bugaenko H8 source, **every all-dotted cap direction with low normalized projected center 11 fails the raw perpendicular attachment**. This covers all three exceptional source corners and every contact scale. Together with the [uniform radial gate](fresh_h9_h8_global_cap_gap_20261014_radial.md), it leaves only cap directions with all twelve normalized projected norm squares at least `sqrt(5)-1`.

The [exact producer](../scripts/fresh_h9_h8_corner_anchors_20261014_polar_gate.py) and [certificate](fresh_h9_h8_corner_anchors_20261014_polar_gate.json) pass generation and saved replay. They check literal physical Gram identities, outward 192-bit rational interval bounds and three symbolic identities. The [independent review](fresh_h9_h8_corner_anchors_20261014_polar_review.md) also passes generation and saved replay using separate 128-bit and sparse-polynomial arithmetic. No sampled cap or optimizer is a proof dependency.

Here the raw attachment consists of the old twelve walls, one all-dotted cap r and the source walls perpendicular to r. Additional nonperpendicular walls can cut its uncovered end and are outside this theorem. Neither this result nor the preceding scalar gates exclude arbitrary larger completions retaining the source face.

## Coverage forces a small wall-4 contact

Normalize r and every old normal to square two. Put `a=(sqrt(5)-1)/4`, `B=4a` and `L²=2(B+2)`. Assume `kappa_11²<B`. The [source gap theorem](fresh_h9_h8_facet_cover_20261014.md) and the [two-anchor theorem](fresh_h9_h8_corner_anchors_20261014.md) imply:

- The projected center of old wall 11 is one of source vertices 0, 32 or 33.
- Its projected norm square is at least a, while every other projected norm square is at least B.
- The other two exceptional source vertices are each at distance at least `D_A=acosh(9+4sqrt(5))` from that anchor, and are mutually at least that far apart.

Old wall 11 is right to every old wall except 4 and 8. The checker proves that the maximum anchored reach of a ball belonging to any of those right neighbors is less than `D_A`. The anchor's own ball also cannot reach either other exceptional vertex. Every high-norm ball has diameter at most `2 asinh(sqrt(2/B))<D_A`, so no single ball covers both required vertices. Therefore balls 4 and 8 must each cover one. In particular, ball 4 must have anchored reach at least `D_A`.

Let `A=-(r,n_11)`, `D=-(r,n_4)` and `C=-(n_11,n_4)`, using normalized n_11. Exact physical pairings give `44/5<C<9`. Since `kappa_11²<B`,

`A²<2(B+2)<169/25`, hence `A<13/5`.

If `D>=4`, then `kappa_4²=D²/2-2>=6`. Anchored reach decreases as either projected norm increases, and increases with C. Even at the relaxed values `kappa_11²=a`, `kappa_4²=6`, `C=9`, the reach is less than `D_A`. Thus coverage forces **D<4**.

The comparison uses cosh expressions only. For projected norm squares x,y and old negative product C, the center-distance cosh is

`(sqrt((x+2)(y+2))+C)/sqrt(xy)`.

It decreases in both x and y. Adding the high-ball radius preserves that monotonicity, and the resulting reach bound is approximately 3.519, below `D_A≈3.580`.

## Four old roots give a norm contradiction

The literal old roots n_0 and n_2 are mutually orthogonal and both orthogonal to n_4 and n_11. Their norms are two. The span of n_4,n_11 has Lorentz Gram

`[[2,-C],[-C,2]]`.

Its orthogonal complement is positive definite. Because centers 0 and 2 are high, their contact squares satisfy `(r,n_i)²>=L²>6`. Orthogonal projection therefore gives

`r² >= L² - 2(A²+D²+C A D)/(C²-4)`.

The magnitude of the negative term increases with A,D and decreases with C>2. The rational bounds above yield

`2(A²+D²+C A D)/(C²-4) < 2857/918`.

Consequently `r²>6-2857/918=2651/918>2`, contradicting the normalization `r²=2`. The certificate checks the positive margin `815/918` exactly.

The argument uses coverage only to bound D. It does not require enumerating E8 contact profiles, optimizing a source rotation, or choosing discrete low-anchor scales. The remaining fixed-source direction search must keep all projected norm squares at least B.
