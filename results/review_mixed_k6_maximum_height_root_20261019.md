# Root review: clique-six maximum-height triangle cap

19 October2026. **PASS, full written-source review.** Source
`right_angled_mixed_k6_maximum_height_triangle_cap_20261019.md`, SHA256
`cf000f2f49c8b11ce20f21678258a860a50b90eee3ed5396177645c5cf93fb1d`.
No graph enumeration, optimizer, scientific program or Lean job was used.

The triangle packing argument and the three outside-edge cases check.
For common-neighbor counts at least3,2,1 the respective bounds52,55,44
follow from the exact core-neighborhood unions and Mantel's nine-edge
bound. In the zero-common-neighbor case, a total at least57 forces a
five-neighbor vertex and two others of degree at least four. The
three-disjoint-triangle construction correctly makes every contributing
outside neighborhood avoid its omitted core vertex. Replacing the
six-clique by that five-clique preserves every triangle count.

The five selected degree patterns follow from the47 pair-count target.
I recomputed each convex objective value and each of its four gradient
prefix sums. They are exactly1,1,30/29,35/33,39/37, with nonnegative
prefix sums and zero total gradient. The last row's only nonconstant
successive-coordinate difference has zero prefix coefficient. Thus
summation by parts proves the global ordered-coordinate lower bound;
the calculation does not assume that an optimizer has equal coordinates
or that the positive parts remain active at an unknown minimizer.

The cap56 therefore holds when the actual maximum-height bad-triangle
graph has clique number six. Together with the previously reviewed K7
and K8 cases, this settles clique numbers at least six. Clique number
at most five remains unresolved. There is still no bound on an ordinary
global triangle average, no required carrier density, and no improvement
of the mixed dimension ceiling12.
