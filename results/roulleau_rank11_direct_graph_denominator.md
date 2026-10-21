# Closure of Roulleau's two open direct-graph rank-11 stocks

The two remaining direct weighted graphs in rank 11 are now transcribed and
certified.  Their full nef chambers are published finite-volume noncompact
`H10` Coxeter polytopes, but their complete diagram-subgroup denominator gives
**no proper fixed section above `H9`**.

| source stock | weighted graph group | complete subgroup result | published full chamber | lineage |
|---|---|---|---|---|
| `U+E8(-1)+A1(-1)`, 12 walls | trivial | one subgroup; no proper fold | `H10/12`, 18 finite vertices, 1 cusp, 99 elliptic edges | explicit wall-by-wall conjugacy with the intrinsic all-root `U+E8(-1)+A1(-1)` Vinberg chamber; full index 1 |
| `U+D8(-1)+A1(-1)`, 14 walls | `C2` | identity and one involution; the involution has fixed rank 10 and gives only `H9` | `H10/14`, 35 finite vertices, 4 cusps, 218 elliptic edges | three incomparable reflection-stable index-2 determinant-two overclosures: two even `U+E8(-1)+A1(-1)` maxima, each containing the source group with exact index 1023, and one odd maximum whose numerical index remains open |

The graph transcriptions are pinned by the published divisor arithmetic.  In
the first stock, `A1,...,A11` have determinant 2, the dependent `A12`
coordinate is integral, and the published ample class has square 848 and the
stated intersections.  In the second, the published nonconsecutive basis
`A1,...,A9,A11,A12` has determinant 8; the coordinates of `A10,A13,A14` are
integral; and the ample class has square 208 and precisely the fourteen stated
degrees.  These checks also disambiguate the overlapping thick segments in
the printed TikZ graph.

The closure above `H9` is structural as well as computational.  A nonidentity
faithful action on a rank-11 lattice has fixed rank at most 10.  The first
weighted graph has no nonidentity symmetry; the second has exactly one, its
integral Néron--Severi action is faithful, and its fixed rank is exactly 10.

Neither source chamber is new: both are explicitly published by Roulleau.
The second `H10/14` relation graph is absent from both pinned CoxIter lanes,
but that bounded absence is not a priority claim.  Roulleau's general-surface
`C2` has tag `t` and acts trivially on `NS`; the nontrivial diagram involution
of the second stock acts nontrivially and therefore is not that surface
automorphism.

Replayable files:

- [`roulleau_rank11_direct_graph_denominator.json`](../data/enumerations/roulleau_rank11_direct_graph_denominator.json)
- [`certify_roulleau_rank11_direct_graph_denominator.py`](../scripts/certify_roulleau_rank11_direct_graph_denominator.py)
- [`test_roulleau_rank11_direct_graph_denominator.py`](../tests/test_roulleau_rank11_direct_graph_denominator.py)

