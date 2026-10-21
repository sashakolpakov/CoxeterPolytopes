# A bounded three-cap roof closes one heptagon but adds open boundary

This is a **numerical diagnostic**, not an exact new chamber or a finite-volume construction. Generation and a fresh saved replay both pass. No roots in this package are promoted into the frozen 171-wall branch.

Starting from candidates0/77 in the [two-cap package](fresh_h20_replacement_outer_patch_pairs_20261017.md), the search tests seven right-T finite endpoint profiles, old-anchor labels2 through6, and both norm-completion branches. It obtains35 endpoint seeds on173/174 and32 on171/173. All1,120 seed pairs pass their mutual native Coxeter gate. For each pair, a middle root is prescribed right to173 and at labels2 through6 to the two seeds. The two exact-form norm-completion formulas are evaluated numerically.

Among56,000 real branches,34,777 fail old finite-vertex strictness and21,215 fail retained native Coxeter gates. The remaining **eight** have all retained native gates, the two closing spherical supports, and a complete numerical target-H3 boundary ledger. Each closes facet173 as a compact heptagon. Seven have26 finite/10 exterior positive-side-polytope vertices and13 open hyperbolic edges; the last has25/11 and14 open edges. The base two-cap branch had21/9 and10 open edges. The maximum squared exterior radius remains1.2186340819392139.

Candidate2 uses only labels2 and5 for the new prescriptions:175 is right toT,173,159 and meets174 atπ/5;176 is right toT,171,173,161;177 is right toT,173 and meets175,176 atπ/5. These are recoverable numerical contact recipes. Their new outer facets remain unfinished, so the larger roof does not close an exterior component.

The native ambient gate test uses the inherited T-projection offsets; all candidates here happen to be right toT. The ledger exhausts positive-side plane triples and clips pair-intersection lines by the unit ball and every retained halfspace, as in the [boundary package](fresh_h20_replacement_boundary_20261017.md). Exterior Euclidean vertices are not claimed to be hyperbolic vertices. Floating tolerances and the bounded profile/label domains are explicit in the [source](../scripts/fresh_h20_replacement_three_cap_roof_20261017.py) and [saved data](fresh_h20_replacement_three_cap_roof_20261017.json).

This branch has176 actual walls (labels through177, with31 and158 absent). It remains outside the exact171 certificate. It is retained to show that closing one new facet can increase the unfinished boundary; it is not the chosen next construction branch.
