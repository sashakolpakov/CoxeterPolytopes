# Releasing old150 opens its triangular cavity

Generation and a fresh saved replay pass. This package is a **bounded numerical boundary diagnostic**, not an exact replacement chamber.

The input is the173-wall two-cap branch0/77 from the [outer-patch pairs](fresh_h20_replacement_outer_patch_pairs_20261017.md). Its target-H3 boundary has21 finite and9 exterior positive-side-polytope vertices,10 open hyperbolic edges, and maximum squared radius1.2186340819392139. Old150 is a compact triangular facet on142,145,152; its only nonright T-contact isπ/8 to151.

Deleting150 loses its three finite vertices: the result has18 finite/10 exterior vertices and12 open edges. In particular142/145 and142/152 become open, while145/152 becomes an edge with two exposed ideal ends. Adding any of the18 continuation caps from the [belt step](fresh_h20_replacement_belt_step_20261017.md) whose sole failed native gate was150 gives20 finite/10 exterior vertices and13 open edges. The best maximum squared radius falls to1.111867538677528, but the triangular cavity remains open and142 loses compactness.

Every one of those18 caps is rechecked against all172 retained native walls before the complete numerical target boundary is recomputed. Deletion changes the chamber: no persistence of vertices involving150 is asserted. These results support replacing the triangular facet, rather than promoting deletion alone.

The [source](../scripts/fresh_h20_replacement150_boundary_20261017.py) and [data](fresh_h20_replacement150_boundary_20261017.json) retain every candidate recipe, native-label list, and full boundary ledger. The inherited ledger distinguishes exterior Euclidean vertices from genuine hyperbolic vertices, and clips feasible pair lines by all retained halfspaces and the unit ball. Neither floating gates nor the18-input set give an exact or all-parameter theorem.
