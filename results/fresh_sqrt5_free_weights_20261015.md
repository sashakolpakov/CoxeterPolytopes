# An optimistic free-weight diagnostic over Q(sqrt5)

The [exact LP certificate](fresh_sqrt5_free_weights_20261015.json) is feasible after allowing individual distance-two weights to decrease. For the 210 included face constraints, its exact optimal common component deficit is \(21/10\): the proposed weights satisfy \(w(C)\le2\operatorname{rank}(C)-21/10\) on the 28 component types represented in this finite relaxation. **This is not a valid global weight rule or a dimension bound.**

The [checker](../scripts/fresh_sqrt5_free_weights_20261015.py) assigns 266 independent nonnegative weights to connected marked-pair orbits. Each weight depends only on the entire connected component and marked pair; disconnected marked pairs have weight zero. Unlike the [additive barrier](fresh_sqrt5_sealed_weights_20261015.md), no distance-two baseline is fixed. The constraints comprise:

- All 112 sealed singular no-short quadrilateral cores, each requiring total weight at least one.
- All 76 three-black markings of the complete minimal Lannér catalogue, each requiring total weight at least two.
- All 22 literal E7/E8 triangle stars supplied by the frozen, independently reviewed quadratic-field one-root catalogue, again requiring total weight at least two.

For the last group, the new checker reconstructs each full Gram directly from its literal contacts. It checks physical Lorentz inertia, positive-definite conjugate Gram, the projected three-normal Gram, all three feasible timelike corner rays and a common strict future center. Thus these additional constraints are backed by actual compact local triangle faces. They do not classify every possible triangle with an E7/E8 corner in an arbitrary larger star.

The optimizer proposes rational primal and dual coefficients. `Fraction` arithmetic verifies every constraint, nonnegativity, dual feasibility and exact equality of the primal and dual objectives. Both generation and saved replay pass. The saved coefficients substantially reduce some distance-two weights, demonstrating why the fixed-baseline obstruction does not exclude more general component-local weights.

Four major coverage gaps remain. First, the complete quadrilateral union enumeration was restricted to faces with no marked distance-two corner. Quads discarded as automatically paid by the old distance-two baseline are **not** automatically paid after those weights decrease; their unions have not been enumerated here. Second, arbitrary finite white extensions of the minimal triangle stars impose further constraints. Third, the nonsingular and dotted quadrilateral frontiers are omitted. Finally, the represented 28 component types do not exhaust all finite connected diagrams up to rank 29, or all embeddings of the present marked subdiagrams into them.

A route to a genuine rank-29 certificate would need all those coverage obligations. One possible finite formulation takes a variable for every marked-pair orbit in every allowed spherical component of rank at most 29, and bounds a partial-star payment by every one of its induced finite extensions. This report records a feasible necessary relaxation and an explicit next target, not a proposed numerical improvement to the dimension ceiling.

```sh
python3 scripts/fresh_sqrt5_free_weights_20261015.py --verify
```
