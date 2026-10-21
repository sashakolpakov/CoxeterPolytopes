# Complete edge star of the 813-wall replacement

The new finite vertex has **22 outgoing edges, all open**. None has another existing finite or ideal endpoint. More strongly, **no new acute-compatible Coxeter wall that cuts the saved strict escape can end any of these 22 edges**, at either a finite or ideal point. This second statement covers every angle label and real normal, not a bounded numerical root bank.

These are local results for the [813-wall changed core](fresh_h22_core_replacement_20261013_review.md). They do not exclude intermediate caps that leave the saved escape untouched, later steps from new vertices, caps elsewhere, or further wall replacements. No finite-volume example is claimed.

## The complete current star

The actual finite vertex has norm `−29716` and exactly 22 incident walls: the 20 retained old cusp walls I and the two new mutually orthogonal caps. Its positive-definite incident Gram is `G_I ⊕ 2I2`, so it is a simple vertex and releasing each incident wall exhausts its outgoing edge directions.

The [exact star checker](../scripts/fresh_h22_replacement_edges_20261013.py) reconstructs this Gram and its rational inverse. If d_j is the dual to the released wall, an edge has parametrization

`v - t d_j`, with `t>=0`, `v²=-29716`, `(v,d_j)=0`.

The exact light-cone parameter is `sqrt(29716/d_j²)`. At every one of the 22 resulting null points, precisely the 21 retained edge walls have zero pairing; all **792** other wall pairings are strictly negative. The checker then increases the parameter by an explicit rational factor and verifies a future spacelike point with the same 21 zero pairings and every other inequality still strict. Thus these boundary points are not bounded ideal vertices: the edge continues into an open escape region.

There is also a structural explanation. Every edge retains at least one new cap. Every old wall outside I is strictly dotted to both new caps, so its hyperplane cannot meet a retained cap hyperplane even at the ideal boundary. The other initially incident wall meets the outgoing edge only at its starting vertex. Hence no existing wall can terminate the edge. Each new cap facet is a cone based at the one finite vertex; its 21 finite-neighbor walls have positive-definite Gram.

The [star certificate](fresh_h22_replacement_edges_20261013.json) records all 22 supports, duals, null points, and exact spacelike continuations. The current counts are zero other finite endpoints, zero other ideal endpoints, and 22 open directions. This is the entire star of this vertex, not an exhaustive enumeration of all chamber vertices.

## Why a direct endpoint cutter cannot remove the saved escape

The saved strict spacelike point E has norm **1586**. The [positive-cone certificate](fresh_h22_replacement_escape_dual_20261013.json) gives two literal rational vector identities of the forms

`E = sum positive coefficients * retained walls - (161/2) W7`,

`E = sum positive coefficients * retained walls - 50 W168`.

For any proposed normal r acute to every retained wall, `(r,E)>0` therefore forces **both** `(r,W7)<0` and `(r,W168)<0`. The [endpoint checker](../scripts/fresh_h22_replacement_edges_20261013_endpoint_gate.py) independently reconstructs all vectors and replays both identities; it does not trust an LP objective or rounded coefficients.

Normalize the proposed root to norm two. A nonright Coxeter contact with a wall of norm n has pairing magnitude at least `sqrt(n/2)`, because the smallest nonzero Coxeter cosine is `cos(pi/3)=1/2`. Parallel and dotted contacts only increase that lower bound.

For an outgoing edge support with Gram G, a finite or ideal endpoint involving the new wall requires the enlarged Gram to be positive semidefinite. Consequently its edge pairing vector p must satisfy

`p^T G^-1 p <= 2`.

The checker verifies that every G inverse is entrywise nonnegative. It then retains only the mandatory contacts with W7 and W168 that remain in the support, assigning their smallest possible nonright magnitudes and zero to all other entries. Acuteness and the nonnegative inverse make this a lower bound for every possible profile, including additional nonright contacts.

For **all 22 supports**, the resulting projected norm is strictly greater than two. The smallest lower bound is **3**, attained when W7 is released. Releasing W168 gives **15/4**; releasing either new cap gives **18+10sqrt(2)**. The [exact endpoint certificate](fresh_h22_replacement_edges_20261013_endpoint_gate.json) contains every bound and independently verifies the inverse identities.

Thus the enlarged Gram cannot even be positive semidefinite. A cap cutting E cannot meet any current outgoing edge at a finite or ideal endpoint. This includes arbitrary finite Coxeter labels and parallel limits; it has no arithmetic-field assumption.

The conclusion does not say that no wall can end these edges. A wall may end an edge while leaving E feasible. Subsequent walls based at the resulting new vertices are outside the proved obstruction. A finite-volume completion would eventually need some wall cutting E, but this proof does not control every possible sequence of intermediate additions.

## Search that led to the exact test

A small [numerical continuous screen](fresh_h22_replacement_edges_20261013_target_scout.json) first examined 14 relevant old-wall releases, keeping the other 19 old incident walls and allowing nonright contact with a new cap. Every relaxed LP had a negative maximum on E. No finite-label quadratic search was executed after those screens rejected the families. This numerical observation is superseded by the exact two-contact Schur argument above, which covers all 22 edges and additional nonright contacts.

## Replay

```sh
python scripts/fresh_h22_replacement_edges_20261013.py --verify
python scripts/fresh_h22_replacement_edges_20261013_endpoint_gate.py --verify
```

Both pass. The first checks all edge null points and outward continuations with exact sparse-radical arithmetic. The second uses exact rational Gram inverses and radical signs and independently replays the mandatory-contact vector identities. Neither changes the frozen 813-wall construction.
