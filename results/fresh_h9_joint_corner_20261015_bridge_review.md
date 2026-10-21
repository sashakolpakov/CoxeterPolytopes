# Independent review of the octic H9 bridge

The separate exact replay verifies the new bridge cap and the complete finite-vertex and finite-edge census of the resulting 33-wall chamber. It has **92 finite vertices, 253 finite edges and 322 open finite-corner directions**. The new wall joins two old finite-edge components, reducing their number from eight to seven. The inherited strict spacelike escape survives, so this remains an infinite-volume chamber.

The [independent checker](../scripts/fresh_h9_joint_corner_20261015_bridge_review.py) uses a newly constructed degree-eight SymPy primitive-element field and 384-bit outward rational intervals for signs. It imports no producer arithmetic or chamber engine. The [certificate](fresh_h9_joint_corner_20261015_bridge_review.json) binds the frozen [primary construction](fresh_h9_joint_corner_20261015_bridge.json), its source, and the [independently complete 32-wall predecessor](fresh_h9_joint_pair_20261015_two_review.md). Generation and saved replay both pass; the review is final and frozen.

## Exact cap and chamber checks

The checker reconstructs the particular solution and one-dimensional kernel from the nine prescribed physical contact equations. The norm equation independently yields

    delta = 1988340828/1990805 + (4448872824/9954025)sqrt(5).

The primitive-element minimal polynomial has degree eight. Both norm-two branches are checked; exactly one has a common strict center and admissible pairings with all 32 retained walls. It is right to `[0,1,2,5,6,10,29]`, has angle pi/5 at wall 8 and pi/3 at wall 30, and is strictly dotted to the remaining 23 walls. The review checks all 528 pair relations and direct timelike relative-interior witnesses for all 33 facets.

Every old finite ray is reconstructed from its nine incident walls and compared projectively with the saved predecessor. Its spherical Gram and all 33 chamber inequalities are checked. All 90 old finite vertices are strictly retained by the new cap.

## Exhaustive finite-vertex and edge calculation

A finite vertex on the new wall can use only its nine finite-angle neighbors. Testing all nine subsets of eight such neighbors gives exactly two positive-definite nine-wall supports:

    [0,1,2,5,6,8,10,29,32],
    [0,1,2,5,6,10,29,30,32].

Their physical nullspaces give future timelike rays satisfying every chamber inequality, with exactly the displayed incidences. Any new finite vertex lies on the added wall. A finite acute Coxeter corner is spherical and simple; equivalently, some rank-nine positive-definite incident subset must occur, and the two successful supports have strict inequalities at every other wall. The support enumeration therefore exhausts new finite vertices, including possible degeneracies.

Finite edges are independently rebuilt from pairs of finite vertices sharing eight walls, with the common Gram checked positive definite. Convexity ensures the segment is an actual edge. All 250 old finite edges remain, and precisely three new edges occur:

* old vertex 20 to new vertex 90;
* new vertex 90 to new vertex 91;
* new vertex 91 to the formerly isolated joint corner, old vertex 21.

Together with the old edge 19--20, these form the five-vertex path `19--20--90--91--21`. A separate union-find calculation checks the complete component partition. The equality `9*92-2*253=322` counts the remaining open finite-corner directions.

The exact inherited escape has positive square, negative pairing with the common future center and strictly negative pairing with all 33 walls. This proves infinite volume. Neither this review nor the new primary construction claims a complete exterior-ray or ideal-vertex census for the 33-wall chamber. The complete old cone is inherited from its frozen independent review; all geometry affected by the new wall within the finite-vertex/edge scope is rechecked here.

```sh
python scripts/fresh_h9_joint_corner_20261015_bridge_review.py
python scripts/fresh_h9_joint_corner_20261015_bridge_review.py --verify
```
