# Independent exact check of the joint 27-facet H9 chamber

PASS. The [standalone checker](../scripts/fresh_h9_makarov_joint_20261013_exact_independent.py)
imports no producer arithmetic. Its [certificate](fresh_h9_makarov_joint_20261013_exact_independent.json)
binds the frozen [27-wall source](fresh_h9_makarov_joint_20261013_exact.json),
the old seed, and its own script by SHA-256.

Each cap is read as A+B sqrt(D), with A,B,D in Q(sqrt(5)) and D positive
at the geometric embedding. The checker implements its own rational-pair
arithmetic and nested quadratic sign comparisons. Mutual pairings involve
at most two additional radicals. The comparison also handles dependent
radicals and exact cancellations; no floating-point sign test is used.

It independently verifies:

* Lorentz signature (9,1), all 27 wall norms, and a common strict timelike center.
* All 351 wall pairs: 161 right angles, 31 angles pi/3, 12 angles pi/5,
  and 147 dotted pairs.
* Actuality of every facet, by projecting the strict center to its wall:
  mutual acuteness preserves all other strict inequalities and the point
  remains timelike.
* All ten old finite vertices, including their exact incident walls and
  positive-definite rank-nine Gram matrices.
* All 40 distinct new first finite endpoints. For each point the checker
  reconstructs its location on the declared old edge, checks its positive
  edge parameter, negative norm and future orientation, and evaluates all
  27 wall inequalities exactly. Its eight old incident walls have positive
  rank eight. The owner cap is independent because it is strict at the old
  vertex. The resulting rank-nine point is an actual finite vertex.
  Feasibility at both ends of the segment proves firstness by linearity.
* The strict integer spacelike escape and its orientation relative to the
  center. This proves infinite volume.

The independently checked endpoints give at least 50 finite vertices in
total. This is not an exhaustive census of the chamber's new vertices,
and the checker does not independently enumerate the search recipe bank.
The chamber has eleven seed walls plus sixteen new caps; it does not retain
the earlier twelfth wall.

```sh
python3 scripts/fresh_h9_makarov_joint_20261013_exact_independent.py --verify
```

The [separate proof review](fresh_h9_makarov_joint_20261013_exact_review.md)
checks the producer's reconstruction, endpoint and facet arguments.
No new finite-volume or compact example is established.
