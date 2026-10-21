# Complete finite spectator families of five actual local quadrilaterals

This package certifies **682 compact local quadrilateral configurations in rank 13**, covering every finite orthogonal spectator forest for five previously certified cores. It identifies 1,852 corner forests, including 1,660 absent from the initial corner-orbit model.

| Local core | Required spectator rank | Complete spectator forests |
|---|---:|---:|
| 4202 | 7 | 196 |
| 4266 | 7 | 196 |
| 7746 | 5 | 47 |
| 7661 | 5 | 47 |
| B6 first-extension class 39 | 7 | 196 |

Only six contexts had all four corners represented in the [first repaired LP](fresh_sqrt5_rank13_joint_corner_repair1_20261017.md). Three still received zero payment: core 7661 with A1+A4 or A1+D4 spectators, and the B6 core with A2+D5 spectators. The other 676 contexts introduce at least one missing whole vertex type. The saved output includes every corner's exact marked-forest code and its predecessor payment when defined.

The [source](../scripts/fresh_sqrt5_rank13_five_core_spectators_20261017.py) and [saved certificate](fresh_sqrt5_rank13_five_core_spectators_20261017.json) independently recheck the physical and conjugate inertias, integral reflection coefficients, strict center, four finite future corner rays, and all six black-pair intersections for each fixed core. The inherited cores are [4202](fresh_sqrt5_adjacent_discount_witness_20261016.md), [4266 and 7746](fresh_sqrt5_e6_discount_witness_20261016.md), [7661](fresh_sqrt5_multiplicity_geometry_review_20261016.md), and [B6](fresh_sqrt5_global_b6_discount_witness_20261016.md).

Every allowed finite spectator component of rank at most seven is realized directly over K = Q(sqrt(5)). For a tree edge with squared normalized product t², choose a root scale s at one end and s t² at the child. Give each root norm 2s and the parent–child pair product −s t². The Cartan coefficients are then −t² and −1. All scales are totally positive, every normalized bond is checked, and the component's Gram is checked positive definite in both embeddings.

There is also a direct proof of conjugate positivity: after normalization the finite Gram is 2I−A, where A is nonnegative and has spectral radius below two. The nonphysical embedding only decreases the absolute 5- and 10-bond magnitudes; Perron monotonicity preserves strict positive definiteness. A separate narrow mathematical read confirmed this gauge and orthogonal-extension argument; it did not independently replay the census.

For each spectator forest, the full form is the orthogonal sum of its positive form with the core Lorentz form. Its Gram inertias are therefore (13,1,1) physically and (14,0,1) at the other embedding; all 225 Cartan entries are integral. The original four corner rays acquire zero spectator coordinates, and imposing all spectator walls recovers exactly the original compact quadrilateral. A strict timelike center is obtained by adding −epsilon S^(-1)1 in the spectator space. The source saves and exactly checks a positive dyadic epsilon with negative full norm for every context.

Generation and fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_five_core_spectators_20261017.py
python scripts/fresh_sqrt5_rank13_five_core_spectators_20261017.py --verify
```

The spectator enumeration is complete for these five fixed cores and the stated arithmetic finite-type domain. It is not a complete rank-13 face catalogue and does not assert compact finite-volume ambient realizations. The complete batch is used by the [next partial-model repair](fresh_sqrt5_rank13_joint_corner_spectator_repair_20261017.md), with variables and full vertex budgets added for every newly encountered forest.
