# Complete vertex-type coverage audit and two actual unpaid local quadrilaterals

The first rank-13 corner-orbit rule fails two actual compact local quadrilaterals, with **exact payment zero** instead of the required one. Every corner already belongs to one of the rule's 630 represented whole vertex types.

| Frozen local core | Orthogonal spectator forest | Rank-13 payment |
|---|---|---:|
| 4266 | A2+D5 | 0 |
| 7746 | D5 | 0 |

The [source](../scripts/fresh_sqrt5_rank13_joint_corner_separation_20261017.py) and [saved certificate](fresh_sqrt5_rank13_joint_corner_separation_20261017.json) reconstruct the frozen [E6 local quadrilaterals](fresh_sqrt5_e6_discount_witness_20261016.md), add the displayed spherical blocks, and check the extended geometry exactly. Each 15-normal Gram has physical inertia (13,1,1) and conjugate inertia (14,0,1). A 14-dimensional basis form is Lorentzian physically and positive definite at the other embedding. All 225 Cartan coefficients are integral. A strict timelike center, all four finite future corner rays, and both orientations of all six possible black-wall intersections are checked. Exactly the four stipulated corners are feasible.

This proves local arithmetic and geometric admissibility. It does not claim a compact finite-volume ambient polytope. Orthogonal spherical spectators leave the compact quadrilateral in the common white-wall face unchanged.

The whole-vertex coverage gap is also exact. The complete finite-component list over the Q(sqrt(5)) squared-label domain consists of A, B, D, E, F4, H3, H4 and I2(5), I2(6), I2(10), with A2 and B2 representing I2(3) and I2(4). Enumerating all component multisets of rank 13 gives:

| Coverage quantity | Count |
|---|---:|
| All finite forests | 6,789 |
| Forests represented in the partial model | 630 |
| Absent forests | 6,159 |
| Absent forests with baseline at least 12 | 2,794 |
| Absent forests with baseline strictly above 12 | 2,174 |

Among the 21,492 represented pair-orbit variables, only 3,634 occur in a saved face constraint. The rule uses 839 nonzero variables, discounts 5,844 orbits relative to the baseline, and increases 648. Of the discounted orbits, 4,310 occur in no saved face constraint. Every represented whole forest has a discounted orbit. The saved output includes the full finite-forest census and all changed orbit indices, multiplicities and exact payments.

Generation and fresh saved replay passed:

```sh
python scripts/fresh_sqrt5_rank13_joint_corner_separation_20261017.py
python scripts/fresh_sqrt5_rank13_joint_corner_separation_20261017.py --verify
```

The [original LP certificate](fresh_sqrt5_rank13_joint_corner_orbit_lp_20261017.md) remains valid on its stated partial domain. The [first repaired model](fresh_sqrt5_rank13_joint_corner_repair1_20261017.md) adds these two exact face constraints. Neither result establishes complete face coverage or a dimension bound.
