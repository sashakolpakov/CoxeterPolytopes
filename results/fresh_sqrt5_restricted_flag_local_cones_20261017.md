# The new restricted flag dual survives the complete local simplex-packing cone

The new legacy-export restriction still gives an exact lower bound above the rank-13 target 12 after checking **every nonnegative charge depending only on a finite vertex forest and the orbit of a maximal-simplex port block**:

\[
 \frac{267632410797552528000000000000000}
 {21640812474569839372569203000501}
 =12.367022315453632\ldots>12.
\]

This is a method obstruction on the saved partial face domain, not an ambient dimension theorem or a realization of the dual as a polytope. It strengthens the previous finite-cut diagnostic: neither the missing old block features nor further nonnegative simplex-block charges at a single vertex remove this particular obstruction. Joint geometric restrictions involving simplex and prism configurations remain a separate avenue.

## Inputs and scope

The input is `legacy_export_seed` in the frozen [support diagnostic](fresh_sqrt5_flag_local_packing_dual_support_20261017.md). It has 2,219 positive face masses: 1,864 maximal-simplex or maximal-edge rows and 355 selected **actual** dotted-quad provenances. Its exact initial ratio is

\[
2676324107975525280/216106965733328291
=12.3842565597\ldots.
\]

All containing-simplex masses already belong to the declared legacy export; the support restriction restored no simplex row outside it. The new calculation leaves every face mass, including every actual quad provenance, unchanged and raises only vertex masses.

The [full flag model](fresh_sqrt5_flag_local_packing_full_model_20261017.md) is rebuilt and compared with its saved preparation. All 82,806 face obligations and 6,558 vertex budgets remain present. The initial exact dual is replayed against that model. After the local repair, a direct exact load recount checks all 333,829 positive-cost columns and all 94,481 signed columns, including the 77,391 coface-flag transfers. No inherited coefficient or constraint is removed.

The fresh [source](../scripts/fresh_sqrt5_restricted_flag_local_cones_20261017.py) and [certificate](fresh_sqrt5_restricted_flag_local_cones_20261017.json) bind the frozen support, model, block-reconstruction source and local packing routines. The saved model's recursive `prepared.sources` manifest remains an inherited dependency.

## Complete new feature reconstruction and local checks

The block reconstruction is rerun at the new face masses. In particular, **397 positive simplex rows were absent from the earlier 3,574-row full-block export**. They are included here rather than assigned zero features.

| Quantity | Exact count |
|---|---:|
| Positive simplex/edge rows | 1,864 |
| Positive actual quad rows | 355 |
| Contributing finite vertex forests | 593 |
| Occupied simplex-block orbits | 3,485 |
| Literal simplex-corner occurrences | 5,579 |
| Local packing count vectors | 136,200 |
| Positive rational mixture terms | 2,830 |
| Explicit disjoint literal blocks in mixtures | 12,278 |
| Matching exact local primal/dual optima | 524 |
| Exact lower/cover brackets | 69 |
| New violated local charge cuts | 31 |
| Vertex masses increased | 99 |
| Open or resource-limited forests | 0 |

The simplex row dimensions are 1:1,015, 2:209, 3:278 and 4:362. Dimension-one rows represent maximal edges and carry their two singleton port incidences, even though their curvature demand and pair-corner coefficients are zero.

For each forest, every root subset of a relevant size is classified by its exact marked weighted-forest code. This produces the complete literal automorphism orbit of each occupied block. The least-available-port recursion enumerates all disjoint packings, retaining the skip branch and every eligible block through that port. Four-bit orbit counts suffice because there are only 13 ports.

One bounded local LP is run per forest; no global LP is run. The resource guards are 1,000,000 cached count-vector entries, 30 seconds of enumeration, and 20 seconds of local optimization per forest. The largest observed vector list has 5,736 elements and the largest cache 69,667 entries. No guard fired.

Every numerical proposal is replaced by exact rational charges and explicit rational mixtures of disjoint literal packings. The 524 equal lower/upper values certify exact local optima. The remaining 69 cases are explicitly **brackets**, not claimed optima; their total gap is

\[
633608025823045753637/79200000000000000000000000000
<8.1\times10^{-9},
\]

and the largest individual gap is below \(6.1\times10^{-10}\). The upper cover certificates are sufficient for the stated obstruction regardless of these gaps. No certificate from the older fixed-alpha 740-forest calculation is reused.

## Why the certificate covers the whole nonnegative local charge family

Fix a finite forest and let \(w_j\) be the exact aggregate mass of simplex-corner blocks in orbit \(j\). For a disjoint literal packing \(P\), let \(a_j(P)\) count its blocks in orbit \(j\). The saved nonnegative rational mixture satisfies

\[
 \sum_P\lambda_P a_j(P)\ge w_j\quad\text{for every occupied orbit }j.
\]

For any nonnegative orbit charge \(q_j\), define its exact vertex cost by

\[
 c(q)=\max_{P\text{ disjoint}}\sum_j a_j(P)q_j.
\]

Then

\[
 \sum_jw_jq_j
 \le\sum_P\lambda_P\sum_ja_j(P)q_j
 \le\Bigl(\sum_P\lambda_P\Bigr)c(q).
\]

The repaired vertex mass is at least the mixture mass. Thus it satisfies **all** nonnegative local simplex-block charge inequalities simultaneously, not only the 31 saved cuts or the 175 finite cuts in the inherited model. Unoccupied block orbits have zero dual marginal; allowing charges on them cannot lower the maximum packing cost, so they do not weaken this conclusion.

For a literal charge that is not constant on an orbit, average each **whole packing** under the automorphism group of the forest. This gives uniform within-orbit marginals and the same argument. Symmetrizing a charge preserves its load against the uniform marginals and cannot increase its maximum packing cost. This is the same orbit-symmetrization argument as in the [complete earlier local-cone report](fresh_sqrt5_all_local_simplex_packing_cones_20261017.md), applied to newly reconstructed features and newly computed mixtures.

The disjoint-packing condition follows from the [maximal-simplex disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md). The calculation deliberately permits every disjoint packing, without asserting that all such packings occur geometrically. Stronger restrictions on their simultaneous realization, or joint simplex/prism conditions, are outside the tested cone.

## A transparent new separating charge

For vertex forest index 947,

\[
3A_1+A_3+A_5+I_2(5),
\]

the canonical port labels give a cost-one charge with the following nonzero block values:

| Literal block family | Charge |
|---|---:|
| \(\{3,4,11,12\},\{4,5,11,12\}\) | \(1/2\) |
| \(\{a,3,4,5\},\ a\in\{0,1,2\}\) | \(1\) |
| \(\{0,1,2\}\) | \(1/2\) |
| \(\{a,b\},\ a\in\{0,1,2\},\ b\in\{3,5\}\) | \(1/2\) |

Its exact maximum over disjoint packings is one. Its load is \(25588547/4000000000\), exceeding the input vertex mass \(58246233/10000000000\) by \(11450269/20000000000\). This is one of 31 genuine improvements, but all improvements together leave the repaired global ratio above 12.

## Verification and artifact status

Generation session 55730 and fresh saved replay session 70710 both terminated with exit 0 and `PASS_RESTRICTED_FLAG_DUAL_COMPLETE_LOCAL_CONES`. Replay reconstructed the full model, every feature and literal orbit, every packing count-vector set, every exact charge and mixture, and every inherited column load; it called no numerical optimizer.

The source hash is `77e346f6f0761ddd5e2e5d14e62f611c9078a24ef7dd13930d5fa84e6d532b2d`; the certificate hash is `43838b04b5ab517f057c8a64f43f9327e18a137f704bd6a99d0bebadb5d48cd5`.

The conditional second test, retaining all simplex rows and the same selected actual quads, was **not run**: the first restricted witness already survives the full local charge family. No new optimizer run or dimension-bound claim follows from this package.
