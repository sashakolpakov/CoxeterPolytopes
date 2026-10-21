# Conditional exclusion of compact dimension23 by farthest forests

20 October 2026. Standalone paper composition, **held pending root's
acceptance of this final assembly**. The new component sources below now
have full root and outside-author written-proof PASS. No new bound is
promoted by this candidate.
No scientific, Lean or TeX execution or canonical manuscript edit.

## 1. Exact mathematical inputs

All paths below are under `results/`. The accepted baseline
`compact_dimension23_farthest_forest_composition_20261020.md`, SHA256
`12d9112ded91d51a209051a5d4a10122349a158f06cc5cf080bc6724e74bcad4`,
already excludes every dimension at least24, using its explicit published
compact ceiling29 and internally proved finite-range tail. We only need
exclude dimension23; no face-inheritance argument is used.

The component inputs concern projections of the SAME actual vector x at
a vertex farthest from a strict interior point, with `||x||^2<1` and
additive squared norms over orthogonal whole finite components.

| Internal input | Source and SHA256 |
| --- | --- |
| Complete ordinary small-type rates, short-A1/rank-two exceptions, B/D>=9 energy>9/14 and component caps A<=20, B<=22, D<=23 | `compact_forest_ordinary_a_parameter_and_dominant_classical_frontier_20261020.md`, `dd6f8b0de2d7d0bc835be3eb4e320013aa40792407fb2a747f35bcf04385d5ad` |
| A_r>=7 energy>1/2 | `compact_large_an_actual_half_uniform_transfer_20261020.md`, `7351b471ba3b0edaf0a524ccf082010e458b7fedd8d063491a75632223bd70e7` |
| A_r>=17 energy>5/7 | `compact_large_an_actual_projection_five_sevenths_20261020.md`, `7a57b66e1b9b68948c6fe14c493527ff3ac0b963144f904ff512ead723e192b0` |
| Localized B/D prefix floors with arbitrary spectators | `compact_classical_generic_prefix_energy_20261020.md`, `f637d7294d17e37581330ca508c914c23a08ce259845614b055e45815d2790ef` |
| A_r>=8 energy>4/7 | `compact_large_an_actual_four_sevenths_uniform_transfer_candidate_20261020.md`, `390bbf35828e490fe05918af404c14ba946d1cddd0deed325c93d0f36522d1a5` |
| A19/A20 energy>6/7, arbitrary spectators | `compact_a19_a20_central_window_actual_six_sevenths_20261020.md`, `8075f92f4766d681cdef42040925586bcba76ebdd99c059ce5d75a706dc3ee0c` |
| Whole D23 exclusion; B22/D22 exclusion with isolated short I | `compact_last_d23_center_and_isolated_i_frontier_candidate_20261020.md`, `a77f88f5bd2c2c02588dffae7be4bc9ad96e83421bea2d525a3911c4659d996c` |

These are accepted internal proved lemmas, not external assumptions of
the desired dimension conclusion. The last D23 source
contains the rational-premise clarification at lines244–246. Its original
53986c86 edition and outside review remain preserved.

## 2. Dominance and the complete spectator budget

Suppose a compact23-polytope exists. Call A_r with r>=7 or B_r/D_r with
r>=9 dominant. Without a dominant, finite spherical classification and
the ordinary rates give total energy greater than22/14 after discarding
the possible one short A1; alternatively one isolated rank-two7 defect
gives energy greater than9/100+21/14. Both exceed1. All finite types,
including E/F/H and arbitrary I2(m), are covered by the small-type input.

Thus a dominant exists. Every dominant pays strictly more than1/2, so
it is unique. Only now may all its spectators be declared ordinary
apart from the following incompatible alternatives. Write S=23-r.

| Spectator case | Sufficient lower budget |
| --- | --- |
| No exception | S/14 |
| One short A1, no rank-two defect | (S-1)/14, discarding that one rank |
| One rank-two7 defect, no short A1 | 9/100+(S-2)/14 |

Both rank-two defects are included: old A2/new7 and old I2(7)/new3.
Their isolated actual walls and the native product gap permit at most
one defect and forbid coexistence with either short-A1 type I or J.
J needs two A1 components; its light direction is nonshort. It may be
handled conservatively by the single short-rank discard above, without
charging its light rank twice. For S>=2 the rank-two budget is larger
than the short-rank budget, since9/100-1/14=13/700>0.

## 3. Every dominant A rank is excluded

The half bound and the spectator table first force S<=7, hence r>=16;
at S>=8 even the weakest short-rank budget gives1/2+7/14=1, strictly
exceeded. The A cap then leaves exactly r16–20.

| Rank | Component bound | S | Lower total after short-rank discard |
| --- | --- | --- | --- |
| 16 | >4/7 | 7 | >4/7+6/14=1 |
| 17 | >5/7 | 6 | >5/7+5/14=15/14 |
| 18 | >5/7 | 5 | >5/7+4/14=1 |
| 19 | >6/7 | 4 | >6/7+3/14=15/14 |
| 20 | >6/7 | 3 | >6/7+2/14=1 |

The no-exception and rank-two-defect budgets are stronger in every row.
Every lower A dominant rank was excluded before this table; every higher
one is excluded by the actual arbitrary-spectator A cap. No A rank remains.

## 4. Every dominant B/D rank is excluded

The9/14 bound and the spectator budgets first force S<=5, hence r>=18.
The complete generic prefix table applies with M=r-6 for both types:
its D rank requirement is r>=M+6 and its B requirement is weaker.

| Rank | M | Actual component energy above | S | Weakest coarse spectator budget |
| --- | --- | --- | --- | --- |
| 18 | 12 | .8689525 | 5 | 4/14 |
| 19 | 13 | .88864 | 4 | 3/14 |
| 20 | 14 | .90715 | 3 | 2/14 |

Every row has total energy strictly above1, even after either debit.
All decimals are exact rationals from the accepted prefix theorem.

At r21, M15 gives energy>.9240775 and S2 remains. A rank-two defect
costs>.09, giving total>1.0140775. Otherwise either an ordinary rank-two
component pays>1/7, or the two components are A1. At most one A1 is short,
so the other pays>1/10 and the total exceeds1.0240775. This also handles
J; its two-rank support cannot create a second short direction.

The stronger ordinary-A1 bound is essential here. Its threshold is
`L_*=asinh(1/sqrt(2))`; even equality gives half-tangent square
`(sqrt(3)-sqrt(2))^2=5-2sqrt(6)>1/10`, since `sqrt(6)<49/20`.

At r22, M16 gives energy>.9401975 and S1 remains. A nonshort A1 pays
>1/10, contradicting the norm. A rank-two defect or J cannot fit in one
rank. The only remaining spectator is exactly isolated short I, excluded
for both B22 and D22 by the new actual coupling theorem.

At r23 there are no spectators. B23 is excluded by the accepted B cap;
the new whole-D23 theorem has exactly the required no-spectator scope.
Higher B/D ranks cannot fit dimension23. Thus no dominant B/D remains.

## 5. Conditional conclusion and limits

Every finite farthest forest has been excluded: either there is no
dominant, or the unique dominant is A, B or D, with every possible rank
and every permitted debit covered above. Conditional on root's acceptance
of the bound component sources and this composition, dimension23 is
impossible. Together with the accepted exclusion of all dimensions>=24,
this would give general compact dimension at most22.

This note does not itself promote that bound, edit the canonical compact23
paper, change arithmetic compact20, or claim a new example or sharpness.
No arbitrary metric face is assumed Coxeter and no new formal verification
or finite-data execution is claimed.
