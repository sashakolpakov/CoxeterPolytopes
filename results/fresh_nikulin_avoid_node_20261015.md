# Avoiding the narrow-facet node at a cost of 2/9

**Uniform local lemma.** In every irreducible affine Coxeter factor of rank r, and for every forbidden node s, one can prescribe two omissions i,j different from s, before either terminating contact row is known, such that

    σ τ E / M + 2C ≥ 2(r−1)/9.                    (1)

Here the roots have square two, C is the sum of the two retained-support Schur costs, σ,τ are the positive affine null heights, M=max(σ/τ,τ/σ), and E is the normalized-row energy. The rows satisfy the actual endpoint hypotheses from the frozen product proof: they are nonnegative, every nonzero entry is at least one, and each endpoint's total Schur cost is at most two. No upper contact or Coxeter-label bound is imposed.

Together with the reviewed non-narrow estimate, this gives a **research proof candidate n≤181** for finite-sided finite-volume real hyperbolic Coxeter polytopes. The new lemma is not yet independently reviewed. The [189 checkpoint](fresh_nikulin_dimension_189_candidate_20261015.md) stays frozen and remains the preceding reviewed result. No novelty, publication-priority, or current-best claim is made.

The [checker](../scripts/fresh_nikulin_avoid_node_20261015.py) and [data](fresh_nikulin_avoid_node_20261015.json) rebuild all 34 exceptional forbidden-node cases, the 16 modified classical cases in ranks 2 through 9, the three E8 inexpensive-contact identities, and the all-rank scalar comparisons. The unbounded-rank localization and geometric framing inputs are explicitly inherited from the frozen proofs and their reviews; finite computations are not substituted for those arguments.

## Definitions and fixed choices

For an affine Gram matrix G with positive kernel δ, use the normalization min δ=1. For the two prescribed omissions i,j let S_i,S_j be the retained supports and put

    C = p_Siᵀ G(S_i)⁻¹ p_Si + q_Sjᵀ G(S_j)⁻¹ q_Sj,
    σ = δ·p,   τ = δ·q,   M=max(σ/τ,τ/σ).

For normalized point rows define

    R_ab = (e_a/δ_a−e_b/δ_b)ᵀ G⁺ (e_a/δ_a−e_b/δ_b).

The normalized rows p/σ and q/τ have difference of energy E. In particular,

    σ τ E/M = min(σ,τ)² E.                       (2)

The inverse retained Grams are entrywise nonnegative. A retained nonzero entry therefore costs at least the corresponding inverse diagonal. A terminating row cannot vanish on a whole affine factor, since all affine null sums are positive multiples of the initial cusp ray and its null height is positive.

Use these fixed rules in the single factor containing the forbidden narrow-facet node:

* A_r, r≥2: rotate an opposite pair of cycle nodes to avoid s. The existing estimate is rotation invariant.
* D_r, r≥4: use leaves from opposite fork ends, switching to the other leaf if necessary to avoid s.
* B_r: use a left fork leaf and the right short end. Switch left leaves if necessary. If the short end is forbidden, replace that omission by its adjacent interior node.
* C_r: use both ends, replacing a forbidden end by its adjacent interior node. The two endpoint cases are diagram reversals.
* A_1: both omissions are the only node other than s; coincident omissions are allowed in this factor.
* E6: choose the first two of leaves 2,4,6 not equal to s.
* E7: choose 4,7, except when s is 4 or 7, when the pair is 2,5.
* E8: choose 3,7, except when s is 3 or 7, when the pair is 1,6.
* F4: choose 1,3, except when s is 1 or 3, when the pair is 0,4.
* G2: choose the two nodes other than s.

The exact exceptional node numbering and complete Gram matrices are serialized. For E6/E7/E8 they are stars with arm lengths respectively (2,2,2), (1,3,3), and (1,2,5), numbered consecutively from the center along each arm. The rules do not depend on the subsequent contact rows, their zero sets, labels, or heights.

## Unmodified classical cases

The [reviewed contact-cost lemma](fresh_nikulin_energy_charge_20261015_review.md) already gives the stronger inequality

    σ τ E/M + 2C ≥ 2r/9

for the ordinary A/B/C/D choices. A-cycle rotation and D-fork switching preserve those diagrams and estimates. Thus only the modified B and C endpoint choices require new work. This includes the separately checked A9 localization in the earlier review; no unstated large-rank threshold is being applied to rank nine.

## Modified B and C choices

Consider B_r when the right short node is forbidden, and C_r when the left end is forbidden. The other C case is obtained by reversing its entire diagram. In the modified support the forbidden end becomes an isolated A1. Its inverse diagonal is 1/2.

When C=0 the rows are point contacts at their omitted nodes. The weighted point energy is (r−1)/4 for modified B_r with r≥3, and (r−1)/2 for modified C_r. Modified B2 instead has point energy 1. All these values are at least 2(r−1)/9.

For ranks 2 through 5, every retained inverse diagonal is at least 1/2. Consequently C>0 implies 2C≥1≥2(r−1)/9, proving (1) without using energy.

For ranks 6 through 9, the isolated A1 is the only inverse diagonal below one in either deletion support. If C≥(r−1)/9, its charge already pays (1). Otherwise 0<C<(r−1)/9<1, and exactly one retained entry can be nonzero across the two rows. It is the isolated A1 entry, of magnitude b≥1 and cost b²/2. The other row is supported only at its omitted node. With i the modified omission, l the isolated node and j the opposite omission, write

    p=a e_i+b e_l,   q=c e_j,   u=δ_l b/(δ_i a+δ_l b).

Here a≥0, c≥1 and 0≤u≤1. The complete inverse-diagonal lists and the following three resistances are checked exactly for all eight cases:

| Type | R_ij | R_lj | R_il | Resulting energy |
|---|---:|---:|---:|---:|
| B_r | (r−1)/4 | r/4 | 1/4 | (r−1+u²)/4 |
| C_r | (r−1)/2 | r/2 | 1/2 | (r−1+u²)/2 |

The energy formula is the squared-distance interpolation identity

    E=(1−u)R_ij+uR_lj−u(1−u)R_il.

Since min δ=1, every nonzero row has height at least one. Equation (2) therefore bounds the weighted energy below by (r−1)/4 for B, or (r−1)/2 for C. Either is at least 2(r−1)/9, even without adding its positive contact cost.

For every r≥10, the inherited fixed-choice localization gives weighted energy at least (r−5)/4 in modified B and (r−5)/2 in modified C. The closed finite B/D inverse formulas show that the isolated A1 still has the smallest diagonal, 1/2; all other retained inverse diagonals are at least one. Thus C>0 implies 2C≥1. At r=10+t, t≥0, the exact margins are

    (r−5)/4 +1−2(r−1)/9 = (t+9)/36 >0,
    (r−5)/2 +1−2(r−1)/9 = (5t+27)/18 >0.

The C=0 point-contact calculation already covers every rank. This proves (1) for all modified classical cases, without enumerating contact magnitudes or imposing a rank cutoff.

## Direct exceptional cases

For A1,E6,E7,F4,G2 the prescribed pair has both point energy at least the target 2(r−1)/9 and twice the minimum retained inverse diagonal at least that target. The C=0 and C>0 cases are therefore immediate. The records cover every forbidden node, not only diagram-orbit representatives.

The E7 replacement is useful here: pair 4,7 has point energy 3/2 and minimum inverse diagonal 3/2; pair 2,5 has the same point energy and minimum diagonal 2/3. Both meet the new target 4/3. No cheap-contact subcase is needed in E7.

For A1 the target is zero. Its two omissions coincide; two rows supported at this common node have E=C=0. Hence the old target 2r/9 cannot hold uniformly under an avoid-node constraint, even at the level of admissible local contact rows. The additive loss 2/9 in (1) is therefore sharp for this uniform local accounting. This observation does not preclude a different global argument from improving the eventual dimension bound.

## E8 and its replacement omissions

Here the target is 14/9. If C≥7/9 the charge alone suffices. In the remaining case 0<C<7/9, the exact inverse diagonals force a single cheap retained contact.

For pair 3,7 this is node 8 after deleting 7, with cost b²/2. Its null coefficients are δ_i=δ_j=2, δ_l=1 and its resistances are (1/2,1,1/2). This is the previously reviewed E8 argument: if a=0, σ τ E/M≥1 and 2C≥1; if a≥1, both heights are at least two and E=(1+u²)/2≥1/2, so σ τ E/M≥2. Both cases exceed 14/9.

When s is 3 or 7 the new omissions are 1,6. Their null coefficients both equal three. The only inverse diagonals below 7/9 are nodes 7 and 8 after deleting 6, each equal to 2/3. Every other retained diagonal is at least 8/9. Two nonzero retained contacts would cost at least 4/3, so only one can occur. Its row is p=a e_6+b e_l, while the other is q=c e_1, with c≥1.

| Cheap node l | δ_l | R_61 | R_l1 | R_6l | E |
|---|---:|---:|---:|---:|---:|
| 7 | 2 | 2/9 | 7/18 | 1/6 | (4+3u²)/18 |
| 8 | 1 | 2/9 | 8/9 | 2/3 | 2(1+3u²)/9 |

If a≥1, both heights are at least three and E≥2/9. Equation (2) gives σ τ E/M≥2.

If a=0 and l=7, the heights are at least two and three. The weighted energy is at least 4·7/18=14/9, as required. If a=0 and l=8, the energy is at least 8/9 and the cost contributes 2C≥4/3, giving at least 20/9. These are all cases. The only discreteness used is that a is either zero or at least one.

## Consequence for the narrow-facet dichotomy

Use the unchanged point/contact-cost omissions in every non-narrow factor and the above avoid-node choices in the single narrow factor. The [frozen product argument](fresh_nikulin_narrow_20261015_products.md) then supplies actual terminating framing walls x,y. In particular, both outgoing edges retain the narrow node; at an ideal terminal point its positive-definite edge support omits one node per affine component, and a suitable missing wall has a finite angle with the narrow wall. The strict-center projection argument makes their intersection an actual ridge. No merely ideal-point incidence is substituted for Nikulin's framing condition.

Writing N=n−1, summing (1) and the stronger non-narrow estimate gives

    Σ_a σ_a τ_a E_a ≥ M(2N/9−2Σ_a C_a−2/9)
                       ≥ M(2N/9−8−2/9).

The two endpoint Schur budgets still give Σ_a C_a≤4, now including the narrow factor's own cost. The common height ratio and the Lorentz identity give

    −(x,y) = (1/2)Σ_a σ_a τ_a E_a−q−1/q
            ≥ M(N−55)/9.

For n≥182, N≥181 and this is at least 14, contradicting the strict narrow-facet distance bound arccosh 7 for square-two normals. Thus the ideal-vertex branch gives n≤181. If the narrow facet has no ideal vertex, it is compact and the independently reviewed compact-facet bound n≤143 applies. The resulting candidate bound is n≤181.

Only the avoid-node accounting is new. The published narrow-facet theorem, affine classification, common-height identity, endpoint existence and actual-ridge bridge, and compact-facet proof remain the explicitly bound dependencies of the [reviewed 189 checkpoint](fresh_nikulin_dimension_189_candidate_20261015.md). The new lemma requires its own independent review before 181 replaces that checkpoint.

```sh
python scripts/fresh_nikulin_avoid_node_20261015.py
python scripts/fresh_nikulin_avoid_node_20261015.py --verify
```
