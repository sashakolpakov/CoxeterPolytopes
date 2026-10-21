# Unit charge fails; a uniform charge of 4/3 gives a 169 candidate

For the affine-factor energy estimate, coefficient one is impossible even after changing the prescribed omission pair. A coefficient of **4/3** does work uniformly, with a different non-narrow E8 pair and the existing avoid-node choices. It is optimal for this local lemma with slope 2/9 and omissions prescribed before the contact rows. This is a statement about that particular accounting scheme, not an optimal dimension bound.

The resulting **research proof candidate is n≤169** for finite-sided finite-volume real hyperbolic Coxeter polytopes. The [reviewed 181 candidate](fresh_nikulin_dimension_181_candidate_20261015.md) remains frozen. The new argument requires independent review before it becomes the reviewed baseline. No priority or current-best claim is made.

The [checker](../scripts/fresh_nikulin_unit_charge_20261015.py) and [certificate](fresh_nikulin_unit_charge_20261015.json) check all exceptional cases, complete E6/E7/E8 omission-pair obstructions, and symbolic all-rank classical formulas. They use no finite-label cutoff and import no producer arithmetic.

## Precise local statement

Retain the notation and geometric hypotheses of the [product proof](fresh_nikulin_narrow_20261015_products.md) and [avoid-node proof](fresh_nikulin_avoid_node_20261015.md). The affine roots have square two; the positive kernel δ has minimum entry one. For prescribed omissions i,j, contact rows p,q are nonnegative and nonzero, and every nonzero entry is at least one. Each endpoint's full Schur cost is at most two. Write

    σ=δ·p,  τ=δ·q,  M=max(σ/τ,τ/σ),
    C=p_SiᵀG(S_i)⁻¹p_Si + q_SjᵀG(S_j)⁻¹q_Sj,
    W=σ τ E/M=min(σ,τ)² E.

The new inequalities are

    W+(4/3)C ≥ 2r/9                         (ordinary factor),
    W+(4/3)C ≥ 2(r−1)/9                     (factor avoiding a fixed node).

Use the frozen ordinary choices except that E8 now uses nodes **1,6**. Use exactly the frozen avoid-node choices in the factor containing the narrow-facet normal. Thus every choice remains independent of the eventual contact rows. In particular, the new ordinary E8 choice does not replace the separate rule needed when node 1 or 6 must be avoided.

## Exact obstruction to coefficient one

With the old ordinary E8 omissions 3,7, take the row at omission 7 to be p=e8 and the other row q=e3. Their heights are 1,2, the total cost is 1/2, and W=1. Consequently W+C=3/2<16/9; this pair would require a charge at least 14/9.

The certificate checks every unordered omission pair, not just the old choice. A pair whose zero-cost point rows have energy below 2r/9 is ineligible for every coefficient. Among the remaining pairs it tests all unit point-contact rows whose two individual costs are at most two. These are legal rows under the stated local hypotheses and use only contact magnitudes zero and one.

| Type | Eligible omission pairs | Least necessary coefficient among them |
|---|---:|---:|
| E6 | 6 of 21 | 1 |
| E7 | 1 of 28 | 10/9 |
| E8 | 4 of 36 | 4/3 |

The four E8 pairs are `(1,6)`, `(2,6)`, `(2,5)` and `(3,7)`, with respective necessary coefficients `4/3`, `65/48`, `14/9` and `14/9`. For the best pair `(1,6)`, choose p=e1 and q=e8. The heights are 3,1, the costs are 0 and 2/3, and W=8/9. Hence `(16/9−W)/C=4/3` exactly. The uniform upper bound below therefore attains the best possible coefficient within this fixed-slope, preset-two-omission local statement. No claim that these local obstructions alone impose a global dimension lower bound is made.

## A linear loss identity

Choose Euclidean point potentials v_k for the normalized point rows e_k/δ_k, defined up to a common translation. Put R_ab=|v_a−v_b|². The contact-row potentials are convex combinations

    z_p=Σ_k (δ_k p_k/σ)v_k,   z_q=Σ_k (δ_k q_k/τ)v_k.

Set a=z_p−v_i, b=z_q−v_j and

    ℓ_ik=R_ij+R_ik−R_kj.

Direct expansion gives

    E=R_ij−Σ_k (δ_k p_k/σ)ℓ_ik
           −Σ_k (δ_k q_k/τ)ℓ_jk+|a−b|².          (1)

The inverse diagonal for a retained node k after omission i is

    c_ik=δ_k² R_ik.

All inverse entries are nonnegative, so the row cost is at least Σ_k c_ik p_k². Every convex weight δ_k p_k/σ is at most one; when it is nonzero, p_k²≥1. Therefore, if every possible retained contact satisfies `ℓ_ik≤α c_ik` (and similarly for the other row), each linear loss in (1) is paid by α times its diagonal cost. Negative ℓ terms only improve the bound. It follows that

    E+αC≥R_ij,  hence W+αC≥R_ij,                (2)

because σ,τ≥1 and E≥0. This elementary identity uses the energy and cost together, instead of discarding all energy when a contact occurs.

For each factor with target T, the case C≥T/α is immediate. In the remaining case every retained contact has diagonal `c_ik<T/α` and `c_ik≤2`. These restrictions eliminate the exceptional distant contacts in the following coordinate tables.

## All classical ranks

For A_r let m=r+1, take opposite cycle nodes, and let ε=m mod 2. A node at cyclic distance d from an omission has

    c=d(m−d)/m,    ℓ≤(m+ε)d/m.

For r≤5 the minimum retained diagonal r/(r+1) is already at least the small-cost threshold r/6, so a nonzero retained contact is impossible below that threshold. For r=6,7,8, the threshold excludes distance two and farther; only d=1 remains. For r≥9, the individual endpoint budget two excludes d≥3, since `3(m−3)/m>2` for m≥10 and the diagonal increases up to the antipode. Thus d≤2.

In all these remaining cases, multiplying `(4/3)c−ℓ` by 3m gives the lower bound

    d(m−4d−3ε)≥0.

For odd m≥11 the smallest case is m=11,d=2, where equality holds. The cases m=7,8,9,d=1 and even m≥10 are also nonnegative. Opposite-cycle point energy is at least 2r/9. Equation (2) proves the ordinary and rotated avoid-node estimates for every rank. A1 is handled separately below.

For B,C,D use orthonormal coordinates e1,…,er. The following point potentials give all required resistances directly by squared distances; their adjacent differences solve the corresponding normalized affine-root equations.

* B_r: `v0=0`, `v1=e1`, and `vk=(e1+…+ek)/2` for 2≤k≤r. Kernel entries are 1 at nodes 0,1; 2 at nodes 2,…,r−1; and √2 at r.
* C_r: `v0=0`, `vk=(e1+…+ek)/√2` for 1≤k≤r. Kernel entries are 1 at both ends and √2 at the interior nodes.
* D_r: `v0=0`, `v1=e1`, `vk=(e1+…+ek)/2` for 2≤k≤r−2, `v_(r−1)=(e1+…+e_(r−1)−er)/2`, and `vr=(e1+…+er)/2`. Kernel entries are 1 at the four fork leaves and 2 elsewhere.

The ordinary omissions are 0,r. The modified B omission pair is r−1,0; the modified C pair is 1,r. Leaf switches and reversal cover the other forbidden-node cases. In all of these cases every retained contact has `ℓ≤c`, except the following distant contacts, whose cost is already strictly above the threshold T/(4/3):

| Case | Distant contacts | Their inverse diagonal | Threshold |
|---|---|---:|---:|
| Ordinary B | opposite left fork leaves, from omission r | r/4 | r/6 |
| Ordinary D | opposite fork leaves | r/4 | r/6 |
| Ordinary C | opposite end | r/2 | r/6 |
| Modified B, r≥3 | left fork leaves, from omission r−1 | (r−1)/4 | (r−1)/6 |
| Modified C | opposite end from omission 1; node 0 from omission r | (r−1)/2; r/2 | (r−1)/6 |

For example, in ordinary B the interior-node costs and losses are `c=k, ℓ=k/2` from omission 0, and `c=r−k, ℓ=(r−k)/2` from omission r. A same-end fork leaf has c=ℓ=1. In modified B the isolated short end has c=1/2 and ℓ=0, so it causes no linear loss at all. In C the interior costs equal their linear losses. The certificate records the full 24-row index-class table and verifies its symbolic margins with independent nonnegative parameters for the node index and the excess rank; this is not an inference from a finite rank scan.

The point resistances are r/4 for ordinary B/D, r/2 for ordinary C, (r−1)/4 for modified B, and (r−1)/2 for modified C. Each is at least its required target. Equation (2), even with coefficient one on the remaining contacts, proves all these cases. Modified B2 is separate: its point energy is 1 and its minimum positive cost is 1/2, already sufficient for target 2/9. Ordinary B2/B3 are covered by the displayed coordinate formulas without assuming their inverse minima equal one.

## Exceptional factors and every forbidden node

The checker reconstructs the six ordinary exceptional cases and all 34 forbidden-node cases from their literal Gram matrices. Let T be the target and discard the immediate case C≥3T/4. Nonnegative inverse entries then restrict retained contacts to nodes with diagonal strictly below 3T/4.

There are only 46 single-contact occurrences across these records, and one possible two-contact pattern. Every other support is excluded by its diagonal costs, including cross terms when two contacts lie in the same row. There is no restriction on the magnitude of any remaining contact beyond the endpoint cost and the gap between zero and one.

For one retained contact the rows are

    p=a e_i+b e_l,  q=c e_j,
    b,c≥1,  a=0 or a≥1.

Writing `u=δ_l b/(δ_i a+δ_l b)`, the normalized-row energy is the quadratic

    E(u)=(1−u)R_ij+uR_lj−u(1−u)R_il.

If a=0, a valid lower bound for W+(4/3)C is

    min(δ_l,δ_j)² R_lj +(4/3)c_il.

If a≥1, a valid lower bound is

    min(δ_i,δ_j)² min_[0,1] E(u) +(4/3)c_il.

Both bounds meet T in every serialized occurrence. The quadratic minima are evaluated exactly at the two endpoints and, when applicable, the interior stationary point. The critical E8 ordinary case is omission 6, cheap node 8, opposite omission 1, with a=0: it gives equality 8/9+(4/3)(2/3)=16/9.

The only surviving two-contact pattern is ordinary E7: one isolated A1 contact at node 4 after deleting 3 and one at node 7 after deleting 6. Its exact energy is `(1+u²+v²)/2≥1/2`; its total cost is at least one. Thus W+(4/3)C≥11/6>14/9. Three contacts already exceed the cost threshold. This accounts for all exceptional rows, not merely single-contact proposals.

The avoid-node rules remain those of the frozen 181 proof, including E7 pairs 4,7 or 2,5, and E8 pairs 3,7 or 1,6. A1 has ordinary point energy 1/2 and positive cost at least 1/2; its avoid-node target is zero. It remains the exact obstruction to eliminating the additive avoid-node loss 2/9 in this local scheme.

## Global consequence and scope

The same geometric argument supplies actual terminating walls framing the chosen narrow facet, with a common height ratio q and total cost ΣC≤4. Only the single factor containing the facet normal uses the avoid-node target. Consequently

    Σ σ_a τ_a E_a ≥ M(2N/9−(4/3)·4−2/9),
    −(x,y) ≥ M(N−43)/9,  N=n−1.

At n≥170, N≥169, so the lower bound is at least 14, contradicting the strict narrow-facet distance bound arccosh 7 for square-two normals. The ideal-vertex branch gives n≤169. The unchanged reviewed compact-facet branch gives n≤143. Their combination is the stated 169 research proof candidate.

For clarity about constants, coefficient one with the same total cost four and avoid-node penalty 2/9 would give 163, not 145. The value 145 would correspond to eliminating the Schur-cost loss entirely. The exact E7/E8 examples explain why that unit-charge route does not work.

The published narrow-facet theorem, endpoint existence, common-height product identity, actual-ridge incidence bridge and compact-facet proof are inherited mathematical dependencies, bound through the frozen 181 checkpoint. The new results are the local obstruction, linear loss accounting and coefficient 4/3. Computation checks the exact identities and finite exceptional case coverage; the all-rank and geometric arguments remain mathematical proofs requiring independent scrutiny, not a formal proof-assistant verification.

```sh
python scripts/fresh_nikulin_unit_charge_20261015.py
python scripts/fresh_nikulin_unit_charge_20261015.py --verify
```
