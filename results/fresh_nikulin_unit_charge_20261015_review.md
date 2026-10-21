# Independent review of the four-thirds contact-cost estimate

The new local estimate and its dimension-169 consequence pass independent arithmetic and mathematical review. This is an internal research review. External validation, publication priority and current-best status remain unestablished.

The [review checker](../scripts/fresh_nikulin_unit_charge_20261015_review.py) and [certificate](fresh_nikulin_unit_charge_20261015_review.json) import no arithmetic from the [primary proof](fresh_nikulin_unit_charge_20261015.md). They reuse only elementary Fraction matrix operations from an earlier independent checker. Exceptional energies are recovered from grounded rational linear systems, with every equation checked, rather than the producer's regularized inverse. The nonsimply-laced matrices are treated by rational diagonal congruences; the literal physical Gram entries and the normalized kernel squares are compared with the primary data.

## Coverage of the local lemma

For each of the 40 ordinary and forbidden-node exceptional cases, the review rebuilds every retained inverse diagonal. It enumerates **every subset** of contacts whose diagonal costs are below the trivial-charge threshold, adding the inverse cross terms for contacts in the same row. All inverse entries are nonnegative, so setting every nonzero magnitude to one gives a valid support-cost lower bound. This establishes support coverage for arbitrary real contact magnitudes at least one, with no finite-label cutoff.

There are 46 single-contact occurrences and one remaining two-contact support. In the single-contact case, the energy is a convex quadratic in the normalized mixture weight. The review computes its three coefficients and exact minimum separately, then checks both possibilities for the omitted-node contact: zero, or at least one. The latter distinction supplies the correct lower bound on the two null heights. The sole two-contact support is ordinary E7; polarization of the independently reconstructed resistance matrix gives `(1+u²+v²)/2`. Its cost is at least one. Both cases meet their stated targets.

The complete E6/E7/E8 omission-pair census is independently reconstructed. Zero-cost rows eliminate the ineligible pairs; all unit point-contact rows satisfying each endpoint budget at most two give the remaining necessary coefficients. Their best values are 1, 10/9 and 4/3. In E8, the omissions `(1,6)` and contact rows at nodes 1 and 8 attain the 4/3 obstruction. This proves optimality for the **specified local lemma**, whose slope is fixed at 2/9 and whose omissions must precede the contact rows. These local rows need not be realizable together in a whole hyperbolic chamber, and no global optimality follows.

## The all-rank step

The linear identity used by the producer has the correct sign. For point potentials `v_i`, write `z_p=v_i+a` and `z_q=v_j+b`. Expanding their squared distance gives

    E=R_ij−Σw_k(R_ij+R_ik−R_kj)
           −Σz_k(R_ij+R_jk−R_ki)+|a−b|².

Each coefficient is a convex weight. When its contact is nonzero, that weight is at most one and the squared contact magnitude is at least one. Since the retained inverse has nonnegative entries, its cost dominates the sum of diagonal contributions. Thus a linear loss bounded by `(4/3)c_ik` is paid by the corresponding cost. Nonpositive losses require no payment. Finally `W=min(σ,τ)²E≥E`, since the normalized affine kernel has minimum entry one and both rows are nonzero. These steps use neither arithmeticity nor a bound on the contact labels.

The classical point potentials have an independent direct verification. For B and D use roots `α0=−e1−e2`, `αk=e_k−e_(k+1)`, with final root `√2 e_r` in B and final two roots `e_(r−1)−e_r`, `e_(r−1)+e_r` in D. In C use `α0=−√2 e1`, the same difference roots, and `αr=√2 e_r`. The prefix vectors in the primary proof satisfy, for every index,

    α_l·v_k = 1_(l=k)/δ_k − 1_(l=0)/δ_0.

This follows by subtracting adjacent coordinates; the two fork/end equations check the boundary cases. It identifies the potentials with normalized point contacts, proves their squared-distance formula, and gives `c_ik=δ_k²R_ik` by deleting node i. It is an all-rank identity, not an extrapolation from small diagrams.

The reviewer computes distances from these coordinate vectors by symbolic interval intersections. It checks 24 index classes partitioning the required B/C/D rows, with independent nonnegative parameters for the node index and excess rank. Every class either satisfies the four-thirds loss inequality or has cost strictly above the trivial-charge threshold. Path reversal and fork automorphisms cover the stated symmetric rows and forbidden-node choices. The exceptional modified B2 pair is checked directly. No producer resistance table is imported.

For A, direct substitution in `R(d)=d(m−d)/m` gives linear losses `d(m±ε)/m` for the two directions from an opposite omission, where `m=r+1` and `ε=m mod 2`. Below the trivial-charge threshold, ranks 2–5 allow no retained contact, and ranks 6–8 allow only distance one. For `m≥10`, the endpoint budget excludes distances at least three because `3(m−3)/m>2`; monotonicity up to the antipode completes this restriction. For the remaining distances, the cleared margin is `d(m−4d−3ε)`, nonnegative in both parity classes. A1 and its coincident avoid-node omissions are included among the exact exceptional cases. Cycle rotation avoids any prescribed node in all other ranks.

## Global hypothesis and threshold check

The [frozen 181 binder](fresh_nikulin_dimension_181_candidate_20261015.md) supplies the geometric dependencies: actual ridge neighbors of the narrow facet, terminating walls at finite or ideal endpoints, a common null-height ratio across affine factors, the Lorentz pairing identity, and the compact-facet branch. The new omission in ordinary E8 is chosen before the contact rows and does not alter the separate avoid-node rule. Thus the framing argument still applies.

There are two endpoint budgets, each at most two. Exactly one affine factor contains the narrow-facet normal and incurs the avoid-node deficit 2/9. With `N=n−1`, the new local estimate therefore gives

    Σ σ_a τ_a E_a ≥ M(2N/9−16/3−2/9),
    −(x,y) ≥ M(N−43)/9,    M≥1.

At `n=170`, this is at least 14. It contradicts the strict narrow-facet distance bound `arccosh 7` for square-two normals. At `n=169`, the scalar lower bound alone is smaller than 14, so no further one-unit improvement has been silently inferred. The ideal-vertex branch gives 169; the inherited compact-facet branch gives 143. These two alternatives exhaust the chosen narrow facet.

The checker certifies the exact arithmetic and finite coverage supporting these arguments. The coordinate identities, geometric bridge and published classifications remain mathematical proof dependencies, not a formal proof-assistant development.

```sh
python scripts/fresh_nikulin_unit_charge_20261015_review.py
python scripts/fresh_nikulin_unit_charge_20261015_review.py --verify
```
