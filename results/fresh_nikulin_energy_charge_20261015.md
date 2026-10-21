# Charging affine-factor losses to contact cost: a bound of 189

**Research proof candidate.** Every finite-sided finite-volume real hyperbolic Coxeter polytope in H^n satisfies **n≤189**.

This improves the separately frozen [225 candidate](fresh_nikulin_dimension_225_candidate_20261015.md) by one change: the loss in a spoiled affine factor is charged to its actual Schur-complement contact cost. The narrow-facet component retains the earlier loss allowance of two. No further choice optimization or finite-label enumeration enters this report. The published inputs, actual-ridge endpoint selection, affine product scaling, and compact-facet branch are unchanged. Publication priority and current-best status are not asserted.

The [checker](../scripts/fresh_nikulin_energy_charge_20261015.py) and [certificate](fresh_nikulin_energy_charge_20261015.json) reconstruct the exceptional inverse diagonals and three small resistance identities exactly, and check the all-rank inequalities below. They bind the earlier geometric proofs rather than importing producer arithmetic. Independent review is separate from this primary package.

## Setup retained from the 225 proof

Choose a narrow facet F. If it is compact, the frozen [compact-facet lemma](fresh_nikulin_narrow_20261015_compact_facet.md) gives n≤143. Otherwise F contains an ideal vertex. Its link is a product of affine Coxeter simplices, of ranks r_a summing to N=n−1.

Use exactly the fixed omission choices from the [product-cusp proof](fresh_nikulin_narrow_20261015_products.md): opposite cycle nodes for A, the designated end nodes for B/C/D, and the recorded pairs for the exceptional types. Only the component containing F's normal uses the old avoid-node modification. These choices precede the selection of terminating contact rows.

Write x,y for the resulting two norm-two terminating walls, chosen to share actual ridges with F. In a component let p,q be their nonnegative negated contact rows, δ its positive affine kernel, and σ=δ·p, τ=δ·q. The ratio σ/τ is common to all components. Put M=max(σ/τ,τ/σ)≥1. Let E_a be the corresponding normalized-row energy in the Euclidean cusp quotient. The exact identity is

    −(x,y) = (1/2) Σ_a σ_a τ_a E_a − σ_a/τ_a − τ_a/σ_a.

For the two deletion supports S_i,S_j in this component define its **contact cost**

    C_a = p_Siᵀ G(S_i)⁻¹ p_Si + q_Sjᵀ G(S_j)⁻¹ q_Sj.

Each endpoint's total Schur cost is at most two, so Σ_a C_a≤4. The inverses are entrywise nonnegative. Every nonzero contact is at least one, including contacts with an omitted node; no upper bound on its label or magnitude is used.

We prove, for every component other than the one containing F's normal,

    σ_a τ_a E_a / M ≥ 2r_a/9 − 2 C_a.                 (1)

If C_a=0, both rows are supported only on their omitted nodes. The frozen point-contact bound already gives σ_aτ_aE_a/M≥2r_a/9. It remains to consider C_a>0.

## Classical and direct exceptional factors

A positive contact in a retained row costs at least the corresponding inverse diagonal, because its magnitude is at least one and all cross terms are nonnegative. Thus C_a is at least the smallest inverse diagonal among the two retained supports.

For affine A_r, r≥2, this minimum is r/(r+1). For r≤8, nonnegativity of energy bounds the deficit by 2r/9, and

    2r/(r+1) − 2r/9 = 2r(8−r)/(9(r+1)) ≥ 0.

For r≥9 the unchanged localization argument gives energy at least (r−8)/4 after division by M, so the deficit is at most 2−r/36. Here

    2r/(r+1) − (2−r/36) = (r−8)(r+9)/(36(r+1)) ≥ 0.

Both ranges prove (1). Affine A_1 is separate: its retained inverse diagonal is 1/2 and its baseline is only 2/9, so the zero-energy bound suffices.

For the unmodified B/D end deletions, the retained inverse minimum is at least one in ranks r≥4. The B_2 and B_3 presentations have minimum respectively 1/2 and 3/4; their zero-energy baseline is already smaller than twice that value. For r≤9 the bound 2r/9≤2 suffices. For r≥10, the old localization gives energy at least (r−5)/4 after division by M, so the deficit is at most 5/4−r/36<2. This proves (1) for B and D.

For C_r, r≥2, either end deletion leaves finite B_r, with inverse minimum one. The zero-energy argument covers r≤9. In larger rank the retained bound (r−5)/2 already has deficit at most 5/2−5r/18≤2. Thus (1) holds here as well. These claims use the closed inverse formulas in the frozen classical proof; the checker separately verifies the low B ranks where the generic minimum-one statement would be false.

The direct exceptional cases are also immediate from their exact inverse diagonals:

| Affine type | Rank | Fixed omitted nodes | Point energy | Minimum retained inverse diagonal |
|---|---:|---|---:|---:|
| A_1 | 1 | 0,1 | 1/2 | 1/2 |
| E_6 | 6 | 2,4 | 4/3 | 4/3 |
| F_4 | 4 | 0,4 | 1 | 1 |
| G_2 | 2 | 0,2 | 2/3 | 2/3 |

In each row the point energy is at least 2r/9, and twice the last entry is also at least 2r/9. Zero energy therefore suffices whenever a retained contact occurs. The node numbering is fixed by the full matrices in the certificate.

## The inexpensive E_7 and E_8 contacts

A minimum-diagonal argument alone would fail in these two types. The fixed E_7 omissions 3,6 each isolate an A_1 node, respectively 4,7. The fixed E_8 omissions 3,7 have just one inverse diagonal below one: node 8 after omitting 7. Each of these inexpensive diagonals is 1/2; every other retained inverse diagonal is at least one. The checker verifies the entire two inverse-diagonal lists, rather than only these three selected entries.

If C_a≥1, inequality (1) follows from nonnegative energy since 2r/9<2 for r=7,8. Suppose 0<C_a<1. Nonnegativity of inverse entries implies that there is exactly **one** retained nonzero contact across the two endpoint rows: it is one of the inexpensive A_1 contacts. Two such contacts would already cost at least one, and any other retained contact would cost at least one by itself.

Orient the two rows so that this contact, of magnitude b≥1, belongs to x. Let i be x's omitted node, l the inexpensive retained node, and j the other endpoint's omitted node. The kernel is normalized so that

    δ_i=δ_j=2,   δ_l=1.

The contact rows consequently have the exact forms

    p = a e_i + b e_l,       q = c e_j,
    σ = 2a+b,               τ = 2c,
    C_a = b²/2.

Here c≥1 because the second row cannot vanish, and a is either zero or at least one. This last alternative is the only discreteness used: a nonright Coxeter contact cannot have magnitude strictly between zero and one.

For nodes u,v let R_uv be the squared energy of the normalized point-row difference e_u/δ_u−e_v/δ_v. In each of the three inexpensive-contact cases, exact inversion gives

    R_ij=1/2,   R_lj=1,   R_il=1/2.                  (2)

If a=0, the normalized first row is the point row at l. Hence E_a=1. Since both σ,τ≥1,

    σ τ E_a/M = min(σ,τ)² E_a ≥ 1,
    2C_a=b²≥1.

Their sum is at least two, which exceeds the required baseline 2r/9.

If a≥1, both σ and τ are at least two. Put u=b/(2a+b). The first normalized row is the convex combination (1−u)e_i/δ_i+u e_l/δ_l. The standard squared-distance interpolation identity and (2) give

    E_a = (1−u)R_ij+uR_lj−u(1−u)R_il
        = (1+u²)/2 ≥ 1/2.

Therefore σ τ E_a/M = min(σ,τ)² E_a≥2, so (1) follows even without adding its contact cost. The same argument applies when the inexpensive contact belongs to y. This closes all E_7/E_8 possibilities, including arbitrary large finite labels and dotted omitted contacts.

## Summing the costs and completing the dichotomy

The one component containing F's normal keeps the frozen estimate

    σ_a τ_a E_a/M ≥ 2r_a/9 − 2.

There is no need to charge that component's cost, or to optimize its omission choices. Summing (1) over the other components and using Σ_a C_a≤4 gives

    Σ_a σ_a τ_a E_a ≥ M(2N/9 − 8 − 2).

Since σ_a/τ_a+τ_a/σ_a≤2M, the Lorentz identity implies

    −(x,y) ≥ M(N/9−7).

For n≥190, N=n−1≥189 and the right side is at least 14. Both x,y are actual framing walls of the narrow facet, so this contradicts its strict distance bound arccosh 7. Thus the ideal-vertex branch now gives n≤189.

If the chosen narrow facet has only finite vertices, finite-sidedness and finite volume make it compact, and the unchanged compact-facet proof gives n≤143. The two cases therefore yield the stated candidate bound n≤189.

This argument changes the energy accounting only. It does not change the 225 baseline artifacts, the endpoint-existence argument, cusp rank decomposition, positive common null-height ratio, or the proof that the selected ending walls meet the exact narrow-facet framing condition. Its classification and primary-source scope are those documented in the frozen [combined candidate and audit](fresh_nikulin_dimension_225_candidate_20261015.md). Independent scrutiny of the new contact-cost step is required before using this sharper value as a checked checkpoint.

```sh
python scripts/fresh_nikulin_energy_charge_20261015.py
python scripts/fresh_nikulin_energy_charge_20261015.py --verify
```
