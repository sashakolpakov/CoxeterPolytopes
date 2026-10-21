# Independent hostile review of the contact-cost bound

**PASS.** No gap was found in the new charging argument giving the research candidate **n≤189** for finite-sided finite-volume real hyperbolic Coxeter polytopes, using the explicitly inherited geometric and published inputs. The review makes no publication-priority or current-best claim. The previously frozen 225 proof and review are unchanged.

This note reviews [the new argument](fresh_nikulin_energy_charge_20261015.md). Its [independent checker](../scripts/fresh_nikulin_energy_charge_20261015_review.py) and [saved certificate](fresh_nikulin_energy_charge_20261015_review.json) use only standard-library `Fraction` arithmetic. They import no producer code or algebra engine. Exceptional energies are calculated by grounding a coordinate in rationally congruent affine Gram matrices, rather than using the producer's regularized inverse. All complete physical Gram entries, retained inverse diagonals, and relevant point energies agree.

## Scope of the new proof

The proposed improvement changes only the loss assigned to a non-narrow affine factor. The prior [hostile 225 audit](fresh_dimension225_hostile_20261015.md) supplies the finite/ideal narrow-facet dichotomy, actual-ridge terminating-wall construction, positive common null-height ratio, relative factor scales, and compact-facet bound n≤143. The classical localization and fixed omission choices are inherited from their bound primary proofs. No geometric hypothesis is silently strengthened, and the new calculation does not re-prove these inputs through a finite matrix census.

Let C be the sum of the two endpoint Schur costs in one factor. Their global sum is at most four, because each norm-two endpoint contributes at most two. Each retained inverse is entrywise nonnegative, and every nonzero norm-two Coxeter contact has magnitude at least one. Thus a retained contact costs at least its inverse diagonal, regardless of how large its label or dotted product is. An omitted contact is either zero or at least one as well; it has not been assumed finite.

The quotient energy is nonnegative. Writing M for the common maximum of the two reciprocal null-height ratios gives the exact identity

    sigma*tau/M = min(sigma,tau)^2.

This identity is unchanged by a separate positive normalization of each affine kernel. Therefore the scale conventions in the exceptional calculations do not introduce a hidden common-scale assumption.

## Exceptional factors and the low-cost cases

The independent checker reconstructs A1, E6, E7, E8, F4 and G2, with the fixed omissions from the earlier product proof. For F4 and G2 it clears radical edges by positive diagonal congruence and tracks the diagonal scaling exactly. The actual inverse diagonals, rather than the diagonals of the rescaled matrices, are compared. Positive definiteness of every retained block and nonnegativity of every inverse entry are checked.

The direct A1/E6/F4/G2 factors satisfy both needed inequalities: the zero-contact point energy is at least 2r/9, and twice the minimum retained inverse diagonal is at least 2r/9. Nonnegative energy then handles every positive cost. The low B2 and B3 presentations are independently reconstructed; their minimum retained inverse diagonals are 1/2 and 3/4, respectively. They are not incorrectly assigned the generic lower bound one.

For E7 and E8, costs at least one are harmless because 2r/9<2. If 0<C<1, the complete inverse lists leave precisely one nonzero retained contact across both rows. It is an isolated A1 entry of inverse diagonal 1/2. Two inexpensive contacts would already cost at least one; any other retained contact would do so alone. Thus, after possibly exchanging endpoints, the rows really are

    p = a*e_i + b*e_l,    q = c*e_j,
    sigma = 2a+b,         tau = 2c,       C=b^2/2,

with b,c≥1 and a either zero or at least one. The zero row is excluded by the inherited positive null-height argument. The independent matrices give delta_i=delta_j=2 and delta_l=1 in all three cheap cases: E7 has two, E8 one.

The grounded calculations independently yield R_ij=1/2, R_lj=1 and R_il=1/2. Hence the interpolation energy is exactly (1+u²)/2, where u=b/(2a+b). If a=0, the energy is one and min(sigma,tau)≥1, so weighted energy plus twice the cost is at least 1+1=2. If a≥1, min(sigma,tau)≥2 and energy≥1/2, so weighted energy alone is at least two. Both exceed the required baseline 2r/9. The inequalities remain valid for arbitrary allowed contact magnitudes; no angle truncation is present.

## Classical factors and the A9 boundary

The all-rank charging comparisons have positive, explicit expressions. For A_r with 2≤r≤8 the zero-energy margin is 2r(8−r)/(9(r+1)). For r=9+t, t≥0, the large-rank margin is

    (t+1)(t+18) / (36(t+10)) > 0.

The corresponding B/D and C margins at r=10+t are 37/36+t/36 and 41/18+5t/18. The inverse-minimum and energy assumptions here are exactly those in the bound classical proof; no rank sampling is substituted for their closed formulas.

There is one range extension requiring an explicit check: the old product report stated its localization estimates for r≥10, while the new A argument starts at r=9. In A9 the inverse diagonal at distance three from an end is 21/10>2; all five middle entries exceed two. Thus only the first two and last two retained nodes can be contacted. Opposite deletions in the affine 10-cycle leave two intervening arcs of length one. Their effective conductance is 1/2, which is stronger than the needed (9−8)/4=1/4. The new boundary rank is therefore covered directly, rather than inferred from the old stated range.

## Global conclusion and limitations

Every non-narrow factor now satisfies weighted energy divided by M at least 2r/9−2C. Summing its cost losses gives at most eight. The one narrow factor retains the old loss allowance two, even if it itself also uses contact cost. This double allowance is conservative; no cost is omitted from the global upper bound.

Consequently total energy is at least M(2N/9−10). In the Lorentz identity the reciprocal-height term is at most 2M, giving −(x,y)≥M(N/9−7). At n≥190, N=n−1≥189, so this is at least 14. The inherited actual-ridge result puts both selected walls in the framing set of the strict arccosh(7) narrow facet, contradicting that bound for square-two normals. The compact-facet branch remains n≤143. The two branches give n≤189 with the same generality as the earlier candidate.

This conclusion does not rely on arithmeticity, simplicity at ideal vertices, a single affine cusp factor, a finite label bound, or a common relative scale for distinct cusp factors. Its published-source and geometric dependencies are explicit. This is an independent research proof audit, not a claim of external peer review or a literature-priority determination.

Generation and saved replay, syntax, source manifests and local links pass.

```sh
python scripts/fresh_nikulin_energy_charge_20261015_review.py
python scripts/fresh_nikulin_energy_charge_20261015_review.py --verify
```
