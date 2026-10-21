# Independent review of the common-height charge and dimension 166

The [independent checker](../scripts/fresh_nikulin_coupled_charge_20261015_review.py) passes for the [new local estimate](fresh_nikulin_coupled_charge_20261015.md). It reconstructs exceptional energies using grounded rational linear algebra and rechecks the complete support census, rather than importing the producer's arithmetic. This supports an **internally reviewed research proof candidate**, not external validation or a priority claim.

For the same prescribed omissions as the frozen [four-thirds proof](fresh_nikulin_unit_charge_20261015_review.md), let `M` be the common height ratio and `C` the retained-support cost. The new statement is

    M ≤ 21/20  ⇒  W + (7/6) C ≥ 2r/9,

with target `2(r−1)/9` in the single factor containing the narrow-facet normal. If the individual height lower bounds have squares `a,b`, the smaller actual height has square at least `max(min(a,b),max(a,b)/(21/20)²)`. This elementary improvement is valid for every nonnegative contact magnitude and is the only additional height hypothesis.

The independent reconstruction covers all 40 exceptional omission records and all 54 single-contact occurrences. It minimizes each normalized energy quadratic exactly. It also enumerates every subset of cheap contacts, including the nonnegative inverse-Gram cross costs. Only one multiple-contact support survives: ordinary E7, whose independently established energy is `(1+u²+v²)/2`; its cost is at least one, so its lower bound `5/3` exceeds the target `14/9`. No larger support survives. The review preserves every forbidden-node choice and verifies the literal normalized Grams.

For affine A, inverses of independently grounded path matrices recover the cycle resistances for all sizes 2 through 19. Each possible retained contact is checked against both its endpoint budget and the new strict small-cost threshold. Four symbolic distance/parity classes cover every larger cycle. For B, C and D, the review reuses the explicitly identified, independently proved coordinate identities in the frozen four-thirds review, and checks all 24 classes against the new threshold. These inherited identities and the published affine classification remain mathematical dependencies; this is not a new formalization of them.

The global arguments concerning actual framing ridges, terminating endpoints, one common height ratio, and the compact-facet alternative are inherited from the [reviewed dimension-169 binder](fresh_nikulin_unit_charge_20261015_dimension169.md). They are not inferred from these local arithmetic checks. Writing `N=n−1`, the two bounds are

    M ≤ 21/20:  −(x,y) ≥ M(N−40)/9,
    M ≥ 21/20:  −(x,y) ≥ M(N−43)/9.

For `N=166+t`, `t≥0`, these are at least `14+t/9` and `287/20+7t/60`, respectively. Both contradict the strict narrow-facet threshold 14. The compact-facet branch still gives 143. Thus the combined candidate is **n≤166**, without an arithmetic, compactness or simple-cusp assumption. At `N=165,M=1`, the new scalar lower bound is below 14; this particular estimate does not establish 165.

The [saved certificate](fresh_nikulin_coupled_charge_20261015_review.json) binds the exact producer files and inherited review dependencies. The producer's frozen pending-review wording records its earlier status.

```sh
python scripts/fresh_nikulin_coupled_charge_20261015_review.py --verify
```
