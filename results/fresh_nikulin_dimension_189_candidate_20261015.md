# Reviewed research proof candidate: dimension at most 189

**Candidate theorem.** Every finite-sided finite-volume real hyperbolic Coxeter polytope in H^n satisfies **n≤189**. There is no arithmeticity assumption, no simple-ideal-vertex restriction, and no facet-count or finite-label cutoff.

This package combines the frozen [225 baseline](fresh_nikulin_dimension_225_candidate_20261015.md), the new [contact-cost estimate](fresh_nikulin_energy_charge_20261015.md), and its [independent proof and arithmetic review](fresh_nikulin_energy_charge_20261015_review.md). It preserves every earlier artifact. The [binder](../scripts/fresh_nikulin_dimension_189_candidate_20261015.py) and [certificate](fresh_nikulin_dimension_189_candidate_20261015.json) record their source, data, and report hashes and validate the recursive dependency manifests.

The conclusion is presented as a reviewed research proof candidate, not as an externally peer-reviewed theorem, a publication-priority assertion, or a determination of the current strongest known bound. Nikulin's ICM report already announces the unrestricted inequality n<300. That historical claim remains distinct from the explicit estimates reconstructed here. [Primary report, printed pp.657–659](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

## Exhaustive case split

Nikulin's narrow-facet theorem supplies a facet F whose framing hyperplanes have pairwise distance strictly below arccosh 7. For outward normals of square two, the selected framing walls x,y therefore cannot satisfy −(x,y)≥14. The construction uses walls sharing actual ridges with F, so they belong to the precise framing set in the theorem. [Theorem A, printed p.655](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf).

If F contains an ideal vertex, its horospherical link is a product of affine Coxeter simplices with total rank N=n−1. The frozen product argument selects two outgoing edges, and terminating walls x,y that meet F in actual ridges. Their normalized affine null-height ratio q is common across all factors. With M=max(q,1/q)≥1, their Lorentz pairing is

    −(x,y) = (1/2) Σ_a σ_a τ_a E_a − q − 1/q.

At each terminating endpoint the total retained Schur cost is at most two. If C_a is the sum of the two endpoint costs in factor a, then Σ_a C_a≤4. The new estimate proves

    σ_a τ_a E_a / M ≥ 2r_a/9 − 2C_a

for every factor except the one containing F's normal. That one factor retains the previous uniform loss allowance of two. Consequently

    Σ_a σ_a τ_a E_a ≥ M(2N/9−10),
    −(x,y) ≥ M(N/9−7).

For n≥190, N≥189 and the final lower bound is at least 14, contradicting narrowness. This branch therefore gives n≤189.

If F has no ideal vertex, finite-sidedness and finite volume make F compact: its vertices are among the finitely many finite and ideal vertices of the ambient projective polytope. The already reviewed [compact-facet lemma](fresh_nikulin_narrow_20261015_compact_facet.md) gives n≤143. It uses ambient vertex diagrams and quasi-Lannér extraction; it does not assume that F is itself a Coxeter polytope or that its ambient face stars have no parabolic subdiagrams.

These cases exhaust the chosen narrow facet. Thus n≤max(189,143)=189.

## What the sharper estimate adds

The 225 proof charged a fixed loss to each of at most eight spoiled non-narrow affine factors. The new proof charges their loss to the actual total Schur cost, bounded by four independently of how many factors occur.

For classical factors and affine A_1,E_6,F_4,G_2, the cost of one retained nonright contact, together with the inherited localization estimate, directly pays the required loss. The two exceptional issues are E_7 and E_8, whose chosen deletions can leave an inexpensive A_1 contact. Exact inverse-diagonal checks show that cost below one permits only a single such contact. Three exact resistance values, 1/2,1,1/2, then give the interpolation identity E=(1+u²)/2. Splitting according to whether the omitted-node contact is zero or at least one proves the same charge inequality. This is an all-label argument, not a search with a maximum Coxeter label.

The independent review reconstructs the exceptional arithmetic and checks the geometric and inequality cases. The earlier independent reviews continue to cover common null-height scaling, the actual-ridge endpoint selection, the prescribed omission choices, and the compact-facet branch. The primary energy-charge report's final pending-review sentence is historical: that report was frozen before the separate review finished; its completed review is bound in this package.

The binder checks source closure, the exhaustive case combination, and the exact scalar threshold. It imports no producer arithmetic. Published classifications and the geometric proof arguments remain explicitly identified mathematical inputs; these computational checks are not a formal proof-assistant development.

```sh
python scripts/fresh_nikulin_dimension_189_candidate_20261015.py
python scripts/fresh_nikulin_dimension_189_candidate_20261015.py --verify
```
