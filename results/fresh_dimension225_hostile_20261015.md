# Hostile audit of the narrow-facet dimension argument

**The repaired argument passes this proof audit.** The product-cusp branch gives n≤225 when the chosen narrow facet contains an ideal vertex. The separate coarse compact-facet lemma gives n≤143 if any facet is compact. Together with Nikulin's narrow-facet theorem, these establish n≤225 for a finite-sided, finite-volume real hyperbolic Coxeter polytope, under the standard published affine and quasi-Lannér classification inputs stated below. No arithmeticity, simplicity at ideal vertices, or finite-label cutoff is assumed. This is an independent proof review and small exact certificate, not a novelty or current-best literature claim.

The initially suggested compact-facet bound n≤31 did **not** follow from its proposed import. Vinberg's finer compact-chamber weights use global absence of parabolic subdiagrams. We found an exact counterexample to that localization premise. The replacement proof uses graph distance at most17 and permits affine diagrams; it does not invoke the unsupported finer transfer.

The reviewed producers are the [product-cusp proof](fresh_nikulin_narrow_20261015_products.md) and [compact-facet proof](fresh_nikulin_narrow_20261015_compact_facet.md). The [independent checker](../scripts/fresh_dimension225_hostile_20261015.py) and [saved data](fresh_dimension225_hostile_20261015.json) import no producer arithmetic. They independently reconstruct six exceptional affine energies by grounding a potential, verify the scale and rank identities, check the parity thresholds and distance counts, and certify the rational counterexample below. The geometric audit is given explicitly here; these finite calculations are not presented as a machine verification of every geometric statement or of the classification inputs.

## Published hypotheses checked

[Nikulin, ICM1986, printed pp.654–655](https://www.mathunion.org/fileadmin/ICM/Proceedings/ICM1986.1/ICM1986.1.ocr.pdf) defines the relevant reflection groups using finite-volume chambers, including nonarithmetic groups. Theorem A supplies a facet whose framing hyperplanes have pairwise distance strictly below arccosh7. His subsequent discussion explicitly treats ideal vertices on that facet. Thus compactness of the ambient chamber is not an input to this narrow-facet theorem. For outward normals of square2, the strict hyperplane-distance bound contradicts a pairing ≤−14.

[Prokhorov, printed pp.407–408, Propositions1–2](https://www.mathnet.ru/php/getFT.phtml?jrnid=im&paperid=1484&what=fullteng) supplies minimal indefinite/quasi-Lannér extraction and maximal diagram order10. His Proposition3 also records the spherical graph-distance pair count. We read these primary statements directly, rather than inferring them from a secondary dimension bound.

[Vinberg, printed pp.70–71,76–77,81–88](https://www.mathnet.ru/php/getFT.phtml?jrnid=mmo&option_lang=rus&paperid=444&what=fullt) supplies the simple-polytope counting input and spherical pair bound. His original fine-weight face proof invokes Proposition9, whose no-parabolic assertion assumes a compact ambient chamber. The new compact-facet proof replaces that step.

## Product-cusp audit

The following points address the possible rank, scale and endpoint gaps.

1. **Full cusp rank and scales.** A finite-volume ideal link is a product of Euclidean Coxeter simplices. If the factor ranks are r_a, their sum is N=n−1. The positive affine combinations c_a are positive multiples k_a c of one future null vector. For terminating walls x,y set σ_a=−x·c_a and τ_a=−y·c_a. The ratio q=σ_a/τ_a is independent of a. In the quotient, dividing the factor difference by k_a and multiplying its energy by σ_aτ_a cancels k_a². Consequently the identity

   −x·y = (1/2)Σ_a σ_aτ_a E_a − q −1/q

   uses no equal-scale or arithmetic assumption. The checker verifies this cancellation in independent null coordinates with separate symbolic factor scales.

2. **Actual outgoing edges and endpoints.** Omitting one node per affine factor gives a positive-definite support S of N walls and a vertex of the horospherical link. It determines an actual edge. Finite volume gives a finite or ideal second endpoint. An additional terminating wall is not another original cusp wall: a deleted original cusp wall is strict along the corresponding outgoing edge away from the original cusp. Acuteness and the full cusp span then give σ_a,τ_a>0. Indeed, vanishing contact with one factor would make the common cusp pairing zero; nonnegative contact rows force all factor rows to vanish, leaving a norm-two vector in the null orthogonal complement of the full cusp span, which is impossible.

3. **The framing-wall choice covers ideal endpoints.** The retained N-node support lies in the full affine diagram at a second ideal endpoint. Each component of rank r contains at most r nodes of this positive-definite support. Equality of the totals forces exactly one missing node in each component. If the narrow-facet root lies in a component of rank at least2, its missing node has a finite angle with it. If its component is affine A1, choose a missing node of another component; the resulting pair is orthogonal. A second component exists in the relevant n≥226 range. Thus the chosen terminating wall has a positive-definite pair with the narrow-facet root. Projection of a strict interior point proves an actual ridge, so this is precisely a wall controlled by Nikulin's theorem.

4. **The four-contact budget is valid at both endpoint types.** The Schur complement of S in S+x is nonnegative, so p_SᵀG_S⁻¹p_S≤2. A positive-definite diagonal-two Z-matrix has a nonnegative inverse and inverse diagonal at least1/2. Every nonzero normalized Coxeter contact has magnitude at least1, including any parallel or dotted contact. Hence at most four retained entries are nonzero. The two walls spoil at most eight factors; marking the narrow-facet factor separately gives at most nine. No upper bound on the deleted-node pairing is needed.

5. **Omissions are selected before the contact rows.** In a non-narrow factor the fixed point-contact choices are retained even if the factor later proves spoiled. Low ranks r≤9 use only energy nonnegativity. In high ranks, only affine A/B/C/D occur. The closed inverse formulas give the asserted localization. For example, in modified B_r the left row is supported in nodes0,1,2, while the right row is in r−3,r−2,r−1,r; normalized coordinates3 through r−3 have fixed difference1/2, giving r−5 coordinates. The analogous C modification leaves coordinates4 through r−2 with fixed difference1. D has at least r−4 such coordinates, and the weaker r−5 estimate suffices. Positive normalized contact rows are convex combinations of the localized point rows, so the coordinates stay fixed for arbitrary allowed contact magnitudes.

6. **The cycle bound is uniform.** For A_r, write m=r+1. The two separating arcs have lengths floor(m/2)−4 and ceil(m/2)−4. Minimizing their constant-flow energy gives their product divided by their sum. For r≥10 this exceeds (r−8)/4: the surplus is1/4 when m is even and (k−4)/(2(2k−7)) when m=2k+1. The checker verifies the identities. Thus no sampled-rank inference is needed.

7. **Final normalization.** Unspoiled energy is at least M(2r/9), where M=max(q,1/q). Each marked factor loses at most2M. The independently reconstructed exceptional scale-invariant costs are A1:1/2, E6:4/3, E7:2, E8:2, F4:1, G2:2/3. Summing and using q+1/q≤2M gives −x·y≥M(N/9−11). For n≥226, N≥225 and this is at least14. It contradicts the strict narrow-facet bound. This proves exactly the conditional product-cusp result.

## A concrete obstruction to the original fine-weight import

Use the following Gram matrix, with each normal of square2:

    G = [[ 2,-1, 0, 0],
         [-1, 2,-1,-1],
         [ 0,-1, 2,-1],
         [ 0,-1,-1, 2]].

It has determinant−3 and signature(3,1). Its inverse has diagonal (0,−2,−4/3,−4/3), and every entry is nonpositive. The four vectors v_i=−G⁻¹e_i are future causal, satisfy v_i·n_j=−δ_ij, and give one ideal and three finite simplex vertices. Their sum is a strict center of square−22. Facet0 is the compact triangle with vertices v_1,v_2,v_3. Nevertheless, its three black nodes1,2,3 form affine A2. All proper principal diagrams are positive semidefinite; the whole four-node diagram is indefinite but not Lannér. Hence this compact triangle's star contains **no Lannér subdiagram**.

This explicitly disproves the attempted no-affine/Lannér-extraction localization. It does not disprove the numerical face-weight inequality: in this example the original three corner weights are each2. A new argument is necessary, and the following coarse one supplies it.

## Hostile audit of the coarse compact-facet repair

**Minimal indefinite diagrams have the required classification scope.** Such a diagram is connected and every proper principal subdiagram is positive semidefinite. If a proper connected affine block A occurs, choose an adjacent outside node v. Its pairing with A's positive kernel is strictly negative. Restricting to that kernel and v gives an indefinite two-dimensional form. Minimality forces the whole diagram to be A+v; its coupling to the radical makes the full form nonsingular of signature(|L|−1,1). Any other proper affine block must likewise have codimension one. Thus neither zero eigenvalues nor an affine block with disconnected additional nodes evade the Lannér/quasi-Lannér classification. If all proper blocks are positive definite, the usual one-node Schur complement gives the same nonsingularity. The published order bound10 therefore applies.

**The opposite-pair block really is indefinite.** Let W be the n−2 white normals of a compact quadrilateral Q, and U a pair of opposite black sides. If W+U were positive definite, acute projection would give an actual finite common point of the opposite sides in Q. If it were singular positive semidefinite, a positive kernel of an affine component gives a nonzero future null combination c. Positivity of the whole block implies c is orthogonal to all of W+U; acuteness gives c·n≤0 for every other chamber wall. Thus c is an ideal point on Q, contradicting its compactness. Both alternatives are excluded. This reasoning uses global acuteness and compactness of Q, and permits cusps elsewhere in P.

**The graph paths belong to actual corner diagrams.** For a compact triangle, a minimal indefinite diagram contains all three black nodes and at most seven white nodes. Two pairs of black nodes have connecting paths avoiding the third; each path has at most eight edges. For a compact quadrilateral, minimal indefinite diagrams L⊂W+U and M⊂W+V each contain their opposite black pair and at most eight white nodes. They cannot be orthogonal in an index-one space. An overlap or connecting edge supplies a path between the two black sets; a shortest such path has only white internal nodes and length at most8+1+8=17. Its endpoints are consecutive sides of Q. These paths lie in the full positive-definite ambient vertex diagrams where the weights are counted.

Weight1 on ambient node pairs at distance at most17 therefore gives triangle sum≥2 and quadrilateral sum≥1. For k≥5, nonnegativity gives the required sum≥5−k. The spherical pair bound is17n at an ambient finite vertex. It also follows elementarily: components of order≤35 have at most binom(r,2)≤17r pairs; higher-rank spherical diagrams are paths or D forks, with at most r pairs at any specified distance1,…,17.

A compact facet F has dimension d=n−1 and is simple, since its vertices are finite vertices of P. Its two-faces are compact two-faces of P. Nonnegative weights on its angles count only a subset of all ambient vertex pairs, so the upper bound remains17n=17(d+1), rather than17d. The symmetric nonnegative h-vector gives average lower bounds d(d−6)/8 for even d and (d−1)(d−5)/8 for odd d, strictly. For d=143+2t and144+2t, subtracting17(d+1) gives respectively

    (1/2)t²+36t+3/2,
    (1/2)t²+(73/2)t+19.

Both are positive for t≥0. Thus a compact facet forces n≤143. All constants and the rank shift have been independently checked.

## Exhaustiveness of the dichotomy and conclusion

A finite-volume polytope is the convex hull of finitely many finite or ideal vertices. A facet without ideal vertices is compact. Apply Nikulin's theorem to choose a narrow facet. If it has an ideal vertex, the audited product-cusp argument gives n≤225. Otherwise the audited coarse compact-facet argument gives n≤143. The cases exhaust the chosen facet. Therefore their combination yields n≤225 under the stated finite-sided, finite-volume Coxeter hypotheses.

The audit does not promote a fixed-lattice or fixed-source obstruction into a dimension theorem: neither proof fixes a field, lattice, source polytope, facet count, or cusp factor count. It does rely on the published narrow-facet theorem, affine Coxeter classification, minimal-indefinite order bound, and simple-polytope h-vector facts. Broader literature priority and publication status are separate from correctness of this derivation.

```sh
python scripts/fresh_dimension225_hostile_20261015.py
python scripts/fresh_dimension225_hostile_20261015.py --verify
```
