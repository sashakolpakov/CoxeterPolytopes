# Independent review of selective distance12 and dimension136

**PASS.** The [producer proof](fresh_selective_distance12_dimension136_20261017.md) gives a complete internally checked computer-assisted proof of dimension at most136 for finite-sided, finite-volume real hyperbolic Coxeter polytopes, without arithmeticity or cusp restrictions. Its separate compact-facet conclusion is ambient dimension at most103. No unresolved geometric case or arithmetic discrepancy was identified in this review.

The [independent checker](../scripts/fresh_selective_distance12_dimension136_20261017_review.py) and [saved certificate](fresh_selective_distance12_dimension136_20261017_review.json) use standard-library Fraction arithmetic, literal graphs and an independently organized matching enumeration. They import no producer arithmetic, SymPy, NumPy or optimizer. The mathematical review below is separate from those calculations. This is internal review, not external refereeing or proof-assistant formalization; no publication-priority claim is made.

## Coverage of the minimal-core split

The inherited actual-edge argument makes both marked deletions of each minimal indefinite opposite-side core positive definite, even when some quadrilateral vertices are ideal. The opposite closed sides remain disjoint. The complete order6–10 quasi-Lanner catalogue and the order-at-most5 elementary radius bound therefore have the hypotheses required by the producer.

The checker reconstructs all384 finite-black-pair contexts directly from the frozen literal diagrams. It recomputes distances and75 ports at distance at least5, obtaining18 allowed ports. For diagrams with a4 bond, a diagonal sqrt2 congruence gives a rational Gram; the inverse diagonal is rescaled back to the norm-two gauge. Every tested black deletion is independently checked positive definite. The bound excludes inverse diagonal equal to2 as well as larger values: the relevant block including the external white is contained in an actual edge star, which is positive definite, not merely semidefinite.

Under the assumption of no marked white path of length at most12, a disjoint joining bond must have endpoint radii summing to at least12, hence both at least5. The18 allowed ports yield exactly115 ordered joins of six underlying types. The checker verifies that both rejected affine blocks are wholly white in every relevant marked context, not just in an unmarked representative.

For overlapping cores, every shared root has its two radii summing to at least13. Thus each radius is at least6. The16 relevant marked contexts each have at most two such roots; intersections of three or more roots are impossible. Rather than using the producer's filtered combinations of candidate identifications, the reviewer enumerates left and right subsets and all bijections between them. It recovers208 candidate identifications and162 consistent ones:

| Outcome | Count |
|---|---:|
| External-neighbor inverse obstruction |120|
| Wholly white affine E8 |6|
| L25 copies sharing their far node |27|
| L25 copies sharing the reversed far edge |9|

An external neighbor used by an inverse rejection is a white root: its distance from the other core's black pair is too large for it to be black. The enumeration treats the matching as the complete intersection, so that this neighbor really is outside the first core. The independent checks reproduce every rejection category and surviving identification.

## Disjoint joins and unspecified cross contacts

The far-near and near-near white blocks have strictly positive affine kernels. Increasing a bond magnitude or adding another acute contact decreases their quadratic form on that positive kernel. Thus their obstruction applies beyond the minimal simple-bond representative.

After those types are removed, every marked context pair permits only one remaining cross edge. The checker verifies this uniqueness and a wholly white A3 tail at each endpoint. Eliminating those tails gives endpoint residuals4/3, forcing squared cross contact below16/9<2; consequently the only possible nonright Coxeter contact is simple. No finite-label cutoff is used, and parallel or dotted contacts are excluded by the same bound.

Independent rational Gaussian elimination reproduces the four two-center Schur matrices and positive determinants:

| Right core/port | Schur matrix | Determinant |
|---|---|---:|
|L18/7|[[-1/15,-1/10],[-1/10,-7/30]]|1/180|
|L22/8|[[-5/66,-1/11],[-1/11,-7/44]]|1/264|
|L22/4|[[-1/15,-1/10],[-1/10,-1/5]]|1/300|
|L25/9|[[-7/78,-1/13],[-1/13,-7/78]]|1/468|

Each eliminated block is independently positive definite. The residual matrices are negative definite, proving negative index at least2. This excludes all disjoint possibilities in the actual Lorentzian ambient space.

For each surviving overlap context the reviewer also tests every unspecified cross edge. It either creates a marked path of length at most12 or, in the shared-vertex case only, creates the wholly white affine triangle across the shared node. Hence the two surviving core Grams are fixed; their isolation proofs do not silently assume other cross contacts vanish.

## All36 supporting-white isolation cases

The19-root shared-vertex core has the stated one-dimensional Gram radical and a timelike vector in its span. Its Gram-null combination is therefore an actual linear relation: it is null and orthogonal to a timelike vector, whose orthogonal complement is positive definite. This justifies pairing an arbitrary additional supporting white with that relation.

The independent checker reconstructs all27 preferred corners, including every finite/affine block, inverse diagonal, positive affine kernel and allowed extra contact. Their marked components have order16 or17. Every allowed contact has a nonzero radical coordinate; the negative and positive coordinates lie respectively within distance2 of the two opposite marked pairs. Acuteness and the actual radical relation force contacts on both sides whenever any contact exists, producing a white path of length at most6. Thus every additional supporting white is orthogonal to the core.

The18-root shared-edge core is nonsingular. Its eliminated positive block and residual

`[[-5/66,-1/11],[-1/11,-5/66]]`

give signature(17,1); its independently recomputed determinant is−1. The nine preferred corners have order16 and marked distance13. Their PSD tests leave possible contacts only in the two end sets{1,2,3} and{11,12,13}. A contact at both ends gives a marked path of length at most4.

The independently calculated full inverse block at either end is

`[[8,10,5],[10,14,7],[5,7,4]]`.

For a nonzero one-end contact vector b, all entries have the same nonpositive sign and every nonzero magnitude is at least1. Hence b^T G^(-1)b>=4. Because the core already contains the unique negative direction, the root's orthogonal residual must have nonnegative square; instead its square is2−b^T G^(-1)b<=−2. This contradiction holds in every ambient dimension and for arbitrary contact magnitudes. It is not a positive-definite Schur argument applied incorrectly to a Lorentzian block.

At a finite preferred corner, simplicity identifies its isolated component with the actual ambient corner component. At an ideal preferred corner, the inherited actual factor-fiber lemma places the entire varying affine factor in W plus the two side roots. Thus hidden incident cusp walls cannot enlarge an isolated finite or affine marked component. This is the essential justification for the small-component payment in a nonsimple cusp.

## Global payment and exact dimension calculation

The selective rule pays a marked pair if its distance is at most12 or its connected component has order at most17. The core analysis pays every quadrilateral. The inherited triangle argument supplies two paths of length at most8. Components of order at most17 have at most8 times their order in pairs, while larger classical finite or affine diagrams satisfy the inherited all-rank distance counts. Therefore finite vertices cost at most12n and a cusp costs at most

`V sum_i min(12,r_i/2)`.

The convex-order transfer applies with cap24 exactly as with the earlier caps: transfers between two ranks below24, or two above24, preserve the charge and do not increase the curvature. The two rank regions must remain distinct; a transfer across the cap is not being assumed. The resulting normal forms contain copies of24, at most one smaller part and at most one larger part.

The reviewer independently enumerates these normal forms for d=136,...,145. Integer convolution and a separate inclusion-exclusion moment calculation agree on all1,175 profiles. Every margin is positive; the least is172096/1953125. The finite-vertex slack beta_n−12n is positive for every n=137,...,146, and the shape inequality retains its strict positive5(n−1) term. Together with the completed146 upper tail, this proves136 without a new hypothesis on cusp rank or factor count.

At n=136, the profile(24,24,24,24,24,15) has margin−59157/3906250. The checker confirms this boundary; no135 conclusion follows from the same scalar argument. For the compact-facet branch it independently obtains the parity polynomials with coefficients(12,208,4) and(112,212,4), proving ambient dimension at most103.

## Validation and limits

Final generation and a fresh saved replay both pass. Syntax, all ten direct bindings, all19 producer input bindings and local report links pass. The producer source/data/report are bound at their final frozen hashes. The computations are independent exact rechecks of the displayed finite cases; the completeness of the26-type quasi-Lanner catalogue, the actual-face incidence theorem, all-rank graph budgets, convex-order theorem and shape inequality remain explicit inherited mathematical inputs. They are not replaced by bounded numerical testing. Source, data and report are final/frozen.
