# Independent review of the two finite-edge caps

**PASS; the 90-wall chamber is still infinite-volume.** The [separate checker](../scripts/fresh_h20_e7_finite_edgewalk_20261014_review.py) and [certificate](fresh_h20_e7_finite_edgewalk_20261014_review.json) verify the [construction](fresh_h20_e7_finite_edgewalk_20261014.md) and its finite-edge frontier.

The checker uses sparse squarefree radical arithmetic from an earlier independent review, with outward rational enclosures for base-field signs and exact relative-norm comparisons for the bridge's positive quadratic radical. It imports neither the construction nor its arithmetic implementation. It reconstructs both caps from their shared-wall projections and E7 fundamental weight, then checks all 177 new wall pairs, the common center, and 90 actual-facet witnesses.

All 84 specified finite points have exactly their advertised 20 incident facets and are strict at every other wall. Their distinct active sets prove they are different vertices. Each of the 31 forward endpoints has a positive edge parameter, negative norm, the correct time orientation, and all 90 inequalities. Convexity then certifies its segment from the old vertex as an actual edge.

The independent finite-support enumeration uses all 19-subsets of each new cap's complete finite-neighbor set, rather than the producer's recursive clique traversal. FLINT principal determinants implement Sylvester's criterion independently of its rational elimination routine. This reproduces 3,080 candidates: 40 positive-definite supports and 3,040 failures. The 40 accepted supports exactly match the saved new vertices, including nine incident to both caps. The list of 145 finite edges also agrees.

The complete finite-endpoint classification of the 800 germs at those 40 new vertices is independently checked: 259 have a second finite vertex and 541 do not. A germ retaining a new cap is covered by the complete new-cap census. Every germ releasing its sole new cap has the certified old endpoint. This argument does not require a complete census of the old 88-wall chamber and does not classify the ideal endpoints of the remaining germs.

The new parallel pair `(9,88)` has a future null sum satisfying all 90 inequalities, with exactly 20 incident facets. Its scaled Gram matrix has rank 18. A separately checked vector of square `98/27` is orthogonal to the null ray and every incident normal. Modulo the null direction it gives a nonzero Euclidean translation preserving every inequality of the horospherical link. The link is therefore unbounded. Its component matrices and positive affine kernels are also checked.

Finally the inherited rational escape remains strictly feasible and has positive square

```text
259854917001622621 / 62500000000000000.
```

Together with the strict timelike center, it proves infinite volume. The frozen 88-wall Coxeter proof and completeness of the two previously filled cusp links are inherited; their preservation is checked directly. No new finite-volume example is claimed.

Generation and saved replay pass with the source hashes in the certificate.
