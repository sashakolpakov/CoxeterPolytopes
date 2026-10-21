# Independent review of the direct B2 coupling

**PASS.** The direct B2 coupling theorem and the coherent42-wall example have no issue in the stated scope. The common strict spacelike escape applies to every assignment of the allowed inducing profiles, including assignments not mutually compatible. Extra walls that do not induce one of the forty fixed source facets remain outside the theorem.

The [independent checker](../scripts/fresh_h20_b2_h18_coupling_20261014_review.py) and [review certificate](fresh_h20_b2_h18_coupling_20261014_review.json) bind the frozen [primary construction](fresh_h20_b2_h18_coupling_20261014.md). They use a new standalone Fraction implementation in the basis `1,√2`, without producer or project arithmetic imports.

The replay reconstructs the nineteen-dimensional source coordinates from the literal40-root Gram and a verified invertible19-root submatrix. It checks all1600 source Gram entries, reconstructs the orthogonal B2 summand and its weights, and independently obtains ambient inertia `(20,1,0)` by exact symmetric elimination. It then checks all74 profile norms, the complete sets of24 forced ordinary and16 potentially special facets, every bank-center and bank-escape sign, all861 coherent-example pairings,42 actual facet points with1764 signs, and the actual source-face point. The nonrational two-cycle and rational indefinite two-wall witness are also reproduced.

The all-label profile argument is sound. For normalized ambient norm2, the B2 inverse has diagonal1 and positive off-diagonal `1/√2`. Two nonzero support contacts already cost at least `2+√2`; a single contact of label at least4 costs at least2. Parallel or dotted support contacts cost still more. A facet-inducing wall must have positive projected norm, so only the ordinary profile and the two single-π/3 profiles remain.

A source π/3 edge prevents either endpoint from being special: one special endpoint gives a nonzero product strictly between−1 and0, while two special endpoints give a positive product. The literal source therefore leaves exactly the saved16 potentially special facets. The independent replay reconstructs the entire74-normal union and verifies the escape against every member.

The strict-center/strict-spacelike argument proves infinite volume without a vertex census: the segment crosses the future null cone strictly inside every wall halfspace, giving an open ideal-boundary region. It applies uniformly to every coherent direct assignment. The actual-facet projection formula and the source-face restriction are valid because the retained example is acute and has the checked common timelike center.

The final general cone criterion is correctly scoped to a finite direct coupling. In transverse coordinates `−aω0−bω1`, `a,b≥0`, either special color bounds both variables above on a bounded source slice. An exact full extreme-ray census then provides a global test; future-causal rays yield finite volume, while a feasible spacelike ray can be perturbed toward a strict center to exhibit infinite volume. This criterion does not assert that the present coherent local contacts close the chamber.

The original finite-volume H18 source and its gluing proof remain inherited dependencies. This review does not repeat that source's vertex census or claim an independent classification of H18 sources. The conditional non-quasi-arithmetic conclusion concerns any hypothetical finite-volume completion retaining the witness walls; the present group has infinite covolume.

Generation and saved independent replay pass. A fresh primary saved replay, source hashes, Python compilation and local links pass.
