# Five coherent additions to the H5 face

The [previous 118-wall chamber](fresh_h9_higher_face_20261013.md) extends to 123 actual H9 facets. Five new walls retain every previous wall; two have prescribed finite angles with another new wall. One addition closes an actual open H5 edge and creates a finite endpoint. The H5 face now has 24 facets, 26 finite vertices, 42 finite edges and 46 open rays, in one finite-vertex component. **Neither the H5 face nor the H9 chamber has finite volume:** the exact certificate contains strict spacelike escapes for both.

This is a coherent metric construction outside the preceding roof recipes, but it does not complete the H5 face. The first four cuts reduce remote gaps without creating finite vertices; only the fifth creates the new endpoint. No general obstruction to further completion is claimed.

## The five exact normals

All five normals have norm 2 and are perpendicular to the fixed H5 walls

`seed2, compact1, compact2, compact3`.

Their other five prescribed contacts are:

| New wall | Prescribed walls, in order | Angle labels |
| --- | --- | --- |
| closure0 | seed6, nextH4, roof3, adaptive0, roof0 | 2,6,5,5,4 |
| closure1 | wall17, nextH4, roof0, adaptive1, crosspair0 | 2,2,6,8,2 |
| closure2 | wall17, nextH4, roof0, adaptive1, crosspair0 | 8,8,8,8,16 |
| closure3 | seed6, crosspair0, secondAdaptive0, closure0, adaptive0 | 3,2,6,5,3 |
| closureEndpoint | seed6, nextH4, facePrism, roof5, closure0 | 2,2,2,2,4 |

Every unlisted pairing with an earlier wall is strictly dotted. In particular, closure3 meets closure0 at π/5, and closureEndpoint meets closure0 at π/4. These are exact finite equalities, not conclusions drawn from approximate cosine comparisons.

The [recipes](../data/enumerations/fresh_h9_h5_closure_20261013_recipes.json) determine the normals algebraically. Let E contain the nine prescribed normals, let `M=EGEᵀ`, and set `b_i=−2cos(π/m_i)`. Define

`p=bM⁻¹E`,

`z=c−projection_span(E)(c)`,

`r=p−sqrt((2−p²)/z²) z`,

where c is the retained strict timelike center. Certified nonzero elimination pivots establish that M is invertible; both `z²` and the radicand are positive. The defining identities `p·E=b`, `z·E=0` and `p·z=0` prove `r²=2` and all prescribed pairings exactly. Outward-rounded rational intervals prove the strict signs, dotted inequalities and branch orientation. The finite-angle interval overlaps in the implementation are sanity checks on these identities, not equality proofs.

A first endpoint trial used π/3 at closure0. The exact affine test rejected it: closureEndpoint, nextH4 and closure0 would form the Euclidean triangle with labels 2,3,6. The retained π/4 choice avoids that circuit. Numerical pair compatibility alone would not have caught the compactness obstruction.

## Actual face and volume checks

The [certificate](../data/enumerations/fresh_h9_h5_closure_20261013_verify.json) verifies all 600 new pairs: 45 finite and 555 dotted. The unchanged strict timelike center and nonpositive pairings imply that all 123 inequalities define actual facets. Indeed, its orthogonal projection onto any wall remains timelike and satisfies every other inequality strictly.

Each inserted wall has exactly nine finite neighbors. All 512 subsets containing it are tested. Across the five insertions the 2,560 cases comprise 1,024 exclusions by a dotted pair, 896 positive-definite subsets and 640 indefinite subsets; there is no affine subset. Together with the inherited 118-wall certificate, this excludes positive affine subdiagrams in this literal configuration. It does not imply finite volume.

The four fixed face walls are unchanged and all five new walls are orthogonal to them. Thus the H5 face has exactly the old 19 facets plus these five, with their physical Coxeter angles. An old wall already disjoint from the face cannot acquire an intersection when further inequalities are imposed. The complete finite-vertex census tests all `binomial(24,5)=42,504` face-facet subsets: 42,465 contain a dotted pair, 26 are positive definite and 13 are indefinite. Every positive case has an explicit timelike, future-oriented point strict on every other H9 wall. This proves the 26 vertices are actual and gives the full finite-vertex census.

The new endpoint has the five incident H5 facets

`seed6, nextH4, facePrism, roof5, closureEndpoint`.

Its predecessor has

`seed6, nextH4, facePrism, roof5, crosspair0`.

Their four common walls define the previously open edge that is now closed. Sharing four H5 facets gives 42 finite edges and 46 remaining open rays, with `2·42+46=5·26`. All 26 vertices belong to one component. The affine census excludes ideal endpoints supported by retained facets; no claim is made that the remaining light-cone ends resist new walls.

The H9 strict escape has exact norm

`(42139596315036540421702734738 − 26043379655490313420562864501 φ) / 8000000000000000000000000`

`≈0.06535631696078>0`,

where `φ=(1+sqrt(5))/2`. Its future-center pairing and every one of 123 strict inequalities are checked. Its projection onto the four fixed walls' orthogonal complement is also spacelike and is strict on every nonfixed inequality; its norm is enclosed near 0.06535622187390. This independently establishes infinite volume of the H5 face. The first new cap is also checked to cut the preceding 118-wall strict escape.

The proposal for this remaining escape came from the numerical gap supported by `wall17, nextH4, roof6, secondAdaptive0, closure1`. That support is a next construction target, not an exact vertex classification; the strict rational witness and its infinite-volume implication are exact.

## Search scope and replay

The numerical proposals came from continuous angular feasibility on positive two-dimensional complements after imposing four face contacts, followed by a fifth finite anchor. Three different remote supports and coherent compatible subsets were considered. The finite profile labels were bounded, nonpositive residual complements were skipped in that scout, and the subset selection was heuristic. A later spherical-edge search supplied the retained endpoint. These searches do not exhaust all cap profiles, labels, replacements or continuous completions.

The [verifier](../scripts/fresh_h9_h5_closure_20261013_verify.py) imports no numerical scout or optimizer:

```sh
python scripts/fresh_h9_h5_closure_20261013_verify.py --verify
```

It reconstructs the earlier 106-wall core, uses source-bound rigorous enclosures and exact relation tables for the preceding 12 roof walls, then defines and verifies these five new normals. Its dependency is the already certified 118-wall construction; it does not claim to be a second implementation of that predecessor. The frozen predecessor enclosures have denominator `10^80`; the loader preserves those same rational bounds if a reviewer increases the working interval precision. The norm and projection identities refer to the exact algebraic normals specified by those predecessor certificates.

The new 123-wall certificate and its saved-result replay pass. The inherited 118 files are unchanged. This checkpoint supplies a verified additional finite endpoint and two new metric couplings, while preserving explicit evidence that substantial closure remains necessary.

The parent's [source and higher-precision review](fresh_h9_h5_closure_20261013_review.md) also passes. It checks the projection formulas and census logic and replays the new construction at160 decimal places using the same interval engine, while retaining the correctly rescaled80-digit predecessor enclosures. This is not a separate implementation of all123 normals.
