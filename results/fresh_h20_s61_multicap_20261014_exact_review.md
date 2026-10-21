# Independent review of the S61 right-chain completion

**PASS.** Independent generation and saved replay both pass. The review source, certificate and this note are frozen; no mathematical issue was found.

The review covers the new 26-wall addition to the frozen 107-wall chamber. It does not repeat the predecessor's complete geometry census. The [separate checker](../scripts/fresh_h20_s61_multicap_20261014_exact_review.py) imports only earlier independent field arithmetic, not the new producer's arithmetic. Its [certificate](fresh_h20_s61_multicap_20261014_exact_review.json) binds the final source and input JSON files.

The physical Lorentz frame and all 107 normalized old-wall covectors are reconstructed directly. The raw old normals 96 and 97 are tied coefficientwise to the saved orthogonal decomposition, including their transverse weight of square 2/3. This justifies their projected S61 norms 4/3 and the signs used in the endpoint equation.

The radius parameter is defined by an exact rational interval. The checker evaluates the chain using complex powers by repeated squaring: the spatial coordinate is `(-sqrt(a²+1)+i)*(C-iS)^j`, with `C=a²/(a²+2)` and `S=2sqrt(a²+1)/(a²+2)`. Its derivative is evaluated by the closed power rule. This differs from the producer's repeated rotation and automatic differentiation. Outward 512-bit rational intervals certify the opposite endpoint signs and a strictly positive derivative throughout the interval. Existence and uniqueness therefore follow for the exact old algebraic coefficients.

The equalities are proved symbolically, rather than inferred from intervals containing zero. The identities `C²+S²=1`, spatial norm `a²+2`, and `(a²+2)C=a²` prove every norm-two and consecutive-right relation. The central reflected pair is right because `-a²+(a²+1)-1=0`. The defining endpoint equation and reflection give the two exact old-side right contacts.

All 3,107 new pairs are replayed: 495 right and 2,612 strictly dotted. The common center is checked directly. It gives actual facets by the acute orthogonal-projection argument. All three bounded old cusps are strictly retained.

For the intrinsic face, the checker independently evaluates all 406 pairs of boundary lines. Exactly 29 consecutive intersections are feasible, future timelike, and satisfy every old ambient inequality; each of the other 377 intersections violates a polygon inequality. The two central chain normals and old normal 95 have spatial parts spanning R² with a strictly positive zero-sum. Their positive time coefficients make the spatial feasible region bounded. Thus the complete face is a compact 29-gon. Its two inherited angles need not be Coxeter angles; no intrinsic Coxeter-polygon claim is used.

The complete new finite-neighbor lists justify exactly 27 new ambient finite vertices. Each new wall has the eighteen S61 normals and only its two cyclic neighbors as finite neighbors. Those two neighbors are mutually dotted, so any finite vertex must use the entire S61 support and one neighbor. Both endpoint types are explicitly realized. Together with acute preservation, this gives at least 168 known finite vertices overall, not a complete ambient census.

Every new normal has a strictly future-timelike projection to the already compact S75 plane, so that face remains unchanged. Finally, the same rational spacelike escape has square `38017800175279/500000000000`, remains strict for all new walls, and has the correct orientation relative to the new center. **The 133-wall chamber still has infinite volume**, despite the compact completion of both named planes.
