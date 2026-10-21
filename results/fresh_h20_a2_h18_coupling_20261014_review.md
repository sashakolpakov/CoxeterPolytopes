# Independent verification of the direct A2 obstruction

The direct A2 coupling of this fixed forty-facet H18 source always has infinite volume. The independent [checker](../scripts/fresh_h20_a2_h18_coupling_20261014_review.py) and [certificate](fresh_h20_a2_h18_coupling_20261014_review.json) reconstruct the full 138-normal relaxed profile bank and exhibit a common strict spacelike escape. This independently confirms the [primary scalar-bound proof](fresh_h20_a2_h18_coupling_20261014.md).

The review uses SymPy arithmetic in `Q(√2,√3,√5)`, with no producer or project arithmetic imports. Every nonzero sign is certified by a rational interval obtained from integer square-root bounds; no floating-point sign is used. It reconstructs all 1,600 source Gram entries, checks ambient inertia `(20,1,0)`, all 25 source-π/3 transition cases, all 138 normal squares, and all 276 strict escape and center pairings.

For normalized new-wall square two, the A2 inverse Gram has diagonal `2/3` and off-diagonal `1/3`. Two nonright contacts already use all available positive square. A single label at least six does the same. Thus the complete inducing-profile list is the ordinary profile and either color at labels three, four or five. No upper label cutoff is assumed.

An ordinary/special pair above a source π/3 edge would have inner product strictly between minus one and zero, forbidden for norm-two Coxeter walls. Ordinary/special status therefore propagates along the source's connected 24-node π/3 component. That component contains the path `0–1–2`, whose endpoints are orthogonal. Two special lifts of those endpoints have positive inner product for every pair of colors. Hence all 24 component facets are forced ordinary, leaving sixteen possible special facets. The review checks all 27 edges and the literal induced-path witness.

The independent witness uses source scale `3207/1000`, transverse coordinates `(-1,-1)` and source-center perturbation `1/10^9`, different from the producer's `1/10^7`. Its exact squared norm is positive, while all 138 wall pairings are negative. A common timelike center is also checked, with negative pairing against the escape. The segment between these points crosses the future null cone strictly inside every inequality, giving an open ideal-boundary region and therefore infinite volume for every coherent direct assignment.

The original finite-volume H18 source theorem remains inherited. This result excludes only direct couplings using the two A2 support walls and walls inducing the fixed source facets. Additional walls, changed source metrics and other H20 constructions remain open. No full vertex census or assignment enumeration is claimed.

Generation and saved independent replay pass.
