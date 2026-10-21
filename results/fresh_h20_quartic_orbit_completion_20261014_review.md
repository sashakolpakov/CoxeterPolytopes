# Review of the quartic changed-chamber scout

The reflection-bank construction, finite-angle alphabet, and necessary cusp-edge test are sound within the stated **numerical** scope. No exact endpoint exclusion or finite-volume example follows from this scout.

The [separate checker](../scripts/fresh_h20_quartic_orbit_completion_20261014_review.py) imports neither producer. Using Python integers, it reconstructs the complete set of original roots and their images under one original reflection, verifies every archived provenance recipe, and recomputes every norm. The bank has exactly 29,043 distinct oriented vectors: 2,163 of norm two and 26,880 of norm four. This check also avoids relying on the producer's fixed-width integer arithmetic. The [review certificate](fresh_h20_quartic_orbit_completion_20261014_review.json) binds both scouts, both numerical archives, the bank, and their source data by hash.

For a rational ambient root projected along a normal over a quartic field, its squared cosine with a projected seed root belongs to that field. For `m >= 3`, the degree of `cos²(pi/m)` is `phi(m)/2`; hence it must divide four. This gives exactly the producer's necessary superset `{2,3,4,5,6,8,10,12,15,16,20,24,30}`, with `m=2` added separately. The independent checker verifies the irreducibility of all 44 stored quartic polynomials and reconstructs the label list. The argument does not require the quartic field to be Galois or every listed cosine to occur in it.

The seed link is a product of affine simplices, as certified in the upstream local census. Omitting one wall from each component gives all its vertices and therefore all outgoing hyperbolic edges. The remaining 19 normals have positive-definite Gram `A`. A cap of projected norm `N` and pairing vector `b` can supply a finite or ideal other endpoint only if `N - bᵀ A⁻¹ b >= 0`. Strict positivity is necessary for a finite endpoint; equality is retained for a possible ideal endpoint and does not certify a bounded cusp link. Splitting the quadratic cost by the mutually orthogonal affine factors is valid. A wall perpendicular to the starting cusp either contains such an edge or never closes its other end, so the restriction to strictly negative cusp pairings is appropriate for closing caps.

The checker exhausts the archived omission tuples independently and verifies their counts: 220 sections, 42,420 edges, 8,136 numerically accepted cap proposals, and 96 near-angle proposals. The archived endpoint lists are all empty. This is an exact check of the **records and their combinatorial coverage**, not an independent recomputation or certification of the floating-point Schur signs.

The producer evaluates its algebraic normal at high precision but then uses double-precision geometry and fixed tolerances. Those gates can neither prove absence at borderline parameters nor establish exact Coxeter angles. The added caps are compared with the retained cusp walls, not with one another. Only five shifts and one ambient reflection step are searched. Both source docstrings, result scopes, and the main report preserve these limitations. No correction is needed.

Replay the independent integer and archive checks with:

```sh
OPENBLAS_NUM_THREADS=1 python scripts/fresh_h20_quartic_orbit_completion_20261014_review.py --verify
```

This command does not overwrite or rerun the numerical producers.
