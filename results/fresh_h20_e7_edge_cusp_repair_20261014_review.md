# Independent review of the five-wall cusp repair

**PASS: 95 actual facets, three specified bounded cusp links, and still infinite volume.** The [independent checker](../scripts/fresh_h20_e7_edge_cusp_repair_20261014_review.py) constructs the five new normals directly from the frozen 90-wall chamber and its null ray. Its [construction certificate](fresh_h20_e7_edge_cusp_repair_20261014_review.json) does not depend on the repair producer or numerical scout. The [endpoint certificate](fresh_h20_e7_edge_cusp_repair_20261014_endpoints_review.json) additionally compares the primary construction and checks its accepted finite endpoints and ridge witnesses.

Write the new cusp ray as `ν=n9+n88`, and let `L` be the saved tangent vector, of square `98/27`. The quotient roots are the negatives of the highest roots of the finite A5 and A1 factors, the negative highest short root of the B3 factor, and `±(3√3/7)L`. Their shifts by multiples of `ν` are, respectively,

```text
1, 2, 1, 1 + 4√3/7, 10√3/7.
```

Separate sparse radical arithmetic verifies the norm-two identities and all 460 new unordered pairs: 107 right, ten π/3, one π/4, two π/6, five parallel and 335 dotted. The common strict timelike center proves all 95 facets actual by orthogonal projection and acuteness. Every previous finite vertex is strict at each new wall; the two previous bounded cusps are likewise preserved.

The repaired cusp is incident to exactly 25 of the 95 walls. Its six orthogonal Gram components have ranks `8,5,3,1,1,1`, summing to 19. Each component is positive semidefinite of corank one, with an explicitly positive kernel and a positive-definite principal submatrix. These components bound the full 19-dimensional horospherical link. This conclusion is independently computed, rather than inferred from a numerical rank.

The endpoint mode checks exact equality with all 95 primary normals and all new pair records. It then verifies twelve distinct new finite vertices: six on wall 92 and six on wall 93. Each point is future timelike, has exactly its stated 20 incident walls, and is strict on all other inequalities. Its exact positive edge parameter and the preserved old endpoint certify a first intersection along the stated edge. Eight additional finite-angle ridge witnesses are also checked.

The six wall-93 endpoints each share 19 independent incident walls with the repaired cusp. Their open segments to that null ray are strict on every other wall, so they give six actual finite-to-proper-ideal edges. The independently checked construction therefore includes 96 specified finite vertices and these six paths into the bounded cusp. This is not a complete vertex or edge census; the producer's full bounded sweep of 541 earlier germs is not independently re-enumerated here.

The same rational strict escape survives all five walls, with square

```text
259854917001622621 / 62500000000000000 > 0.
```

Thus the chamber remains infinite-volume. Local cusp completeness and additional verified paths do not establish a finite-volume completion. Generation and saved replays pass with the source hashes in the two certificates.
