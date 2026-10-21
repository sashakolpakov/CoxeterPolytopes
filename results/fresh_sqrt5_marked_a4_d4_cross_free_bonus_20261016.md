# Same-type A4/D4 cross payments remove one barrier but leave another

Adding the two missing root-pair orbits from quadrilateral 64 removes its exact two-term barrier at 36. The enlarged model nevertheless has the exact lower bound

**256032020701238/7254500362725 ≈ 35.29285380 > 35.**

This is an exact obstruction to the specified enlarged payment family. The finite vertex relaxation reached numerical lower bound 35.29943836 after four solves and 14,011 vertex cuts; neither value is asserted to be the fully converged optimum.

The new payments are confined to endpoints in distinct A4 components and leaves in distinct D4 components. If there are h copies and m eligible roots per component, there are m²h(h−1)/2 eligible unordered pairs. Each pays 2γ/[m²h(h−1)], so their total is γ. A marked eligible root removes fraction 2/(mh); an ineligible root removes none. The model treats h<2 as zero. The certificate checks all counts h=2–5 and every root literally. A separate agent checked this normalization algebra, without claiming an independent review of the implementation or face coverage.

There are 206 new cross-context coefficients, distinguished by actual same-type count, exact A1 count, capped I2(5) count, and H3 presence. Together with all earlier free coefficients, the model has 30,598 variables. It carries exact common counts for A1, I2(5), H3, H4, A4, and D4, with 2,238 rank-21 count contexts and complete marked finite-forest separation. All 1,979 inherited quadrilateral rows remain, now in 191,218 common contexts. The 163 selected triangle families give 340,502 marked contexts. The combined system has 48,170 distinct face profiles.

The exact dual contains 337 vertex terms, 59 triangle terms, and 847 quadrilateral terms. Every dual vertex is recounted by literal root pairs. The largest vertex masses are 2A1+H4+3D5 marked in A1 (about 43.51%), 3A1+2D4+2D5 marked in A1 (23.53%), and A1+4D5 marked in A1 (16.36%).

The [certificate](fresh_sqrt5_marked_a4_d4_cross_free_bonus_20261016.json) also gives an explicit local reassignment for the [preceding two-term barrier](fresh_sqrt5_marked_h4_remainder_free_bonus_20261016.md). Set its D5 weight to 3/8, and take A4 totals at counts 2,3,4,5 to be 1,2,4,5; take D4 totals to be 9/8,9/4,9/2,5. Every possible common A4/D4 spectator distribution for that quadrilateral still pays at least 1, while A1+4D5 and every marked-A1 rank-21 forest made only from A4 and D4 have budget at most 35. This assignment is deliberately scoped to that local subsystem.

The new dominant rows are prior 172, 1260, and 1246. Literal inspection identifies three missing mixed orbits: A4 endpoint with the D5 degree-two root, A4 endpoint with the D5 long-arm endpoint, and A4 endpoint with a D4 leaf. A corner whose forest is D4+D5 is not automatically a mixed pair: in row 1260, both marked roots in that corner lie inside D5. The certificate stores the complete three rows and their actual root positions to prevent this confusion.

The [source](../scripts/fresh_sqrt5_marked_a4_d4_cross_free_bonus_20261016.py) passed generation, final literal-witness refresh, and a fresh saved replay. All preceding frozen artifacts are retained.

**Scope:** restricted-model obstruction only. Newly affected triangles and paid-corner quadrilaterals still need complete coverage; no dimension-bound promotion or compact arithmetic example follows.
