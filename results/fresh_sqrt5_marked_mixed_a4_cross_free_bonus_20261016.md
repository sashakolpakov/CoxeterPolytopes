# Three mixed root-pair payments leave an exact gap above 35

The enlarged model has the exact lower bound

**7811987653852521769865/223193384043020776594 ≈ 35.0009821633 > 35.**

The gap is certified by rational weak duality. It is an obstruction to this payment family, not a dimension bound or a numerical rounding error. The finite vertex relaxation reached 35.0010871171 after seven solves and 19,836 cuts; its fully converged optimum is not claimed.

The three added mixed orbits are the literal missing corners identified by the [same-type A4/D4 diagnostic](fresh_sqrt5_marked_a4_d4_cross_free_bonus_20261016.md): an A4 endpoint with the D5 degree-two root, an A4 endpoint with the D5 long-arm endpoint, and an A4 endpoint with a D4 leaf. Their eligible root multiplicities are respectively 2×1, 2×1, and 2×3.

For a A4 copies and b target copies, each eligible pair pays γ/(2mab), giving total γ. An eligible marked A4 root removes γ/(2a); an eligible marked target root removes γ/(mb). Both counts are fixed exactly. The dynamic program assigns the unmarked total to the A4 copies and subtracts the marked-root loss; a signed component debit is bookkeeping, while every complete-forest coefficient is nonnegative. Either count zero gives zero payment. The certificate includes 206 literal pair-count checks. Before generation, 98 full-profile comparisons against literal root-pair summation passed across every feasible pair of copy counts and all facet-root orbits. A separate agent checked the normalization algebra only.

The model has 30,995 variables, including 397 new mixed coefficients. Its complete marked finite-forest search uses 3,165 exact contexts for A1, I2(5), H3, H4, A4, D4, and D5. It retains all 1,979 inherited baseline-zero quadrilateral rows in 227,408 common contexts, and 163 selected triangle families in 400,044 marked contexts. There are 48,686 distinct face profiles.

The exact dual has 380 vertex terms, 53 triangle terms, and 1,108 quadrilateral terms. Every dual vertex is separately recounted by literal root pairs in the producer replay. The largest masses now lie at 4H4+A5 marked at an A5 endpoint (about 32.30%) and A1+4H4+F4 marked in A1 (32.13%). The second shares H4 compensation with triangle corners whose facet mark lies on an A-type tail. Separating an external A1 facet mark from other external non-H4 marks is a concrete next refinement; it has not been applied in this package.

The [certificate](fresh_sqrt5_marked_mixed_a4_cross_free_bonus_20261016.json) also records a bounded coverage scout: attach one finite component of rank at most 17 to node 0 of triangle core 51 by a 3-edge, testing every rooted orbit from the finite catalogue. Fifty extensions have all three corners spherical. These rows were **not** added to the present LP. Among them, a white D5 attached at its long endpoint gives corners H4+D5, A1+D8, and H3+D6. Attaching at a short leaf instead gives H4+D5, A1+E8, and H3+E6. Three common external H4 components raise these corners to rank 21. Such branched extensions must be checked before interpreting any future budget-35 candidate more strongly.

The [source](../scripts/fresh_sqrt5_marked_mixed_a4_cross_free_bonus_20261016.py) passed generation, final diagnostic refresh, and a fresh saved replay. An interrupted initial run produced no certificate; the final run removed redundant canonicalization without changing the model.

**Scope:** restricted payment-model obstruction. Selected triangle coverage is incomplete, including the documented new extension rows; paid-corner quadrilaterals also require complete coverage. No dimension theorem or compact arithmetic example follows.
