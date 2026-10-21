# Independent check of the complete local packing cover

**FINAL/FROZEN.** Generation and a fresh saved replay both pass for the
[independent checker](../scripts/review_sqrt5_all_local_packing_cover_20261017.py)
and its [exact certificate](review_sqrt5_all_local_packing_cover_20261017.json).
The checker imports no producer module, optimizer or packing enumerator.
It uses standard-library rational arithmetic and bit-mask intersections.

All 740 contributing finite forests are covered. The check reconstructs
the 6,011 required block-orbit marginals from the unchanged face masses,
then verifies all 5,255 explicit packing-mixture terms. Each of the 23,324
literal blocks belongs to its asserted orbit, and blocks in a packing are
pairwise disjoint. All mixture masses are positive rationals. Their covered
marginals dominate every required marginal exactly.

This independently establishes the part of the
[complete local-cone result](fresh_sqrt5_all_local_simplex_packing_cones_20261017.md)
needed for its universal obstruction. Write w_j for the required mass of
orbit j and a_j(P) for its multiplicity in a disjoint packing P. The saved
mixtures satisfy sum_P lambda_P a_j(P) >= w_j. Hence every nonnegative
charge q_j with packing capacity at most one satisfies

\[
 \sum_j w_jq_j
 \le\sum_P\lambda_P\sum_j a_j(P)q_j
 \le\sum_P\lambda_P.
\]

Replacing the vertex mass by the larger of its previous mass and this
exact cover mass preserves all old constraints and pays every such charge
simultaneously. The independent check recomputes the entire repaired mass
vector and obtains

\[
 \frac{69774201720926976000000000000000}
      {5448173768070467817598921440139}
 =12.80689726342\ldots>12.
\]

The underlying finite-forest automorphism orbits and positive simplex
features are taken from the frozen all-subset package. Averaging each
literal packing under the full forest automorphism group gives uniform
marginals in each orbit; no symmetry of an ambient chamber is assumed.
This review does not repeat that orbit reconstruction, the enumeration
proving the local lower optima, or the inherited signed-transfer checks.
Only the explicit upper cover is needed for the obstruction above.

The scope is nonnegative charges depending on the finite forest and
maximal-simplex block orbit. Additional nonincident walls, coface flags,
geometric restrictions on which packings can occur, and mixed simplex/prism
constraints can carry further information. The mixtures are abstract
packings, not asserted configurations in one hyperbolic polytope. This
is an internal review, not external refereeing or a dimension theorem.

```sh
python3 scripts/review_sqrt5_all_local_packing_cover_20261017.py --verify
```
