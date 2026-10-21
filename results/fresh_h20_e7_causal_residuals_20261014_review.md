Review **PASS** for the [86-wall replacement](fresh_h20_e7_causal_residuals_20261014_exact.json) and its [87-wall contact refinement](fresh_h20_e7_causal_residuals_20261014_contact.json). Both fresh producer `--verify` runs passed. No source/proof issue was found.

The [separate checker](../scripts/fresh_h20_e7_causal_residuals_20261014_review.py) reconstructs the E7 and A11+A2 fundamental weights directly from the frozen predecessor. It uses the earlier independent nested-quadratic arithmetic, extended by √43, rather than either construction producer or their sparse interval arithmetic. Its [saved certificate](fresh_h20_e7_causal_residuals_20261014_review.json) binds the producer sources and JSONs, not the prose report.

With `w=ω2+ω4+ω6`, the exact identities are `w²=95/2` and

`y=−√2 λ0−x λ12+tT`, `y²=−91/2`, `r=−w+y`, `r²=2`.

The first construction has `(x,t)=(2,5)` and deletes old walls 61 and 75, whose nonzero cosines lie strictly between 0 and 1/2. The refinement has `x=√(86/11)`, `t=17/√11`, retains old wall 75 with pairing −1 and deletes only old wall 61. The relation `3t²=71+x²` explains the preserved norm throughout this change.

The separate replay verifies all **171 new pairs** across the two constructions: 85 in the first and 86 in the second. It also checks each new facet point against every retained wall, both old-escape cuts, all 88 new-cap signs at the 44 old finite vertices, and the strict signs preserving both completed cusps. For the refined cap it independently reconstructs the projection of the strict center onto the two-root span with old wall 75. The Gram matrix is `[[2,−1],[−1,2]]`; the resulting future timelike point has exactly those two walls active. This proves an actual π/3 ridge, not merely an allowed pair value.

The full H13 face remains unchanged because all fifteen source walls remain, every deleted wall had future timelike projection to that face, and each new cap has future timelike intrinsic projection. Thus the deleted and added inequalities are redundant within the face. The retained cusp incident sets are unchanged, so their previous positive affine components and quotient rank 19 continue to prove bounded links. These points were reviewed in the source and checked by the fresh primary replays.

The finite-neighbor census is appropriately scoped and complete for finite vertices on the new cap. Any such vertex would supply nineteen other finite-contact walls with a positive-definite twenty-root Gram matrix. The first cap has just nineteen finite neighbors; the refinement has twenty. The scripts test every corresponding nineteen-wall subset, rejecting a nonfinite old pair, a nonpositive old Gram pivot, or a negative final Schur complement. This excludes finite vertices on the new cap in each saved chamber. The separate review checker does not repeat that subset census; it was source-reviewed and rerun through the primary verifiers.

Finally, the root's new strict spacelike escape was independently checked against all 86 and 87 normals. Both saved chambers therefore still have infinite volume. The contact refinement gives a genuine finite-angle ridge while preserving the complete H13 face and its two bounded ambient cusps; it does not supply a new finite vertex or global closure.

Separate generation and saved replay pass:

```sh
python scripts/fresh_h20_e7_causal_residuals_20261014_review.py
python scripts/fresh_h20_e7_causal_residuals_20261014_review.py --verify
```
