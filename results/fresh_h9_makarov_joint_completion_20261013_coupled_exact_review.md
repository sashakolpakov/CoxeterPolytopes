# Review of the two-cap residual-face construction

**PASS.** The root reviewer read the [producer](../scripts/fresh_h9_makarov_joint_completion_20261013_coupled_exact.py), checked the geometric implications, and ran its saved exact replay after the helper-source manifests were hardened. The [review manifest](fresh_h9_makarov_joint_completion_20261013_coupled_exact_review.json) binds that final state. This review is distinct from the separate arithmetic implementation maintained by the independent checker.

Both new normals are projectively defined over `Q(sqrt(5))`. Their positive squares and normalization factors multiply exactly to two. Each is right to the same seven original walls and strictly dotted to all twenty others; the two new normals are mutually right. In particular, the old polar wall's unequal square is used in the dotted comparison. Thus the 55 added pairs are all valid, and the prior 351 pairs are inherited from the frozen verified source.

The common center is strict for both new inequalities. The seven-wall support has positive-definite Gram. Adding the two orthogonal positive directions gives a positive-definite nine-wall support. Orthogonal projection of the center produces the archived future timelike point; checking all 29 signs gives precisely the nine expected equalities. This proves an actual finite vertex. The standard acute projection argument also proves that every old and new inequality remains an actual facet.

Both new walls strictly cut the previous exact common-bank escape. The new archived integer vector nevertheless has positive square, negative center pairing, and strictly negative pairing with every one of the 29 normals. A segment from the strict timelike center to this strict spacelike vector crosses the ideal boundary through an open feasible patch. Therefore this chamber still has infinite volume.

```sh
python3 scripts/fresh_h9_makarov_joint_completion_20261013_coupled_exact.py --verify
```

The exact replay passed. The new construction is a valid continuation outside the excluded old dictionary, not a compact or finite-volume example.
