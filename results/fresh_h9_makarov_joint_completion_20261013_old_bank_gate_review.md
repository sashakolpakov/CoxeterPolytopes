# Review of the retained-chamber cap-bank obstruction

**PASS, with the stated finite-bank scope.** The root reviewer read the reconstruction, outward interval primitives and all-label angle gate, and ran the saved producer replay successfully. This is a proof/code review with an exact replay, not a second independent implementation of all 25,101 reconstructions. The [review manifest](fresh_h9_makarov_joint_completion_20261013_old_bank_gate_review.json) binds the reviewed sources.

Nine prescribed pairings have rank nine over the base field. Their affine solution line intersects the norm-two quadric in two distinct real roots. Exactly one lies in the archived proposal's rational coordinate box; the other is rigorously separated. Thus the norm equation is an algebraic identity and the floating proposal only selects the exact branch. The unequal original polar-wall length is correctly included in the prescribed pairings.

The escape is strictly spacelike, has negative pairing with the center, and is strict on all 27 retained walls. Every recipe either remains strict there (24,152 entries), pairs positively with a retained cap (553), or lies in a forbidden finite-angle interval (396). The angle gate encloses labels 2 through 120 and conservatively retains the entire larger-label tail, as well as parallel and dotted relations. Hence no finite label was discarded by a numerical tolerance or a field-degree assumption. A positive-proportional duplicate of a retained wall already misses the escape and is harmless.

The conclusion follows for any acute-compatible finite selection retaining all 27 facets: the same strict spacelike point survives, giving an open ideal escape and infinite volume. The proof does not exclude replacement of a retained facet, a different quadratic branch outside the certified box, different recipes, or arbitrary completing normals.

Replay:

```sh
python3 scripts/fresh_h9_makarov_joint_completion_20261013_old_bank_gate.py --verify
```
