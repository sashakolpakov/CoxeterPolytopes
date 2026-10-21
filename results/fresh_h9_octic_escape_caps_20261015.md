# Two separate octic additions, both of infinite volume

Two exact norm-two walls separately extend the previous 29-wall chamber. Each cuts its old escape, satisfies all 29 Coxeter pair conditions, preserves all 30 actual facets, and creates a finite endpoint. The two new walls have positive mutual pairing and therefore cannot be retained together in this acute-angled chamber. A single exact strict spacelike escape survives either addition. Neither extension has finite volume.

The [producer](../scripts/fresh_h9_octic_escape_caps_20261015.py) and [certificate](fresh_h9_octic_escape_caps_20261015.json) use quadratic extensions of the previous quartic field. For support
`[0,1,3,5,6,7,10,27]`, the extra anchor is wall 2 with label 3; the only nonright support contact is wall 27 with label 3 or 5. Nine independent linear contacts leave a quadratic norm equation. Its radical squares are respectively

\[
50+22\sqrt5,\qquad (245+111\sqrt5)/4.
\]

Both are nonsquares in the inherited field. This does not obstruct geometric walls: the exact extension arithmetic certifies them. Each has seven right pairings, two finite nonright pairings and twenty dotted pairings with the retained walls. Their positive mutual pairing is certified in the compositum.

The [independent review](fresh_h9_octic_escape_caps_20261015_review.md) reconstructs the contacts in reverse order using separate algebraic arithmetic. It verifies all 58 old-wall pairings, 1,800 facet witness signs, 60 endpoint signs, both old-escape cuts, the common new escape and mutual incompatibility. Generation and saved replay pass.

The [field classification](fresh_h9_octic_angle_20261015.md) supplies the complete finite-angle alphabet for these extensions and their compositum. A separate [152-candidate selection diagnostic](fresh_h9_octic_angle_20261015_selection.md) finds no numerically compatible pair in the saved finite bank, and the same escape survives all 152. That bank diagnostic is numerical and bounded; it does not exclude every wall in either field.

```sh
python scripts/fresh_h9_octic_escape_caps_20261015.py --verify
python scripts/fresh_h9_octic_escape_caps_20261015_review.py --verify
```
