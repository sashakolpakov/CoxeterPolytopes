# Review of the four-anchor H22 bridge extension

The [four-anchor verifier](../scripts/fresh_h22_contact_bridge_20261013_big_verify.py) passes mathematical source review and a fresh saved-result replay. Its [certificate](fresh_h22_contact_bridge_20261013_big_verify.json) completes the direct-contact test for the certified 305-cap component: ten of the 88 edge families are already empty, and all **312** additional contact intervals in the other 78 families are strictly empty. No issue was found.

The four positive-height anchors have original type indices `0,2,6,22`, physical indices `812,816,824,854`, and squared norm two. These indices are reconstructed from the frozen 860-wall ordering. The separately source-bound component certificate identifies its other 301 members as precisely the positive parent null caps. Thus the extension tests every member missing from the [earlier 301-anchor proof](fresh_h22_contact_bridge_20261013_review.md).

The verifier reconstructs the four positive-definite rank-22 incident Gram matrices, checks their entrywise nonnegative inverses, and rebuilds every Lorentzian edge frame. It rechecks the lower bounds and the ten empty families against the frozen primary certificate. The inherited normalization is `r=(t(v-e)-t^-1(v+e))/sqrt(2)`, with `t>1`, and the positive-height restriction remains explicit.

For an anchor of norm n, with `a=(v-e,w)`, `b=(v+e,w)` and `T=t²`, finite or parallel contact requires `(aT-b)²<=4nT`. For every new test, the shared exact helper proves `a<0` and a positive discriminant before enclosing the larger quadratic root. The verifier then checks that this outward upper bound is **strictly below** the necessary lower bound. No exceptional sign case is silently skipped. Equality, singleton intervals, and parallel endpoints are retained until this strict comparison excludes them; no numerical width threshold enters the certificate.

Together, the two certificates exclude a positive-height cap that is right to a 21-wall edge at one of the four specified vertices, remains acute with all 860 retained walls, and directly meets this 305-cap component at a finite angle or in parallel. They do not exclude indirect chains of new caps, nonright edge profiles, other vertices or faces, lower-height caps at these vertices, or replacement of retained walls. This is not a general obstruction to completing the H22 chamber.

This audit is source and theorem review plus an external replay of the author's exact verifier, not a second implementation of interval arithmetic. The earlier 301-anchor proof was not rerun or modified. [Review metadata](fresh_h22_contact_bridge_20261013_big_review.json) records the reviewed source and certificate hashes.

```sh
PYTHONPATH=src:. python scripts/fresh_h22_contact_bridge_20261013_big_verify.py --verify
```
