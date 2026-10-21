#27-wall fixed check: generation and separate replay passed

19 October 2026. This note records actual execution separately from the
unchanged source-preparation and pinned mathematical completion reports.
Both runs are now complete; the completion addendum records the separate
replay. The initial RAM deferrals below are historical events.

The coordinator's guarded generation passed with source
`scripts/check_h19_twenty_seven_wall_seed_20261019.py`, SHA256
`4fed900c65b747dbf99226d16fd0e99b7e1f3106fd3269e0ebce9882455f9cf1`.
Its actual result checks all351 wall pairs,27 actual-facet projections,
three finite-vertex witnesses, the positive alternate-root product and
the explicit feasible spacelike escape. The nine initial p directions
have eight finite first-hit endpoints and one spacelike endpoint.
This is not a full vertex or cone census.

The generation took0.234 seconds and approximately56.74MB process-tree
RSS. It passed start gates at187.49GiB disk and6.19GiB available RAM.
The single-lock wrapper retained its59-second/1GiB limits and parent-death
handling. An earlier generation attempt had been deferred before launch
at5.8GiB RAM; that attempt ran no checker body and made no certificate.

|Actual output|SHA256|
|---|---|
|`results/h19_twenty_seven_wall_seed_20261019.json`|`b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d`|
|`results/h19_twenty_seven_wall_seed_generation_20261019.json`|`4710fe2a23b43852cddc97e434aa6a2d7d5595bfdef9ca1327ad16993b2fa8b8`|

**Historical status before recovery: replay was pending.** Its attempted launch was deferred at
5.8GiB available RAM, below the6GiB start gate. No replay child, retained
lock or replay ledger exists from that attempt. A future actual replay
had to receive its own fresh resource gate; elapsed time did not permit
bypassing it. That deferred attempt supplied no replay PASS.

The executed algebra supports the unchanged written conclusion: the
new q wall preserves the two inherited full41-vertex H8 facets, p and
p_s, and adds p_q, giving at least85 actual finite vertices. The explicit
feasible spacelike point and strict centre imply an open ideal-boundary
patch, so this27-wall seed still has infinite volume. The complete H8
facet enumeration remains inherited. No finite-volume example, complete
mirror census, nonreflectivity theorem or dimension improvement follows.

## Completion addendum,19 October 2026

After available RAM recovered, the coordinator's **separate replay
actually passed**, in0.234 seconds at approximately57.11MB process-tree
RSS, starting at187.49GiB disk and6.19GiB available RAM. The source and
certificate bytes above remained unchanged. Both generation and replay
are now complete; the earlier prelaunch deferrals remain recorded as
historical events and were not treated as executed checks.

Replay ledger: `results/h19_twenty_seven_wall_seed_replay_20261019.json`,
SHA256 `ea2f724679ca6cc3f95423a495f8a18913257e0c906b872d371755e92628fa45`.
