# Guarded D16 node3 verification ledger

The complete proof and exact rule are in
[the composition report](fresh_d16_node3_rank16_large_d_guard_transfer_20261018.md).
All new arithmetic ran in one foreground process at a time under the
shared exclusive lock,59s alarm, start/end resource gates and1GiB RSS checks.
All processes exited and the lock is free. No frozen predecessor changed;
no Lean or giant result JSON load was used.

|Check|Result|Evidence|
|---|---|---|
|Geometry and guarded positive-family generation|PASS|chunk6bd08b|
|Complete finite maximum generation|PASS:1487/4|chunk5a020c|
|Current retained-output generation|PASS:216 literals/564 later corners|chunk e3d6e3|
|Staged map-count label correction and three final regenerations|PASS|chunk b1591e|
|Fresh separate-process replay of all three final certificates|PASS|chunk141104|
|Independent new-rule checker, including budget|PASS|chunk d1e4ad|

The count-label correction distinguishes24 changed maps from48 total
saved maps in the same fourteen contexts. All24 changed maps were already
handled by the geometry filter. The rule, local proof and finite values
did not change. Only these new, not-yet-frozen D16 staging files were
regenerated; the completed A16 packet and all earlier evidence stayed fixed.

[The independent checker](../scripts/check_d16_node3_guard_independent_20261018.py)
imports the frozen, parent-reviewed A16 baseline for unchanged costs and
independently writes the new D16 delta. It imports no new D16 producer.
It passes65 streamed source-hash checks and70 positive-definite principal
checks, reconstructs full contact domains and native paying roles, verifies
all exact inverse/radical identities and forced completion entries, and
checks the absence of changed outputs in the recorded current interfaces.

For the three complete old top classes it independently verifies all new
table weights,179,130 Bellman inequalities,7,482 attaining predecessors
and2,116 explicit guard options. Their complete new maxima are370,371,371.
It separately verifies2,472 Bellman inequalities and114 predecessors for
the complete new-coupon cover, with maximum11513/32. The latest-increase
classification and streamed original38-class maximum summary complete
the global upper1487/4. The actual formal witnessA17+D16+D24 attains it.

The complete historical census and selected-map enumeration are inherited,
not rerun. Matrix/certificate execution does not replace the written
arbitrary-support proof or an independent mathematical read. Parent review
of the completed new composition is separately pending at freezing.
The dimension bound remains57 because371.75 exceeds curvature364.
