# Coupled A16 verification ledger

The new local theorem and exact costs are in
[the complete composition report](fresh_a16_double_large_d_transfer_research_20261018.md).
All arithmetic below used the shared exclusive lock,59s alarm, start/end
resource gates and1GiB RSS checks. All processes exited; no job remains.
Frozen predecessors were read only. No Lean was compiled.

|Check|Result|Evidence|
|---|---|---|
|Complete upper generation|PASS:16 classes,371|chunk d26066|
|Upper separate-process replay|PASS:16 classes,371|chunk2341f0|
|Physical generation|PASS:17 full seals;8 current affected literal cases|chunk e1d6c4|
|Physical separate-process replay|Same complete certificate|chunk e28069|
|Independent exact checker|PASS in0.44s|chunk ecdaf2|

The checker is
[check_a16_double_large_d_independent_20261018.py](../scripts/check_a16_double_large_d_independent_20261018.py).
It imports only standard-library arithmetic and the resource gate;
no mathematical producer is imported. It checks48 streamed source hashes,
153 positive-definite principals,17 exact core inertias and radical/inverse
relations, freshly reconstructed contact domains and native marked roles,
432 literal payment sums, all eight original whole D13 alternatives,
and the complete16-class/77-reserved-D arithmetic cover. It verifies
38,656 Bellman inequalities and1,824 attaining predecessors independently.

The arithmetic upper is a relaxation: it discards the A16/D17 negative
changes while including all positive costs. It proves that every newly
increasing forest costs at most371. The unchanged forest2D16+D25 costs372,
so372 remains the exact global formal finite maximum. The particular
forestA1+2A16+D24 costs363.5, and2A16+D25 costs365.5. The dimension bound
remains57. Saved-map completeness and the earlier universal routing
theorems are inherited inputs, and no historical census was rerun.

This ledger records execution and a separate exact checker. The universal
guard-persistence and routing proof is written in the linked report;
finite arithmetic alone is not presented as that proof. Parent mathematical
review of the proposed1010 geometry preceded the final report; review of
the completed composition is separately pending at the time of freezing.
