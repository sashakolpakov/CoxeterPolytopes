# Independent source review of the disappeared-PID monitoring repair

25 October 2026. **Source clearance accepted for the exact repair and
focused regression source below.** Root owns any resource-admitted test
execution. No test, process probe or proof job was run by this reviewer.

Root authored the repair and tests. I read the complete actual diff,
the full focused test source and the surrounding process-identity,
ownership, group-signal, pre-release probe and resource-runner interfaces.
I did not edit the monitor, tests, frozen sources or execution records.

## Exact bindings

|Object|SHA-256|
|---|---|
|Preserved original monitor|`517c0c78347c0aa81cd16489afcf7f11b7bce064490919b9cf5ff94c88e59cda`|
|`scripts/monitor_owned_proof_tree_20261024.py`, corrected successor|`b07f84fbcc5b35e484696f0515b8e233ea975f13c2da6ff75c0d3b54b99ed5db`|
|`scripts/test_owned_tree_disappeared_eperm_20261025.py`|`51dad3f76aa9b127724a83ee8acc7502d1676f8070a41c2f3785ef1bde2ffa71`|

The original and unaccepted first drafts are preserved under
`results/source_archive/owned_tree_monitor_before_disappearance_repair_20261025/`.
The first monitor draft4d86545c used stripped empty output, which could
accept arbitrarily large whitespace; it also let subprocess failure
replace the original identity diagnostic. I requested strict output
semantics and timeout/error regression coverage. The accepted successor
addresses both. The first draft and its nine-case test receive no source
clearance from this report.

## Exact accepted branch

Only Darwin `proc_pidinfo` returning length0 with errnoEPERM triggers
the new query. Existing length0 ESRCH/ENOENT absence handling is
unchanged. The new command is an argv call for the EXACT requested PID,
with no shell, no broad match and a three-second timeout:

`ps -p PID -o pid=`.

Only return code1 AND stdout exactly the empty string AND stderr
exactly the empty string yields `None`. Whitespace, newlines, oversized
whitespace, a printed PID, nonzero diagnostics, return0 or another
exit status all retain the original refusal. The code captures the
original errno before spawning ps, so subprocess activity cannot corrupt
that diagnostic. OSError and SubprocessError, including timeouts and
permission failures, fall through to that same identity refusal.

This is a bounded-time confirmation using the operating system ps
single-PID semantics. Accepted output is strictly empty; it is not an
enforced heap-allocation cap on `capture_output`. The selected command
does not request a whole process table. No alternative process identity
or birth is invented from its output, and it cannot grant ownership of
another process. A present denied process still refuses monitoring.

The exact diff does not change valid136-byte identity decoding, PID/birth
checks, Linux behavior, pre-child monitoring admission, owned-tree
discovery, trusted-group rules, signals or cleanup ownership. It changes
neither proof resource thresholds nor the bounded inner runners. The
existing documented polling and check-to-signal races remain; this is
not an atomic PID facility or a proof that every Darwin EPERM means
disappearance.

## Fourteen intended focused cases

The tests instantiate the reader without calling its real OS setup and
mock only its requested Darwin proc result and subprocess query. They
make no actual PID query and launch no mathematical worker. The changed
exception mocks use `side_effect`, so they exercise the exception path
rather than returning exception objects.

The14 INTENDED cases cover confirmed absence; present PID; PID output
despite return1; permission stderr; return2; empty return0; stdout
whitespace; stderr newline;65-space output; timeout; PermissionError;
the unchanged ESRCH/ENOENT cases; and malformed length135. Each EPERM
case asserts the exact query/timeout; the unchanged branches assert no
query. The source's normal-Python check protects its assertions.

These cases have not run in this review. They test the branch protocol,
not the occurrence or root cause of the previously observed kernel
failure. Root's retained failed executions and subsequent process
inspection remain separate evidence. No retry or mathematical PASS is
credited here, and successful scientific batches need no rerun merely
because the monitoring source changed.
