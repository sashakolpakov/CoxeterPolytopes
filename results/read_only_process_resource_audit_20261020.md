# Read-only process resource audit

20 October 2026; `/root/resource_cleanup`, at root's request after the
resource gate deferred a scientific child. This is a sampled process
observation, not a proof execution or a claim of recovered memory.

The sandbox denied the initial `ps`. Approved process-monitor access was
then used for `ps` fields **PID, PPID, RSS and executable name only**, and
`lsof -a -p <identified PIDs> -d cwd -Fn` for working directories. No full
command arguments, environment variables or unrelated file descriptors
were inspected. The two reads occurred close together; RSS is a snapshot.

| Process | PID | PPID | RSS, KiB | Working-directory attribution |
|---|---:|---:|---:|---|
| Lean4.32.2 | 22644 | 79589 | 5,578,800 | `/Users/sasha/math_gkm/.scratch/concept-driver-20261020-general-0001/deployment-preflight/positive/source` |
| Python parent of that Lean | 79589 | 1 | 17,792 | `/Users/sasha/math_gkm` |
| Repository-directory Python | 23071 | 1279 | 22,864 | `/Users/sasha/CoxeterPolytopes` |
| Codex app-server parent | 1279 | 1277 | 897,024 | Executable name identifies the Codex app server; no cwd inspection was needed |

The large Lean RSS is approximately5.32GiB and belongs to the separate
`math_gkm` workspace, rather than CoxeterPolytopes. The other matched
Python processes had `/Users/sasha/bayesilisk` working directories. No
`lake` process or repository-directory Lean process was present in the
snapshot. The small repository Python was approximately22.3MiB, parented
by the app server; name/cwd data alone do not establish a reclaimable
scientific orphan. Root subsequently reported its admission watch
terminal without starting a proof child and considered that small Python
likely related to the watch; that identification is root's assessment,
not a command-argument inspection by this reviewer.

**No large Coxeter-owned proof process was established as reclaimable.**
No process was signalled, stopped or restarted; no foreign workspace or
application was changed. No cleanup, proof worker, Lean, TeX or test was
run. No resource limit was raised. The audit created only this note.
