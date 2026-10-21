# Retained numerical dual proposal and conservative repair

This packet records the authorized one-thread rerun of the small compact-rank-21
support successor. It retains numerical multipliers so subsequent rational
repairs need no additional optimizer. The final useful exact certificate is in
`fresh_compact21_priority_support_dual_repaired_20261018.json`.

The numerical objective was 40.004479162454544 on 2,805 rows, 6,121 columns and
17,387 nonzero coefficients. The initial repair rounded face masses down and
forest masses up, then scaled all face masses by one common factor. One tiny
numerical face multiplier demanded a column with zero forest capacity. That
common factor therefore became zero, producing only the valid lower bound
29287327459/1000000247. This is a deliberately conservative failed recovery,
not evidence that the support admits a value at most 40.

The retained numerical dual has 2,717 nonzero multipliers. No numerical
multiplier is used as a proof. The later packet clips the single offending
face coefficient locally and checks rational weak duality on every support
column.

The guarded run took 1.399 seconds, used at most 135,020,544 bytes of combined
driver/worker RSS, and emitted 100,750 bytes. Its hard limits were one thread,
59 seconds, 1 GiB combined RSS and 512 KiB output. Initial free disk was
197.88 GiB and available RAM was 8.82 GiB. The shared lock was released and no
temporary output was retained. The execution guard hash is recorded separately
from mathematical source pins, because changing resource-control tooling does
not change the stated rational inequalities.

Verification: internally generated and all 2,805 rows replayed with rational
arithmetic; no Lean translation or new dimension bound is claimed.

Commands, run separately under the shared lock:

```sh
python3 -B scripts/fresh_compact21_priority_support_exact_dual_20261018.py --verify
python3 -B scripts/fresh_compact21_priority_support_dual_repaired_20261018.py --verify
```
