# Source review and higher-precision replay of the H9 H4 cap

The new86th wall in [the cap archive](../data/enumerations/fresh_h9_remote_next_20261013_caps.json) passes a source-level proof review and a100-digit replay. This reuses the constructor's exact and outward-rounded interval implementation; it is **not a separate arithmetic implementation**.

The frame calculation gives norm2 and the prescribed right, `pi/3` and `pi/10` relations. The replay certifies all76 remaining new/old relations as strictly dotted, preserves the strict center, verifies an actual `H4+5A1` finite vertex and the closed edge joining it to the prior finite vertex, and checks the remaining strict future spacelike escape. All269 new interval enclosures are contained in the corresponding archived80-digit enclosures. The86-wall chamber is still infinite-volume.

[Replay result](fresh_h9_remote_next_independent_20261013.json) and [checker](../scripts/fresh_h9_remote_next_independent_20261013.py):

```sh
python scripts/fresh_h9_remote_next_independent_20261013.py --verify
```
