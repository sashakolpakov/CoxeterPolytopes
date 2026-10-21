# Complete forest separation rejects the first full-face candidate

The first exact full-native-face vector is not a full feasible rule. Its
complete finite-forest maximum is

    722951420029/10000000000 = 72.2951420029 > 40.

An attaining forest is 5A1+D4+2D5+I2(5), with exact hub counts (5,1,0,0).
The old baseline there is 23. The certificate gives its complete literal
bonus row and right-hand side, and retains the 64 strongest distinct violating
maximizers. Another omitted forest, 9A1+3A4, has budget exactly 72. These are
failures of the candidate on omitted forest contexts, not new obstructions to
the funded H4 rule or to polytopes.

The calculation includes all 68 finite component types through rank 21 and
every admissible exact count of entire A1, I2(5), H3 and A4 components: 729 exact
contexts in 261 capped coefficient contexts. Connected bonuses, normalized
A1/I2(5)/H3 crosses, both A4 source-root orbits and A4-to-A4 pair bonuses are
included. The latter contribute one eighth of all sixteen ordered source-root
pairs per A4 component. This is the coefficient convention of the old producer.

For each remaining rank/context the exact DP retains the three best distinct
old forests. Only H4+A1+2E8 and D5+2E8 have decreased budgets, so a best unchanged
forest survives after excluding both. All four increased source forests are
evaluated explicitly. The D5 budget uses the worst of its active and inactive
states, namely the uniform saving min(1,2 beta). All source capacities in this
packet use the separately proved strengthened value one.

Generation passed in terminal chunk `92744c`. The guard recorded 1.163 seconds,
83,492,864 bytes combined peak RSS and 22,674 bytes output; worker peak was
51,789,824 bytes. Limits were 59 seconds, 1 GiB RSS and 1 MiB output. Initial
free disk was 197.93 GiB and available memory was 8.79 GiB. No optimizer ran,
and the process exited, released the lock and removed temporary output.

```sh
python3 -B scripts/fresh_compact21_full_candidate_forest_separation_20261018.py --verify
```

This packet records the first negative separator result. The next refinement
must add a maximizer for each violated exact context, rather than only these
64 globally largest rows. A new context-complete packet will preserve this
one unchanged and cross-check its coefficient extractor against the existing
literal forest support. No dimension improvement, Lean verification, external
review or full feasible weight rule is claimed.
