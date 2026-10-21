# The two priority rows remove the latest saved obstruction above forty

18 October 2026. **Exact saved-dual repair, not a primal feasibility proof.**
Replacing the D5+2E8 budget by the two geometrically justified
active/inactive priority rows lowers the best repair of the latest
saved A4-source dual, with every other multiplier held fixed, from
40.02199019 to

  7158707133073/179000000000 = 40-1292866927/179000000000,

approximately 39.99277728, strictly below forty. Thus that particular
saved obstruction does not survive this successor. A different dual
above forty may still exist. No optimizer, complete model rebuild,
new primal certificate or dimension theorem is claimed.

The proposed rows are proved in the complete extreme-E8 consumer
packet, together with the no-half-unit strengthening of the active
D5 rule. They preserve the old unpaid quadrilateral rows and all
triangle and long-face routes. The old nonnegative bonus variables
remain the variables of the successor. No fixed half-units are
added to the golden-field baseline.

## Literal exact calculation

The source is the frozen
[A4-source model](fresh_sqrt5_dimension21_a4_source_orbits_20261016.md).
The [streaming checker](../scripts/fresh_compact21_extreme_e8_priority_saved_dual_20261018.py)
reads its connected orbits, variable ordering and 2,804 sparse
dual terms one bounded object at a time, never loading the 26-MiB
file. It identifies canonical pair (0,6), reconstructs all literal
column loads and checks weak duality on every touched column.
There are 93,045 bonus variables and 6,121 touched columns including
the upper-bound coordinate. The extreme-pair bonus has column 10350;
the D5+2E8 forest has term index 2374.

Write b for that bonus. Its old forest coefficient is two and its
old baseline is 36. The new rows are

  old_bonus_total-z<=-35,
  old_bonus_total-2b-z<=-36.

The old forest mass is mu=1462087/50000000. Its total old column
load on b, including all face and other forest terms, is

  L=-4386337/76075000000.

Split the mass into lambda on the active row and mu-lambda on the
inactive row. Every other column load stays unchanged. The b load
increases by 2(mu-lambda), so the exact least admissible active
mass is

  lambda=max(0,mu+L/2)=1111186101/38037500000.

The inactive mass is 4386337/152150000000. Both are nonnegative
and sum to mu. The active row has right-hand side one larger
than the old row, so the dual objective decreases by lambda.
All 6,121 literal column inequalities remain valid after this
change. Direct summation of the changed right-hand sides agrees
with old_objective-lambda and gives the value stated above.

This split is best among repairs which hold all other old
multipliers fixed and preserve total target-forest mass. It is
not a dual optimum for the full successor. The old face and
forest coefficient interpretation remains inherited from the
frozen model; this checker independently verifies the literal
rational arithmetic and the exact two-row replacement, without
rebuilding the geometric input catalogue.

## Reproduction and next bounded calculation

```sh
python3 -B scripts/fresh_compact21_extreme_e8_priority_saved_dual_20261018.py --verify
```

The [small certificate](fresh_compact21_extreme_e8_priority_saved_dual_20261018.json)
pins the frozen old model and new geometry coefficient packets.
Generation passed in chunk `b29a3d`, in 0.219 seconds, peak RSS
29,032,448 bytes. Fresh replay passed in `226fb8`, in 0.169
seconds, peak RSS 28,983,296 bytes. Both held the shared lock,
used a 59-second alarm and checked RSS below one GiB. Each input
object is capped at two MiB; source hashes are streamed. No
temporary files, bytecode cache, numerical solver or Lean process
was used. Independent review of this diagnostic is pending.

The next useful scout is the sparse active system: the 2,804 old
positive-support rows plus the second target-forest row, on the
6,121 touched coordinates. This is a small proposal problem.
Any resulting candidate must subsequently be checked exactly
against every full-model face case and complete finite-forest
separation. Passing only the old support is not global feasibility.
No nontrivial optimizer has been launched.
