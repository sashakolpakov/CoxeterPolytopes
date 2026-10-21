# Independent check of the D19 fork-height corollary

21 October2026. **Written mathematical PASS.** This checks only the
elementary corollary from the accepted complete centre theorem; it does
not reopen or replace its full review. No jobs or source edits were used.

Exact bindings:

| Source under `results/` | SHA-256 |
|---|---|
| `arithmetic_d19_centre_complete_actual_profile_floor_20261021.md` | `3d88932839dbd25efa3d2f4ce0d433d0b779ca15a79721021579f6669abcdb14` |
| `review_arithmetic_d19_complete_centre_floor_and_fork_corollary_root_20261021.md` | `c16a19d7f6f10644ee6da4792182770793e887fdeb02af42c48acb0004a09cd6` |
| `review_arithmetic_d19_complete_centre_actual_profile_floor_outside_author_20261021.md` | `66cd099665d3c1b4a6032ff50ff606d37b26583255b619b44e4d7d069ed5ceae` |

All three current hashes match. The reviewer previously supplied the full
outside-author centre review. The derivation below is independent of the
root's corollary calculation, conditional on that accepted centre theorem
and its actual compact arithmetic first-type whole-D19 setting.

For the SAME original chamber vector x, positivity against the two unit
fork roots (e18-e19)/sqrt2 and(e18+e19)/sqrt2 gives
x18>|x19|. With S17=sum(x1,...,x17), the standard D19 duals yield

    h17=sqrt2*S17,
    h18=(S17+x18-x19)/sqrt2 > h17/2,
    h19=(S17+x18+x19)/sqrt2 > h17/2.

The accepted strict centre floor therefore gives either fork f

    h_f > sqrt(17/2)-21/80.

Both sqrt(19/2) and sqrt(17/2) exceed5/2, so rationalization gives

    sqrt(19/2)-sqrt(17/2)
      =1/[sqrt(19/2)+sqrt(17/2)] <1/5=16/80.

Consequently BOTH forks satisfy the claimed strict inequality

    h_f > sqrt(19/2)-37/80.

Since37/80<12/25, each also satisfies the loss-.48 floor. Every step
uses the same x, and no separately optimized vector or added energy
estimate enters. This is a height consequence only: it excludes neither
a literal cap satisfying the floor nor the whole D19 forest, and supplies
no additional long-prefix floor or dimension18 conclusion. The frozen
centre source and its previous reviews remain untouched.
