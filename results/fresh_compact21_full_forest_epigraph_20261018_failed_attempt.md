# First bounded Bellman solve reached its numerical time limit

The first complete Bellman-epigraph solve did not produce a primal. HiGHS
returned status 1, "Time limit reached", with primal status None. This is a
numerical timeout, not an infeasibility proof or a dimension result. The final
tool chunk was `eb82ba`; the live session was 45739 and exited with code one.

The implementation reconstructed all 510 existing literal forest coefficient
rows and right-hand sides successfully before calling the solver. Structural
counting gave 243 capped context groups, 729 exact tuples, 1,659 auxiliary rank
states, 26,469 Bellman rows, 115 off-target rows, 727 ordinary budgets and nine
explicit budgets. Including native faces and existing support, there are
79,274 rows and 94,706 global column indices. The exact Bellman formulation
and its independent mathematical review remain valid; no numerical result was
inferred from the timeout.

The solve used one thread, an at-most-45-second solver cap with time reserved
inside a 59-second total alarm, and a live 1 GiB combined RSS watchdog. Cleanup
used the signal-safe final process-group kill after reaping the leader. The
process exited and released the shared lock, and no candidate or temporary
artifact remains. This first source did not persist its live peak RSS or exact
pre-solve build time on failure, so those are not reported as measured facts.

The next authorized strategy is distinct: eliminate bonus columns absent from
all negative native-face coefficients, use a fixed z=39.99999 and zero
feasibility objective, and call the IPM solver. It must prove the elimination
by checking all non-face bonus coefficients are nonnegative, record progress
before solving, and retain a compact failure ledger if needed. That successor
uses a new source stem. This failed attempt remains byte-identical.
