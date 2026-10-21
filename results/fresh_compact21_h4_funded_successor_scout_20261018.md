# A below-40 exact primal on the funded-H4 support successor

The proved H4 angle transfer clears the previous selected-support obstruction.
At transfer value t=1/4, the selected successor has the verified rational primal

\[
z=\frac{49994751497}{1250000000}=39.9958011976<40.
\]

The exact margin below 40 is 5248503/1250000000. This is **not yet a dimension
bound**: full quadrilateral and finite-forest separation remains outstanding.
No full-model feasibility is claimed.

The new model retains the 2,805-row D5 active/inactive support successor. It
adds the three previously absent source forests A5+2E8, A1+B4+2E8 and B5+2E8,
and the row t<=1. The H4+A1+2E8 budget gains coefficient -t, and the A4 source
gains +t. The other source coefficients are +2t, +t and +t. The separately
frozen H4 geometric packet proves the underlying face payments, including all
simultaneous H4 corner modifications and compatibility with the prior D5 rule.

The sparse support has 2,809 rows, 6,131 occurring columns and 17,439 nonzero
coefficients. The compact JSON retains 2,669 positive rational bonus values,
the rational transfer value, 2,720 exact negative dual multipliers, and all
2,722 nonzero numerical dual multipliers for later recovery without another
optimizer. Every selected row is checked literally against the rational primal.
All literal weak-duality columns are also checked. The exact dual lower bound
is 13331933563155/333333333416, close to the numerical objective
39.99580119525319. No exact optimality equality is asserted.

The numerical solver only proposes a vector. Positive bonuses are rounded up,
all negative face constraints are made feasible by a common nonnegative bonus
scaling, t is recovered rationally in [0,1], and z is recomputed as the largest
exact selected forest budget. The independent rational replay then checks all
2,809 constraints. The dual repair rounds capacities conservatively, clips only
overloaded face demands, and adds any needed t-cap multiplier before its exact
column replay. A numerical feasibility tolerance is never used as proof.

Generation completed in terminal chunk `0958af`, exit zero. The guard recorded
3.81 seconds, one thread, 136,167,424 bytes peak combined driver/worker RSS,
243,556 bytes worker output, and 104,185,856 bytes worker peak RSS. Hard limits
were 59 seconds, 1 GiB combined RSS and 512 KiB output. Initial free disk was
198.52 GiB and available RAM was 7.17 GiB. The process exited, the global lock
was released, and temporary output was removed. The final certificate is
244,917 bytes. No Lean process or unrestricted proof scheduler was started.

```sh
python3 -B scripts/fresh_compact21_h4_funded_successor_scout_20261018.py --verify
```

This command only replays rational constraints and source hashes; it never
imports NumPy/SciPy or reruns the solver. Mathematical input pins are separate
from the mutable execution-guard provenance.

Verification status: internally generated; primal and dual checked in both the
worker and driver; fresh separate-process replay can use the command above.
The H4 geometry has also received the root agent's mathematical review.
No complete Lean translation or external peer review is claimed. The immediate
next task is a bounded, streaming separation of every old native face case and
all finite forests, with the new t and D5 rules composed exactly.
