# Existing-proof LP recovery: 57, 17 and 15

For the self-contained mathematics, notation and input/consequence map for
all three cases, read [the uniform LP paradigm](UNIFORM_LP_PARADIGM.md).
This note records the completed executions and their scope.

The user's clarified task is to recover the existing proofs through LP,
retaining their proved local inequalities and complete domains. In that
scope the finite-volume recovery was already complete. The compact and
compact arithmetic global LPs now also have actual solver proposals,
complete exact rational certificate checks and fresh byte-for-byte model
regeneration. No smaller bound or primitive Gram relaxation is claimed.

## Results and actual execution

| Case | Excluded dimensions | Certificate status |
| --- | --- | --- |
| General finite volume | Every n >= 58 | Existing exact minimax value 377 and accepted all-dimension incidence/tail composition; retained without rerunning unchanged work |
| General compact | Every n >= 18 | New global capacity-dual certificate, with all actual ranks 18–29 represented and the explicit published ceiling 29 handling the tail |
| Compact arithmetic | Every n >= 16 | New joined capacity-dual certificate: arithmetic ranks 16/17 plus the unchanged compact 18–29 catalogue and the same external tail |

The [finite-volume completion](../../audit/joint58_complete_lp_execution_20261110.json)
binds matching exact feasible and dual certificates for its specified family.
There was no numerical optimizer run for that older certificate. Its
dimension-58 calculation is accompanied by separate valid higher-dimensional
and ideal-vertex arguments; larger-dimensional exclusion is not inferred
from dimension 58 alone.

The new jobs use the same [global assembler](../../scripts/run_existing_proof_energy_lp_20261118.py)
and the unchanged exact sparse LP transport/checker. Their actual sizes are:

| Master | Branch domains | Local contradiction leaves | Primal variables / rows | Dual variables / rows |
| --- | ---: | ---: | ---: | ---: |
| Compact | 829 | 151 | 3,471 / 4,876 | 3,157 / 2,224 |
| Compact arithmetic, including compact tail | 917 | 167 | 3,869 / 5,147 | 3,324 / 2,488 |

`EXACT_LP_FEASIBLE` in the execution files refers to the **capacity dual**.
Its checked point proves the normalized strict-budget primal infeasible.
Every row, bound and rational coefficient is checked; a native solver status
alone is not used. A second execution regenerates both complete LP files
and the complete catalogue and independently rechecks the rational point.

| Actual job | Watchdog elapsed seconds | Peak process-tree RSS bytes |
| --- | ---: | ---: |
| Compact export, solve and exact check | 0.454 | 131,645,440 |
| Compact fresh regeneration and exact replay | 0.222 | 74,137,600 |
| Arithmetic joined export, solve and exact check | 0.469 | 134,086,656 |
| Arithmetic fresh regeneration and exact replay | 0.231 | 75,497,472 |

All four actual watchdog records are PASS. They use the existing 59-second,
1-GiB, one-thread runner and shared exclusive lock. Their exact byte-ratio
admission samples pass the unchanged 2%-disk/15%-RAM policy. Temporary replay
models are removed by their owning context. Completed models and certificates
remain evidence, totalling only a few megabytes.

The complete compact evidence is under
[`results/existing_proof_compact17_energy_master01_20261118/`](../../results/existing_proof_compact17_energy_master01_20261118/),
including `primal.jsonl`, `dual.jsonl`, `catalogue.json`, `candidate.json`,
`completed.json` and `verified.json`. The arithmetic evidence is under
[`results/existing_proof_arithmetic15_energy_master01_20261118/`](../../results/existing_proof_arithmetic15_energy_master01_20261118/).
The [source/root acceptance](../../audit/existing_proof_energy_lp_recovery_source_root_acceptance_20261118.json)
binds the independent reviews, exact sources and execution contract.
The [independent execution review](../../audit/existing_proof_energy_lp_execution_independent_review_20261118.json)
checks all four monitored jobs and opaque artifact bindings. The final
[root completion acceptance](../../audit/existing_proof_lp_recovery_completion_root_acceptance_20261118.json)
accepts all three existing-proof LP recoveries at their retained inputs.

## The common energy LP for the two compact cases

Each branch b is a complete admitted local configuration at its actual
dimension. For compact polytopes one original farthest vector has
nonnegative orthogonal projection energies e_j with sum strictly below one.
A proved local estimate supplies a rational row

    sum_(j in S_l) e_j >= c_l.

The support S_l identifies the actual components whose energy it uses.
Joint inequalities share those components' variables. Several estimates on
the same component are simultaneous lower bounds, whose effect is their
maximum; they are not separately added costs.

For every branch use nonnegative mass tau_b, scaled energies E_bj and
deficit delta_b, with

    sum_(j in S_l) E_bj >= c_l tau_b,
    delta_b = tau_b - sum_j E_bj >= 0,
    sum_b delta_b = 1.

A hypothetical actual configuration lifts to this LP by dividing its
energies and its branch mass by its positive energy deficit. Conversely a
positive deficit in the LP identifies a branch whose energy consequences
can hold below one. This normalization avoids an artificial fixed epsilon.

The solver finds nonnegative local weights alpha_bl and a common eta >= 1:

    sum_(l with j in S_l) alpha_bl <= 1,
    sum_l c_l alpha_bl >= eta.

The first row is the capacity of each original energy component. It prevents
repeated use of a defect or any other shared support. These inequalities
give sum E_bj >= tau_b for every branch, contradicting the normalized
positive total deficit. The numerical objective maximizes available local
payments while preserving this feasibility problem; no numerical optimality
claim is needed for the exact contradiction.

Some local proofs end with a null-vector, native-contact, arithmetic or
acuteness contradiction rather than an energy floor. Their complete matched
domains receive tau_b = 0. They are disclosed imported local theorems; no
energy cost is invented for them. In particular the 151 compact and 16
additional arithmetic leaves are not 167 new LP-derived geometric proofs.

## Inputs, generated structure and consequences

The following remain **proved inputs**: original-vector geometry and
orthogonal energy accounting; spherical classification and actual defect
compatibility; the named local cost inequalities and their complete
applicability domains; explicitly matched local contradiction theorems;
and the published compact ceiling 29. Arithmetic additionally retains its
field, computed-output and applicable lattice/group premises.

The compact cover is **generated uniformly** from component families.
It separates no dominant component, at least two, and exactly one. A
dominant is A_r with r >= 7 or B_r/D_r with r >= 9. Its energy is at least
one half. Every nondominant spectator of total rank S has aggregate lower
bound (S-1)/14, allowing the unique permitted defect debit. Coarse branches
keep these rows explicitly. Any context needing finer information retains
every spherical spectator partition and every actual defect profile; all
ordinary dihedral parameters remain symbolic over their complete domains.
No historical eleven-case forest list is assumed. The precise rules and
source anchors are in the [compact formulation](compact17_complete_energy_lp_recovery_20261118.md).

The resulting **LP consequences** include exclusion of zero or multiple
dominants, exclusion of sufficiently large spectator remainders, and the
simultaneous closure of all remaining branch energy accounts. We do not
insert n <= 17, n <= 15 or an entire exact-rank exclusion as a local row.
The large-component and whole-component contradiction leaves remain explicit
local hypotheses in their original domains.

Arithmetic uses the **already proved pruned rank-16/rank-17 covers as inputs**:
19 and 14 source groups expand to 49 and 39 complete parameter-family
branches. Their pruning is not derived by this initial LP. The [arithmetic
formulation](compact_arithmetic15_global_energy_lp_recovery_20261118.md)
lists every row, alternative, computed attachment and local exit. The same
assembler handles those branches and the compact higher ranks. Thus the
global linear combination is uniform, while the retained condition families
have their own geometric and arithmetic source theorems.

Some supplied rows are already visibly redundant: when a stronger floor is
available on the identical support and domain, its weaker predecessor adds
no energy restriction. In particular the prefix-plus-gain floor dominates
its corresponding prefix floor. No general minimal-constraint claim is made;
a zero weight in one certificate would not establish logical redundancy.

## Threshold monotonicity and stopping point

Freeze the complete catalogue and its local rows. At threshold N retain
exactly the branches whose actual dimension is at least N. Increasing N
only removes branches. A feasible point at the larger threshold extends by
zero masses to the smaller one. Therefore infeasibility at N implies
infeasibility at every larger threshold in the catalogue. The independent
published ceiling excludes every dimension above 29. This proves the full
dimension bound without assuming that faces inherit Coxeter or arithmetic
properties.

The new uniform global LP recovery is complete in the user's selected
**existing-proof-input scope**. The task stops here. Deriving the local
theorems or arithmetic pruning anew from a smaller primitive condition set,
minimizing the input families, changing the objective for new bounds and SDP
strengthening await the user's next instruction.
