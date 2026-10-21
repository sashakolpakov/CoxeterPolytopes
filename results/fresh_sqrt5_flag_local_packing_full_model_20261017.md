# Preserved rank-13 model with vertex–coface flags and finite packing cuts

**FINAL/FROZEN.** The enlarged partial model still cannot attain budget 12.
One bounded solve and exact rational recovery give

\[
 {345782670398463480\over27027032875363841}
 \le T_*\le
 {3198492461931\over250000000000},
\]

or

\[
 12.793956036279\ldots\le T_*\le12.793969847724.
\]

This successor adds vertex–coface incidence information while preserving every
inherited face obligation, variable, and budget. The flag refinement breaks
the old fixed witness, but the optimizer finds a different obstruction.
Neither outcome is a dimension conclusion.

## Exact model

| Quantity | Count |
| --- | ---: |
| Inherited face obligations, all preserved | 82,806 |
| Complete vertex budgets retained | 6,558 |
| Inherited variables, all preserved | 350,744 |
| Added vertex–coface flag transfer columns | 77,391 |
| Added positive packing columns | 175 |
| Total variables | 428,310 |
| Active variables with some positive face entry | 162,484 |
| Total signed transfer columns | 94,481 |
| Exact face-matrix nonzeros | 801,967 |

The 175 positive columns consist of 148 previously maximizing
[all-subset cuts](fresh_sqrt5_all_subset_port_packing_20261017.md), four
[intersecting-orbit clique cuts](fresh_sqrt5_intersecting_block_cliques_20261017.md),
and 23 new cuts from the
[complete local packing-cone check](fresh_sqrt5_all_local_simplex_packing_cones_20261017.md).
This is an explicitly finite cut family: it does not insert all subset
inequalities or the entire local charge cone.

The [flag transfers](fresh_sqrt5_flag_refined_base_capacity_20261017.md)
use the complete maximal-simplex star, its vertex, and an immediate coface.
Each column has a negative simplex coefficient and positive actual-prism
coefficients. Its positive entries are checked against the frozen map from
all 1,858 selected provenances to their expanded quad rows. There are 13,345
nonzero positive flag entries. The previous shared, forest-labelled, and
other signed columns remain intact.

Every added rational positive column is multiplied by a positive common
denominator in both its face coefficients and its vertex cost. The exact
scale is saved. New positive packing features are zero outside their exported
simplex classes and zero on quadrilaterals. Flag positive coefficients are
zero outside the selected actual quad provenances; their negative simplex
terms remain. These are valid weaker restricted payments, not an assertion
that the omitted geometric coefficients vanish.

The preparation checks every inherited row coefficient, vertex-budget
coefficient, variable record, and old signed-column record for exact
preservation. It also checks the actual provenance map and all integer
scalings. The resulting model digest is

`742bf404f3095d33583d5b42dd3ffca9160bcede90896ba1b1e3639b3767680c`.

## One solve and exact recovery

The active numerical matrix has 89,364 constraints, 162,485 columns including
the maximum-budget variable, and 876,528 nonzeros. The single `highs-ipm` call
terminated optimally in 52.146 seconds, within its 300-second cap, after 12
iterations. Its numerical objective was `12.793969847152075`. The saved output
records `solver_invocations=1`; the entrypoint refuses a second solve. No
repeat solve was run.

Rational primal recovery rounds the numerical vector and adds a common
nonnegative correction only to inherited positive pair/port variables.
Every face row has a checked positive such coefficient, so the maximum exact
deficit determines a valid repair. All rows and vertex budgets are then
recounted rationally. This yields the exact upper bound above, which is
strictly greater than 12.

Rational dual recovery rounds nonnegative face masses and monotonically
raises each containing-simplex mass to satisfy every signed transfer.
Every signed column has positive terms only on quad rows and one negative
simplex term, so these repairs cannot invalidate another signed column.
Recomputing every positive-cost column yields valid vertex masses and an
exact lower bound. The standalone verifier additionally checks multiplier
nonnegativity, indexing, uniqueness, and both target-status flags. This yields
the exact lower bound above, also strictly greater than 12. The new dual has
4,428 positive face multipliers and 2,104 positive vertex masses. These are
the masses actually certified; the previous fixed-alpha local covers are
not imported as certificates for this changed witness.

## Sources, reviews, and remaining scope

- [Prepared source and guarded solve/replay entrypoint](../scripts/fresh_sqrt5_flag_local_packing_full_model_20261017.py)
- [Saved prepared model and both exact certificates](fresh_sqrt5_flag_local_packing_full_model_20261017.json)
- [Independent source and mathematical review](fresh_sqrt5_flag_local_packing_full_model_20261017_review.md)
- [Frozen predecessor model](fresh_sqrt5_labelled_packing_full_model_20261017.md)

The root and frontier agent independently read the preservation, sign,
scaling, restriction, and recovery arguments and found no gap. Those reviews
do not duplicate the optimizer. A separate fresh saved replay rebuilt the
whole exact model, matched its digest, and verified both rational certificates.
Syntax, all fifteen direct prepared-source hashes, and local report links
also passed.

Reproduce the saved verification without optimization using
`python3 scripts/fresh_sqrt5_flag_local_packing_full_model_20261017.py --verify`.

This remains the compact arithmetic `Q(sqrt(5))`, rank-13 partial-model lane.
The retained 6,558 vertex forests are not the full finite-forest universe,
and the selected face domain omits other geometric cases and fixed
discriminant classes. Even an exact primal budget at most 12 would require
further coverage before a dimension conclusion. No fixed-alpha packing
cover is silently reused after changing the face masses.
