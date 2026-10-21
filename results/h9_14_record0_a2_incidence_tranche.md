# Record 0, compact schedule `a=2`: an eight-source incidence no-go tranche

## Result

Record 0 has exactly **36** source orbits of two missing four-sets after the
complete local-link, forced-shadow, strict-Macaulay, and full missing-graph
automorphism reduction.  The following eight source orbits admit no compact
Coxeter-nerve incidence completion:

```text
0, 8, 9, 10, 12, 23, 31, 34.
```

Each fixed-source exclusion is terminal in two independently constructed
models: one-worker OR-Tools CP-SAT returns `INFEASIBLE`, and Z3 returns
`unsat`.  The models agree exactly on every ground-set and constraint count.

This is a completed tranche, not a theorem for all of `a=2`.  The other 28
source orbits remain outside the certificate.  In particular, this result
does not exclude Record 0, does not change the global `134 + 8` compact-shadow
partition, and says nothing about finite-volume noncompact chambers.

The terminal scout ledger at the freeze is:

| compact-incidence status | source-orbit indices | count |
|---|---|---:|
| two-solver `UNSAT` | `0,8,9,10,12,23,31,34` | 8 |
| explicit `SAT` incidence witness found | `1,2,3,4,5,6,7,11,13,14,15,16,17,18,19,20,21,22,25,26,27,28,33,35` | 24 |
| `OPEN` after both recorded time limits | `24,29,30,32` | 4 |

Here `SAT` means only that the incidence overcover has a solution.  It does
not mean that Coxeter labels, a Lorentzian Gram matrix, compactness, or even a
hyperbolic realization survives.  For example, the first source-18 incidence
witness dies at the complete compact-label gate on its first A5 extension.
The `OPEN` rows are timeouts, never exclusions.

## Exhaustive source reduction

For the graph

```text
M?????????????o??
```

the exact source computation gives

| stage | count |
|---|---:|
| raw unordered pairs of locally admissible A4 sets | 386,760 |
| pair-capacity feasible | 384,780 |
| strict forced-shadow/aggregate viable | 384,780 |
| full automorphism orbits | 36 |

The orbit-size multiset sums to 384,780, so the quotient covers the entire
viable labelled ground.  The full group is the missing-edge graph
automorphism group `S_11 x C_2`: `S_11` permutes vertices `2,...,12`, while
`C_2` interchanges `0` and `1`.

For a fixed unordered A4 pair, the symmetry slice is also lossless.  The
isolate vertices split into four membership-column classes.  Adjacent swaps
generate the symmetric group on every class, and all valid lifts of active
swap and A4-row swap generate the quotient.  A global lexicographic minimum
in every finite orbit satisfies every installed generator comparison.

## Exact compact-incidence model

The schedule has 176 rank-nine facets and 38 independent missing five-sets.
The dual model selects complements of facets, which are five-subsets of the
14-node ground set.  It imposes all of the following simultaneously:

1. exactly 176 selected facet complements;
2. a facet witness for every required face of order at most four;
3. exactly 38 independent missing five-sets, including the forced A4 shadow;
4. closure of every nonface of order at least six by a missing set of order
   at most five; and
5. degree zero or two over every six-set, the simplicial ridge law.

The fourth condition is the compact Lanner-order gate.  Every minimal
nonface in the nerve of a compact hyperbolic Coxeter polytope supports a
classical Lanner diagram, hence has order at most five.  The model is an
overcover: it does not yet impose Coxeter labels or a Lorentzian Gram matrix.
Consequently, infeasibility is already a valid compact no-go for the fixed
source.

## Frozen evidence

The machine-readable certificate is
[`h9_14_record0_a2_incidence_tranche.json`](../data/enumerations/h9_14_record0_a2_incidence_tranche.json).
It records, for every excluded source, the literal A4 pair, all exact model
counts, the source stabilizer order and generator movements, both terminal
solver statuses, and strict theorem/non-theorem flags.

Production and independent replay:

- [`certify_h9_14_record0_a2_incidence_tranche.py`](../scripts/certify_h9_14_record0_a2_incidence_tranche.py)
- [`independent_h9_14_record0_a2_tranche_z3.py`](../scripts/independent_h9_14_record0_a2_tranche_z3.py)
- [`explore_h9_14_record0_a2_dual_facets.py`](../scripts/explore_h9_14_record0_a2_dual_facets.py)

## Scope ledger

| statement | status |
|---|---|
| eight fixed `a=2` sources are compact-incidence impossible | **proved** |
| all 36 `a=2` sources are excluded | **false / not claimed** |
| Record 0 is compact-Coxeter impossible | **not proved** |
| finite-volume noncompact Record 0 is excluded | **not tested** |
| an exact Gram survivor was produced | **no** |
| a new hyperbolic Coxeter polytope was produced | **no** |

The next exact step is to close or produce Coxeter-label/Gram survivors in
the remaining source orbits.  Incidence feasibility alone is not promoted to
a hyperbolic example.
