# Status of the dimension-163 proof

The repository contains a **complete, internally checked computer-assisted
proof** that every finite-sided, finite-volume real hyperbolic Coxeter polytope
has dimension at most163. No unresolved mathematical step is currently
identified in its final reviewed chain. Arithmeticity, compactness, simple
cusps and a bound on the number of facets are not additional hypotheses.

The authoritative [final case split](fresh_nikulin_two_budget_frontier_20261015_dimension163.md)
binds the [ordered-budget argument and independent arithmetic review](fresh_nikulin_two_budget_frontier_20261015_review.md),
the all-rank estimates and the inherited geometric proof. The latter has an
explicit [internal geometric audit](fresh_dimension225_hostile_20261015.md),
including actual framing-wall incidence, common cusp-height ratios, endpoint
Schur budgets and the compact-facet alternative. Standard published narrowness,
classification and simple-polytope results are stated inputs.

The computational portions have exact certificates and separately implemented
checks. The geometric portions have written internal reviews. These reviews
were performed within the agent workflow; they are not external mathematician
refereeing. The proof is not formalized in Lean or another proof assistant.
External validation, publication priority and current-best literature status
remain separate unestablished matters.

The word *candidate* in frozen historical reports expresses that validation
status; it does not identify a remaining open case or unproved local inequality.
It should not be read as reducing the result to numerical evidence. Those
frozen reports and their hashes are preserved.

On2026-10-17 the final case/dependency checker was rerun successfully, checking
30 recursive manifests and83 source files. This rechecks that package's
integrity, exact final inequalities and case split; it does not rerun every
historical calculation or formally verify the geometric arguments.

```sh
python3 scripts/fresh_nikulin_two_budget_frontier_20261015_dimension163.py --verify
```

The gap between163 and the known examples motivates a stronger theorem. It
is not a missing step in the existing argument.
