# Independent review of the ordered-budget dimension-163 candidate

**PASS.** The new ordered endpoint-budget argument is valid within the
explicitly inherited geometric framework. The
[review checker](../scripts/fresh_nikulin_two_budget_frontier_20261015_review.py)
and [certificate](fresh_nikulin_two_budget_frontier_20261015_review.json)
independently reconstruct the new arithmetic, without producer, project, or
prior-review arithmetic imports. The checked
[primary source](../scripts/fresh_nikulin_two_budget_frontier_20261015.py) and
[data](fresh_nikulin_two_budget_frontier_20261015.json) remain unchanged.

This is an internally reviewed research proof candidate, not a claim of
publication priority, current-best literature status, or formal verification.
The primary report's pending-review statement is historical; this separate
artifact supplies that review without changing the frozen report.

## Independent arithmetic and domain coverage

The checker builds literal norm-two affine Grams for all relevant types,
computes their positive kernels, and grounds the rational weighted Laplacian
\(\operatorname{diag}(\delta)G\operatorname{diag}(\delta)\). A separate
Fraction Gaussian elimination gives the Euclidean point-potential Gram.
For every omitted node in the finite census, its inverse physical spherical
Gram is independently checked against

\[
J_{kl}=\frac{\delta_k\delta_l}{2}
       \bigl(R(i,k)+R(i,l)-R(k,l)\bigr).
\]

All inverse entries are nonnegative and every diagonal is at least one half.
Since every nonzero contact is at least one, a support with endpoint cost at
most 3/2 must be a subset of the nodes whose inverse diagonal is at most 3/2.
The review enumerates every such subset and tests its **exact physical** cost,
including all cross terms. This separately confirms the primary's radical
lower-bound exclusion logic. No feasible support is lost by an interval test,
and the retained costs are exactly rational.

The finite domain contains 91 models, 1,116 retained-support pairs, and 3,630
continuous full-support domains. Both omitted-entry alternatives, zero and
at least one, are included. Empty full rows are excluded only because the
inherited geometry gives positive affine height. The zero-retained-cost
cases satisfy the independently checked point-energy bound W≥T.

For each continuous domain, the saved barycentric vector and equality
multipliers are checked against a separately reconstructed Gram. The
nonnegativity, simplex constraints, gradient slack and complementary
slackness prove the global convex minimum. This verifies all 3,630 exact
KKT certificates, including singular minimizer sets. It does not depend on
trusting a numerical optimizer or finding a new minimizer independently.

If the full support forces heights at least s and t and the actual common
height ratio is at most B=21/20, then

\[
\min(\sigma,\tau)^2\ge
\max\{\min(s,t)^2,\max(s,t)^2/B^2\}.
\]

This bound is increasing in both support-height lower bounds. The exact
rational enclosures for their square roots therefore have the correct
direction. Multiplying by the nonnegative convex energy gives a valid lower
bound for W over the entire contact domain, not only at its unit-contact
point.

The review derives the 28-record Pareto frontier independently from all 1,395
positive-loss domains and verifies every saved domination map. A separate
forward reachability calculation, which permits all record sequences rather
than using the producer's ordered-multiset recursion, proves:

| Ordered remaining budgets | Maximum relaxed loss |
|---|---:|
| (10/9, 4/3) | 20/9 |
| (1, 3/2) | 22/9 |
| (3/2, 1) | 22/9 |

The recursion terminates because each positive-loss record costs at least one
half in total. It permits at most one modified factor. Negative-loss factors
can be dropped when seeking an upper bound. Replacing avoid-node A/D cases,
or already-avoiding B/C cases, by the larger ordinary target enlarges the
domain and does not violate the one-modified-factor restriction.

## All-rank tails

The finite census is supplemented by an analytic, not sampled, tail proof.
The literal prefix coordinates give the following complete inverse-diagonal
reductions. Interior distances below mean integers at least two.

* For ordinary B, omission 0 leaves node1 of cost1; the other node classes
  have costs k≥2 or r/2. Omission r leaves r−1 of cost1; the other classes
  have costs r/4 or r−k≥2. For modified B, omission r−1 leaves r−2 and r
  of costs1 and1/2, while its other classes have costs (r−1)/4 or
  r−1−k≥2. The opposite omission0 again leaves only node1.
* For ordinary C the only nodes left at the two ends are1 and r−1, each
  of cost1. The other costs are at least two or r/2. Modified omission1
  leaves0 and2 with costs1/2 and1; the remaining end cost is (r−1)/2.
* For D, omission0 leaves only1, and omission r leaves only r−1, each of
  cost1. Other classes have cost at least two or r/4.

These exclusions hold for all B/C ranks r≥8 and D ranks r≥18. Inserting the
coordinate formulas shows that arbitrary convex means in the surviving full
supports differ in fixed coordinates. Their energies are at least
(r−2)/4, (r−3)/4, (r−2)/2, (r−3)/2 and (r−2)/4, respectively. In both B
cases one full support has all kernel coefficients at least √2, so the
height-square factor is at least 2/B². The resulting W−T expressions are
affine in r−8 or r−18 with nonnegative constant and positive slope. The
checker verifies these symbolic expressions and separately checks their
physical coordinate identities at the boundary ranks.

For A cycles of size m≥18 even or m≥19 odd, the resistance is
\(R(d)=d(m-d)/m\). Its concavity shows that any retained distance at least
two exceeds the cost cap. Two distinct retained entries also exceed the cap,
so each full support is its omission plus at most one adjacent node. Two
singleton full supports have separation at least floor(m/2)−2. If a full
support has two nodes, its height is at least two. Expanding the difference
of the two means leaves a nonnegative square and a linear loss at most
\(2(m+(m\bmod2))/m\). The independent symbolic calculation reproduces the
primary's two positive parity polynomials and proves W≥T at every remaining
rank. Rotations and fork swaps implement the stated avoid-node choices;
avoid-A1 has target zero. Thus no unexamined large rank or dihedral-label
cutoff remains.

## Failing-factor split and dimension threshold

Completeness of the old 15 possible unit-charge defects is inherited from
the frozen [independent 164 review](fresh_nikulin_defect_frontier_20261015_review.md).
The new review independently reconstructs their ordered cost floors. The two
same-row E7 cases are eliminated by four fresh convex certificates; these
give W+3/2>14/9 even when the omitted entry is allowed to vary. The remaining
13 support occurrences have exactly the stated E8 or E7-type floors.

If no factor actually fails unit charge, the total loss is at most the sum
of the two endpoint budgets, namely four. Otherwise choose any actually
failing factor. Bounding its own loss by T and applying the full residual
envelope to all other factors yields either
16/9+20/9=4 or 14/9+22/9=4. This includes additional failing factors in the
residual product; it does not require uniqueness or assume that every
possible defect is actually realized. Allowing a modified residual factor
even when the selected factor is modified only weakens the bound.

Consequently the near-balanced pairing bound is M(N−37)/9. The unchanged
four-thirds estimate, using the same omissions, gives M(N−43)/9 outside
the window. At N=163+t these are at least 14+t/9 and 14+7t/60. They contradict
the inherited strict narrow-facet bound at every t≥0. With the independently
reviewed compact-facet branch at143, the combined candidate is n≤163.

The loss-four local witness and its rank163 padding are independently
rechecked. They show sharpness of this relaxed loss constant, not existence
of a hyperbolic chamber or optimality of the dimension bound.

The narrow-facet theorem, actual terminating-wall incidence, common height
ratio, geometric Schur budgets, and compact-facet branch are inherited with
their frozen review chain. The new review does not repeat or weaken those
hypotheses. Generation and saved replay pass, as do recursive source
manifests, syntax and report links.

```sh
python3 scripts/fresh_nikulin_two_budget_frontier_20261015_review.py
python3 scripts/fresh_nikulin_two_budget_frontier_20261015_review.py --verify
```
