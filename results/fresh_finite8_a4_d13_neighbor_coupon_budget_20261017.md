# Finite-eight cost with the whole A4/D13-neighbor half coupon

For every finite Coxeter vertex of dimension `n<=80`, the complete frozen
[finite8/A17/D11 payment rule](fresh_finite8_a17_d11_coupon_budget_20261017.md)
remains within **8n** after adding the following half payment:

> For marks in distinct whole components, pay one half when one mark is an
> endpoint of a whole A4 and the other is the unique non-leaf neighbor of
> the fork center of a whole D13.

Count every unordered marked pair once. The new D13 root is the first root
on the long arm adjacent to the fork center. It is neither the center nor a
short leaf, so these pairs do not duplicate the inherited D13 coupon. The
whole-component hypotheses are essential. Applying this counting allowance
to a particular quadrilateral requires a separate local transfer proof.
This package changes no ideal payment and proves no new ambient bound.

The [source](../scripts/fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.py)
and [certificate](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.json)
are separate successors to the frozen
[spectrum/C3 checkpoint](spectator_spectrum_and_c3_frontier_20261017.md).
That checkpoint is not modified or extended by this package.

## Counting proof

Retain the notation of the frozen budget. Let `a` be the total rank of whole
components of rank at most5. Let `t` count whole components D11 through D16,
A13 and A15; let `u` count whole A17 components. Let `b` count whole D13
components, so `0<=b<=t`.

All same-component unit payments, actual leaf-nine payments and inherited
half-coupons remain unchanged. The predecessor proves that small components
have base slack at least `6a`, every component counted by `t` has base slack
at least8 and at most three old eligible marks, and each A17 has two eligible
marks and zero base slack. All other components have nonnegative base slack.
These estimates use its all-rank finite-type base proof, not an enumeration
of forests up to dimension80.

A whole D13 has all `binom(13,2)=78` same-component pairs paid and hence slack
`8*13-78=26`. Relative to the old lower bound8, each D13 supplies18 additional
units of slack and uses two extra ranks relative to the old minimum11.
Thus

```
a+11t+2b+17u <= n <= 80.
```

Each whole A4 has exactly two endpoints and the D13 has one new eligible
root. The two new pairs each cost one half, giving one unit per A4/D13
component pair. The number of whole A4 components is at most `a/4`.
Consequently the added payment is at most `a*b/4`, while the refined slack
is at least `6a+8t+18b`. The excess above8n is therefore at most

```
E_new = a*(3t/2+u-6+b/4)-8t-18b
      = E_old + b*(a/4-18),
E_old = a*(3t/2+u-6)-8t.
```

The old rank inequality `a+11t+17u<=80` still holds, so the frozen theorem
gives `E_old<=0`. If `b=0`, there is no new term. If `b>=1`, then `t>=1` and
the refined rank inequality gives `a<=80-11-2=67`. Hence

```
b*(a/4-18) <= -5b/4 < 0.
```

This proves the enlarged rule in every dimension through80. No new table
of partitions or geometric existence assumption is needed. In particular,
the estimate deliberately ignores the additional slack available from a
whole A4, so it is a conservative allowance.

## Exact certificate and scope

The independent-of-optimizer Fraction checker verifies the component
costs and the analytic dominance identity. It also checks both endpoints
of the affine function in `a` for every allowed integer `n,t,b,u`, for every
`1<=n<=80`. There are3,465 distinct endpoints in total and138 at dimension80.
Counting both endpoints separately when they coincide at `a=0` gives140
at dimension80; this explains the difference from the exploratory count.
All38 positive-coefficient triples at dimension80 have strictly negative
maximum excess. The1,487 inherited budget endpoint records are checked as
well, and their frozen source/data/report hashes are bound.

**FINAL/FROZEN.** Primary generation and a separate fresh saved replay both
terminated exit0/PASS. Root read the complete source and mathematical proof,
including the distinct new root role, component-pair count, refined rank/slack
bounds and analytic dominance; no gap was identified. That review did not
duplicate the3,465 endpoint arithmetic run. All4 direct source bindings,4
local report links and source syntax pass. These checks establish the finite
budget only. They do not assert that a proposed D17 pair can be removed,
that its replacement corners contain the required whole A4 and D13, or that
a new global curvature comparison succeeds.

```
python3 scripts/fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.py
python3 scripts/fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.py --verify
```
