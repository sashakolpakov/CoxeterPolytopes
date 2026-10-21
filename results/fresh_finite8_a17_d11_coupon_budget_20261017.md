# Finite-eight cost with whole A17 and D11 half-coupons

**FINAL/FROZEN: generation, fresh saved replay and separate full mathematical/source review PASS.**

For a finite Coxeter vertex of dimension n<=80, the following rule has total
cost at most **8n**. Retain the unit payments of the
[finite-eight budget](fresh_finite8_leaf9_half_coupon_budget_20261017.md): two
marks in the same whole component pay one if their distance is at most8,
their distance is9 with an actual whole-component leaf mark, or that component
has rank at most17. For marks in distinct whole components, add one half when
one component has rank at most5 and the other is one of:

* D11,...,D16, marked at the fork center or either short leaf;
* A13 or A15, marked at its unique middle;
* A17, marked at path position5 or11, with positions numbered0,...,16.

Count each unordered marked pair once. All earlier half-coupons are retained.
The A17 rule covers both sealed cycle18 coupons from the
[ideal-five local theorem](fresh_finite8_ideal5_local_replacement_20261017.md).
That theorem needs only an A1 on the small side. Allowing all ranks at most5,
and also allowing D11, is a proved counting enlargement; this statement does
not claim a new geometric equality configuration uses either enlargement.

## Counting proof

The predecessor proves the all-rank base costs: at rank r>=18, types A and B
cost8r-34 and type D costs8r-32. Every component of rank r<=17 costs at most
binom(r,2)<=8r. These bounds include all finite labels through the finite-type
classification; their validity does not come from extrapolating a graph sample.

Let a be the total rank of small components of rank at most5. Their base
slack is at least6a. Let t count whole components D11,...,D16,A13,A15.
Each uses at least11 ranks, has at most three eligible marks, and has base
slack at least8. Let u count whole A17 components. Each uses17 ranks and has
two eligible marks. Its base slack is zero, so it cannot be included in t.
All remaining components have nonnegative base slack.

There are at most a*(3t+2u) eligible cross-component pairs, at half a unit each.
Consequently the excess above8n is bounded by

    E = a*(3t/2+u-6)-8t,
    a+11t+17u<=n<=80.

If the coefficient of a is nonpositive, E<=0. Otherwise maximize a at
80-11t-17u. The complete possible (t,u) pairs and maxima are:

| t | u | Maximum a | Coefficient of a | Maximum E |
|---:|---:|---:|---:|---:|
|3|2|13|1/2|-35/2|
|4|1|19|1|-13|
|4|2|2|2|-28|
|5|0|25|3/2|-5/2|
|5|1|8|5/2|-20|
|6|0|14|3|-6|
|7|0|3|9/2|-85/2|

Completeness is elementary: t<=7; for t<=2, the allowed u makes the
coefficient nonpositive, while the rank inequality bounds u by2 for t=3,4,
by1 for t=5 and by0 for t=6,7. When t=3 only u=2 has positive coefficient.
All displayed excesses are negative. This proves the
counting theorem in every dimension through80 without a rank-partition
enumeration.

The [source](../scripts/fresh_finite8_a17_d11_coupon_budget_20261017.py) and
[certificate](fresh_finite8_a17_d11_coupon_budget_20261017.json) additionally
check both endpoints of the affine function E for every allowable integer
n,t,u. The whole-component and actual-leaf hypotheses are essential. The
geometric assertion that every triangle and quadrilateral receives enough
payment is a separate obligation; **no ambient dimension improvement follows
from this counting theorem alone**.

## Verification scope

Generation and a separate fresh saved replay both terminated exit0/PASS for
all1,487 endpoint checks. A second agent read the complete proof and source,
including every component slack, the cross-component pair count and all seven
positive-coefficient cases, and found no gap. That review did not duplicate
the1,487-case arithmetic run. This package reuses the predecessor's all-rank
finite base proof and verifies the new cross-component allowance. No external
refereeing or proof-assistant formalization is claimed.
