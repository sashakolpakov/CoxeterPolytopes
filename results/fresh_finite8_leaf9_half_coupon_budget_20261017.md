# Finite distance8 and marked-leaf distance9 payments cost at most8n

For a finite Coxeter vertex of rank n<=80, give a unit to a marked pair in
the same whole connected component when at least one condition holds:

* the marked distance is at most8;
* the distance is9 and at least one mark is a leaf of the whole component;
* the whole component has rank at most17.

Add the [previous broad finite half-payments](fresh_broad_finite_half_coupons_20261017.md):
one half for a pair in different whole components, one of rank at most4
with arbitrary mark, and the other D12,...,D16 marked at its fork center or
short leaf, or A13/A15 marked at its middle. Each eligible pair is counted
once. **The total cost is at most8n.**

The [source](../scripts/fresh_finite8_leaf9_half_coupon_budget_20261017.py) and
[certificate](fresh_finite8_leaf9_half_coupon_budget_20261017.json) check the
identities and inequalities below. This is a vertex-budget theorem only.
Whether these payments cover every required face remains a separate
geometric obligation, including the actual-leaf and whole-component premises.
No new ambient dimension bound follows from this budget alone.

## All-rank finite-component cost

For rank r<=17 the base cost is binom(r,2)<=8r. For r>=18 the finite Coxeter
classification leaves only A_r, B_r and D_r. The A/B graph is a path; its
number of pairs at distance at most8 is

    sum_{k=1}^8(r-k)=8r-36.

The two endpoints each contribute one further distance9 pair. These are
distinct at r>=18, giving total8r-34.

Construct D_r from an A_(r-1) path by duplicating one terminal leaf at its
neighbor. The added leaf has eight spine vertices and its mate within
distance8, contributing nine pairs. Hence the short-pair count is

    8(r-1)-36+9=8r-35.

Each of its three leaves has exactly one vertex at distance9. No two leaves
are at distance9 in this rank range, so the three additional pairs are
distinct. The resulting base cost is8r-32. These formulas hold for every
r>=18; checking ranks18,...,80 merely corroborates the identities in the
range used below. No graph sample is used to extrapolate the classification.

## Half-payments fit the remaining slack

Let a be the sum of ranks of components of rank at most4. Their total base
slack below8 times rank is at least13a/2, since

    8r-binom(r,2)>=13r/2 for1<=r<=4.

Let t count the listed large coupon components. They have rank at least12,
base slack at least8 (the minimum is D16) and at most three eligible marked
nodes each. The new cross-component cost is at most3at/2. Every other
component has nonnegative slack by the preceding all-rank formulas. Thus

    total cost<=8n+a*(3t/2-13/2)-8t,
    a+12t<=n<=80.

Here t<=6. For t<=4 the coefficient of a is nonpositive. For t=5 and6,
the largest possible a is20 and8, respectively, giving maximum excesses
-20 and-28. This proves the8n budget throughout the stated dimension range.
The checker also verifies both possible endpoint values of a for every
admissible t and n=1,...,80.

## Verification scope

Generation and a separate fresh saved replay both terminated exit0/PASS,
including126 graph-formula checks and622 budget endpoints. A separate agent
read the complete source and mathematical argument, checking the all-rank
path/fork counts, component slacks, cross-payment count and seven cases for
t; that review does not repeat the622-case computation. This package is
final/frozen. The finite-type classification
is an inherited mathematical input. No triangle or quadrilateral coverage,
new construction, external refereeing or formalization is asserted.

    python3 scripts/fresh_finite8_leaf9_half_coupon_budget_20261017.py
    python3 scripts/fresh_finite8_leaf9_half_coupon_budget_20261017.py --verify
