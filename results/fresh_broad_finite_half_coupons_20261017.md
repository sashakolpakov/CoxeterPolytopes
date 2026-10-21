# Broad finite-component half-payments retain the9n budget

For every finite Coxeter vertex of rank n<=108, the following payments have
total at most9n. Retain the previous unit payment when the two marked roots
belong to the same finite component and either their distance is at most9
or that whole component has order at most17. Add a half-payment for a pair
in different whole components if:

* one component has rank at most4, with its marked node arbitrary;
* the other is D12,D13,D14,D15,D16, with its mark the fork center or one of
  the two short-arm leaves; or is A13 or A15, marked at its unique middle.

Each eligible cross-component pair receives one half, once. This includes
all the preceding A1/D16, A2/D14, A1/D15, A2/D13 and A1/A15-center coupons.
It also allows A13-center and the D12/D14 finite corners arising in the
new local analysis. Eligibility always concerns whole actual components.

The [exact certificate](fresh_broad_finite_half_coupons_20261017.json) and
[checker](../scripts/fresh_broad_finite_half_coupons_20261017.py) verify the
following elementary budget. The inherited
[finite distance9 estimate](fresh_selective_distance9_dimension108_20261017.md)
bounds the base cost of every component of rank r by9r. For r<=17 the
base cost is exactly binom(r,2), since every pair is paid.

Let a be the sum of the ranks of all components of rank at most4. Their
base-cost slack is at least15a/2, because

    9r-binom(r,2) >= 15r/2       for1<=r<=4.

Let t count the listed large components D12,...,D16,A13,A15. Each has
rank at least12, base-cost slack at least24, and at most three eligible
marked nodes. The total new cross-component cost is therefore at most
3at/2. All other components keep their original nonnegative slack. Hence

    cost <= 9n + a*(3t/2-15/2)-24t,       a+12t<=n.

For n<=108, t<=9. If t<=5 the coefficient of a is nonpositive. At
t=6,7,8,9 the respective maximum values of a are36,24,12,0, giving
excess bounds-90,-96,-138,-216. Thus the displayed surplus is always
nonpositive. This proves the claim for every finite vertex in that range,
independently of its other components or finite labels.

The checker uses exact rational arithmetic for all ten count cases and
the two endpoint choices of a in every dimension1,...,108. The endpoint
checks corroborate the linear inequality; they are not a cutoff in the
finite-diagram classification. Generation and a separate saved replay both
terminated exit0/PASS. A separate agent read the complete source and proof
and checked the component slack, cross-pair count and all count cases. That
read does not repeat the inherited finite classification or prove any face
payment. This source/data/report are final/frozen.

This is a budget theorem. A separate geometric argument must establish
that a face receives these payments, including that the components remain
whole after all supporting whites are added. No ambient dimension bound
or new polytope is claimed by this package alone.

    python3 scripts/fresh_broad_finite_half_coupons_20261017.py
    python3 scripts/fresh_broad_finite_half_coupons_20261017.py --verify
