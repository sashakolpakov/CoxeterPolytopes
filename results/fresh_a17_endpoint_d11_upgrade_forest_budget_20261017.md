# Exact finite budgets after the A17 endpoint deletion

The finite-corner rule with the whole-A17 endpoint unit removed has exact
maximum544 in rank69. The same maximum holds after upgrading the two
whole-A1/whole-D11-short-leaf halves to full units. Thus either explicit
rule fits the rank69 finite curvature coefficient544 with zero slack.
The D11 upgrade's geometric face coverage remains a separate obligation;
this package proves only its finite-forest cost. No ambient dimension
bound follows from this budget alone.

The [source](../scripts/fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.py)
and [certificate](fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.json)
compute both variants separately in every rank1 through75. They retain
the [complete inherited coupons](fresh_finite8_d12_b10_coupon_budget_20261017.md)
and use exactly the [single proposed A17 endpoint deletion](fresh_finite_a17_endpoint_payment_support_20261017.md).
The additional D11 rule is precise: for each marked pair in distinct whole
components A1 and D11, upgrade weight1/2 to1 when the D11 mark is one
of its two short fork leaves. Its fork-center weight remains1/2. All pairs
are counted once.

## Why rank17 requires two representatives

The [previous exact budget](fresh_exact_finite8_forest_budget_20261017.md)
proved componentwise type domination for the unchanged rule. Its rank17
representative A17 previously had the same base136 as B17/D17 and also
two eligible cross-coupon positions. After removing one A17 pair, that
argument no longer applies at rank17:

|Type|Same-component cost|Broad cross-coupon positions|
|---|---:|---:|
|A17|135|2|
|D17|136|0|

Keep **both** types. B17 is dominated by D17 with equal base cost and
no coupons. All other rank representatives remain as before: A_r for
rank1 through9, B10 at10, D_r at11 through16 and at least18. The D11
upgrade only strengthens D11's domination of other rank11 types. The
small-type representatives A1/A2/A4 retain every special coupon.
Noncrystallographic and exceptional finite types remain covered by the
published finite-type classification and inherited same-rank domination.

For rank at most17 the unchanged base is binomial(r,2), except A17 now
has135. For ranks at least18, the dominating D base is8r-32. There is
no assertion that a maximizing abstract forest occurs at an actual
finite vertex. Every actual finite vertex belongs to the counted class,
which suffices for the upper bound.

## Exact cost and complete maximization

Let a_i count A_i components for1<=i<=5, put `a=sum(i*a_i)`, let t count
D11 through D16, and let u,b10,b11,b12,b13 count A17,B10,D11,D12,D13.
The exact cost of a representative forest in the upgraded variant is

    sum_components base(type)
      + a*(3t/2+u)
      + a1*b12 + a2*b10 + a4*b13
      + a1*b11.

Omit only the last term for the first variant. The broad three-mark D
half-coupons and two-mark A17 half-coupons give the second term. The
inherited D12 leaf upgrade, A2/B10 halves and A4/D13 neighbor halves give
the next three terms. The new D11 addition has two leaves, each increased
by1/2 against an A1, giving exactly `a1*b11`. This is an increment to the
existing halves, so the same marked pair is never counted twice.

Enumerate every count vector of the eight coupon-bearing types
`B10,D11,D12,D13,D14,D15,D16,A17` with total rank at most75. There are
1,608 vectors, as in the inherited complete enumeration. **D17 is not
included in its A17 coordinate.** The residual additive components have
ranks1 through9 or17 through75, represented by A_r in the first range
and D_r in the second.

For each fixed large vector, double the costs. A residual rank r has
base `r*(r-1)` at r<=17 and `16*r-64` at larger r. At r<=5 add
`r*(3t+2u)`, and additionally `2*(b12+b11)` for r=1, `2*b10` for r=2,
or `2*b13` for r=4. In the first variant set the upgrade coordinate b11
to zero while retaining D11 in t.

Only `(t,u,b10,b12,b13,upgraded_b11)` matters. There are377 tuples in
the first variant and675 in the upgraded variant. For every tuple k,
compute all76 Bellman values

    F_k(0)=0,
    F_k(s)=max_(allowed r<=s) (F_k(s-r)+weight_k(r)).

This is a complete induction over residual rank: delete any component
to obtain a predecessor, and conversely attach the selected component
to an attaining predecessor. Then maximize the large base plus the
appropriate residual value over all feasible large vectors. Each saved
maximum includes an attaining **typed** forest, checked by a separate
direct evaluation of the full displayed cost formula. The type labels
ensure an A17 coupon cannot be accidentally assigned to D17.

## Results near the frontier

Both variants have identical maxima in every rank through75. In
particular:

|Rank|Maximum|One maximizing forest|
|---:|---:|---|
|66|513|D15+3D17|
|67|528|D16+3D17|
|68|544|4D17|
|69|544|A1+4D17|
|70|549|A2+4A17|
|71|555|A3+4A17|
|72|562|A4+4A17|
|73|570|A5+4A17|
|74|574|A1+A5+4A17|
|75|579|A2+A5+4A17|

At rank69, A1+4A17 also attains544: its A17 base is540 and its retained
cross-coupons total4. A1+4D17 attains544 without any cross-coupons.
Keeping only the reduced A17 base at rank17 would have incorrectly
underestimated the rank68 maximum, among other configurations.

At69 the finite curvature coefficient is544, so this local refinement
removes the earlier finite deficit4. The existing ideal rank-profile
inequalities remain an independent obligation. In a future global
composition, zero finite slack can be sufficient together with positive
ideal margins and the inherited strict positive global constant; this
package does not establish those ideal margins.

## Verification scope

Final generation and a separate fresh saved replay both terminated
exit0/PASS after binding the frozen A17 local package. All22 direct
bindings, syntax and5 local report links pass. The two runs
exhaust the complete large-vector and Bellman domains,
verify every attaining forest, compare the deletion-only result against
the frozen old maximum, and compare the two new variants in all ranks.
The finite-type domination proof is inherited except for the explicit
two-representative rank17 correction and the nonnegative D11 increment.
There is no geometric realizability, new face coverage, external refereeing,
formalization or new ambient dimension claim in this budget calculation.
Root completed a full primary source and mathematical report read. It
passed the A17/D17 distinction, unchanged domination at all other ranks,
the precise `a1*b11` increment, all1,608-vector enumeration, complete
residual rank domain and Bellman recurrence, typed attaining forests,
and zero rank69 finite slack. The review did not infer any ideal or
ambient dimension conclusion from the finite budget. Source, certificate
and report are final/frozen.

    python3 scripts/fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.py
    python3 scripts/fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.py --verify
