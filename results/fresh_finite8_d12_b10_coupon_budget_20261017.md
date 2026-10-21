# Finite cost eight with two order-13 transfer coupons

At every finite vertex in ambient dimension at most75, the following enlarged
corner rule has total cost at most8n. Retain the entire
[previous rule](fresh_finite8_a4_d13_neighbor_coupon_budget_20261017.md), including
the A4/D13-neighbor half-payments, and make exactly these changes:

* Upgrade the existing half-payment to one for a pair in distinct whole A1
  and D12 components, when the D12 mark is a short fork leaf. Its fork center
  is not upgraded.
* Add a half-payment for a pair in distinct whole A2 and B10 components,
  when the B10 mark is its terminal node incident to the doubled bond.
  Either node of A2 is eligible.

Each unordered marked pair is counted once. These two additions have disjoint
component types, and neither duplicates the A4/D13-neighbor rule. The
[source](../scripts/fresh_finite8_d12_b10_coupon_budget_20261017.py) and
[certificate](fresh_finite8_d12_b10_coupon_budget_20261017.json) verify the exact
finite integer inequality below. Local geometric transfers are separate
obligations; this counting statement alone removes no ideal payment.

Let a be the total rank of components of rank at most5. Let t count the
inherited large types D11 through D16, A13 and A15; let b and c count the
D12 and D13 components among them. Let u count A17 and v count B10.
Then

    b+c<=t,  a+11t+b+2c+17u+10v<=n<=75.

The inherited same-component cost is at most8 times component rank.
Small components supply at least6a of slack. Each of the t large components
supplies at least8; D12 supplies30, hence22 extra, and D13 supplies26,
hence18 extra. Whole B10 has all45 unordered pairs paid and supplies35
slack. A17 supplies no slack. All other components have nonnegative slack.

The old small/large half-payments cost at most a(3t/2+u). The retained
A4/D13-neighbor addition costs at most ac/4. For the D12 upgrade, each A1
has two eligible D12 pairs whose weights each rise by1/2. Its added cost
is thus one per A1/D12 component pair, at most ab. Each A2/B10 component
pair has two eligible pairs of weight1/2, giving cost one; the number of
A2 components is at most a/2, so this costs at most av/2.

Consequently the surplus above8n is at most

    E=a*(3t/2+u-6+b+c/4+v/2)-8t-22b-18c-35v.

For every fixed nonnegative integer t,b,c,u,v satisfying the rank constraint,
E is affine in a. Its maximum on the permitted interval occurs at a=0 or
a=n-11t-b-2c-17u-10v. The checker exhausts all those integer tuples for
every n from1 through75 and both distinct endpoints. All13,034 values are
nonpositive exact rational numbers. This proves E<=0 throughout the whole
permitted domain. It is a finite proof with an explicit rank bound, not a
sampling of finite Coxeter forests.

The restriction to75 is intentional. A previously completed dimension
bound may handle higher ambient dimensions in a subsequent global argument;
this report makes no assertion that these additional coupons fit through80.

## Verification status

Generation and a fresh separate-process saved replay both terminated
exit0/PASS on all13,034 exact endpoints. A separate full internal source and
mathematical review passed the component slacks, disjoint coupon increments,
rank relaxation, affine endpoint reduction and complete integer loop ranges.
That review did not duplicate the13,034-row calculation. No local geometric
transfer was inferred from the budget check. Source, data and report are
final/frozen. No new ambient bound, external refereeing or proof-assistant
formalization is claimed here.

    python3 scripts/fresh_finite8_d12_b10_coupon_budget_20261017.py
    python3 scripts/fresh_finite8_d12_b10_coupon_budget_20261017.py --verify
