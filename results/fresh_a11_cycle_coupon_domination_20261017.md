# The cycle12 A11 coupons fit the revised exact finite budget

Adding the two [cycle12 coupons](fresh_cycle12_ideal5_transfer_20261017.md)
does not increase the [revised finite-forest maxima](fresh_a17_endpoint_d11_upgrade_forest_budget_20261017.md).
In particular the exact maximum at ambient rank69 remains544, equal to
the finite curvature coefficient. The rule includes the whole-A17
endpoint deletion, the full A1/D11-short-leaf upgrade and every inherited
coupon. This is a finite-budget theorem, not an ambient dimension bound.

The new pairs join whole A1 to the unique middle node of whole A11,
at unit weight, or either node of whole A2 to that A11 middle, at weight
one half. No old coupon uses A11 as its large component. Their total
interaction costs with one A11 are therefore1 against A1,1 against A2,
and0 against every other component.

Replace each A11 by D11. Both same-component costs are55. Under the
existing upgraded D11 rule, interactions cost5/2 against A1,3 against
A2, and3r/2 against any other component of rank r<=5. Interactions
with larger components vanish. Each total is at least the corresponding
new A11 total. Interactions between two large components vanish before
and after replacement, so the replacement works successively throughout
an arbitrary forest without double counting. All pair weights replace
or supplement the stated distinct roles; an inherited half is never
counted twice when upgraded to one.

The [source](../scripts/fresh_a11_cycle_coupon_domination_20261017.py)
and [certificate](fresh_a11_cycle_coupon_domination_20261017.json) check
the interaction comparison against all228 representative finite types
through rank75, including exceptional types and a symbolic I2 family.
The arbitrary dihedral label is irrelevant to these cross-coupon roles.
The rank bound75 is the domain of the inherited exact finite theorem.

The resulting forest is one counted in that exact theorem and has at
least the enlarged-rule cost. Its maxima are therefore valid upper bounds.
Conversely every saved maximizing forest contains no A11, so its cost
is unchanged and still attains the same maximum. The checker directly
recounts all76 typed witnesses, including rank0, using the separate
graph-pair implementation from the finite-budget reviewer.

At rank69 the forests A1+4D17 and A1+4A17 still attain544. At rank68,
4D17 still attains544, exceeding the curvature coefficient527. The new
coupons do not remove that lower-rank finite obstruction.

Generation and a fresh separate-process saved replay both terminated
exit0/PASS on228 partner comparisons and76 attaining forests. A separate
full source and mathematical read, including the imported direct graph
pair-count routines, passed the complete role comparison, successive
rank-preserving replacements and unchanged attainment argument. That
read did not repeat the228/76 arithmetic checks or inherited Bellman
calculations. Source, certificate and report are final/frozen. The
finite-type classification and inherited exhaustive Bellman certificates
retain their stated verification scopes.
