# Two exact A12-source transplants

The [I12,1 source](fresh_h12_a12_source_20261014.md) supplies two missing facets in the [E8-face frame](fresh_h20_contact_closure_20261014_e8_face.md):

\[
R=-\lambda_3+2T,\qquad S=-\lambda_{11}+3T.
\]

They have norms two and one and mutual pairing zero. Chain reversal replaces the weights by λ10 and λ2.

The [checker](../scripts/fresh_h20_contact_closure_20261014_h12_transfer.py) and [certificate](fresh_h20_contact_closure_20261014_h12_transfer.json) test both orientations against all 76 walls. Both caps are acute to every wall, but the norm-two cap has forbidden finite angles:

| Orientation | Norm-two cap: forbidden walls | Norm-one cap |
|---|---|---|
| Forward | 0, 18, 20, 31 | 19 right, one π/4, 56 dotted |
| Reversed | 0 | 19 right, one π/4, 56 dotted |

All offending walls belong to the original 60-wall core. Their squared forbidden cosines are nonrational in Q(√79). Such a value cannot be cos²(π/m): degree two restricts m to 5,8,10,12, with fields Q(√5), Q(√2), or Q(√3). The rational possibilities are checked directly.

For example, the reversed norm-two cap pairs with wall0 by

\[
-\frac{202+46\sqrt{79}}{351}.
\]

Its normalized squared cosine lies strictly between 3/4 and 4/5. The forward wall0 value is already strictly between zero and 1/4.

Neither individually admissible norm-one cap cuts the saved 76-wall escape. Their pairings, forward then reversed, are

\[
\frac{7787240030-297185372131\sqrt{79}}{234000000000}<0,\qquad
\frac{947220461579-297185372131\sqrt{79}}{234000000000}<0.
\]

This excludes only the two fixed transplants retaining all 76 walls, not other H12 sources or changed cores. In fact, [deleting wall0 makes the reversed pair coherent](fresh_h20_contact_closure_20261014_h12_replacement.md).

    python scripts/fresh_h20_contact_closure_20261014_h12_transfer.py --verify
