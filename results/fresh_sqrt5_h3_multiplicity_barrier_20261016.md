# A small exact H3 multiplicity barrier

The [source](../scripts/fresh_sqrt5_h3_multiplicity_barrier_20261016.py) and
[data](fresh_sqrt5_h3_multiplicity_barrier_20261016.json) prove that the current
presence-only H3 cross-payment family satisfies

\[
T\ge4997/139\approx35.94964028777.
\]

Let β be its H3×A4-end coefficient with actual A1 count 2, no I2(5), and H3
present. Quad 7661 with common external counts (1,0,4,0,0,0,0) and white facet
root 4 has payment 1−δ+β/15, so β≥15δ. The denominator reflects five actual H3
copies at that corner. The marked-A1 forest 2A1+H3+4A4 has baseline 23 and
eight A4 endpoints, giving T≥23+8β≥23+120δ. The marked-A1 forest A1+D20 gives
T≥38−19δ when the distance-two discount is zero.

Combining the latter two inequalities with coefficients 19/139 and 120/139
cancels δ and yields 4997/139. The checker reconstructs every sparse profile
from literal root pairs, records the exact four-inequality Farkas combination
(including the distance-two-zero restriction), and verifies that all residual
bonus coefficients are nonnegative.

Generation and fresh saved replay pass. This is a payment-model obstruction,
not a geometric nonexistence theorem. The same β is used at H3 counts 5 and 1;
giving these counts separate coefficients removes this specific algebraic
coupling. It does not by itself certify successor feasibility or complete
affected-face coverage.
