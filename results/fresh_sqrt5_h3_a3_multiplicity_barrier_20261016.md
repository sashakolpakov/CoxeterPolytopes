# The analogous A3 endpoint multiplicity barrier

The [source](../scripts/fresh_sqrt5_h3_a3_multiplicity_barrier_20261016.py) and
[exact data](fresh_sqrt5_h3_a3_multiplicity_barrier_20261016.json) prove that the
A4-only count-split family still satisfies

\[
T\ge6042/169\approx35.75147928994.
\]

Let β be the unchanged H3×A3-end coefficient with actual A1 count 3, no I2(5),
and H3 present. Quad 7668 with common external counts (1,0,4,0,0,0,0) and white
facet root 4 has payment 1−δ+β/15, so β≥15δ. The marked-A1 forest
3A1+H3+5A3 has baseline 18 and ten A3 endpoints, giving
T≥18+10β≥18+150δ. The marked-A1 forest A1+D20 gives T≥38−19δ.

Taking weights 19/169 and 150/169 on the latter two inequalities cancels δ and
gives 6042/169. The checker reconstructs all literal profiles and records the
exact four-inequality Farkas combination, including the distance-two-zero
restriction. Remaining bonus coefficients are nonnegative. Generation and
fresh saved replay pass.

This is a restricted payment-model obstruction. It shows that splitting the
H3 count for A4 endpoints alone leaves the same defect at A3 endpoints. It
does not prove geometric nonexistence or complete affected-face coverage.
