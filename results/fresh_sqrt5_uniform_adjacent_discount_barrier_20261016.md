# A small exact obstruction to the uniform adjacent discount

The [source](../scripts/fresh_sqrt5_uniform_adjacent_discount_barrier_20261016.py)
and [data](fresh_sqrt5_uniform_adjacent_discount_barrier_20261016.json) isolate
an exact **T≥38** obstruction in the present payment family, assuming the
distance-two discount is zero.

Sealed quad 4266 has payment exactly 1−δ at zero reserved spectator counts,
for every white facet mark. Sealed quad 7746 has the same property and uses
only bonds 3 and 5. Neither profile has any available bonus variable in this
family. The face inequality therefore forces δ≤0; nonnegativity forces δ=0.

Both quads have the same discounted E6 adjacent-pair orbit: an edge from the
trivalent root to a degree-two neighbor on either length-two arm. There are
two such edges in E6. Their literal corner graphs and colored pair codes are
saved. Quad 4266 has an E6 paid corner; quad 7746 has E6+I2(5).

The marked-A1 vertex A1+D20 has literal payment

\[
38-19\delta_{\rm distance2}-19\delta_{\rm adjacent}
  +\text{nonnegative bonus terms}.
\]

Adding its budget inequality to 19 copies of the zero-distance-two restriction
and 19 copies of either quad inequality yields T≥38. The script saves the exact
combined columns and checks that all remaining bonus coefficients are
nonnegative. This is a three-inequality Farkas argument, independent of the
large numerical optimizer.

The bounded scan of the 35 sealed stopping rows also finds bonus-free profiles
2−2δ in classes 5830, 7504, 7506, and 7514, with paid A6, D6, or A7 corners.
Those impose only δ≤1/2 and do not alone force δ=0.

Generation and fresh saved replay pass. This is a payment-family obstruction,
not a geometric nonexistence theorem. Local arithmetic admissibility of the
literal quads is an additional geometric question addressed separately.
Exempting the identified E6 orbit from the discount, or paying another corner
of these faces, breaks this specific algebraic obstruction; it does not
establish feasibility of the resulting model.
