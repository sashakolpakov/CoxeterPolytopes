# Narrow independent review of the compact-dimension-29 reduction

**PASS for the necessary-condition reduction.** The [independent checker](../scripts/fresh_compact_dimension_bound_20261015_review.py) and [certificate](fresh_compact_dimension_bound_20261015_review.json) reproduce the h-vector coefficients and nine rank partitions, with 2,126 typed forests and 1,354 after removing (H_4). Generation and saved replay pass. No compact dimension is excluded.

This review inherits the connected-component weight table and global planar-face inequality from the [primary package](fresh_compact_dimension_bound_20261015.md), which explicitly cites Vinberg's original proof. It does not independently repeat the 66 classical graph counts, exceptional (E_8) picture, or geometric face proof. The arithmetic implementation imports no producer code and enumerates bounded component multiplicities rather than recursively enumerating partitions.

For a simple 29-polytope, each vertex belongs to (\binom{29}{2}=406) two-faces. The strict odd-dimensional Nikulin bound gives (f_2>98f_0), hence total weight greater than (84f_0). Pairing the symmetric dual h-vector terms independently gives

\[
f_2-98f_0=\sum_{i=0}^{14}(14-i)(15-i)h_i.
\]

Double-counting face-angle weights therefore proves the primary's exact balance identity. Its lower bound 1050 uses only (h_0=1) and nonnegativity. Since each vertex has weight at most 87, at least 350 vertices must have weight greater than 84; negatively weighted excesses at other vertices can only increase that requirement.

From the stated component table, the deficit (3r-w(r)) is nonnegative, is zero exactly at ranks 4,7,8,9, and is one at rank 10. Every other rank has deficit at least three. An extremal rank-29 forest consequently has only those five ranks and at most two rank-10 components. The independent bounded-multiplicity count reconstructs all nine possibilities and their unordered typed multiplicities. In particular the unbounded (I_2(m)) family cannot enter: its rank-two deficit is four, independent of (m). There is no label cutoff.

The separate arithmetic field reduction is valid: an (H_4) component has the cyclic product (4\cos^2(\pi/5)=(3+\sqrt5)/2), forcing (\sqrt5) into the cyclic-product field. Thus it is absent in either specified field not containing (\sqrt5). This does not exclude label 4 or B/F components.

The proposed next step needs an occurrence argument, which the primary author supplied during review and which passes independently. Set every weight at a face with at least five sides to zero, keeping the original triangle/quadrilateral weights. The face lower bound remains valid because (0\ge5-k) there. Modified vertex totals only decrease, so they remain at most 87, and the same balance still forces at least 350 vertices with modified subtotal greater than 84. Each such vertex has original weight greater than 84, belongs to the listed types, and is incident to a triangle or quadrilateral. Consequently a complete rejection of all such small-face stars through all listed types would exclude compact29. No such rejection has been performed; the current result is the necessary-condition reduction.

```sh
python scripts/fresh_compact_dimension_bound_20261015_review.py --verify
```
