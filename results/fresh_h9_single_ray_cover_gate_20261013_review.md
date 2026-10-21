# Independent review of the four separate single-ray obstructions

The theorem is sound with its stated scope: **fix one H4 fundamental-weight direction in the eleven-wall seed**. Any finite Coxeter-compatible collection of caps in that direction leaves infinite volume. The proof does not apply to mixtures of directions, mixed transverse profiles, deletion of seed walls, or arbitrary additions to the twelve-wall chamber.

The reviewed [proof](fresh_h9_single_ray_cover_gate_20261013.md), [producer](../scripts/fresh_h9_single_ray_cover_gate_20261013.py), and [certificate](fresh_h9_single_ray_cover_gate_20261013.json) have a separate [sparse-radical checker](../scripts/fresh_h9_single_ray_cover_gate_20261013_review.py) and [audit artifact](fresh_h9_single_ray_cover_gate_20261013_review.json). The checker imports no producer arithmetic. It re-enumerates all 127 support relaxations, obtaining 34 empty cases and 122 vertex records, and independently verifies every saved support minimum and covering margin.

The potentially delicate continuous step is valid. Acuteness against the complete compact old prism puts every old cap component in its future timelike cone. For a fixed zero set, the relaxed polyhedron has a pointed recession cone contained in a face of that cone. A nonempty such polyhedron has vertices and equals their convex hull plus its recession cone. The positive Lorentz norm \(\sqrt{-u^2}\) is concave, and the reverse triangle inequality shows that adding a future recession vector cannot lower it. Its minimum is therefore attained at a vertex. The finite support census consequently bounds the entire continuous polyhedron, not just the enumerated points.

All seven six-constraint systems are nonsingular. Thus selecting six active constraints, including every forced equality, misses no vertex of any relaxation. The polar lower bound \(2/5\) is strictly weaker than the actual first nonzero Coxeter pairing. Every support with at most four zero pairings has square at least
\[
(-38+64\sqrt5)/25>1+\sqrt5.
\]
For the relevant (j=3,m=3) norm, only four of the ten old vertex zero sets survive this relaxation. Checking the gap against all ten old vertices is therefore a valid stronger check.

The ball comparison is also correct. With (A=B_{jj}-2/t^2), (A'=B_{jj}-2/s^2), and (P=B_{ji}^2/B_{ii}), acuteness bounds the center distance by (B_{jj}/\sqrt{AA'}). The sum-of-radii expression is
\[
\frac{P+\sqrt{(2/t^2-D)(2/s^2-D)}}{\sqrt{AA'}},
\qquad D=B_{jj}-P.
\]
The strict separation condition simplifies exactly to (D(t^2+s^2)>2). All denominators are positive, and the square-root factors are nonnegative whenever the coverage balls exist. The saved diameter inequalities exclude a single isolated ball covering the old prism.

For (j=2), all relevant balls are separated. For (j=3), every label-4 or label-5 ball is separated from every other ball, while the label-3 balls have the certified common gap. A finite union of these separated closed components cannot cover the connected old prism unless one component covers it; the diameter and vertex-ball checks rule out those possibilities. This argument does not require an unjustified bound on the number of centers or on their fields.

Finally, an uncovered point can be chosen in the interior of the old prism because the complement of finitely many closed balls is relatively open. All added cap inequalities are strict at its corresponding ideal vector. Perturbing the transverse unit vector into the interior of the H4 sector preserves those strict inequalities and produces an open patch of ideal boundary. Infinite volume follows. Finiteness of the cap collection is used here and in the separation argument.

Independent checker generation and saved replay **PASS**. The primary saved replay also passes. No mathematical issue was found in the reduction, arithmetic, or finite-volume inference.

```sh
python scripts/fresh_h9_single_ray_cover_gate_20261013_review.py --verify
python scripts/fresh_h9_single_ray_cover_gate_20261013.py --verify
```
