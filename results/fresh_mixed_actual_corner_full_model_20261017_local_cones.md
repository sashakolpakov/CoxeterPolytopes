# Complete local simplex-packing covers at the new mixed-model masses

Reconstructing the entire local simplex-block packing family at the newly optimized face masses lowers the exact witness from **12.50758862113165 to 12.405107765824441**, still above the target12. The 4,661 face masses remain fixed. There is no new global optimization, target12 feasibility claim or dimension theorem.

The [source](../scripts/fresh_mixed_actual_corner_full_model_20261017_local_cones.py) and [certificate](fresh_mixed_actual_corner_full_model_20261017_local_cones.json) bind the frozen [mixed finite model](fresh_mixed_actual_corner_full_model_20261017.md). Generation and a fresh saved replay both finished with exit0 and `PASS_NEW_MIXED_ACTUAL_CORNER_DUAL_COMPLETE_LOCAL_PACKING_COVERS`. Replay invokes no optimizer.

The new positive simplex support consists of2,753 rows:460 four-simplices,500 three-simplices,597 triangles and1,196 maximal edges. The source reconstructs8,483 actual corner occurrences,5,173 complete literal port-block orbits and all760 contributing rank13 vertex forests. It reads the current optimized alpha; no old fixed-alpha cover is reused. Progress records, if reused within this run, are bound to both the exact model digest and the complete current face-mass list.

Every forest receives a verified rational covering mixture of disjoint literal blocks. There are4,363 mixture terms in total. Matching exact local charges certify675 local optima; the remaining85 forests have exact covering brackets. There are no resource-open cases. The implementation retains the prior bounded local enumeration/solve safeguards:at most1,000,000 count-vector entries,30seconds enumeration and20seconds for an individual local solve. These are small per-forest packing problems, not a rerun of the global model.

Each saved mixture is checked directly:every block belongs to its full literal orbit, every term has pairwise-disjoint port masks, every mass is positive, and the resulting marginal for every orbit covers its newly required alpha-mass. Let w_j be those required orbit masses and let a_P be the saved nonnegative mixture masses. The certificate proves

\[
\sum_P a_P n_j(P)\ge w_j,
\qquad M=\sum_Pa_P.
\]

For any nonnegative orbit-local charge whose cost on every disjoint packing is at most1, its load is at mostM. Simultaneously averaging each whole packing over the vertex-forest automorphism group gives uniform literal marginals; it does not choose independent incompatible orientations for separate blocks. Thus these mixtures cover the entire nonnegative local simplex-block charge family at once, not just the finite charges recovered by the local solves. The [disjoint-port theorem](fresh_simplex_disjoint_ports_20261017.md) and the [earlier complete-cone implementation](fresh_actual_corner_enlarged_model_20261017_local_cones.md) provide the inherited combinatorial framework. Additional geometric compatibility among different faces is a separate restriction.

There are132 newly violated exact local charges and205 increased vertex allowances. For example, at3A1+A3+H3+H4 (vertex982), an exact cost-one charge has load52311003/2000000000, exceeding the old allowance2742255213/110000000000. At A1+A2+A3+H3+H4 (vertex2559), another has load4143329/500000000, exceeding408115897/55000000000. Their literal orbit charges, complete simplex-row coefficients and disjoint packing certificates are saved. The132 charge cuts can be exported to a future finite model, but no new global solve was run here.

After setting each vertex allowance to at least its covering-mixture mass, the source independently recounts **all432,021 exported model columns**, including96,176 zero-budget signed transfers and52 signed vertex-cost incidence columns. All1,908 quadrilateral masses remain fixed. The numerator is unchanged, while the new normalizer is

\[
D=19983545698211739397486458908069/
19819800000000000000000000000000.
\]

The exact resulting ratio is

\[
\frac{247898037929594052000000000000000}
{19983545698211739397486458908069}
=12.405107765824441\ldots>12.
\]

This certifies an obstruction for the inherited finite model together with the entire represented nonnegative local simplex-packing family. It does not settle selected simplex/prism joint charges at the changed masses; those require a separate scan. Unexpanded quadrilateral provenances, missing face and forest classes, other fields, and simultaneous ambient geometric compatibility remain outside this conclusion.
