# Support of the enlarged actual-corner dual

The new finite-model lower certificate is12.52181135130475… . Dropping every positive Q row outside the common actual-provenance map and recomputing all vertex loads gives

\[
\frac{1286851235806609200}{108107450811568633}
=11.903446304081223\ldots.
\]

Thus the currently optimized obstruction materially uses the unexpanded merged quadrilateral support. The restricted value below12 is not primal feasibility and does not show that these geometric classes are dispensable.

The [source](../scripts/fresh_actual_corner_enlarged_model_20261017_support.py) rebuilds the [enlarged finite model](fresh_actual_corner_enlarged_model_20261017.md), verifies its exact dual, and saves the [support and exact restricted certificate](fresh_actual_corner_enlarged_model_20261017_support.json). It performs no optimization.

The full dual has4,932 positive face masses:

- 1,465 maximal-simplex aggregates and1,173 maximal-edge incidence rows, totaling2,638 S rows.
- 397 explicitly split actual Q rows.
- 1,897 unselected merged Q rows, to which no new common provenance is assigned.

There are1,724 positive vertex allowances. Of the2,638 S rows,997 were absent from the preceding fixed-witness local-cover feature set. Their old packing coverage cannot be inferred from matching some other rows.

The restriction retains every positive S mass and all397 selected actual Q masses, removes the1,897 other Q masses, and computes the exact finite-model dual allowances again. No new S row is restored and no retained face mass changes. Every finite exported column is checked, including zero-budget transfers, the52 signed vertex-cost incidence columns, and the294 actual-corner columns.

All57 formerly positive raw23 classes remain positive after reoptimization, but their total Q mass falls from0.1164465876 to0.0610840531. Their objective contribution is0.610840531. Among the294 new raw23 columns,291 have positive load and22 are exactly tight in the recovered certificate. The52 older forced-incidence columns are all slack. These facts describe the new certificate; they do not establish a unique cause of the optimum.

Final generation54869 and fresh saved replay64600 terminated PASS. Source SHA256 is `fa709d98a2058fdef8fcc871be39bcbfda5f84d85863fb4189737ed28584c4d8`; JSON SHA256 is `c5455f404e0ce65dd4bb4a5a2516157f980d2f1a33ff1741e9cbda689c449a46`. Replay reconstructs the model and both load computations without an optimizer. The report asserts no complete joint-family coverage at these changed masses.

## Next actual-provenance cohort

The largest single-raw-source cohort among the unselected support consists of36 raw25 rows, all canonical prism core21. Each row has exactly one saved source context in the [combined catalogue](fresh_sqrt5_rank13_combined_disjoint_ports_20261017.json), whose core map is bound by the [quadrilateral separation catalogue](fresh_sqrt5_rank13_packing_quad_separation_20261017.json). Thus expanding this cohort requires no choice among alternatives and no duplicated face rows. Its total objective contribution is0.24393453; the first five rows contribute0.1747449582. This is a ranking of the current exact dual, not a new incompatibility theorem.

All these rows have demand6 and represent4-prisms. With the19 already selected positive raw25 rows, the proposed joint raw25 cohort has55 classes and440 finite corner occurrences. Complete actual-corner reconstruction against the existing raw23 family is the next targeted geometric calculation; it has not been run here.

| Face row | Saved source context | Exact face mass |
|---|---|---|
| 26798 | 28187 | 13977849/2000000000 |
| 26862 | 28251 | 69299801/10000000000 |
| 26832 | 28221 | 28713807/5000000000 |
| 26872 | 28261 | 5197193/1000000000 |
| 26864 | 28253 | 42653007/10000000000 |
| 26911 | 28300 | 16609773/10000000000 |
| 26796 | 28185 | 13975167/10000000000 |
| 26834 | 28223 | 13523059/10000000000 |
| 26851 | 28240 | 2641753/2000000000 |
| 26801 | 28190 | 12767123/10000000000 |
| 26927 | 28336 | 11314387/10000000000 |
| 26870 | 28259 | 136227/156250000 |
| 26789 | 28178 | 3661817/10000000000 |
| 26751 | 28140 | 3597231/10000000000 |
| 26743 | 28132 | 2832967/10000000000 |
| 26760 | 28149 | 24171/100000000 |
| 26765 | 28154 | 543479/2500000000 |
| 26827 | 28216 | 1072373/5000000000 |
| 26848 | 28237 | 397057/2000000000 |
| 26793 | 28182 | 1928917/10000000000 |
| 26749 | 28138 | 166539/1000000000 |
| 26787 | 28176 | 495217/5000000000 |
| 26741 | 28130 | 543889/10000000000 |
| 26803 | 28192 | 475531/10000000000 |
| 26758 | 28147 | 25503/625000000 |
| 26962 | 28513 | 182081/10000000000 |
| 26846 | 28235 | 53949/5000000000 |
| 26730 | 28119 | 52973/10000000000 |
| 26729 | 28118 | 14913/10000000000 |
| 26738 | 28127 | 2979/2000000000 |
| 26762 | 28151 | 23/400000000 |
| 26740 | 28129 | 13/312500000 |
| 26779 | 28168 | 37/5000000000 |
| 26778 | 28167 | 1/250000000 |
| 26850 | 28239 | 3/2500000000 |
| 26959 | 28481 | 3/10000000000 |

For comparison, the entire unselected support comprises1,897 rows and2,426 saved provenances, including434 rows with multiple provenances. The unique-provenance raw25 cohort avoids that unresolved allocation issue.
