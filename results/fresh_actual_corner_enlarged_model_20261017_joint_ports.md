# Selected-provenance joint-port separation of the new dual

After the [complete local simplex-packing repair](fresh_actual_corner_enlarged_model_20261017_local_cones.md), all17 selected joint-port families have been checked at the **new** face masses. They give57 additional violated columns and16 further vertex-allowance increases. The final exact repaired ratio is

\[
\frac{67685455440946980000000000000000}{5444966920826040346699633746793}
=12.43082950275823\ldots>12.
\]

This is a fixed dual satisfying the enlarged finite model, the full local simplex-block cone, and the stated selected-provenance joint families. It is not an optimized value of all these families together, not full geometric coverage of the model's merged Q rows, and not a dimension theorem.

The [source](../scripts/fresh_actual_corner_enlarged_model_20261017_joint_ports.py) and [saved exact scan](fresh_actual_corner_enlarged_model_20261017_joint_ports.json) retain all4,932 face masses. Every one of the2,638 positive S rows contributes its complete literal block data. The397 Q rows in the common actual-provenance map contribute all2,652 actual vertex-corner footprints. The other1,897 positive merged Q rows have **zero coefficients in these new joint charges**. No representative of those rows is chosen or inferred.

Zero coefficients on unselected full-star classes are valid restricted nonnegative charging: dropping their payments weakens the charge family. They do not establish that the unselected geometry is absent. The [support restriction](fresh_actual_corner_enlarged_model_20261017_support.md), whose value is11.90345, shows why that distinction matters here.

For each of815 contributing finite-forest types, the source enumerates every3-,4-, and5-port subset orbit. It averages the S-block intersection histogram and actual Q occupied-star histogram under the same forest action. For S blocks, the exact intersection sizes2–5 are retained; for Q stars, occupied edge counts1–4 are retained. Binomial sum identities check every block/corner multiplicity before any facet is evaluated. In particular, nonmonotone coefficients in the five-port facets are applied to the exact histograms rather than a cumulative approximation.

The families and counts are:

| Port subset size | Proved facets | Exact orbit inequalities | Violated columns |
|---|---:|---:|---:|
| 3 | strengthened triangle/square inequality | 66,345 | 10 |
| 4 | the selected dotted-prism facet | 138,032 | 13 |
| 5 | all15 frozen necessary-domain facets | 3,315,360 | 34 |
| Total | 17 | 3,519,737 | 57 |

Complete domain and histogram-stream digests, binomial-identity digests, all57 exact violated columns, and every required vertex mass are saved. Per-vertex maxima combine the families. Only vertex allowances increase; all face masses remain fixed. Therefore the previously rechecked finite-model columns—including52 signed vertex-cost incidence columns—and the new complete local simplex-packing covers remain valid by monotonicity.

Generation19509 and fresh saved replay88006 terminated PASS. The replay reruns the exact scan without an optimizer. Source SHA256 is `4e370da5cff28fb3c9caded5fcb80bb73cfd253666d42f81b04fbe26c52ff5fa`; JSON SHA256 is `109d41daab1b95ed4b1aeb889fb82faa4fabd6ade9297acc99d95fe734a99e36`.

Construction-audit's independent narrow source read passes the changed-alpha loading, explicit unselected-Q zeros, common group averaging, binomial totals, whole-face multiplicities, and monotone allowance repair. It is a source/accounting review, not a duplicate3,519,737-inequality arithmetic replay. The underlying finite-family theorems remain frozen inputs.
