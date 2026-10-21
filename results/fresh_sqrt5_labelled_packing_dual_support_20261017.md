# The new obstruction persists on already exported face classes

The exact dual from the [one-shot enlarged model](fresh_sqrt5_labelled_packing_full_model_20261017.md) has5,638 positive face masses and758 positive vertex-forest masses. A read-only exact restriction shows that the obstruction is **not explained by the zero new payments on unexported face classes**.

Keeping only masses in the previously exported classes, and checking every original, mixed, integer-packing and signed transfer constraint, gives

\[
 \frac{129790338384790720}{10289212336025873}
 =12.61421517469832\ldots>12.
\]

This certificate has4,037 positive faces and687 positive vertex forests. The signed repair adds **no** outside simplex row. Therefore merely extending the current indicator features to previously unexported classes cannot invalidate this witness: it assigns those classes zero mass. New incidence constraints or other structural refinements can still change the conclusion. This remains a restricted-model lower bound, not a dimension theorem or a globally realizable polytope distribution.

## Full dual support

| Face class | Positive rows | Objective contribution |
| --- | ---: | ---: |
| Exported simplex stars |3,135 |10.04100842392202 |
| Unexported simplex stars |439 |0.104816748689697 |
| Exported quad provenance rows |902 |2.5648907976 |
| Unexported quad rows |1,162 |0.1975242564 |

There are55 positive vertex forests absent from the earlier797-forest fixed witness. The largest vertex masses nevertheless remain at familiar exceptional forests:

- `2A1 + A4 + H3 + H4`;
- `A1 + A3 + H3 + H4 + I2(5)`;
- `2A1 + A3 + A4 + H4`;
- `A1 + A5 + H3 + H4`;
- `2A1 + H3 + 2H4`.

The largest objective row is the already exported4-simplex face599, with mass `106996309/5000000000` and demand20. Its corner forests have model indices1844,2254,2902,3137,3324. The next highest rows576,241,520 and235 are also already exported4-simplex stars.

For comparison, retaining only the outside-class masses and restoring any containing-simplex masses forced by signed constraints gives the weaker exact lower bound

`65645797137782400/25317360834447833 = 2.5929162825084853...`.

That repair introduces504 simplex masses. Neither restricted certificate is asserted optimal. The exported-only certificate is sufficient for the targeted conclusion above.

## Exact data and scope

The [source](../scripts/fresh_sqrt5_labelled_packing_dual_support_20261017.py) reconstructs the complete preserved model and verifies its digest against the saved one-shot artifact. It rechecks the full saved dual, identifies each actual provenance copy with its inherited class, and exports every positive face's kind, required demand, exact mass, objective contribution and touched vertex forests. It also stores the30 largest face contributions and vertex masses and both exact restricted duals in the [JSON](fresh_sqrt5_labelled_packing_dual_support_20261017.json).

No optimization, new face enumeration or parameter search occurs here. Restricting a dual is followed by the necessary signed-column repair and a full exact load recount; face masses are not simply discarded without checking transfer obligations.

Generation78560 and fresh saved replay91574 both terminated with `PASS_EXACT_ENLARGED_MODEL_DUAL_SUPPORT`. The replay reconstructs and compares these exact results. The script and all dependencies are hash-bound. Predecessor files were not modified.

```sh
python3 scripts/fresh_sqrt5_labelled_packing_dual_support_20261017.py
python3 scripts/fresh_sqrt5_labelled_packing_dual_support_20261017.py --verify
```
