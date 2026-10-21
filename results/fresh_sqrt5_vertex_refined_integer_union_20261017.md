# Integer port packing and vertex-labelled base capacity together

The exact union test gives

\[
 \frac{135976481921606560}{11375364259453521}
 =11.953593645021346\ldots<12.
\]

Thus the combined constraints remove this particular fixed-alpha obstruction. This is **not** a feasible primal payment vector, an optimum, complete face coverage or a dimension bound. No optimizer is used here.

The [vertex-labelled capacity certificate](fresh_sqrt5_vertex_refined_base_capacity_20261017.md) separately retained a ratio12.04399. The [integer simplex-port packing certificate](fresh_sqrt5_integer_simplex_port_packing_20261017.md), evaluated at the earlier actual-provenance alpha, separately retained12.07057. One cannot combine their separate normalizers. The [new source](../scripts/fresh_sqrt5_vertex_refined_integer_union_20261017.py) instead recomputes every packing load at the **vertex-refined alpha**, and raises each vertex mass to the maximum required by both families.

The packing export has62,353 columns:59,962 three-port matching orbit columns and2,391 full13-port integer capacities for divisors2,3,4. All5,462 positive simplex rows have complete exported features. The1,582 quad rows receive zero new simplex-packing coefficient. The positive support remains7,044 rows and797 vertex forests.

There are57 newly violated packing inequalities at36 forests. Repairing them adds exactly `23883145927/2880000000000` to the normalizer. The numerator remains `65373308616157/4950000000000`; the new normalizer is `3791788086484507/3432000000000000`.

The [saved exact result](fresh_sqrt5_vertex_refined_integer_union_20261017.json) recounts all211,339 old positive-cost columns, all59,962 mixed-three-port columns, all62,353 packing columns, all11,397 labelled signed transfers, all4,969 actual unlabelled transfers, and all5,693 inherited signed transfers. The1,858 explicit quad provenance masses and all merged mass equalities are unchanged and checked again. All these families therefore share one consistent certificate.

A successor optimization still needs a concrete combined model. A valid restricted-feature extension may assign zero new indicator payment on face classes outside the exported domain, while preserving every old row and charge. The1,582 affected merged quad rows must be replaced by all1,858 actual provenance rows when using the labelled transfer coefficients. Missing face classes and the231 untouched finite forests remain outside the original model's coverage. A single bounded solve is a separate task; this diagnostic contains none.

```sh
python3 scripts/fresh_sqrt5_vertex_refined_integer_union_20261017.py
python3 scripts/fresh_sqrt5_vertex_refined_integer_union_20261017.py --verify
```

Generation54664 and fresh saved replay66967 both terminated with `PASS_EXACT_LABELLED_BASE_INTEGER_PACKING_UNION`. The replay reconstructs the entire result and checks exact equality. Source and all inherited inputs are hash-bound; predecessor files were not edited.

SHA256:

- source: `1076e6c568141637832530e519170186a4dca39df6fc3b7f2c6f1ef4a7b93751`;
- JSON: `d5a4f6c1ca6b75ab84abb7bc1d241489d4b678505ceef129a85b1588a46e9f25`.
