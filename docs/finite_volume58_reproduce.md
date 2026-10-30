# Reproducing the completed general finite-volume bound58

For current reading, see the [current dimension-bounds manuscript](../arxiv/hyperbolic_coxeter_dimension_bounds/main.tex) and the
[submission guide](arxiv_submission_editions_20261028.md). The supporting
arguments and recorded verification scopes below remain distinct from those
current submission builds. Standalone draft typesetting routes have been retired.

Run commands from the repository root with Python3. The completed endpoint
is [fresh_a15_d15_dimension58_20261018](../results/fresh_a15_d15_dimension58_20261018.md).
It proves n≤58 for finite-sided, finite-volume real hyperbolic Coxeter
polytopes. Later ambient58 frontiers concern a proposed57 bound and are
not substitutes for this frozen certificate.

## Levels of verification

1. **Identity and saved composition.** The final binder checks source/report
   hashes, complete saved finite/ideal domains, exact polynomial moment
   identities, local/count interfaces and final curvature arithmetic.
2. **Finite arithmetic replay.** The finite closure verifies all-last-item
   bounds and attaining predecessors of the complete saved tables, every
   item identity, all type comparisons and every contextual option. It does
   not merely trust a saved optimum or an optimizer success flag.
3. **Local geometry replay.** The two final local producers rederive their
   physical core calculations, native-map domains and current literal
   interfaces. Their written arbitrary-support and route-completeness
   arguments remain part of the mathematical proof.
4. **Historical dependency replay.** The transitive manifest identifies
   older geometric, classification and count packets. The endpoint does
   not execute all of them. Repeating the entire history is a distinct
   validation project; neither a successful hash audit nor216 literal
   faces replaces a universal geometric theorem.
5. **Proof-assistant verification.** The plain arithmetic export supports
   kernel checks. The internal geometric lemmas and higher tail must also
   be translated. The current Lean project is partial; this guide does not
   claim an end-to-end Lean proof of58.

## Frozen endpoint and direct inputs

|Artifact|JSON SHA256|
|---|---|
|[dimension58 endpoint](../results/fresh_a15_d15_dimension58_20261018.json)|`60a6dbedc03b3c0690c80ab86f276e5767de786469b46468d0990e4a75d25985`|
|finite closure|`692e782f03d3dbc7ad230545819ccf37f6d9c69a8ed002bdeadbe8703c96c81e`|
|simultaneous baseline|`b192b2c6b9c27fa07a14d7cbf731a62c71c99da2de2460ea589e5d4a6d87a97d`|
|[completed59 higher tail](../results/fresh_depth6_dimension59_20261018.json)|`2b09dac21f742a216bd10229fe25fba3f1e4114b8f43fae1206dcb1840169bab`|

Each JSON's `sources` dictionary pins the source/data/report dependencies
used at generation. Preserve those files unchanged. Use `--verify` on a
frozen producer; do not regenerate over its certificate to hide a mismatch.

## Commands

The smallest endpoint replay is:

```sh
python3 scripts/fresh_a15_d15_dimension58_20261018.py --verify
```

Its expected status starts with
`PASS_COMPLETE_A15_D15_DIMENSION58_COMPOSITION`; the headline is dimension58,
finite upper783/2, ideal surplus91/2 and constant4350.

The complete finite arithmetic replay is:

```sh
python3 scripts/fresh_a15_d15_ambient59_forest_closure_20261018.py --verify
```

The final two geometric replay commands are:

```sh
python3 scripts/fresh_finite_a15_three_coupon_zero_transfer_20261018.py --verify
python3 scripts/fresh_finite_d15_a4_a11_half_zero_transfer_20261018.py --verify
```

These are separate checks. A failure should be investigated at its first
assertion and against the pinned inputs; do not change a frozen expected
result merely to obtain `PASS`. They were already run successfully when the
frozen proof was completed. The manuscript-translation task does not claim
to have repeated all these historical computations.

## Plain data for an independent arithmetic checker

The [exporter](../scripts/export_finite_volume58_arithmetic_20261018.py)
produces finite_volume58_arithmetic_20261018.json.
It reads the frozen count, computes only the exposed item/base/cross formulas,
and copies the saved certificates. It does not optimize, reconvolve or prove
the geometric interpretation. Its schema label is
`finite-volume58-arithmetic-v1` and its status explicitly says
`EXPORTED_FROZEN_DATA_NOT_A_FORMAL_PROOF`.

```sh
python3 scripts/export_finite_volume58_arithmetic_20261018.py --verify
```

For a clean checkout without the exported file, omit `--verify` once.
Generation refuses to overwrite an existing export. JSON round-trip
verification checks both the copied values and all source hashes. This
packaging check is separate from a mathematical checker of the values.

The representation is deliberately simple:

* finite costs are integer units1/16; the cap is59 and upper is6264;
* rational ideal values are `[numerator, positive_denominator]`;
* `actual_types`, `representative_types` and `representative_map` encode
  all180 actual finite type representatives and their simultaneous map;
* base arrays and actual/upper pair matrices are explicit;
* each of5148 `tables` has a key, fixed rank5 count, excluded-rank2 flag,
  item ranks and weights,60 values and60 predecessor ranks;
* each of2574 `expanded_rows` records its fixed representative indices,
  fixed cost, residual rank and complete contextual options;
* each of1890 `ideal_profiles` records ranks,59 polynomial coefficients,
  V,G, normalized cost and both exact margins;
* all58 normalized affine envelopes and final constants are explicit.

A sound independent checker should establish these assertions, not just
compare the final numbers:

1. Reconstruct every961 fixed vector and every2574 rank5 expansion and
   prove the residual split, including fixed A17 versus residual D17.
2. Verify simultaneous base/cross domination for every actual type and
   pair, and all item identities against the fixed components.
3. Check every last-item upper inequality. Induction then bounds arbitrary
   residual multisets. Predecessor equality supplies additive attainment
   but is not needed for an upper bound.
4. Prove coverage by the no-context, reserve-one and reserve-two branches,
   fixing all A2 multiplicities. Verify the two context formulas and every
   winning full-forest equality. Witness counts are presence guards.
5. Check every resulting upper≤6264. Do not assert actual attainment:
   the representative maximizer's actual cost is355.5.
6. Independently enumerate all rank58 terminal profiles with at most one
   positive rank below9. Recompute each polynomial product, not merely
   its symmetry, degree and coefficient sum.
7. Verify the exact curvature and cost margins, minimum normalized1/340,
   absolute minimum91/2, and the arithmetic783/2 and4350.
8. Connect these finite facts to the separately proved geometry,
   concentration, deficit-sum and higher-tail lemmas.

The export includes an `exclusions` list naming the geometric and all-rank
statements not justified by the exported arithmetic. They are internal
proof obligations, not permissions to add axioms for the final theorem.

## Current manuscript

The current dimension-bounds submission replaces the standalone bound58 draft.
Its [submission guide](arxiv_submission_editions_20261028.md) supplies the
source and typesetting command. The mathematical interfaces above retain
their original scope; the current ceiling57 has its own
[proof and reproduction map](finite_volume57_reproduce.md).

## Final arithmetic to keep visible

The upper and lower estimates concern the same deficit:

```text
D5 ≤ (783/2) vf + sum C_actual(v)
D5 ≥ (783/2) vf + sum kappa(v) + 4350
kappa(v) − C_actual(v) ≥ 91/2.
```

Their difference is `0 ≥ 4350 + (91/2) vi > 0`. The finite equality is
sufficient; a strict finite saving is unnecessary. The constant is added
once, and a negative terminal ideal margin would not justify multiplying
the concentration inequalities. The completed59 tail supplies dimensions
at least60, completing the bound58.
