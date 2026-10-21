# Outside-cusp capacity: conditional Lean verification

## Current successor: ideal5 facet propagation

The [new facet-propagation theorem](../results/right_angled_ideal5_facet_cusp_propagation_20261019.md)
has two independent full written reviews:
[first](../results/review_ideal_ra5_facet_cusp_defect_batch_guard_20261019.md),
[second](../results/review_ideal5_facet_cusp_propagation_finite57_20261019.md).
It gives3V(P5)>=13V(Q4)+6F(Q4) for every four-facet and strengthens
the current necessary counts to62/206 in ideal5 and269/5875 in ideal6.
The corresponding entirely ideal six-face fraction in mixed11/12 is
<15/242. The five-face fraction remains<2/23.

The [successor Lean module](../formal/CoxeterProofs/RightAngledCuspPropagation.lean)
checks the propagation and these numerical consequences from the exact
outside-vertex identity, polygon-corner counts and9G<=2N0 actual-opposite
incidence. These geometric interfaces remain explicit hypotheses. The
new geometric inequality needs no published ideal4 classification;
its numerical applications retain the earlier carrier and transport inputs.

The source hash is
`0952426c063a438d5ab9cc6028b3dc51b5487d44f616c4d3afaae0f3d20b4baa`.
[Compilation](../formal/.proof-batches/CoxeterProofs_RightAngledCuspPropagation.json)
took1.02s, and the [separate axiom audit](../formal/.proof-batches/CoxeterProofs_RightAngledCuspPropagation_axiom_audit.json)
took0.34s. Only one new module ran; its two dependencies had capped cached
objects. Gates recorded190.24GiB disk and7.93/7.86GiB available RAM,
one thread,2048MiB. No warning or `sorryAx`; only standard foundational
axioms. Temporary outputs were removed.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledCuspPropagation --audit-axioms
```

The earlier stages below retain their exact verification histories.

## Earlier successor: positive cusp defect

The subsequent [positive cusp-defect theorem](../results/right_angled_ideal4_cusp_defect_research_20261019.md)
has both a [root written review](../results/review_right_angled_ideal4_cusp_defect_root_20261019.md)
and an [independent full review](../results/review_right_angled_ideal4_cusp_defect_finite57_20261019.md).
Its then-stronger necessary counts are62 facets/189 vertices in ideal5
and252 facets/5504 vertices in ideal6; the density fractions in mixed
dimensions11/12 improve to<2/23 and<1/15 respectively. No dimension
ceiling changes.

The successor [Lean module](../formal/CoxeterProofs/RightAngledCuspDefect.lean)
proves the following logic from individual geometric count interfaces:
small polygon corners force good cusps; good opposite facets force
2N0+q>=6 and N0>=1; exact outgoing-edge/triangle counts then give
V(Q)>=3V(R)+3. It combines that result with the explicit carrier and
facet lists to derive62/189 and252/5504, including integer rounding.
The spinor recurrence, actual-wall geometry and the interpretation of
the count interfaces remain written premises, not kernel-formalized
hyperbolic objects. No final size bound is assumed as an axiom.

The single-module batch and printed-axiom audit both passed:
[compile ledger](../formal/.proof-batches/CoxeterProofs_RightAngledCuspDefect.json),
[audit ledger](../formal/.proof-batches/CoxeterProofs_RightAngledCuspDefect_axiom_audit.json).
Source SHA256: `928fe9775396600704d610041f2142fe88432722d276a4471800023e2cc32cd8`.
Compilation took0.59s and the audit0.43s, with one thread and2048MiB
cap, after fresh gates at190.89GiB disk and7.31/7.29GiB available RAM.
Only `propext`, `Classical.choice` and `Quot.sound` were printed; no
`sorryAx` or warning. Temporary outputs were removed.

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledCuspDefect --audit-axioms
```

The earlier capacity-only verification below is preserved as a predecessor.

## Earlier capacity-two counts

The [written geometric proof](../results/right_angled_outside_cusp_capacity_20261019.md)
and its [independent full review](../results/review_right_angled_outside_cusp_capacity_finite57_20261019.md)
give the necessary counts59 facets/180 vertices in ideal dimension5 and
240 facets/4980 vertices in ideal dimension6. The new
[Lean module](../formal/CoxeterProofs/RightAngledCuspCapacity.lean)
now verifies their integer logic and the updated mixed-density consequences.
It does not formalize the hyperbolic geometry or exclude a dimension.

|Lean result|Proved internally|Explicit premises still requiring their written geometric/published interpretation|
|---|---|---|
|`dimension4_capacity_deficit`|2V=6V(F)+2N0+N1, hence V≥3V(F).|Outside vertices partition into zero/one/two-neighbor classes; outgoing edges total4V(F). The geometric at-most-two lemma is not assumed as a final numerical bound.|
|`ideal5_counts`|The carrier has≥33 vertices and≥25 facets; ambient facets≥59; V≥3F+3≥180.|The≥11-vertex three-facet bridge, capacity, actual facet injections, published minimum24 for every four-facet, complete coface incidence and box/Euler identities. The positive carrier is explicitly one of the listed actual four-facets.|
|`ideal6_counts`|Derives every five-facet's59/180 bounds, then F≥240,5643F≤272V and V≥4980.|A list of the actual five-facet interfaces, actual facet injection, complete coface incidence, and the previously proved universal transport inequality99f4≤99F+136V. No final size bound is assumed.|
|Density and mixed-ratio theorems|Ideal-five fractions<4/43 in dimensions11/12 and<14/43 in9/10; ideal-six fraction<5/71 in11/12;13312vi<20669vf and6656vi<5637vf in11 and12 respectively.|The cleared published average bounds and the written incidence/finite-cover inequalities, with the source-access qualifications of the original reports.|

The one-module batch and the separate axiom audit both passed, one thread
with a2048MiB cap. Final source SHA256:
`66d885ff5d8b6a61bbe076ce45dd4b71b5186a59ef462bd8a665209c673e8c78`.
Compilation took0.45seconds and the audit0.43seconds. The audit prints
only `propext`, `Classical.choice` and `Quot.sound`; exact rounding uses
no axioms. There is no `sorryAx` or source warning. The
[compilation ledger](../formal/.proof-batches/CoxeterProofs_RightAngledCuspCapacity.json)
and [axiom-audit ledger](../formal/.proof-batches/CoxeterProofs_RightAngledCuspCapacity_axiom_audit.json)
bind the source, compiler and resource settings. Temporary outputs were
removed by the batch runner.

Reproduce from the repository root:

```sh
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledCuspCapacity --audit-axioms
```

This module is a separate successor interface. The earlier
`ManuscriptRightAngled` source and its assumption register retain their
original scope and status; they are not retroactively described as
compiled by this batch. The existing ideal5/6 TeX still needs the new
capacity and Figure7 arguments integrated as described in the
[integration plan](right_angled_capacity_integration_plan_20261019.md).
