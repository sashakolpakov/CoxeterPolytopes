# Current ideal5/6 conditional Lean successor: source handoff

19 October 2026. **SOURCE_PREPARED_NOT_COMPILED.** No Lean command, proof
batch, axiom audit or scientific program ran in the authoring lane. The
coordinator owns the final resource-gated execution. Existing sources,
catalogue entries, object files and verification ledgers were not changed.

New import-style module:
`CoxeterProofs.RightAngledElevenVertexCounts`.
Its [source](../formal/CoxeterProofs/RightAngledElevenVertexCounts.lean) is
157 lines/7,278 bytes, SHA256
`c1a0f052a7e68feef3adac9e94775b9f2018ffa54f7578db0e645c2be4f42d5d`.
It imports only `CoxeterProofs.RightAngledCuspPropagation` and adds no
custom axiom, proof placeholder, enumeration or data producer.

The primary manuscript is
[ideal_right_angled_six_face_constraints.tex](../manuscripts/ideal_right_angled_six_face_constraints.tex),
SHA256 `bf05cb7b634bc4e5c89822774d32f2c18bbefe24338737d41cef0aa67003184f`.
This module follows its current37/25 carrier,63/211 ideal5,275/6106 ideal6,
`272V>=6039F`, and strict mixed fractions4/47 and15/248. No superseded
final count theorem or density theorem is used.

## Exact implication and premise map

All theorem names below are in `CoxeterProofs.RightAngledElevenVertexCounts`.

|Theorem|Derived conclusion|Explicit inputs / manuscript location|
|---|---|---|
|`carrier_bounds`|Selected actual carrier has at least37 vertices and25 facets|The inherited `+3` defect interfaces, internal bridge giving a selected three-facet with at least11 vertices, independent injection `f(Q)>=2V(R)+3`, and NEW `M=11 -> V(Q) != 36`; `thm:plus3`, `thm:eleven`, `lem:carrier`.|
|`ideal5_counts`|At least63 facets, `3V>=631`, hence at least211 vertices|Actual carrier-to-polytope facet injection and the earlier primitive count implication `3V>=13m+6f`; `thm:fiveprop`, `thm:fivebounds`.|
|`sum_ideal5_facet_minima`|Sum of ALL listed five-facet facet counts is at least63 times their number|One explicit ideal5 interface per listed facet; no final bound is assumed.|
|`ideal6_counts`|At least275 facets, `6039F<=272V`, at least6106 vertices|Complete five-facet list, actual facet injection, exact two-coface incidence, and internal reciprocal transport `99a4<=99F+136V`; `thm:sixbounds`.|
|`exact_rounding`|Exact631 and1660725 products and their211/6106 integer roundings|No geometric hypothesis; closed natural-number arithmetic.|
|`sum_ideal6_facet_minima`|Ideal six-face facet-count sum is at least275 times the number of such faces|The proved ideal6 implication applied to each explicit interface.|
|`mixed_density`|Positive total face counts and `47I5<4N5`, `248I6<15N6`|Complete ideal/nonideal face partitions, nonideal minima16/27, and strict averages20/42; manuscript `sec:mixed`.|

The last two inequalities are exactly the cleared forms of
`I5/N5<4/47` and `I6/N6<15/248`. The theorem also proves `N5,N6>0`, so
the fractional interpretation has no missing zero-denominator condition.
It does not formalize real or rational division, nor a dimension exclusion.

## What the structures assume geometrically

`Ideal5Data.propagation` reuses the earlier explicit actual-count data.
Its lower layers include the positive-carrier/small-facet bridge, outside
vertex partition, outgoing-edge and triangular correspondence, polygon
small-corner inequality, paired-corner upper bound and good-opposite
incidence. Their existing logic derives the `+3` defect. The propagation
layer supplies the separate five-dimensional outside identity and its
corner/incidence counts. The new proof invokes only
`RightAngledCuspDefect.stronger_carrier_bound` and
`RightAngledCuspPropagation.facet_vertex_propagation` from these modules,
not their older final numerical bounds. Extra published-minimum/Euler
fields in the inherited structure remain visible but are not used to
assume the new constants.

`Ideal5Data.elevenVertexEqualityExcluded` is explicitly an INTERNAL
geometric premise: a carrier containing the selected11-vertex three-facet
cannot have36 vertices. It is not described as published. Its human
source is `results/right_angled_ideal4_eleven_vertex_equality_exclusion_20261019.md`,
SHA256 `b9af9b497c69da69a25faced7067b854f0e3b8332f89da57daa12a5a968c2121`.
Together with the previously derived `3M+3` inequality it gives37 if M=11;
for M>=12 that inequality already gives at least39.

`Ideal6Data` requires the actual complete list of five-facets, actual
injection, exact coface incidence and universal reciprocal transport.
The latter represents the internal good-edge transport and sharp103/132
all-moment argument, not a sampled frame result or a published final bound.

`MixedFaceData` represents complete face partitions for an actual mixed
right-angled polytope in ambient dimension11 or12. Each ideal face carries
the applicable internal structure above. The general minima16/27 and the
strict Nikulin--Khovanskii averages20/42 are the separately stated published
premises, including their dimension range and face-inheritance hypotheses.
These are the only newly published numerical inputs; the ideal63/275
minima are proved by this source. Geometric realization, completeness and
incidence of the lists are not proved by merely providing Lean lists.

None of the Lorentz models, actual wall feasibility, polygon spinors,
eleven-vertex equality exclusion, small-facet classification bridge,
reciprocal transport or hyperbolic realization is fully formalized here.
Kernel verification would certify the implications from these explicit
premises. Axiom-printer output is NOT an assumption inventory: structure
fields and theorem parameters do not appear there.

## Coordinator-only capped compile and real axiom audit

Prepared command, **not executed here**:

```text
python3 -B scripts/run_lean_batches_signal_safe_20261018.py CoxeterProofs.RightAngledElevenVertexCounts --audit-axioms
```

The official runner owns the shared nonblocking lock, one-module/one-thread
2048MiB Lean batches, start/abort resource gates, input snapshots and cleanup.
It compiles the new module and runs a separate actual axiom-printer pass.
Its dependency planner may require earlier sequential batches if their
objects or ledgers are stale; it never uses an unrestricted `lake build`.
The three predecessor `.olean` paths exist by read-only metadata inspection,
but their current capped provenance is for the runner to validate, not a
new check claimed here. Preserve and report any legacy-cache qualification.

The new source contains seven `#print axioms` commands, for the seven
mapped principal results above. Thus no extra audit-only source is needed.
Require actual successful compilation and complete bounded axiom output,
with no `sorryAx`, before updating its status to conditional kernel verified.
Record compile and audit ledgers separately from the paper's mathematical
review and TeX build. Standard foundational axioms alone would not remove
the explicit geometric assumptions described above.

Read-only predecessor source hashes at handoff:

|Source|SHA256|
|---|---|
|`RightAngledCuspCapacity.lean`|`66d885ff5d8b6a61bbe076ce45dd4b71b5186a59ef462bd8a665209c673e8c78`|
|`RightAngledCuspDefect.lean`|`928fe9775396600704d610041f2142fe88432722d276a4471800023e2cc32cd8`|
|`RightAngledCuspPropagation.lean`|`0952426c063a438d5ab9cc6028b3dc51b5487d44f616c4d3afaae0f3d20b4baa`|
