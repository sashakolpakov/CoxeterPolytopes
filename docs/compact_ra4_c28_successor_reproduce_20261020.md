# Compact RA4: the C28 boundary-facet successor

20 October2026. Internally reviewed written proof and fixed exact replay
PASS. No external refereeing or end-to-end formal geometry is claimed.

The result is universal f>=65 and Vol>=49pi^2/6 for compact right-angled
hyperbolic four-polytopes. It strengthens the preceding C26 bounds
f>=61 and Vol>=47pi^2/6. A positive facet adjacent to a regular
dodecahedral facet must have at least30 vertices. The new local theorem
excludes an entire regular pentagonal face in every compact RA3 with28
vertices: the Loebell L7 type and both C28 fullerene types.

The [complete proof](../results/compact_ra3_twenty_eight_vertices_no_regular_pentagon_20261020.md)
starts with an actual polyhedron and realizes combinatorial symmetries
using metric uniqueness. In the tetrahedral type every pentagon has an
edge of cosh length4/3, below the regular value phi. In the other type,
regular-pentagon propagation and exact normal equations force three
distinct actual normals into a two-element unit-normal fiber.
The only new finite-classification premise is the published fullerene
count2 at28 vertices: [Brinkmann--Goedgebeur--McKay, Table2, page15](https://arxiv.org/pdf/1207.7010v2).
The two graphs themselves and all metric obstructions are explicit.

Full [root](../results/review_compact_ra3_c28_proof_root_20261020.md)
and [outside-author](../results/review_compact_ra3_c28_finite57_20261020.md)
written reviews pass. The [execution and charge review](../results/review_compact_ra3_c28_execution_root_20261020.md)
records the complete checker review and the boundary-excess-five charge.
It guarantees a uniform dodecahedral facet for W<=37. At volume no larger
than the120-cell, a non120-cell then requires a disconnected D graph;
the connected-D alternative requires W>=38 and f+W>=103.
Neither the unrestricted120-cell minimum nor its second volume is proved.

The [worker](../scripts/check_compact_ra3_c28_graph_metric_20261020.py)
checks two fixed dual spheres, their indicated symmetries and propagation,
and the rational/Q(sqrt5) identities. It does not enumerate fullerenes,
prove metric uniqueness or infer geometric existence from Gram data.
Its certificate states the external and geometric premises explicitly.

Reproduce serially from the repository root. These commands require the
usual resource gates: at least40GiB disk and6GiB available memory before
starting, one scientific process, a59-second/1GiB watchdog and one thread.
The repository runner owns the exclusive lock and cleans its child group.

    python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_compact_ra3_c28_graph_metric_20261020.py --verify

That command recomputes all checks and compares the complete certificate;
it does not overwrite frozen evidence. The original generation used
--write, which refuses an existing certificate. The optional --verify-code
mode is provenance only and is not a substitute for this full replay.

Generation and a separate fresh replay both passed:0.135s total,
56,754,176B peak process-tree RSS,29,156-byte certificate. The two compact
ledgers and source pins are in
`results/compact_ra3_c28_graph_metric_generation_20261020.json` and
`results/compact_ra3_c28_graph_metric_replay_20261020.json`.

The immutable successor manifest is
`results/compact_ra3_c28_successor_20261020_freeze.json`. It includes the
unchanged C26 packet's files and its manifest. The14-page primary compact
RA4 manuscript remains unchanged; these are supporting results in its
existing research packet. No new TeX or Lean build was needed.
