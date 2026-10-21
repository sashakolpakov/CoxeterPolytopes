# Compact RA4: the C26 successor

The [complete written proof](../results/compact_ra3_twenty_six_vertices_no_regular_pentagon_20261020.md)
shows that a compact right-angled three-polyhedron with26 vertices has
no regular pentagonal face. It uses the published count-one C26 fullerene
classification and actual compact3D metric uniqueness as named inputs;
the explicit dual sphere, symmetry normalization and edge formulas are
derived here. Its [full root proof/checker review](../results/review_compact_ra3_c26_proof_and_checker_root_20261020.md)
passes.

Consequently every compact right-angled four-polytope has at least61
facets and volume at least47pi^2/6. The [execution review](../results/review_compact_ra3_c26_execution_root_20261020.md)
also derives the stronger uniform-facet charge. These results are
supporting successors; the frozen14-page compact volume manuscript is
unchanged. The unrestricted120-cell minimum and second volume remain open.

From the repository root, replay one small bounded job:

```sh
python3 -B scripts/run_bounded_research_check_20261018.py --parent-lock -- scripts/check_compact_ra3_c26_graph_metric_20261020.py --verify
```

The watchdog requires process-list access. It enforces the shared lock,
59-second wall limit,1GiB process-tree cap, one thread and the documented
resource gates. Do not bypass it. The stdlib-only fixed checker recomputes
the [9626-byte certificate](../results/compact_ra3_c26_graph_metric_20261020.json),
checking graph links/incidences/symmetries and exact Q(sqrt17) identities.
It does not re-enumerate fullerenes or prove the actual-Q symmetry bridge.

[Generation](../results/compact_ra3_c26_graph_metric_generation_20261020.json)
and [fresh replay](../results/compact_ra3_c26_graph_metric_replay_20261020.json)
passed in0.131s total, with peak56,033,280B. All scientific/control/self
hashes passed before and after. No Lean or TeX job ran. The
[immutable source manifest](../results/compact_ra3_c26_successor_20261020_freeze.json)
pins these sources and ledgers without duplicating them.
