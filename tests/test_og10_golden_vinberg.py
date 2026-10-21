import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_vinberg import (
    og10_field_gram,
    run_og10_golden_vinberg_search,
)
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.probe_og10_golden_vinberg import _digest


def test_og10_golden_vinberg_bound_one_prefix():
    search = run_og10_golden_vinberg_search(1)

    assert [str(value) for value in search.candidate_coordinates] == [
        "-1 + 1*phi",
        "1",
    ]
    assert len(search.cone_roots) == 10
    assert [len(level.roots) for level in search.levels] == [1, 2]
    assert len(search.roots) == 13


def test_archived_og10_bound_eight_hyperideal_witness():
    payload = json.loads(
        Path("data/candidates/og10_golden_vinberg_bound8.json").read_text()
    )
    roots_json = payload["simple_roots"]
    roots = tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in roots_json
    )
    gram = og10_field_gram()

    assert payload["verified"]
    assert payload["facet_count"] == 171
    assert payload["visited_node_count"] == 160234
    assert payload["pair_relation_counts"] == {
        "m2": 1509,
        "m3": 490,
        "m5": 800,
        "ultraparallel": 11736,
    }
    assert payload["roots_sha256"] == (
        "c20e58f95c44f7740f5ec36ee3245db811a363b142f6b9c4d84729c1c07e74d4"
    )
    assert _digest(roots_json) == payload["roots_sha256"]
    assert len(roots) == len(set(roots)) == 171
    assert all(_pair(root, gram, root) == 2 for root in roots)

    witness_payload = payload["hyperideal_witness"]
    witness = verify_hyperideal_extreme_ray(
        roots,
        gram,
        witness_payload["active_facets"],
        dimension=10,
    )
    assert witness.norm == Golden.from_json(witness_payload["norm"])
    assert witness.norm.sign(0) > 0
    assert all(value.sign(0) <= 0 for value in witness.wall_pairings)
    assert _digest([value.to_json() for value in witness.wall_pairings]) == (
        witness_payload["wall_pairings_sha256"]
    )
