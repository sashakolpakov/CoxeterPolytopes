import json
from collections import Counter
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from coxeter_polytopes.verify.chamber import _dot, _mat_vec, verify_hyperideal_extreme_ray
from scripts.probe_og10_golden_vinberg_extension import _digest, _relation


FINAL = Path("data/candidates/og10_golden_vinberg_through_3_plus_4phi.json")
CHECKPOINT = Path("data/candidates/og10_golden_vinberg_through_2_plus_4phi.json")
BOUND_EIGHT = Path("data/candidates/og10_golden_vinberg_bound8.json")


def _roots(payload):
    return tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )


def test_archived_og10_bound_ten_extension_and_hyperideal_witness():
    payload = json.loads(FINAL.read_text())
    roots_json = payload["simple_roots"]
    roots = _roots(payload)
    gram = og10_field_gram()

    assert payload["verified"]
    assert all(payload["checks"].values())
    assert payload["candidate_generation_identity_bound"] == 10
    assert payload["candidate_control_coordinates"] == [
        Golden(2, 4).to_json(),
        Golden(3, 4).to_json(),
    ]
    assert payload["completed_through_control_coordinate"] == Golden(3, 4).to_json()
    assert [level["root_count"] for level in payload["levels"]] == [150, 136]
    assert payload["extension_visited_node_count"] == 557014
    assert payload["facet_count"] == len(roots) == 457
    assert payload["new_facet_count"] == 286
    assert payload["immediate_new_facet_count"] == 136
    assert payload["previous_hyperideal_ray_status"] == "survives_strictly"

    assert len(set(roots)) == len(roots)
    assert all(all(entry.is_integral for entry in root) for root in roots)
    assert all(_pair(root, gram, root) == 2 for root in roots)
    assert all(root[0] == Golden(2, 4) for root in roots[171:321])
    assert all(root[0] == Golden(3, 4) for root in roots[321:])
    assert _digest(roots_json) == payload["roots_sha256"] == (
        "e86b65bb5348f5c0e645e482b26eb6006c20245a578249fea8c8707c7bbdcd9b"
    )

    functionals = tuple(_mat_vec(gram, root) for root in roots)
    relation_counts = Counter(
        _relation(_dot(left, functionals[right_index]))
        for left_index, left in enumerate(roots)
        for right_index in range(left_index + 1, len(roots))
    )
    assert dict(sorted(relation_counts.items())) == payload["pair_relation_counts"] == {
        "m2": 4413,
        "m3": 1810,
        "m5": 3266,
        "ultraparallel": 94707,
    }

    witness_payload = payload["hyperideal_witness"]
    witness = verify_hyperideal_extreme_ray(
        roots,
        gram,
        witness_payload["active_facets"],
        dimension=10,
    )
    assert witness.ray == tuple(
        Golden.from_json(entry) for entry in witness_payload["ray"]
    )
    assert witness.active_facets == (0, 1, 2, 3, 4, 5, 6, 7, 11, 225)
    assert witness.norm == Golden(27755, 0) / 24964 - Golden(0, 32489) / 49928
    assert witness.norm.sign(0) > 0
    assert _digest([value.to_json() for value in witness.wall_pairings]) == (
        witness_payload["wall_pairings_sha256"]
    )
    assert all(value.sign(0) < 0 for value in witness.wall_pairings[321:])

    bound_eight = json.loads(BOUND_EIGHT.read_text())
    old_ray = tuple(
        Golden.from_json(entry) for entry in bound_eight["hyperideal_witness"]["ray"]
    )
    signs = Counter(_pair(old_ray, gram, root).sign(0) for root in roots[171:321])
    assert signs == {-1: 147, 0: 2, 1: 1}


def test_archived_og10_first_extension_checkpoint():
    checkpoint = json.loads(CHECKPOINT.read_text())
    final = json.loads(FINAL.read_text())

    assert checkpoint["verified"]
    assert all(checkpoint["checks"].values())
    assert checkpoint["candidate_generation_identity_bound"] == 10
    assert checkpoint["candidate_control_coordinates"] == [
        Golden(2, 4).to_json()
    ]
    assert checkpoint["completed_through_control_coordinate"] == Golden(2, 4).to_json()
    assert [level["root_count"] for level in checkpoint["levels"]] == [150]
    assert checkpoint["extension_visited_node_count"] == 237278
    assert checkpoint["facet_count"] == len(checkpoint["simple_roots"]) == 321
    assert checkpoint["new_facet_count"] == 150
    assert checkpoint["immediate_new_facet_count"] == 150
    assert checkpoint["previous_hyperideal_ray_status"] == "cut"
    assert checkpoint["roots_sha256"] == (
        "b671e733034e112f3bc4c0d4e3b3eadf9b13a7a247925d647c9ba5dbc7f7ddbb"
    )
    assert _digest(checkpoint["simple_roots"]) == checkpoint["roots_sha256"]
    assert checkpoint["simple_roots"] == final["simple_roots"][:321]
