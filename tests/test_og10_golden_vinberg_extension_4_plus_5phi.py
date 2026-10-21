import json
from pathlib import Path

from coxeter_polytopes.exact.golden import Golden
from coxeter_polytopes.search.bugaenko_vinberg import _pair
from coxeter_polytopes.search.og10_golden_vinberg import og10_field_gram
from coxeter_polytopes.verify.chamber import verify_hyperideal_extreme_ray
from scripts.probe_og10_golden_vinberg_extension import _digest
from scripts import probe_og10_golden_chamber_folding as folding


FINAL = Path(
    "data/candidates/og10_golden_vinberg_through_4_plus_5phi.json"
)
BOUND_TEN = Path(
    "data/candidates/og10_golden_vinberg_through_3_plus_4phi.json"
)


def _roots(payload):
    return tuple(
        tuple(Golden.from_json(entry) for entry in root)
        for root in payload["simple_roots"]
    )


def test_archived_og10_prefix_through_4_plus_5phi_and_ray_history() -> None:
    payload = json.loads(FINAL.read_text())
    archived = json.loads(BOUND_TEN.read_text())
    roots = _roots(payload)
    gram = og10_field_gram()

    assert payload["verified"] is True
    assert all(payload["checks"].values())
    assert payload["certificate_kind"] == (
        "bounded_exact_OG10_norm_two_Vinberg_prefix_through_4_plus_5phi"
    )
    assert payload["completed_through_control_coordinate"] == Golden(4, 5).to_json()
    assert payload["next_admissible_control_coordinate"] == Golden(3, 6).to_json()
    assert "exhaustive only through control coordinate 4+5*phi" in payload["scope"]
    assert "candidate-generation bound 13" in payload["scope"]

    assert payload["candidate_control_coordinates"] == [
        Golden(2, 5).to_json(),
        Golden(4, 4).to_json(),
        Golden(3, 5).to_json(),
        Golden(4, 5).to_json(),
    ]
    assert [level["root_count"] for level in payload["levels"]] == [71, 13, 371, 232]
    assert [level["visited_node_count"] for level in payload["levels"]] == [
        131946,
        6732,
        1143573,
        650813,
    ]
    assert payload["extension_visited_node_count"] == 1933064
    assert payload["archived_facet_count"] == 457
    assert payload["new_facet_count"] == 687
    assert payload["facet_count"] == len(roots) == 1144

    assert payload["simple_roots"][:457] == archived["simple_roots"]
    assert len(set(roots)) == len(roots)
    assert all(all(entry.is_integral for entry in root) for root in roots)
    assert all(_pair(root, gram, root) == 2 for root in roots)
    assert _digest(payload["simple_roots"]) == payload["roots_sha256"] == (
        "d4e898a0c4310108c95ac790bc9f7860a239b266b6a9bd051bdca91e944b9987"
    )

    starts = (457, 528, 541, 912)
    ends = (528, 541, 912, 1144)
    level_digests = (
        "b7cb33523501cc78cba2a6ef23933a78d87386c63ecff58644795a44afe1bf8d",
        "89c0432bc56791664d39e1a426a2f10da56c4a8f7f1e7eed2cd96d6d25c250bc",
        "ed2a48f27fe024be5cf6dbeba7d4ff8da9dde6563fd19622ca731ae94659cb55",
        "4b2189fa0c04e321ac4294d11d08da3af236ef31426fdcfd24ca92ac81c49cc4",
    )
    for level, start, end, digest in zip(
        payload["levels"], starts, ends, level_digests, strict=True
    ):
        assert all(
            root[0] == Golden.from_json(level["control_coordinate"])
            for root in roots[start:end]
        )
        assert _digest(payload["simple_roots"][start:end]) == digest == level[
            "roots_sha256"
        ]

    assert payload["pair_relation_counts"] == {
        "m2": 11629,
        "m3": 5459,
        "m5": 10897,
        "ultraparallel": 625811,
    }
    assert sum(payload["pair_relation_counts"].values()) == 1144 * 1143 // 2

    for history in payload["hyperideal_witness_history"]:
        prefix = roots[: history["facet_count"]]
        witness_payload = history["witness"]
        witness = verify_hyperideal_extreme_ray(
            prefix,
            gram,
            witness_payload["active_facets"],
            dimension=10,
        )
        assert witness.norm.sign() > 0
        assert witness.ray == tuple(
            Golden.from_json(entry) for entry in witness_payload["ray"]
        )
        assert witness.norm == Golden.from_json(witness_payload["norm"])
        assert _digest([value.to_json() for value in witness.wall_pairings]) == (
            witness_payload["wall_pairings_sha256"]
        )

    final = payload["hyperideal_witness"]
    assert final["active_facets"] == [0, 1, 2, 3, 4, 5, 6, 7, 11, 912]
    assert final["norm"] == {
        "a": "39643/6724",
        "b": "-48969/13448",
    }
    assert final["wall_pairings_sha256"] == (
        "24ddd39e09c1e07e2b629ed8c31282b9be54dd0f274cc5b9945833bc003c787e"
    )

    transport = payload["complete_chamber_involution_transport"]
    h = folding._compose(
        folding._root_reflection_matrix(
            roots[0], folding._covector(roots[0], gram)
        ),
        folding._pell_seed(),
    )
    root_lookup = set(roots)
    new_roots = tuple(
        dict.fromkeys(
            image
            for root in roots
            for image in (folding._apply(h, root),)
            if image not in root_lookup
        )
    )
    closed_roots = (*roots, *new_roots)
    closed_lookup = set(closed_roots)
    assert len(new_roots) == transport["images_outside_archived_prefix"] == 373
    assert len(closed_roots) == transport["involution_closed_root_count"] == 1517
    assert all(folding._apply(h, root) in closed_lookup for root in closed_roots)
    assert _digest(
        [[entry.to_json() for entry in root] for root in new_roots]
    ) == transport["new_roots_sha256"] == (
        "10cf5a825b4e2700283c63c6f9b65b9d894cf912fabc73f540c98e42fddbd585"
    )
    assert _digest(
        [[entry.to_json() for entry in root] for root in closed_roots]
    ) == transport["closed_root_set_sha256"] == (
        "ecc3ef936d3c8766899420d5917ba239e642608fa6b13e07a97f3a62b0663039"
    )
