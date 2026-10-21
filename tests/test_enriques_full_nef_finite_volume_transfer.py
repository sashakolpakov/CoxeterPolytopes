from __future__ import annotations

import json

from scripts.certify_enriques_full_nef_finite_volume_transfer import OUTPUT, analyze


def test_archive_is_full_replay() -> None:
    assert analyze() == json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_all_algebraic_and_hyperbolic_conditions_are_equivalent() -> None:
    document = analyze()
    equivalence = document["equivalence"]
    nodes = set(equivalence["nodes"])
    assert equivalence["all_nodes_mutually_equivalent_in_scope"] is True
    for reachable in equivalence["transitive_closure"].values():
        assert set(reachable) == nodes


def test_exact_global_consequence() -> None:
    consequence = analyze()["exact_seven_image_consequence"]
    assert consequence == {
        "image_cardinality": 7,
        "facet_count_multiset": [12, 12, 20, 20, 20, 20, 20],
        "finite_vertices": 2214,
        "ideal_vertices": 294,
        "hyperideal_rays": 0,
        "elliptic_rank_8_edges": 13182,
        "full_affine_neighborhood_audits": 584,
        "all_finite_volume": True,
        "all_noncompact": True,
        "no_compact_full_nef_enriques_member": True,
        "no_fourteen_facet_full_nef_enriques_member": True,
    }
