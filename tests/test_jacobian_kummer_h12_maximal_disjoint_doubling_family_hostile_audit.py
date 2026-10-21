from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path

from scripts.audit_jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile import (
    OUTPUT,
    PINNED_SHA256,
    build,
)


@lru_cache(maxsize=1)
def _payload() -> dict:
    return build()


def test_hostile_archive_is_an_independent_exact_rebuild() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert payload["verdict"] == "PASS_FREEZE"
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["pinned_inputs"] == PINNED_SHA256
    independence = payload["independence"]
    assert independence["production_certifier_imported"] is False
    assert independence["production_certificate_read"] is False
    assert independence["raw_172_facet_Gram_parsed_directly"] is True
    assert "maximal-independent-set" in independence["graph_method"]


def test_raw_gram_recovers_exact_doubling_stock() -> None:
    raw = _payload()["raw_relation_replay"]
    assert raw["facet_count"] == 172
    assert raw["doubling_wall_count"] == 124
    assert raw["doubling_walls_one_based"] == [
        31,
        32,
        33,
        36,
        39,
        42,
        48,
        51,
        56,
        *range(58, 173),
    ]
    assert raw["whole_relation_histogram"] == {
        "m2": 2820,
        "m3": 144,
        "m4": 800,
        "parallel": 1750,
        "ultraparallel": 9192,
    }
    assert raw["doubling_stock_relation_histogram"] == {
        "m2": 264,
        "parallel": 1138,
        "ultraparallel": 6224,
    }
    assert all(raw["checks"].values())


def test_independent_conflict_recursion_recovers_both_unique_maxima() -> None:
    payload = _payload()
    disjoint = payload["maximal_disjoint_graph"]
    ultra = payload["maximal_ultraparallel_graph"]

    assert disjoint["edge_count"] == 7362
    assert disjoint["degree_histogram"] == {"99": 12, "120": 16, "121": 96}
    assert disjoint["maximal_clique_count"] == 100
    assert disjoint["maximal_clique_size_histogram"] == {
        "36": 9,
        "55": 36,
        "74": 42,
        "93": 12,
        "112": 1,
    }
    assert disjoint["maximum_size"] == 112
    assert disjoint["unique_maximum_one_based"] == list(range(61, 173))
    assert disjoint["unique_maximum_relation_histogram"] == {
        "parallel": 664,
        "ultraparallel": 5552,
    }
    assert disjoint["conflict_matching_size"] == 12
    assert disjoint["matching_upper_bound"] == 112

    assert ultra["edge_count"] == 6224
    assert ultra["degree_histogram"] == {"56": 12, "64": 16, "112": 96}
    assert ultra["maximal_clique_count"] == 61
    assert ultra["maximal_clique_size_histogram"] == {
        "32": 32,
        "57": 12,
        "61": 16,
        "96": 1,
    }
    assert ultra["maximum_size"] == 96
    assert ultra["unique_maximum_one_based"] == list(range(77, 173))
    assert ultra["unique_maximum_relation_histogram"] == {"ultraparallel": 4560}
    assert ultra["conflict_matching_size"] == 28
    assert ultra["matching_upper_bound"] == 96
    assert all(disjoint["checks"].values())
    assert all(ultra["checks"].values())


def test_full_picard_core_and_allcock_scope_are_exact() -> None:
    payload = _payload()
    carrier = payload["full_Picard_inheritance"]
    assert carrier["lattice"] == "U(2)+D5(-1)^2+<-4>"
    assert carrier["rank"] == 13
    assert carrier["signature"] == [1, 12, 0]
    assert carrier["moduli_dimension"] == 7
    assert carrier["disjoint_112_exterior_core"]["walls_one_based"] == list(
        range(1, 61)
    )
    assert carrier["ultraparallel_96_exterior_core"]["walls_one_based"] == list(
        range(1, 77)
    )
    assert carrier["disjoint_112_exterior_core"]["Smith_diagonal"] == [1] * 13
    assert carrier["ultraparallel_96_exterior_core"]["Smith_diagonal"] == [1] * 13
    assert all(carrier["checks"].values())

    transfer = payload["Allcock_transfer"]
    source = transfer["primary_source"]
    assert source["doi"] == "10.2140/gt.2006.10.737"
    assert "free product" in source["section_4_setup"]
    assert "any subtree" in source["Theorem_4_1"]
    assert "index I" in source["Corollary_4_2"]
    assert "exponentially" in source["Theorem_4_3"]
    application = transfer["verified_application"]
    assert application["dimension"] == 12
    assert application["seed_status"] == "finite_volume_noncompact"
    assert application["generator_counts"] == [112, 96]
    assert "|V(T)|" in application["index_and_volume"]
    assert transfer["seed_gate"]["cusps"] == 6214
    assert transfer["seed_gate"]["hyperideal_rays"] == 0
    assert all(transfer["seed_gate"]["checks"].values())


def test_claim_boundary_and_hostile_report_are_frozen() -> None:
    boundary = _payload()["claim_boundary"]
    assert boundary["maxima_scoped_to_exact_124_seed_doubling_walls"] is True
    assert boundary["arbitrary_descendant_roots_classified"] is False
    assert boundary["worldwide_priority_claimed"] is False
    assert boundary["full_nef_or_effectivity_claimed"] is False
    assert boundary["surface_automorphism_claimed"] is False

    report = (
        Path(__file__).resolve().parents[1]
        / "results/"
        "jacobian_kummer_h12_maximal_disjoint_doubling_family_hostile_audit.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "pass / freeze" in normalized
    assert "maximal independent sets" in normalized
    assert "theorem 4.1" in normalized
    assert "theorem 4.3" in normalized
    assert "smith 1^13" in normalized
    assert "worldwide priority" in normalized
