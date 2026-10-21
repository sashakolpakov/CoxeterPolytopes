from __future__ import annotations

from functools import lru_cache
import hashlib
import json

from scripts.certify_enriques_type_i_downstairs_anti_gram_census import (
    OUTPUT,
    analyze,
)


ARCHIVE_SHA256 = "df864bed625d1411c7db46d7bbb500606ccb011297e420de8272e512ab2e22ea"


@lru_cache(maxsize=1)
def _payload() -> dict:
    return analyze()


def test_archive_replays_and_all_raw_assignments_are_partitioned() -> None:
    payload = _payload()
    assert payload == json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert hashlib.sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert payload["all_checks_pass"] is True
    alphabet = payload["downstairs_alphabet"]
    assert len(alphabet["m1_edges_zero_based"]) == 10
    assert len(alphabet["m2_edges_zero_based"]) == 3
    assert alphabet["m1_values"] == [-2, 2]
    assert alphabet["m2_values"] == [-4, 0, 4]
    assert alphabet["raw_assignment_count"] == 2**10 * 3**3 == 27648


def test_sheet_switch_quotient_has_exactly_eighteen_complete_orbits() -> None:
    quotient = _payload()["sheet_switch_quotient"]
    assert quotient["orbit_count"] == 18
    assert quotient["orbit_size_histogram"] == {
        "512": 2,
        "1024": 6,
        "2048": 10,
    }
    records = _payload()["exact_inertia_census"]["orbit_records"]
    assert len(records) == 18
    assert sum(record["raw_assignment_multiplicity"] for record in records) == 27648
    for record in records:
        components = record["support_component_count"]
        assert record["raw_assignment_multiplicity"] == 2 ** (12 - components)


def test_exact_inertia_and_unique_generic_geometric_orbit() -> None:
    census = _payload()["exact_inertia_census"]
    assert census["raw_assignment_histogram"] == {
        "p0_n10_z2": 512,
        "p0_n9_z3": 1024,
        "p1_n10_z1": 9216,
        "p1_n11_z0": 512,
        "p1_n9_z2": 2048,
        "p2_n10_z0": 14336,
    }
    assert census["sheet_switch_orbit_histogram"] == {
        "p0_n10_z2": 1,
        "p0_n9_z3": 1,
        "p1_n10_z1": 7,
        "p1_n11_z0": 1,
        "p1_n9_z2": 1,
        "p2_n10_z0": 7,
    }
    generic = _payload()["generic_rank_nine_theorem"]
    assert generic["unique_orbit"] is True
    assert generic["raw_assignment_count"] == 1024
    assert generic["canonical_gram_inertia"] == {
        "positive": 0,
        "negative": 9,
        "zero": 3,
    }
    assert generic["equals_explicit_Type_I_package_up_to_sheet_swaps"] is True


def test_rank_ten_signature_residue_is_killed_by_forced_saturation_root() -> None:
    residue = _payload()["rank_ten_residue"]
    assert residue["unique_signature_candidate_orbit"] is True
    assert residue["raw_assignment_count"] == 512
    assert residue["canonical_gram_inertia"] == {
        "positive": 0,
        "negative": 10,
        "zero": 2,
    }
    forced = residue["forced_saturation_root"]
    assert forced == {
        "downstairs_half_sum_integral": True,
        "downstairs_root_sum": [8, 4, 8, 14, 12, 10, 8, 6, 4, 2],
        "e_is_integral_anti_invariant_and_in_saturation": True,
        "e_square": -2,
        "formula_for_arbitrary_sheet_signs": (
            "e=(eps_11*delta_11+eps_12*delta_12)/2="
            "b_11^eps_11+b_12^eps_12-(u_11+u_12)/2"
        ),
        "four_times_e_square": -8,
        "sheet_switch_invariant": True,
        "wall_pair_one_based": [11, 12],
    }
    assert residue["status"] == (
        "ELIMINATED: forced saturated anti-invariant minus-two root violates L2"
    )
    assert "eliminated" in _payload()["scope"]["rank_ten_status"]
    assert len(_payload()["scope"]["remaining_rank_nine_compiler_work"]) == 3
    assert _payload()["scope"]["new_abstract_polytope_claimed"] is False
