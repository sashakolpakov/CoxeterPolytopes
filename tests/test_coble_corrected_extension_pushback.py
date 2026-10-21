from __future__ import annotations

from functools import lru_cache
import json
from pathlib import Path
import subprocess
import sys

from scripts.certify_coble_corrected_extension_pushback import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/coble_corrected_extension_pushback.json"
SCRIPT = ROOT / "scripts/certify_coble_corrected_extension_pushback.py"


@lru_cache(maxsize=1)
def certificate() -> dict:
    return build_certificate()


def test_corrected_gate_is_closed_immersion_not_bpf_separable() -> None:
    payload = certificate()
    language = payload["corrected_extension_class"]
    assert language["semantic_status"] == (
        "P_INDEXED_FINITE_GEOMETRIC_CERTIFICATE_CLASS_NOT_UNQUALIFIED_P_ONLY"
    )
    assert "existence of the coordinate/Hilbert witness is part of membership" in (
        language["hostile_quantifier_warning"]
    )
    moving = language["tower_gate"]["moving_curve_step"]
    assert "relative closed immersion" in moving["required"]
    assert "length-one and length-two" in moving["required"]
    hostile = payload["hostile_audit"]
    assert "O(2,0)" in hostile["rejected_old_clause"]
    assert "fails to separate two points" in hostile["why_corrected_gate_rejects_it"]
    assert "alone does not separate" in hostile["jet_warning"]


def test_parameter_tower_and_final_hilbert_gates_close_nonemptiness() -> None:
    payload = certificate()
    tower = payload["corrected_extension_class"]["tower_gate"]
    assert "constant-rank projective bundle" in tower["irreducible_parameter_induction"]
    assert "Jouanolou" in tower["moving_curve_step"]["good_open"]
    assert "absolute-irreducibility" in tower["explicit_curve_step"]
    final = payload["corrected_extension_class"]["final_Hilbert_gate"]
    assert "h0(-K)=0" in final["minus_K"]
    assert "nullity one" in final["minus_2K"]
    assert payload["corrected_sufficient_theorem"]["not_necessary"] is True
    assert payload["all_checks_pass"] is True


def test_finite_field_statement_is_geometric_not_false_Fp_rationality() -> None:
    warning = certificate()["corrected_extension_class"]["field_semantics"][
        "finite_field_warning"
    ]
    assert "need not contain an F_p-rational point" in warning
    assert "over the algebraic closure" in warning
    assert "some finite extension" in warning
    assert "explicit rational coefficient vector" in warning


def test_type_vii_is_a_fifth_core_not_a_nonzero_extension() -> None:
    fifth = certificate()["core_calibrations"]["fifth_explicit_core"]
    assert fifth["raw_polytope"] == "H10/72"
    assert fifth["role"] == "FIFTH_ZERO_STEP_CORE_CALIBRATION_NOT_A_NONZERO_EXTENSION"
    assert fifth["minus_K_kernel_dimension"] == 0
    assert fifth["minus_2K_kernel_dimension"] == 1
    assert fifth["strong_Coble"] == "PASS"
    assert fifth["same_marked_polytope"] is True


def test_all_154_mobius_plane_frames_have_no_clean_one_centre_core() -> None:
    scan = certificate()["H11_206_direct_extension_scan"]
    assert scan["full_P2_frame_denominator"] == 154
    assert scan["desargues_profile_frames"] == 144
    assert scan["desargues_frames_with_a_common_avoided_eleventh_centre"] == 0
    assert scan["desargues_clean_one_centre_extension_witnesses"] == 0
    assert scan["type_vii_profile_frames"] == 0
    assert scan["verdict"] == "ZERO_CLEAN_ONE_CENTRE_EXTENSION_WITNESSES"


def test_nine_row_census_is_explicitly_a_failed_scoped_search() -> None:
    census = certificate()["nine_row_extension_census"]
    assert census["row_count"] == 9
    assert len(census["rows"]) == 9
    assert census["zero_step_core_calibrations"] == [
        "H10/78",
        "H10/72",
        "H11/206",
        "H14/60",
        "H17/48",
    ]
    assert census["nonzero_extension_witness_count"] == 0
    assert census["search_verdict"] == "FAILED_NONZERO_EXTENSION_SEARCH"
    assert "neither abstract Coxeter-polytope nonexistence" in census["scope_warning"]


def test_kondo_pushback_restricts_every_accepted_raw_P() -> None:
    pushback = certificate()["AG_pushback"]
    consequence = pushback["Kondo_classification_consequence"]
    assert consequence["all_characteristics_not_2"]["boundary_counts"] == [1, 2, 5, 8]
    assert consequence["all_characteristics_not_2"]["hyperbolic_dimensions"] == [
        10,
        11,
        14,
        17,
    ]
    assert consequence["characteristic_5_dimensions"] == [10, 11]
    lengths = pushback["extension_length_restrictions"]
    assert lengths["mobius_S_3_4_10_char3_base_n2"] == [0, 3, 6]
    assert lengths["polar_chord_char3_base_n5"] == [0, 3]
    assert lengths["grid_projectivity_char3_base_n8"] == [0]
    assert lengths["desargues_triangle_cut_char5_base_n1"] == [0, 1]
    assert lengths["type_vii_tetrahedral_char5_base_n1"] == [0, 1]
    assert pushback["compact_scope"]["verdict"] == "NO_COMPACT_OUTPUT_IN_THIS_COBLE_CLASS"


def test_frozen_archive_and_cli_replay() -> None:
    frozen = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert certificate() == frozen
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--check-output", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0


def test_rating_is_not_semantically_inflated() -> None:
    rating = certificate()["rating"]
    assert rating["overall"] == "L3_PLUS_CORRECTED_SUFFICIENT_EXTENSION_THEOREM; NOT_L4"
    assert rating["no_TOOT"] is True
    assert rating["why_not_L4"].startswith("Primary defect: existence")
    assert "zero genuinely nonzero extensions" in rating["why_not_L4"]
