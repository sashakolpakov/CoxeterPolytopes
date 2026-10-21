import importlib.util
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/independent_hostile_audit_silver_lambda22_global_nonreflectivity.py"
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_global_nonreflectivity_independent_hostile_audit.json"


def load_audit():
    spec = importlib.util.spec_from_file_location("lambda22_global_hostile", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def payload():
    return json.loads(ARTIFACT.read_text())


def test_standalone_isometry_and_infinite_order():
    data = payload()
    assert data["verdict"] == "PASS_GLOBAL_NONREFLECTIVITY"
    isometry = data["isometry"]
    assert isometry["full_basis_form_preservation"] is True
    assert isometry["integral_inverse"] is True
    assert isometry["plane_determinant"] == {"a": "1", "b": "0"}
    assert isometry["infinite_order"] is True


def test_exact_separator_bound_and_positive_witness_pairing():
    data = payload()["separator"]
    assert data["cauchy_control_square"] == {"a": "136", "b": "96"}
    assert data["control_bound"] == {"a": "8", "b": "6"}
    assert data["height_bound"] == {"a": "56", "b": "40"}
    assert data["height_bound_less_than_150"] is True
    assert data["z_pairing_lower_coefficient"] == {"a": "5/3", "b": "-7/6"}
    assert data["no_separating_reflective_mirror"] is True


def test_time_zero_and_frontier_completeness():
    data = payload()
    assert data["time_zero_mirrors"]["complete_norm_class_count"] == 8
    assert data["time_zero_mirrors"]["basepoint_regular"] is True
    frontier = data["bounded_frontier_dependency"]
    assert frontier["H150_control_count"] == 58
    assert frontier["H150_retained_type_count"] == 72
    assert frontier["H150_resolved_nonnegative_types"] == 7
    scan = frontier["cleanroom_necessary_superset"]
    assert (
        scan["control_count"],
        scan["retained_level_count"],
        scan["retained_necessary_type_count"],
        scan["post_H14_nonnegative_type_count"],
    ) == (58, 18, 72, 7)
    assert scan["H41_nonnegative_extension_set"] == [[18, 12, 48]]
    assert scan["H56_124_52_80_nonnegative_extension_set"] == [[22, 14, 0]]
    assert scan["direct_matrix_determinant_regression"] is True
    assert frontier["cleanroom_exact_resolutions"]["all_seven_resolved"] is True
    assert frontier["strict_against_every_actual_root_through_150"] is True


def test_cleanroom_gated_shell_minima_and_pairing_spectra():
    frontier = payload()["bounded_frontier_dependency"]
    shells = frontier["cleanroom_gated_shells"]["shells"]
    assert shells["2"]["below_signed_count"] == 0
    assert shells["2"]["norm_12_abs_pairings"] == [0, 2, 4]
    assert shells["3"]["below_signed_count"] == 0
    assert shells["3"]["maximum"] == 24
    assert shells["6"]["below_signed_count"] == 0
    assert shells["6"]["norm_96_abs_pairings"] == [0, 6, 12]
    assert frontier["cleanroom_S2_norm_residues_mod8"] == {
        "class_count": 4,
        "residue_multiset": [0, 4, 4, 4],
    }


def test_global_logic_and_scope():
    data = payload()
    assert data["host"]["admissible_arithmetic_host"] is True
    assert data["host"]["anisotropic_and_cocompact"] is True
    assert data["chamber_and_volume_logic"]["logic_pass"] is True
    conclusion = data["conclusion"]
    assert conclusion["full_reflection_chamber_infinite_volume"] is True
    assert conclusion["full_integral_reflection_subgroup_infinite_index"] is True
    assert conclusion["silver_Lambda22_host_reflective"] is False
    assert conclusion["field_wide_nonexistence"] is False


def test_recompute_without_production_module_import():
    audit = load_audit()
    gram = audit.parse_gram()
    assert audit.independent_isometry_audit(gram)["infinite_order"] is True
    assert audit.separator_audit()["no_separating_reflective_mirror"] is True
    scan = audit.cleanroom_frontier_scan()
    assert scan["retained_necessary_type_count"] == 72
    assert scan["post_H14_nonnegative_type_count"] == 7
