from __future__ import annotations

import hashlib
import json
from pathlib import Path

from scripts.certify_borcherds_h21_outer_routing_counterexample import (
    OUTPUT,
    ROOT,
    SOURCE,
    SOURCE_SHA256,
    build_certificate,
)


REPORT = ROOT / "results/borcherds_h21_outer_routing_counterexample.md"
SCRIPT = ROOT / "scripts/certify_borcherds_h21_outer_routing_counterexample.py"
EXPECTED_CERTIFICATE_SHA256 = (
    "6f60b32bd1a1d067661e990715e60f192f56867e32dfb361c1019013e32c89ed"
)
EXPECTED_FILE_SHA256 = {
    "source": "dc9053322bf54cedde15fcfc5080666ebd79e64feb065b071880dce6ff1552ca",
    "script": "c4e21ce23ea67ecd028437caf181adf51460872aff8c891fab891fcf60424523",
    "archive": "89dd9af466dcd879fd9764787b7fb6da99bd5222d4a0b17b554af9d6896313a9",
    "report": "fa899977f418d6fa93f7d5270e067c28a648c5e55f2636cb2b362e381534821d",
}


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_certificate_passes_and_is_frozen() -> None:
    payload = build_certificate()
    assert payload["all_checks_pass"]
    assert payload["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256


def test_archived_payload_is_exact_replay() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == build_certificate()


def test_all_file_pins() -> None:
    assert _sha256(SOURCE) == SOURCE_SHA256 == EXPECTED_FILE_SHA256["source"]
    assert _sha256(SCRIPT) == EXPECTED_FILE_SHA256["script"]
    assert _sha256(OUTPUT) == EXPECTED_FILE_SHA256["archive"]
    assert _sha256(REPORT) == EXPECTED_FILE_SHA256["report"]


def test_raw_rank_norm_and_divisibility_split() -> None:
    lattice = build_certificate()["raw_lattice"]
    assert lattice["identification"] == "U + D20(-1)"
    assert lattice["rank"] == 22
    assert lattice["signature"] == [1, 21]
    assert lattice["determinant"] == -4
    assert lattice["short_wall_count"] == 42
    assert lattice["long_wall_count"] == 168
    assert lattice["short_wall_rank"] == 22
    assert lattice["long_wall_rank"] == 21
    assert lattice["short_wall_divisibilities"] == {"1": 42}
    assert lattice["long_wall_divisibilities"] == {"2": 168}


def test_short_root_frame_and_discriminant_action_prove_special_maximality() -> None:
    payload = build_certificate()
    lattice = payload["raw_lattice"]
    discriminant = lattice["discriminant_module"]
    assert lattice["short_wall_coordinate_smith_diagonal"] == [1] * 22
    assert discriminant["structure"] == "(Z/2)^2"
    assert discriminant["generated_group_order"] == 6
    assert discriminant["generated_group"] == "GL(2,F2)=S3"
    assert discriminant["all_subgroup_sizes"] == [1, 2, 2, 2, 4]
    assert discriminant["W_stable_subgroup_sizes"] == [1, 4]
    for key in (
        "short_roots_generate_L_integrally",
        "short_odd_coxeter_graph_connected",
        "wall_reflections_induce_full_GL2_F2",
        "only_trivial_and_full_discriminant_subgroups_are_W_stable",
    ):
        assert payload["checks"][key]


def test_maximality_does_not_depend_on_esselmann_classification() -> None:
    maximality = build_certificate()["reflection_maximality"]
    proof = " ".join(
        maximality["level_3_classification_free_absolute_overgroup_argument"]
    )
    assert "L is contained in M" in proof
    assert "M is contained in L*" in proof
    assert "GL(2,F2)=S3" in proof
    assert "Gamma is contained in R(M)=R(L)" in proof
    assert "not used" in maximality["esselmann_context_not_dependency"]["statement"]
    assert "absolutely maximal" in maximality["conclusion"]


def test_common_physical_square_classes_are_exhausted() -> None:
    exhaustion = build_certificate()["rational_physical_square_class_exhaustion"]
    assert exhaustion["common_square_classes"] == [1, 2]
    table = exhaustion["allowed_physical_magnitudes_by_raw_norm_and_square_class"]
    assert table["2"] == {"1": [2], "2": [1, 4]}
    assert table["4"] == {"1": [1, 4], "2": [2]}
    assert "K=0 contradicts" in exhaustion["square_class_1"]["conclusion"]


def test_adjunction_leaves_exactly_two_dual_colourings() -> None:
    branch = build_certificate()["rational_physical_square_class_exhaustion"][
        "square_class_2"
    ]
    assert branch["raw_w_square"] == -6
    assert branch["physical_K_square_formula"] == "K^2=-12*t^2"
    assert branch["rational_surface_K_square"] == -12
    assert branch["canonical_parameters"] == [-1, 1]
    colourings = branch["only_adjunction_colourings"]
    assert len(colourings) == 2
    assert {entry["canonical_parameter"] for entry in colourings} == {-1, 1}
    assert all(
        len(entry["exceptional_minus_1_short_walls_one_based"]) == 21
        and len(entry["boundary_minus_4_short_walls_one_based"]) == 21
        and entry["long_walls_all_minus_2"]
        for entry in colourings
    )


def test_projective_plane_design_and_contraction() -> None:
    obstruction = build_certificate()["projective_plane_obstruction"]
    assert obstruction["design_parameters"] == {
        "points": 21,
        "lines": 21,
        "points_per_line": 5,
        "lines_per_point": 5,
        "incidences": 105,
        "common_line_for_two_points": 1,
        "common_point_for_two_lines": 1,
        "identification": "the point-line plane PG(2,4)",
    }
    contraction = obstruction["contraction"]
    assert contraction["initial_picard_rank"] == 22
    assert contraction["disjoint_minus_1_curves"] == 21
    assert contraction["contracted_picard_rank"] == 1
    assert contraction["rank_one_minimal_rational_surface"] == "P^2"
    assert contraction["image_square"] == 1


def test_literal_fano_subplane_is_not_representable_over_C() -> None:
    fano = build_certificate()["projective_plane_obstruction"]["fano_subplane"]
    assert fano["point_walls_one_based"] == [2, 4, 6, 8, 12, 22, 33]
    assert fano["line_walls_one_based"] == [5, 7, 31, 34, 40, 47, 51]
    assert fano["standard_frame_one_based"] == {
        "A": 2,
        "B": 4,
        "C": 8,
        "D": 6,
        "E": 22,
        "F": 12,
        "G": 33,
    }
    assert fano["determinant"] == -2
    assert "only in characteristic 2, not over C" in fano["conclusion"]


def test_complex_k3_and_direct_lanes_are_excluded_with_literal_scope() -> None:
    obstruction = build_certificate()["complex_k3_and_direct_obstructions"]
    assert obstruction["raw_eigenspace_rank_required"] == 22
    assert obstruction["complex_k3_picard_rank_maximum"] == 20
    assert "rank(raw image) <= rho(X) <= 20" == obstruction["projected_rank_inequality"]
    assert obstruction["direct_complex_rational_bound_imported"] == "dimension <= 11"


def test_characteristic_two_mixed_carrier_is_preserved() -> None:
    escape = build_certificate()["positive_characteristic_escape_clause"]
    assert "supersingular K3" in escape["surface"]
    assert "42 actual smooth rational (-2)-curve classes" == escape["short_orbit"]
    assert "Cremona involutions" in escape["long_orbit"]
    assert "not every possible AG carrier category" in escape["meaning"]


def test_report_states_exact_outer_verdict_and_escape_clause() -> None:
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    required = (
        "absolutely reflection-maximal",
        "classification-free absolute maximality",
        "RLogComp",
        "Fano-plane incidence",
        "full-trace or projected-root **complex K3**",
        "mixed divisor/action",
        "Esselmann's global reflective-lattice classification is not needed",
        "cannot, by itself, prove the desired global",
    )
    for token in required:
        assert token in text


def test_report_avoids_universal_polytope_or_AG_nonexistence_claim() -> None:
    text = REPORT.read_text(encoding="utf-8")
    forbidden = (
        "P_Bor does not exist",
        "all algebraic realizations are impossible",
        "the global n<=25 theorem is false",
        "Esselmann's theorem is false",
    )
    for token in forbidden:
        assert token not in text
