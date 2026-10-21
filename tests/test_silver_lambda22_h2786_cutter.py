import json
import runpy
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/audit_silver_lambda22_h2786_cutter.py"
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h2786_cutter.json"


def test_frozen_h2786_cutter_certificate() -> None:
    payload = json.loads(ARTIFACT.read_text())
    cutter = payload["explicit_cutter"]
    chronology = payload["chronological_acceptance"]
    assert cutter["height"]["display"] == "2786 + 1970*sqrt(2)"
    assert cutter["norm"]["display"] == "2"
    assert cutter["primitive"] and cutter["reflection_integral"]
    assert cutter["old_10d_over_31_ray_pairing"]["signs"][0] == 1
    assert [
        chronology[key]
        for key in (
            "control_count",
            "post_H14_control_count",
            "gate_compatible_post_H14_type_count",
            "nonnegative_row_count",
            "nonnegative_extension_count",
            "positive_extension_count",
        )
    ] == [374, 349, 1997, 23, 26, 8]
    exclusion = chronology["H2404_exact_exclusion"]
    assert exclusion["finite_candidate_count"] == 13_824
    assert exclusion["solution_count"] == 0
    assert "chronologically accepted" in chronology["conclusion"]
    replacement = payload["replacement_hyperideal_ray"]
    assert replacement["scale"] == "28/87"
    assert replacement["norm"]["signs"][0] == 1
    assert replacement["retained_necessary_type_count"] == 2001
    assert replacement["relaxed_nonnegative_row_count"] == 5
    assert replacement["actual_H14_q4_margin"]["signs"][0] == -1
    assert all(
        row["signs"][0] == -1
        for row in replacement["accepted_wall_pairings"].values()
    )


def test_h2786_slice_partition_and_exact_arithmetic() -> None:
    module = runpy.run_path(str(SCRIPT))
    slices = module["PINNED_SLICES"]
    assert slices[0][:2] == (0, 276)
    assert slices[-1][:2] == (370, 374)
    assert all(slices[i][1] == slices[i + 1][0] for i in range(len(slices) - 1))
    assert [sum(row[index] for row in slices) for index in range(2, 6)] == [
        1997,
        23,
        26,
        8,
    ]
    root = module["explicit_root_audit"]()
    assert root["coefficient_gram"] == [[2818, 1970], [1970, 1378]]
    acuteness = module["accepted_wall_acuteness"]()
    assert acuteness["all_strict"]
    replacement = module["replacement_ray_audit"]()
    assert replacement["scale"] == "28/87"
    assert all(row["resolution"] for row in replacement["relaxed_nonnegative_rows"])


def test_h2404_coset_exhaustion_replays() -> None:
    module = runpy.run_path(str(SCRIPT))
    audit = module["h2404_embedding_exclusion"]()
    assert audit["orthogonal_kernel_rank"] == 20
    assert audit["orthogonal_kernel_gram_determinant"] == 432
    assert audit["finite_candidate_count"] == 13_824
    assert audit["solution_count"] == 0
