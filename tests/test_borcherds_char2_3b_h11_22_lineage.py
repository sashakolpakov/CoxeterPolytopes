import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
LINEAGE = (
    ROOT / "data/enumerations/borcherds_char2_3b_h11_22_lineage.json"
)


def _load() -> dict:
    return json.loads(LINEAGE.read_text(encoding="utf-8"))


def test_entry_27_lattice_and_chamber_match() -> None:
    record = _load()
    assert record["all_checks_pass"] is True
    assert record["census_match"]["entry_one_based"] == 27
    assert record["census_match"]["list_entry"] == [
        243,
        "II_{11,1}(3^{-5})",
        [22, 231, 1480, 6320, 18704, 38986,
         57182, 57785, 38140, 14580, 2304],
        72,
    ]
    assert record["lattice"]["candidate_smith_diagonal_absolute"] == [1] * 7 + [3] * 5
    assert record["lattice"]["discriminant_q_histogram"] == {
        "0": 81,
        "2/3": 72,
        "4/3": 90,
    }
    assert all(record["lattice"]["checks"].values())


def test_all_root_vinberg_stopping_certificate() -> None:
    record = _load()
    audit = record["all_root_vinberg_audit"]
    assert audit["control_square"] == 14
    assert audit["control_divisibility"] == 1
    assert audit["simple_root_pairing_histogram"] == {"1": 12, "3": 10}
    assert audit["height_zero"]["positive_kernel_determinant"] == 3402
    assert audit["height_zero"]["qfminim_total_signed_vectors_norm_at_most_6"] == 1226
    assert audit["height_zero"]["representative_pairs_by_norm"] == {
        "2": 0,
        "4": 111,
        "6": 502,
    }
    assert audit["height_zero"]["norm_6_divisibility_3_pairs"] == 0
    assert [
        (
            shell["negative_norm"],
            shell["control_pairing"],
            shell["qfminim_total_signed_vectors_through_bound"],
            shell["oriented_target_shell_candidate_count"],
            shell["reflective_target_root_count"],
        )
        for shell in audit["positive_shells"]
    ] == [(2, 1, 24, 12, 12), (6, 3, 6534, 1720, 10)]
    assert audit["shell_roots_equal_displayed_simple_roots"] is True
    assert audit["vinberg_stopping_applies"] is True
    assert audit["reflection_subgroup_index"] == 1


def test_primitive_reflective_root_classification() -> None:
    record = _load()
    classification = record["reflective_root_classification"]
    assert classification["discriminant_group_exponent"] == 3
    assert classification["formal_primitive_types_before_discriminant_filter"] == [
        [-2, 1], [-2, 3], [-6, 3]
    ]
    assert classification["excluded_type"] == [-2, 3]
    assert classification["excluded_discriminant_value_mod_2Z"] == "16/9"
    assert classification["allowed_primitive_reflective_types_square_divisibility"] == [
        [-2, 1], [-6, 3]
    ]
    assert record["decision"] == {
        "abstract_chamber_lineage": (
            "known Scharlau--Kirschmer List11/Gram11 entry 27"
        ),
        "full_reflection_chamber": True,
        "literal_k3_fixed_section_priority": (
            "not located in bounded source audit"
        ),
        "reflection_subgroup_index": 1,
    }
