"""Regression checks for the clean-room golden S23 hostile audit."""

from __future__ import annotations

import json
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT
    / "data/enumerations/golden_odd_leech_section_h23_first_packets_independent_hostile_audit.json"
)
SOURCE = ROOT / "scripts/independent_hostile_audit_golden_odd_leech_section_h23.py"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def golden_sign(value: dict[str, str]) -> int:
    a, b = Fraction(value["a"]), Fraction(value["b"])
    rational, radical = 2 * a + b, b
    sr = (rational > 0) - (rational < 0)
    sb = (radical > 0) - (radical < 0)
    if not sr:
        return sb
    if not sb or sr == sb:
        return sr
    return sr if rational * rational - 5 * radical * radical > 0 else sb


def test_standalone_section_and_radical_reconstruction() -> None:
    frozen = payload()
    model = frozen["section_model"]
    assert frozen["verdict"] == "PASS_WITH_SCHEMA_FIX"
    assert model["rank"] == 23
    assert model["determinant"] == 3
    assert model["smith_invariants"] == [1] * 22 + [3]
    assert model["gram_rank_mod_three"] == 22
    assert model["pairing_radical_dimension_mod_three"] == 1
    assert model["minimum"] == 3
    assert model["minimal_shell_count"] == 2576
    assert model["minimal_shell_span_rank"] == 23
    assert frozen["radical_minimum"]["nonzero_vectors_of_raw_square_at_most_48"] == 0


def test_full_construction_b_shell_counts_and_height_zero() -> None:
    frozen = payload()
    shells = frozen["section_shells"]
    assert shells["norm_4_by_coset"] == [53_682, 0]
    assert shells["norm_4_total"] == 53_682
    assert shells["norm_5_by_coset"] == [0, 553_104]
    assert shells["norm_5_total"] == 553_104
    assert frozen["height_zero"]["ambient_type_counts"] == [0, 0, 1, 1]
    assert frozen["height_zero"]["gated_type_counts"] == [0, 0, 0, 0]
    assert frozen["height_zero"]["wall_count"] == 0


def test_both_packets_are_complete_primitive_reflective_acute_and_noncompact() -> None:
    frozen = payload()
    specifications = (("a=phi", 7, 553_104), ("a=2phi", 2, 53_682))
    for label, control_count, root_count in specifications:
        prefix = frozen["prefixes"][label]
        assert len(prefix["complete_controls"]) == control_count
        assert sum(bool(row["necessary_gram_triples"]) for row in prefix["complete_controls"]) == 1
        assert prefix["first_packet"]["root_count"] == root_count
        assert prefix["packet_logic"]["primitive_because_control_is_unit"] is True
        assert prefix["packet_logic"]["mutually_acute"] is True
        assert prefix["packet_logic"]["root_square"] == {"a": "2", "b": "0"}
        assert golden_sign(prefix["first_packet"]["pairing_upper_bound"]) < 0
        assert golden_sign(prefix["witness_square"]) > 0
        assert golden_sign(prefix["closest_wall_pairing"]) < 0
        assert prefix["strict_hyperideal"] is True


def test_clean_room_source_and_schema_repair_are_pinned() -> None:
    source = SOURCE.read_text()
    assert "from coxeter_polytopes" not in source
    assert "def section_model" in source
    assert "def construction_b_count" in source
    assert "def enumerate_controls" in source
    defect = payload()["schema_defect"]
    assert defect["status"] == "FIXED_IN_PRODUCTION_GENERATOR"
    assert defect["replacement_keys"] == {
        "gram_rank_mod_three": 22,
        "pairing_radical_dimension_mod_three": 1,
    }
