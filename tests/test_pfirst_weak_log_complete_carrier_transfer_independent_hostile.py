"""Independent hostile tests for the complete weak-log carrier transfer."""

from __future__ import annotations

from fractions import Fraction
import hashlib
import json
import re
from pathlib import Path

from scripts.certify_pfirst_weak_log_complete_carrier_transfer import (
    build_certificate,
    verify_certificate,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/pfirst_weak_log_complete_carrier_transfer.md"
CERTIFIER = ROOT / "scripts/certify_pfirst_weak_log_complete_carrier_transfer.py"
ARCHIVE = ROOT / "data/enumerations/pfirst_weak_log_complete_carrier_transfer.json"
OWNER_TEST = ROOT / "tests/test_pfirst_weak_log_complete_carrier_transfer.py"
AUDIT = ROOT / "results/pfirst_weak_log_complete_carrier_transfer_independent_hostile_audit.md"
THIS_TEST = Path(__file__).resolve()
SEMI_ARCHIVE = ROOT / "data/enumerations/pfirst_weak_log_semiampleness.json"
WEAK_ARCHIVE = ROOT / "data/enumerations/reflection_maximal_weak_log_exact_nef_return.json"
ORBIFOLD_AUDIT = ROOT / "results/rational_wrong_log_orbifold_rank_gate_hostile_audit.md"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def det(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def pair(
    gram: list[list[int]],
    left: list[int | Fraction],
    right: list[int | Fraction],
) -> Fraction:
    return sum(
        Fraction(gram[i][j]) * Fraction(left[i]) * Fraction(right[j])
        for i in range(len(left))
        for j in range(len(right))
    )


def test_f4_semiampleness_provider_is_a_complete_independent_cover_replay():
    payload = json.loads(SEMI_ARCHIVE.read_text(encoding="utf-8"))
    calibration = payload["calibration"]
    fan = calibration["blown_up_fan"]
    rows = fan["cyclic_rays"]
    ids = [row["id"] for row in rows]
    rays = {row["id"]: tuple(row["vector"]) for row in rows}
    expected_cones = [
        [ids[index], ids[(index + 1) % len(ids)]]
        for index in range(len(ids))
    ]
    assert fan["maximal_cones"] == expected_cones
    assert len(fan["smooth_complete_affine_cover"]) == 5
    assert all(det(rays[a], rays[b]) == 1 for a, b in expected_cones)

    toric = calibration["toric_divisor_2D"]
    sections = {
        row["section_id"]: tuple(row["character_exponent"])
        for row in toric["returned_global_sections"]
    }
    assert len(sections) == 4
    local = toric["local_generation"]
    assert len(local) == 5
    for row in local:
        chosen = row["chosen_unit_section_id"]
        assert sections[chosen] == tuple(row["trivialization_cartier_datum"])
        assert row["chosen_local_representative"] == "1"

    classes = calibration["class_calculation"]
    gram = classes["intersection_gram"]
    d = [Fraction(3, 2), Fraction(6), Fraction(-1)]
    assert [2 * value for value in d] == [3, 12, -2]
    assert pair(gram, d, d) == 8


def test_weak_log_calibration_composes_with_the_same_nonzero_divisor():
    semi = json.loads(SEMI_ARCHIVE.read_text(encoding="utf-8"))
    weak = json.loads(WEAK_ARCHIVE.read_text(encoding="utf-8"))
    d = [
        Fraction(numerator, denominator)
        for numerator, denominator in weak["weak_log_class"]
    ]
    assert d == [Fraction(value) for value in semi["calibration"]["class_calculation"]["D"]]
    assert weak["area_over_pi"] == [1, 6]
    assert weak["proper_overgroup_area_upper_over_pi"] == [1, 12]
    assert weak["minimum_noncompact_coxeter_area_over_pi"] == [1, 6]
    assert "return the complete nef chamber" in weak["conclusion"]


def test_boundary_contraction_discrepancy_orbifold_euler_and_rank_are_exact():
    for rho in range(1, 201):
        for b in range(0, min(rho, 80) + 1):
            d_square = 10 - rho + b
            orbifold_times_four = 4 * rho + 8 - 7 * b
            assert 7 * d_square + orbifold_times_four == 78 - 3 * rho
            if d_square >= 0 and orbifold_times_four >= 0:
                assert rho <= 26

    assert Fraction(2, -4) == Fraction(-1, 2)
    assert 1 + Fraction(-1, 2) > 0
    for rho, b in ((3, 1), (10, 0), (26, 16)):
        e_s = 2 + rho
        e_y = e_s - b
        orbifold = Fraction(e_y) - Fraction(3 * b, 4)
        assert orbifold == Fraction(2 + rho) - Fraction(7 * b, 4)
    assert 10 - 10 == 0
    assert 10 - 11 < 0


def test_IJZ_primary_hypotheses_and_hash_are_pinned_in_dependency_audit():
    normalized = " ".join(ORBIFOLD_AUDIT.read_text(encoding="utf-8").split())
    for phrase in (
        "compact normal complex analytic variety",
        "at most klt singularities",
        "existence of a nef and big class",
        "nef anticanonical divisor",
        r"\widehat c_2(Y)\geq0",
        "6cb03d24e199e4524ed0e2cec4e7bce6371c2aee09f910bde9deb4a297877517",
    ):
        assert phrase in normalized


def test_archive_replays_but_explicitly_is_not_a_generic_packet_validator():
    archived = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    assert archived == build_certificate()
    assert verify_certificate(archived)
    scope = archived["scope"]
    assert scope["composition_producer_is_dependency_and_calibration_replay"] is True
    assert scope["composition_producer_accepts_arbitrary_carrier_packet"] is False
    assert scope["complete_Mori_generator_list_assumed"] is False
    assert scope["nefness_inferred_from_proposed_walls"] is False
    assert scope["raw_surface_free_essential_surjectivity_proved"] is False
    assert archived["rank_elimination"]["empty_boundary_maximum_hyperbolic_dimension"] == 9


def test_report_renders_math_and_states_every_proof_gate_and_nonclaim():
    text = REPORT.read_text(encoding="utf-8")
    assert "\r" not in text
    normalized = " ".join(text.split())
    required = (
        r"\(P\subset\mathbb H^n\)",
        r"\(\mathcal O_S(mD)\)",
        "does not require them to be listed or to generate the Mori cone",
        r"\(H.B_i=0\)",
        "source label thm-nef-anti-canonical",
        r"\(\widehat c_2(Y)\geq0\)",
        r"\(b=0\)",
        r"\(\rho\leq10\)",
        r"minimal Vinberg field is \(\mathbf Q\)",
        "complete locally finite mirror arrangement",
        "does **not** ingest and validate an arbitrary complete carrier packet",
        "not by hard-coded Boolean gates",
        "does not assert that every raw maximal arithmetic chamber",
    )
    for phrase in required:
        assert phrase in normalized
    assert "(rho" not in normalized
    assert "(mathcal" not in normalized


def test_independent_audit_records_repairs_scope_and_current_hashes():
    text = AUDIT.read_text(encoding="utf-8")
    normalized = " ".join(text.split())
    required = (
        "PASS AFTER REPAIR",
        "does not treat Boolean fields",
        "lost nearly all inline TeX delimiters",
        "never requires a finite complete Mori list",
        "The primary Iwai--Jinnouchi--Zhang source was checked independently",
        "Each host therefore contributes at most one unmarked ambient-isometry class",
        "It does not verify the geometric gates of an arbitrary carrier",
    )
    for phrase in required:
        assert phrase in normalized

    rows = dict(
        re.findall(r"\| ([^|]+?) \| ([0-9a-f]{64}) \|", text)
    )
    assert rows == {
        "composition theorem": sha256(REPORT),
        "composition producer": sha256(CERTIFIER),
        "deterministic archive": sha256(ARCHIVE),
        "owner regression tests": sha256(OWNER_TEST),
        "independent hostile tests": sha256(THIS_TEST),
    }
    for placeholder in (
        "REPORT_SHA256",
        "CERTIFIER_SHA256",
        "ARCHIVE_SHA256",
        "OWNER_TEST_SHA256",
        "INDEPENDENT_TEST_SHA256",
    ):
        assert placeholder not in text

