"""Tests for the finite P-first weak-log semiampleness certificate."""

from __future__ import annotations

from copy import deepcopy
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import pytest

from scripts.certify_pfirst_weak_log_semiampleness import (
    build_certificate,
    verify_certificate,
    verify_f4_calibration,
    verify_localized_unit_identity,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/pfirst_weak_log_semiampleness_certificate.md"
CERTIFIER = ROOT / "scripts/certify_pfirst_weak_log_semiampleness.py"
SCHEMA = ROOT / "schemas/pfirst-weak-log-semiampleness.schema.json"
ARCHIVE = ROOT / "data/enumerations/pfirst_weak_log_semiampleness.json"
PINS = {
    REPORT: "a96d6fd0332b3acf9e15afeda559d489436084321f8d8e523683c2e06fb73b58",
    CERTIFIER: "ec6fedeae0392a3420f9886b8a14b47200bc65f889b71dad2b522e5080b10108",
    SCHEMA: "07175a7c9d9dda066e67de91e6b1ad0867a76f97aac199fdea9a2203f35132c8",
    ARCHIVE: "fb1cf84400a2c02e709520acc5de798f276d5ebb5844075a9b9926d9b9e9b80c",
}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def archived() -> dict:
    return json.loads(ARCHIVE.read_text(encoding="utf-8"))


def det(left: tuple[int, int], right: tuple[int, int]) -> int:
    return left[0] * right[1] - left[1] * right[0]


def pair(gram: list[list[int]], left, right) -> Fraction:
    return sum(
        Fraction(gram[i][j]) * Fraction(left[i]) * Fraction(right[j])
        for i in range(3)
        for j in range(3)
    )


def test_frozen_artifacts_and_exact_replay() -> None:
    assert {path: sha256(path) for path in PINS} == PINS
    payload = archived()
    assert payload == build_certificate()
    assert verify_certificate(payload)


def test_nonconstant_localized_bezout_identity_and_mutation() -> None:
    record = archived()["checker_control"]["localized_nonconstant_bezout_identity"]
    assert record["variables"] == ["x", "y"]
    assert record["unit_identity"]["saturation_exponent"] == 1
    assert verify_localized_unit_identity(record)

    broken = deepcopy(record)
    broken["unit_identity"]["section_multipliers"][1]["terms"][0][
        "coefficient"
    ] = "2"
    with pytest.raises(ArithmeticError, match="localized Bezout"):
        verify_localized_unit_identity(broken)

    malformed = deepcopy(record)
    malformed["section_numerators"][0]["denominator_localizer_powers"] = []
    with pytest.raises(ValueError, match="denominator powers"):
        verify_localized_unit_identity(malformed)


def test_f4_star_subdivision_and_complete_smooth_cover_independently() -> None:
    calibration = archived()["calibration"]
    base_cone = calibration["base_fan_before_blowup"]["positive_section_fixed_cone"]
    inserted = calibration["star_subdivision"]["inserted_exceptional_ray"]
    assert inserted == [sum(values) for values in zip(*base_cone)] == [-1, -1]

    fan = calibration["blown_up_fan"]
    rows = fan["cyclic_rays"]
    rays = [tuple(row["vector"]) for row in rows]
    ids = [row["id"] for row in rows]
    assert fan["maximal_cones"] == [
        [ids[index], ids[(index + 1) % len(ids)]] for index in range(len(ids))
    ]
    assert [det(rays[index], rays[(index + 1) % len(rays)]) for index in range(5)] == [
        1,
        1,
        1,
        1,
        1,
    ]
    self_intersections = [
        -det(rays[index - 1], rays[(index + 1) % len(rays)])
        for index in range(len(rays))
    ]
    assert dict(zip(ids, self_intersections)) == {
        "B": -4,
        "G": -1,
        "E": -1,
        "P": 3,
        "F": 0,
    }
    assert verify_f4_calibration(calibration)


def test_f4_weak_log_class_is_nonzero_and_has_the_claimed_pairings() -> None:
    row = archived()["calibration"]["class_calculation"]
    gram = row["intersection_gram"]
    canonical = [Fraction(value) for value in row["canonical_class"]]
    boundary = [Fraction(value) for value in row["boundary_class"]]
    divisor = [-(k + b / 2) for k, b in zip(canonical, boundary)]
    assert divisor == [Fraction(3, 2), Fraction(6), Fraction(-1)]
    assert [2 * value for value in divisor] == [3, 12, -2]
    assert pair(gram, divisor, divisor) == 8
    physical = row["physical_curve_classes"]
    assert [
        pair(gram, divisor, physical[name]) for name in ("B", "E", "G=F-E")
    ] == [0, 1, Fraction(1, 2)]


def test_every_cartier_vertex_is_a_returned_global_section_and_local_unit() -> None:
    calibration = archived()["calibration"]
    ray_rows = calibration["blown_up_fan"]["cyclic_rays"]
    rays = {row["id"]: tuple(row["vector"]) for row in ray_rows}
    coefficients = {
        row["id"]: row["divisor_coefficient_in_2D"] for row in ray_rows
    }
    sections = {
        row["section_id"]: tuple(row["character_exponent"])
        for row in calibration["toric_divisor_2D"]["returned_global_sections"]
    }
    assert set(sections.values()) == {(0, -3), (0, -2), (-2, 0), (-12, 0)}
    for exponent in sections.values():
        for ray_id, ray in rays.items():
            assert exponent[0] * ray[0] + exponent[1] * ray[1] >= -coefficients[ray_id]

    local = calibration["toric_divisor_2D"]["local_generation"]
    assert len(local) == 5
    for row in local:
        chosen = row["chosen_unit_section_id"]
        assert sections[chosen] == tuple(row["trivialization_cartier_datum"])
        assert row["chosen_local_representative"] == "1"
        assert verify_localized_unit_identity(row["localized_unit_certificate"])


def test_bad_toric_section_and_bad_unit_certificate_are_rejected() -> None:
    calibration = deepcopy(archived()["calibration"])
    calibration["toric_divisor_2D"]["returned_global_sections"][0][
        "character_exponent"
    ] = [100, 100]
    with pytest.raises(ArithmeticError, match="Cartier vertex|global section"):
        verify_f4_calibration(calibration)

    calibration = deepcopy(archived()["calibration"])
    local = calibration["toric_divisor_2D"]["local_generation"][0]
    local["localized_unit_certificate"]["unit_identity"]["section_multipliers"][0][
        "terms"
    ][0]["coefficient"] = "0"
    with pytest.raises(ArithmeticError, match="localized Bezout"):
        verify_f4_calibration(calibration)


def test_scope_forbids_circular_or_universal_promotions() -> None:
    payload = archived()
    assert payload["theorem"]["no_wall_inequality_used_to_infer_nefness"] is True
    assert all(
        payload["scope"][key] is False
        for key in (
            "incidence_nonemptiness_proved",
            "facet_effectivity_or_irreducibility_proved",
            "boundary_disjointness_or_smoothness_proved",
            "reflection_maximality_proved",
            "raw_surface_free_exhaustivity_proved",
        )
    )
    mutated = deepcopy(payload)
    mutated["scope"]["raw_surface_free_exhaustivity_proved"] = True
    # Recompute only the digest so the scope gate, rather than the digest, fires.
    mutated_without_digest = deepcopy(mutated)
    mutated_without_digest.pop("content_sha256_without_this_field")
    mutated["content_sha256_without_this_field"] = hashlib.sha256(
        json.dumps(
            mutated_without_digest, sort_keys=True, separators=(",", ":")
        ).encode()
    ).hexdigest()
    with pytest.raises(Exception):
        verify_certificate(mutated)


def test_source_contract_pins_pb_atlas_and_hostile_audited_handoff() -> None:
    sources = archived()["sources"]
    for source in sources.values():
        path = ROOT / source["path"]
        assert path.exists()
        assert sha256(path) == source["sha256"]
    assert sources["weak_log_exact_return"]["sha256"] == (
        "5336e6f58377f0f5ab8220b4760d4472a578a05afd0888f3c27cde1281b221e6"
    )
    assert sources["weak_log_exact_return_hostile_audit"]["sha256"] == (
        "3c948a26b1841dc580a6fd0430c5415c47cf7a7f749a09e73340aa97f5f548c3"
    )
    adapter = archived()["pb_atlas_adapter"]
    assert adapter["adapter_contract_specified"] is True
    assert adapter["arbitrary_B1_B2a_adapter_implemented_end_to_end"] is False
    assert adapter["generic_localized_unit_checker_implemented"] is True
    assert adapter["F4_toric_cover_and_section_provider_implemented"] is True
    assert adapter["current_B2a_is_P2_only"] is True
    assert adapter["ruled_packets_require_a_ruled_chart_provider"] is True
    assert adapter["F4_calibration_uses_the_exact_toric_ruled_provider"] is True


def test_report_states_the_finite_identity_and_exact_handoff_boundary() -> None:
    raw_report = REPORT.read_text(encoding="utf-8")
    report = " ".join(raw_report.split())
    required = (
        r"g_\alpha^{N_\alpha}",
        "Nullstellensatz may be used to discover",
        "both opens, or an exact refinement",
        "Expected dimension is never substituted",
        "The current Stage-B2a engine implements",
        "does **not** read an arbitrary B1 packet",
        "versioned contract specification",
        r"L=2D=3B+12F-2E",
        "literal unit identity",
        "for that supplied carrier packet",
        "never accepted as a semiampleness certificate",
    )
    for phrase in required:
        assert phrase in report
    assert ",quad " not in raw_report
    assert ",qquad " not in raw_report
    assert " qquad" not in raw_report


def test_cli_replays_frozen_archive() -> None:
    subprocess.run(
        [sys.executable, str(CERTIFIER), "--check", str(ARCHIVE)],
        cwd=ROOT,
        check=True,
    )
