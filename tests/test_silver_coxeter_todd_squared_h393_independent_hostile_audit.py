import hashlib
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT
    / "scripts/independent_hostile_audit_silver_coxeter_todd_squared_h393_replacement.py"
)
ARTIFACT = (
    ROOT
    / "data/enumerations/silver_coxeter_todd_squared_h393_independent_hostile_audit.json"
)
REPORT = ROOT / "results/silver_coxeter_todd_squared_h393_independent_hostile_audit.md"
CHECKSUMS = (
    ROOT
    / "data/enumerations/silver_coxeter_todd_squared_h393_independent_hostile_audit.sha256"
)


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_clean_room_source_and_frozen_checksums() -> None:
    source = SCRIPT.read_text()
    assert "from coxeter_polytopes" not in source
    assert "import coxeter_polytopes" not in source
    observed = CHECKSUMS.read_text().splitlines()
    expected_paths = (SCRIPT, ARTIFACT, REPORT)
    assert len(observed) == len(expected_paths)
    for line, path in zip(observed, expected_paths, strict=True):
        digest, relative = line.split("  ", 1)
        assert relative == str(path.relative_to(ROOT))
        assert digest == hashlib.sha256(path.read_bytes()).hexdigest()


def test_independent_chronology_cutter_and_replacement() -> None:
    data = payload()
    assert data["verdict"] == "PASS"
    assert data["owner_claims_match"] is True
    chronology = data["chronology"]
    assert chronology["norm_class_count"] == 16
    assert chronology["control_count"] == 128
    assert chronology["distinct_height_count"] == 87
    assert chronology["through_first_packet"] == 50
    assert chronology["strictly_between"] == 76
    assert chronology["at_cutting_height"] == 2
    assert [row["gram_type_count"] for row in chronology["frontier"]] == [
        1, 2, 2, 1, 4, 11
    ]
    assert [row["directional_pair_count"] for row in chronology["frontier"]] == [
        5, 10, 42, 19, 140, 715
    ]

    cutter = data["explicit_cutter"]
    assert cutter["difference_gram"] == [28, 84, 46]
    assert cutter["coefficient_gram"] == [380, 204, 278]
    assert cutter["h_pairings"] == [22, 16]
    assert cutter["norm"]["display"] == "2"
    assert cutter["reflection_integral"] is True
    assert cutter["primitive"] is True
    assert cutter["accepted"] is True
    assert cutter["first_packet_root_pairing_maxima"][0]["display"] == (
        "-18 - 13*sqrt(2)"
    )
    assert cutter["first_packet_root_pairing_maxima"][1]["display"] == (
        "-86 - 60*sqrt(2)"
    )
    assert [row["root_pairing"]["display"] for row in cutter["intermediate_cross_gram"]] == [
        "0", "0", "-7 - 4*sqrt(2)"
    ]
    assert cutter["same_height_q2"]["positive_integral_candidates"] == [
        "1", "1*sqrt(2)", "2"
    ]
    assert cutter["same_height_q2"]["distinct_positive_pairing_impossible"] is True

    replacement = data["replacement"]
    assert replacement["scale"] == "7/4"
    assert replacement["norm"]["display"] == "205/12 - 12*sqrt(2)"
    assert replacement["norm"]["signs"][0] == 1
    assert replacement["tight_margin"]["display"] == "1/2 - 1*sqrt(2)"
    assert replacement["tight_margin"]["signs"][0] == -1
    assert replacement["all_necessary_margins_strict"] is True


def test_clean_room_regeneration() -> None:
    completed = subprocess.run(
        [sys.executable, str(SCRIPT), "--verify"],
        cwd=ROOT,
        text=True,
        capture_output=True,
        timeout=300,
        check=True,
    )
    assert "PASS independent hostile silver K12^2 H393 audit" in completed.stdout
