import json
from pathlib import Path

import pytest

from coxeter_polytopes.combinatorics.h9_14_record17 import (
    EXPECTED_RESIDUAL_ORBITS,
    EXPECTED_STAGE_COUNTS,
    TARGET_GRAPH6,
    audit_record17,
)
from scripts.certify_h9_14_record17 import (
    EXPECTED_CERTIFICATE_SHA256,
    OUTPUT,
    X_ARCHIVE,
    verify,
)
from scripts.independent_h9_14_record17_audit import (
    EXPECTED_SUMMARY_SHA256,
    audit as independent_audit,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data/enumerations/h9_14_nerve_support.json"


def _record() -> dict[str, object]:
    source = json.loads(SOURCE.read_text(encoding="utf-8"))
    survivors = [record for record in source["records"] if record["survives_second_layer"]]
    record = survivors[17]
    assert record["graph6"] == TARGET_GRAPH6
    return record


@pytest.fixture(scope="module")
def proof_audit():
    return audit_record17(_record()["dehn_sommerville_schedule"], X_ARCHIVE)


def test_record17_whole_shadow_census_is_exact(proof_audit) -> None:
    assert proof_audit.verified is True
    enumeration = proof_audit.enumeration
    assert enumeration.coupled_state_count == 11_545
    assert enumeration.stage_counts == EXPECTED_STAGE_COUNTS
    assert len(enumeration.witnesses) == 21
    assert len(enumeration.residual_orbits) == 8
    observed = tuple(
        (
            orbit.representative.missing_four_count,
            orbit.representative.x_graph6,
            orbit.representative.p_local_mask,
            orbit.representation_count,
        )
        for orbit in enumeration.residual_orbits
    )
    assert observed == EXPECTED_RESIDUAL_ORBITS


def test_record17_hostile_replay_checks_every_residual_link() -> None:
    hostile = independent_audit()
    assert hostile["summary_sha256"] == EXPECTED_SUMMARY_SHA256
    assert hostile["residual"]["labelled_representation_count"] == 21
    assert hostile["residual"]["orbit_count"] == 8
    assert hostile["residual"]["orbit_branch_histogram"] == {
        "0": 1,
        "1": 4,
        "2": 3,
    }
    assert all(
        orbit["ridge_degrees_two"]
        and orbit["dual_connected"]
        and orbit["all_links_mod2_homology_spheres"]
        for orbit in hostile["residual"]["orbits"]
    )


def test_frozen_record17_certificate_reproduces() -> None:
    manifest = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert manifest["certificate_sha256"] == EXPECTED_CERTIFICATE_SHA256
    assert manifest["correction_note"].startswith(
        "An exploratory intermediate report incorrectly called the a=0 state unique"
    )
    assert verify(manifest) == {
        "record_index": 17,
        "graph6": TARGET_GRAPH6,
        "theorem_status": "certified_exact_eight_orbit_residual",
        "labelled_representations": 21,
        "incidence_complex_orbits": 8,
        "branch_histogram": {"0": 1, "1": 4, "2": 3},
        "certificate_sha256": EXPECTED_CERTIFICATE_SHA256,
    }
