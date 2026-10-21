from __future__ import annotations

from hashlib import sha256
import json

import pytest

from scripts.audit_coble_nodot_15edge_lane_hostile import OUTPUT, build


ARCHIVE_SHA256 = "ee24a6e8915b5e949b4b32950639f893d54ee7064ec439ff46a8b812cfbfddcd"


@pytest.fixture(scope="module")
def certificate() -> dict:
    return build()


def test_hostile_archive_replays_without_production_helpers(certificate: dict) -> None:
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert sha256(OUTPUT.read_bytes()).hexdigest() == ARCHIVE_SHA256
    assert certificate == archived
    assert certificate["all_checks_pass"] is True
    assert certificate["method_independence"]["imports_production_lane"] is False
    assert certificate["method_independence"]["imports_production_helpers"] is False


def test_direct_support_and_overclosure_denominator(certificate: dict) -> None:
    assert certificate["direct_support_rank_histogram"] == {
        "10": 2,
        "11": 151,
        "12": 5_588,
    }
    records = certificate["class_records"]
    assert [record["graph6"] for record in records] == [
        "K??E@_KK_YWC",
        "K?AA@_ghECHG",
    ]
    assert [
        (record["radical_quotient_determinant"], record["even_unimodular_extension_count"])
        for record in records
    ] == [(-4, 2), (-1, 1)]


def test_two_marked_type_ii_returns_and_one_exact_rejection(certificate: dict) -> None:
    branches = [
        branch
        for record in certificate["class_records"]
        for branch in record["extension_branches"]
    ]
    finite = [branch for branch in branches if branch["classification"] == "finite_volume_noncompact"]
    rejected = [branch for branch in branches if branch["classification"] == "infinite_volume_hyperideal"]
    assert len(finite) == 2
    assert len(rejected) == 1
    for branch in finite:
        assert (
            branch["root_cone"]["finite"],
            branch["root_cone"]["ideal"],
            branch["root_cone"]["hyperideal"],
        ) == (38, 17, 0)
        returned = branch["full_Type_II_return"]
        assert returned["wall_count"] == 30
        assert (
            returned["cone"]["finite"],
            returned["cone"]["ideal"],
            returned["cone"]["hyperideal"],
        ) == (1_333, 207, 0)
        assert all(returned["checks"].values())
        assert returned["marked_isometry_determinant"] in (-1, 1)
    assert (
        rejected[0]["root_cone"]["finite"],
        rejected[0]["root_cone"]["ideal"],
        rejected[0]["root_cone"]["hyperideal"],
    ) == (18, 8, 24)
    assert rejected[0]["root_cone"]["hyperideal_ray_sha256"] == (
        "0c5fb0b79f03096f5d51b55824c62295d47fb9be751aec0932077d85a61897cc"
    )


def test_scope_remains_below_l4(certificate: dict) -> None:
    assert all(certificate["scope_audit"].values())
    assert "scoped inverse calibration" in certificate["conclusion"]
    assert "not strict L4" in certificate["conclusion"]
    assert "a new polytope" in certificate["conclusion"]
