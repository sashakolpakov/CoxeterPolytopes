import hashlib
import json
from pathlib import Path

from scripts.certify_rational_wrong_log_maximal_finiteness import (
    DEFAULT_OUTPUT,
    PINNED_SOURCES,
    binary_subspace_count,
    build_certificate,
    gaussian_binomial,
    rendered,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/rational_wrong_log_reflection_maximal_finiteness.md"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_exact_archive_replay_and_source_pins():
    certificate = build_certificate()
    assert DEFAULT_OUTPUT.read_text(encoding="utf-8") == rendered(certificate)
    assert all(certificate["checks"].values())
    for relative, expected in PINNED_SOURCES.items():
        assert sha256(ROOT / relative) == expected


def test_ag_rank_elimination_is_symbolic_not_horizon_search():
    certificate = build_certificate()
    ag = certificate["ag_elimination"]
    assert ag["derived_inequality"] == "3rho<=78"
    assert ag["maximum_picard_rank"] == 26
    assert ag["maximum_hyperbolic_dimension"] == 25


def test_binary_host_bound_is_exact_and_finite():
    assert gaussian_binomial(4, 2) == 35
    assert [binary_subspace_count(rank) for rank in range(1, 6)] == [
        2,
        5,
        16,
        67,
        374,
    ]
    certificate = build_certificate()
    assert certificate["finite_host_bound"]["coarse_total"] == (
        5515561289142524504797715176506813442521589215285286
    )


def test_adapted_host_and_reflection_scope():
    certificate = build_certificate()
    host = certificate["adapted_host"]
    assert host["interval"] == "L(2)<=M<=L(2)^*=(1/2)L"
    assert host["wall_checks"]["boundary_root_norm"] == -2
    assert host["wall_checks"]["minus_one_root_norm"] == -2
    assert host["wall_checks"]["minus_two_root_norm"] == -4
    assert host["wall_checks"]["all_facet_reflections_preserve_adapted_host"]
    scope = certificate["theorem_scope"]
    assert scope["reflection_maximality_required"]
    assert not scope["raw_essential_surjectivity_claimed"]
    assert not scope["hostile_incidence_claimed"]


def test_report_states_theorem_and_open_boundary_cleanly():
    raw = REPORT.read_bytes()
    assert b"\r" not in raw
    assert not any(byte < 32 and byte not in (9, 10) for byte in raw)
    text = " ".join(raw.decode("utf-8").split())
    for needle in (
        "weak-incidence AG/ABSW theorem",
        "\\boxed{n\\leq25}",
        "strictly includes the wrong-log locus",
        "at most one unmarked chamber isometry class for each host",
        "raw reflection-maximal arithmetic chamber is not yet proved",
        "hostile-incidence complement remain separate gates",
    ):
        assert needle in text
    assert (
        "does **not** assert that every raw reflection-maximal arithmetic "
        "Coxeter chamber has such a rational-surface carrier"
    ) in text


def test_archive_is_valid_json():
    with DEFAULT_OUTPUT.open(encoding="utf-8") as handle:
        assert json.load(handle)["schema"].endswith(".v1")
