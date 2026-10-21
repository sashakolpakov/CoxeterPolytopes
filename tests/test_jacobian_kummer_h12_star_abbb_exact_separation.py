from __future__ import annotations

import hashlib
import json
from pathlib import Path

import scripts.certify_jacobian_kummer_h12_five_vertex_p5_abbb_census as abbb
import scripts.certify_jacobian_kummer_h12_five_vertex_p5_b_le_2_census as p5
import scripts.certify_jacobian_kummer_h12_five_vertex_star_census as star
import scripts.certify_jacobian_kummer_h12_small_subtree_census as seed
import scripts.certify_jacobian_kummer_h12_star_abbb_exact_separation as cert


def _archive() -> dict:
    return json.loads(cert.OUTPUT.read_text(encoding="utf-8"))


def test_frozen_inputs_and_exact_denominators() -> None:
    payload = _archive()
    assert payload["all_checks_pass"] is True
    assert all(payload["checks"].values())
    assert payload["pinned_inputs"] == cert.PINNED_SHA256
    replay = payload["replay"]
    assert replay["star_779_row_count"] == 6_117
    assert replay["star_779_unique_exact_degree_profile_key_count"] == 5_949
    assert replay["ABBB_row_count"] == 36_247
    assert replay["cross_shape_exact_byte_match_count"] == 0
    assert replay["star_779_archive_order_exact_key_stream_sha256"] == (
        "98496abf4ea3c653eb98730be30ebf1d3af7afed44484534d65902ddac364dc1"
    )
    assert replay["ABBB_archive_order_exact_key_stream_sha256"] == (
        "9c4811d7fafb06b343982c72e99836346b2ee9e8a7bd9218cfca1477feb2381e"
    )


def test_exact_cross_shape_and_union_theorem() -> None:
    theorem = _archive()["theorem"]
    assert theorem["overlapping_facet_header"] == 779
    assert theorem["star_779_type_count"] == 6_117
    assert theorem["ABBB_type_count"] == 36_247
    assert theorem["exact_cross_shape_ordinary_Coxeter_isomorphism_count"] == 0
    assert theorem["star_779_plus_ABBB_pairwise_nonisomorphic_type_count"] == 42_364
    assert theorem["complete_star_type_count"] == 16_872
    assert theorem["complete_star_plus_ABBB_pairwise_nonisomorphic_type_count"] == 53_119
    assert theorem["reflection_index_and_volume_multiple"] == 5
    assert theorem["finite_volume"] is True
    assert theorem["compact"] is False


def test_comparison_is_exact_bytes_not_hash_equality() -> None:
    payload = _archive()
    invariant = payload["exact_invariant"]
    assert invariant["profile_shape"] == [779, 5]
    assert invariant["byte_length_per_polytope"] == 7_790
    assert invariant["cross_comparison"] == "literal Python bytes equality/membership"
    assert "not the cross-shape separation test" in invariant["hash_role"]

    source = cert.__file__
    assert source is not None
    text = Path(source).read_text(encoding="utf-8")
    assert "if degree_key in star_degree_keys:" in text
    assert "hashlib.sha256(degree_key).hexdigest()" in text


def test_selected_source_rows_replay_to_their_exact_degree_bytes() -> None:
    star_payload = json.loads(cert.STAR_ARCHIVE.read_text(encoding="utf-8"))
    abbb_payload = json.loads(cert.ABBB_ARCHIVE.read_text(encoding="utf-8"))
    _carrier, gram, roots, interior, _bounded = seed._load_inputs()
    star_engine = star.ExactStarBoundaryEngine(gram, roots, interior)
    p5_engine = p5.ExactP5BoundaryEngine(gram, roots, interior)

    star_row = next(
        row for row in star_payload["orbit_representatives"] if row["facet_count"] == 779
    )
    abbb_row = abbb_payload["orbit_representatives"][0]
    star_key = cert._rebuild_degree_key(star_engine, star_row, expected_shape="K1,4")
    abbb_key = cert._rebuild_degree_key(p5_engine, abbb_row, expected_shape="P5")
    assert len(star_key) == len(abbb_key) == 7_790
    assert hashlib.sha256(star_key).hexdigest() == star_row[
        "ordinary_degree_profile_multiset_sha256"
    ]
    assert hashlib.sha256(abbb_key).hexdigest() == abbb_row[
        "ordinary_degree_profile_multiset_sha256"
    ]
    assert star_key != abbb_key


def test_carrier_direction_and_novelty_boundary() -> None:
    boundary = _archive()["carrier_and_novelty_boundary"]
    assert boundary["direction"].startswith("SURFACE_TO_POLYTOPE")
    assert boundary["same_weak_full_Picard_K3_carrier"] == (
        "U(2)+D5(-1)^2+<-4>"
    )
    assert boundary["full_nef_effectivity_or_surface_action_claim"] is False
    assert boundary["raw_POLYTOPE_TO_SURFACE_recognizer_claim"] is False
    assert boundary["project_cross_shape_accounting_promoted"] is True
    assert boundary["worldwide_priority"] == "open"
    assert boundary["worldwide_novelty_claim"] is False


def test_report_preserves_the_project_scope() -> None:
    report = (
        cert.ROOT
        / "results/jacobian_kummer_h12_star_abbb_exact_separation.md"
    ).read_text(encoding="utf-8")
    normalized = " ".join(report.lower().split())
    assert "53,119" in report
    assert "literal-byte" in normalized
    assert "worldwide priority remains open" in normalized
    assert "not a nef" in normalized
    assert "not a `polytope -> surface`" in normalized
