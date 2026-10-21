import json

from scripts.certify_ag_to_absw_fixed_data_outer_routing_synthesis import (
    OUTPUT,
    SOURCES,
    build_certificate,
)


def test_archive_replays_exactly() -> None:
    result = build_certificate()
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == result
    assert result["all_checks_pass"] is True
    assert all(result["checks"].values())
    assert result["checks"]["field_scope_has_no_local_shorthand_regression"] is True


def test_all_six_frozen_inputs_are_pinned() -> None:
    result = build_certificate()
    assert set(result["source_sha256"]) == set(SOURCES)
    assert len(result["source_sha256"]) == 6
    for name, (_path, expected) in SOURCES.items():
        assert result["source_sha256"][name] == expected


def test_scope_is_literal() -> None:
    scope = build_certificate()["scope"]
    assert scope["projected_ambient_host_finiteness_needed_for_abstract_count"] is False
    assert scope["universal_complex_projected_root_k3"] is False
    assert scope["x3_mixed_divisor_k3_positive"] is True
    assert scope["x3_maximal_only_universality_settled"] is False
    assert scope["h21_absolute_reflection_maximality"] is True
    assert scope["h21_current_complex_stack_empty"] is True
    assert scope["h21_characteristic_two_mixed_curve_cremona_positive"] is True
    assert scope["global_absw_proved"] is False
    assert scope["future_carrier"].startswith("mixed divisor/action")
