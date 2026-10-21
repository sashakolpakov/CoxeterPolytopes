import json
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_s26 import (
    FIRST_PACKET_HEIGHT,
    NORM_CLASSES,
    S26_SHELL_COUNTS,
    enumerate_controls,
    first_packet_audit,
    first_packet_root_is_reflective,
    h1000_minimal_direction_audit,
    height_zero_gated_gram_triples,
    s26_catalogue_audit,
)


ROOT = Path(__file__).resolve().parents[1]
JSON_PATH = ROOT / "data/enumerations/silver_s26_h26_first_packet.json"


def test_catalogue_and_complete_local_alphabet() -> None:
    assert s26_catalogue_audit()
    assert len(NORM_CLASSES) == 6
    assert all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES)


def test_first_packet_chronology_and_hyperideal_certificate() -> None:
    audit = first_packet_audit(shell_four_count=S26_SHELL_COUNTS[4])
    assert FIRST_PACKET_HEIGHT == 5 + 4 * SQRT2
    assert audit.control_count == 12
    assert audit.distinct_height_count == 9
    assert audit.earlier_control_count == 10
    assert audit.earlier_surviving_control_count == 0
    assert audit.first_norm == 2
    assert audit.first_control == 2 + SQRT2
    assert audit.first_target == 12 + 8 * SQRT2
    assert (
        audit.surviving_triple.u_square,
        audit.surviving_triple.v_square,
        audit.surviving_triple.u_dot_v,
    ) == (4, 4, 4)
    assert audit.wall_count == 102_180
    assert audit.wall_margin == 2 - SQRT2
    assert audit.wall_margin.sign(0) > 0
    assert audit.witness_lorentz_norm == 4 - 2 * SQRT2
    assert audit.witness_lorentz_norm.sign(0) > 0
    assert audit.strict_hyperideal
    assert not audit.finite_volume and not audit.compact
    assert len(enumerate_controls()) == 12


def test_explicit_norm_four_root_and_frozen_payload() -> None:
    d = (-1, 0, 0, 0, 1) + (0,) * 21
    assert first_packet_root_is_reflective(d)
    payload = json.loads(JSON_PATH.read_text())
    assert payload["hyperbolic_dimension"] == 26
    assert payload["height_zero"]["reflecting_wall_count"] == 0
    assert payload["first_packet"]["wall_count"] == 102_180
    assert payload["geometry"]["strict_hyperideal_open_patch"] is True
    assert payload["geometry"]["finite_volume_through_prefix"] is False


def test_complete_h1000_necessary_superset_and_strict_replacement() -> None:
    audit = h1000_minimal_direction_audit()
    assert audit.control_count == 134
    assert audit.distinct_height_count == 104
    assert audit.nonempty_control_count == 30
    assert audit.gated_binary_type_count == 256
    assert audit.nonempty_height_count == 24
    worst = audit.worst_row
    assert worst.control.height == 455 + 322 * SQRT2
    assert worst.control.norm_class.norm == 2
    assert worst.control.control == 16 + 11 * SQRT2
    assert (
        worst.triple.u_square,
        worst.triple.v_square,
        worst.triple.u_dot_v,
    ) == (460, 226, 322)
    assert worst.pairing.maximum_pairing == 37 + 26 * SQRT2
    assert (
        worst.pairing.maximizing_u_pairing,
        worst.pairing.maximizing_v_pairing,
    ) == (37, 26)
    assert worst.pairing.retained_extension_count == 147
    assert worst.normalized_wall_bound == 28 + 21 * SQRT2
    assert worst.witness_margin == Silver(-51, 37) / 55
    assert worst.witness_margin.sign(0) > 0
    assert audit.witness_scale == Silver(43) / 55
    assert audit.witness_lorentz_norm == Silver(8572) / 3025 - 2 * SQRT2
    assert audit.witness_lorentz_norm.sign(0) > 0
    assert audit.all_margins_strict and audit.strict_hyperideal
    assert not audit.finite_volume and not audit.compact
