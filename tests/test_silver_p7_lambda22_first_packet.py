import json
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_p7_lambda22 import (
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_WALL_COUNT,
    P7_GENERATOR,
    P7_LAMBDA22_COEFFICIENT,
    P7_LAMBDA22_NORM_CLASSES,
    ambient_gram_triples,
    dual_direction_audit,
    enumerate_p7_lambda22_controls,
    first_packet_pairing_from_dot,
    first_packet_reflection_integral,
    first_packet_root,
    gated_gram_triples,
    height_zero_gram_triples,
    lambda22_basis_audit,
    lambda22_pair_silver,
    lambda22_smith_invariants,
    p7_lambda22_factorization_audit,
    p7_lambda22_lorentz_pair,
    strict_hyperideal_witness,
    strict_witness_norm,
    strict_witness_pairing_from_last_coordinate,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_p7_lambda22_first_packet.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_exact_host_and_complete_sixteen_class_alphabet() -> None:
    assert lambda22_basis_audit()
    assert lambda22_smith_invariants() == (1,) * 20 + (2, 6)
    assert p7_lambda22_factorization_audit()
    assert P7_LAMBDA22_COEFFICIENT == Silver(-1, 2)
    assert P7_GENERATOR == Silver(3, 1)
    assert len(P7_LAMBDA22_NORM_CLASSES) == 16
    assert len({row.norm for row in P7_LAMBDA22_NORM_CLASSES}) == 16
    assert {row.dyadic_valuation for row in P7_LAMBDA22_NORM_CLASSES} == {
        0,
        2,
        3,
        4,
    }
    assert {(row.ternary_gate, row.p7_gate) for row in P7_LAMBDA22_NORM_CLASSES} == {
        (False, False),
        (False, True),
        (True, False),
        (True, True),
    }


def test_complete_controls_and_first_nonempty_chronology() -> None:
    controls = enumerate_p7_lambda22_controls()
    assert len(controls) == 16
    assert len({row.height for row in controls}) == 10
    assert max(row.height for row in controls) == FIRST_PACKET_HEIGHT
    assert sum(row.height < FIRST_PACKET_HEIGHT for row in controls) == 13
    survivors = [row for row in controls if gated_gram_triples(row)]
    assert len(survivors) == 1
    assert survivors[0].height == FIRST_PACKET_HEIGHT
    assert survivors[0].norm_class.norm == 2
    assert ambient_gram_triples(survivors[0]) == ((4, 4, 4),)
    assert gated_gram_triples(survivors[0]) == ((4, 4, 4),)
    simultaneous = [row for row in controls if row.height == FIRST_PACKET_HEIGHT]
    assert [row.norm_class.norm for row in simultaneous] == [Silver(1), Silver(2), Silver(4)]
    assert ambient_gram_triples(simultaneous[2]) == ((16, 4, 8), (8, 8, 8))
    assert gated_gram_triples(simultaneous[2]) == ()


def test_all_height_zero_classes_are_empty() -> None:
    assert all(
        height_zero_gram_triples(norm_class) == ()
        for norm_class in P7_LAMBDA22_NORM_CLASSES
    )


def test_explicit_first_packet_root_and_same_height_acuteness() -> None:
    d = (1,) + (0,) * 21
    assert lambda22_pair_silver(d, d) == 4
    root = first_packet_root(d)
    assert p7_lambda22_lorentz_pair(root, root) == 2
    assert first_packet_reflection_integral(d)
    assert first_packet_pairing_from_dot(2) == Silver(-4, -4)
    assert first_packet_pairing_from_dot(2).sign(0) < 0


def test_strict_hyperideal_witness_is_exact() -> None:
    assert dual_direction_audit()
    witness = strict_hyperideal_witness()
    assert p7_lambda22_lorentz_pair(witness, witness) == Silver(31) / 4 - 5 * SQRT2
    assert strict_witness_norm() == Silver(31) / 4 - 5 * SQRT2
    assert strict_witness_norm().sign(0) > 0
    assert strict_witness_pairing_from_last_coordinate(2) == Silver(-2, -3) / 4
    assert strict_witness_pairing_from_last_coordinate(2).sign(0) < 0
    assert all(
        strict_witness_pairing_from_last_coordinate(last).sign(0) < 0
        for last in range(-2, 3)
    )


def test_frozen_artifact_scope_and_shell_census() -> None:
    data = payload()
    assert data["local_mirror_norms"]["class_count"] == 16
    assert data["height_zero"]["reflecting_wall_count"] == 0
    assert data["positive_prefix"]["control_count"] == 16
    assert data["positive_prefix"]["earlier_control_count"] == 13
    packet = data["first_packet"]
    assert packet["unique_primitive_component"]["wall_count"] == FIRST_PACKET_WALL_COUNT
    assert packet["minimal_shell_audit"] == {
        "algorithm": "PARI qfminim(G,4,,flag=2)",
        "coordinate_abs_0": 13860,
        "coordinate_abs_1": 10752,
        "coordinate_abs_2": 336,
        "coordinate_abs_values": [0, 1, 2],
        "coordinate_note": (
            "Counts by |d_22| are for the one-representative-per-sign matrix; "
            "the absolute spectrum is sign invariant."
        ),
        "maximum": 4,
        "representative_stream_sha256": (
            "9b752381c953619264cf7eb1aba9b73480e6de94553a68db65eec8bfda30f428"
        ),
        "representatives_mod_sign": 24948,
        "signed_count": 49896,
    }
    assert data["strict_hyperideal_witness"]["strict"] is True
    assert data["strict_hyperideal_witness"]["norm"]["signs"][0] > 0
    assert data["strict_hyperideal_witness"]["maximum_wall_pairing"]["signs"][0] < 0
    assert data["discovery"] is None
    assert "later mirrors may cut it" in data["scope"]
