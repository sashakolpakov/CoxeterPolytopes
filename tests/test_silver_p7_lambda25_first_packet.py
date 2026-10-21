import json
from pathlib import Path

from sympy import Matrix

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_lambda25 import (
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_WALL_COUNT,
    FULL_MINIMAL_SHELL_COUNT,
    HEIGHT_ZERO_VECTOR,
    P7_GENERATOR,
    P7_LAMBDA25_COEFFICIENT,
    P7_LAMBDA25_NORM_CLASSES,
    PAIRING_RADICAL_MINIMAL_COUNT,
    PAIRING_RADICAL_MINIMUM,
    enumerate_p7_lambda25_controls,
    first_packet_pairing_from_dot,
    first_packet_reflection_integral,
    first_packet_root,
    gated_gram_triples,
    height_zero_gram_triples,
    height_zero_vector_audit,
    lambda25_catalogue_audit,
    lambda25_pair,
    lambda25_pairing_radical_basis,
    lambda25_smith_invariants,
    p7_lambda25_lorentz_pair,
    strict_hyperideal_witness,
    strict_witness_first_packet_pairing,
    strict_witness_height_zero_pairing,
    strict_witness_norm,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_p7_lambda25_first_packet.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_catalogue_smith_form_and_pairing_radicals() -> None:
    assert lambda25_catalogue_audit()
    assert lambda25_smith_invariants() == (1, 1) + (2,) * 22 + (16,)
    assert PAIRING_RADICAL_MINIMUM == {1: 8, 2: 8, 4: 16, 8: 16, 16: 64}
    assert PAIRING_RADICAL_MINIMAL_COUNT == {
        1: 196_656,
        2: 98_256,
        4: 2,
        8: 2,
        16: 2,
    }
    indices = {
        modulus: abs(int(Matrix(lambda25_pairing_radical_basis(modulus)).det()))
        for modulus in (1, 2, 4, 8, 16)
    }
    assert indices == {1: 1, 2: 2**2, 4: 2**26, 8: 2**50, 16: 2**74}


def test_complete_twenty_class_norm_alphabet() -> None:
    assert P7_LAMBDA25_COEFFICIENT == Silver(-1, 2)
    assert P7_GENERATOR == Silver(3, 1)
    assert len(P7_LAMBDA25_NORM_CLASSES) == 20
    assert len({row.norm for row in P7_LAMBDA25_NORM_CLASSES}) == 20
    assert {row.dyadic_valuation for row in P7_LAMBDA25_NORM_CLASSES} == {
        0, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    }
    assert {row.p7_gate for row in P7_LAMBDA25_NORM_CLASSES} == {False, True}


def test_unique_height_zero_wall() -> None:
    assert height_zero_vector_audit()
    surviving = [
        (row, height_zero_gram_triples(row))
        for row in P7_LAMBDA25_NORM_CLASSES
        if height_zero_gram_triples(row)
    ]
    assert len(surviving) == 1
    assert surviving[0][0].norm == 16
    assert [
        (item.u_square, item.v_square, item.u_dot_v) for item in surviving[0][1]
    ] == [(16, 0, 0)]
    assert lambda25_pair(HEIGHT_ZERO_VECTOR, HEIGHT_ZERO_VECTOR) == 16


def test_complete_control_chronology_and_first_packet() -> None:
    controls = enumerate_p7_lambda25_controls()
    assert len(controls) == 17
    assert len({row.height for row in controls}) == 9
    assert sum(row.height < FIRST_PACKET_HEIGHT for row in controls) == 14
    surviving = [row for row in controls if gated_gram_triples(row)]
    assert len(surviving) == 1
    assert surviving[0].height == FIRST_PACKET_HEIGHT
    assert surviving[0].norm_class.norm == 4
    assert surviving[0].control == Silver(2, 2)
    assert [
        (item.u_square, item.v_square, item.u_dot_v)
        for item in gated_gram_triples(surviving[0])
    ] == [(8, 8, 8)]


def test_explicit_reflecting_minimal_root_and_same_height_acuteness() -> None:
    minimal = (0, 1) + (0,) * 23
    assert lambda25_pair(minimal, minimal) == 8
    assert first_packet_reflection_integral(minimal)
    root = first_packet_root(minimal)
    assert p7_lambda25_lorentz_pair(root, root) == 4
    assert first_packet_pairing_from_dot(4) == Silver(-8, -8)
    assert first_packet_pairing_from_dot(4).sign(0) < 0


def test_strict_hyperideal_witness() -> None:
    witness = strict_hyperideal_witness()
    assert p7_lambda25_lorentz_pair(witness, witness) == Silver(17, -2)
    assert strict_witness_norm() == Silver(17, -2)
    assert strict_witness_norm().sign(0) > 0
    assert strict_witness_height_zero_pairing() == -16
    assert strict_witness_first_packet_pairing() == Silver(-6, -2)
    assert strict_witness_first_packet_pairing().sign(0) < 0


def test_frozen_artifact_records_exact_gate_and_no_go() -> None:
    data = payload()
    assert data["local_mirror_norms"]["class_count"] == 20
    assert data["height_zero"]["wall_count"] == 1
    assert data["positive_prefix"]["control_count"] == 17
    assert data["positive_prefix"]["earlier_control_count"] == 14
    packet = data["first_packet"]
    assert packet["full_minimal_shell_count"] == FULL_MINIMAL_SHELL_COUNT
    assert packet["reflecting_s2_shell_count"] == FIRST_PACKET_WALL_COUNT
    assert packet["excluded_by_reflection_gate"] == 98_400
    shell = packet["shell_audit"]
    assert shell["s2_signed_count"] == 98_256
    assert shell["s2_representatives_mod_sign"] == 49_128
    assert shell["s2_nonorthogonal_count_mod_sign"] == 0
    assert shell["s2_last_coordinate_nonzero_count_mod_sign"] == 0
    assert data["strict_hyperideal_witness"]["strict"] is True
    assert data["strict_hyperideal_witness"]["norm"]["signs"][0] > 0
    assert "reflection gate" in data["discarded_full_shell_polar"]["reason"]
    assert data["discovery"] is None
    assert "later mirrors may cut it" in data["scope"]
