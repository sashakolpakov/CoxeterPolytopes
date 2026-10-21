import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_t26 import (
    FIRST_PACKET_CONTROL,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_NORM,
    FIRST_PACKET_WALL_COUNT,
    NORM_CLASSES,
    POLAR_BOUNDARY_NORM,
    POLAR_STRICT_FIRST_PACKET_MARGIN,
    POLAR_STRICT_NORM,
    POLAR_VERTEX,
    POLAR_VERTEX_CONTACT_COUNT,
    POLAR_VERTEX_SQUARE,
    T26_KISSING_NUMBER,
    T26_SMITH_INVARIANTS,
    enumerate_controls,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    lorentz_pair,
    polar_strict_witness,
    polar_strict_witness_audit,
    polar_vertex_audit,
    t26_catalogue_audit,
    t26_pair,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_t26_first_packet.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_catalogue_and_complete_eight_class_alphabet() -> None:
    assert t26_catalogue_audit()
    assert T26_SMITH_INVARIANTS == (1,) * 25 + (3,)
    assert len(NORM_CLASSES) == 8
    assert len({row.norm for row in NORM_CLASSES}) == 8
    assert {
        (row.dyadic, row.p7, row.discriminant_three) for row in NORM_CLASSES
    } == {
        (dyadic, p7, ternary)
        for dyadic in (False, True)
        for p7 in (False, True)
        for ternary in (False, True)
    }


def test_complete_first_packet_chronology() -> None:
    controls = enumerate_controls()
    assert len(controls) == 10
    assert len({row.height for row in controls}) == 7
    assert sum(row.height < FIRST_PACKET_HEIGHT for row in controls) == 8
    surviving = [row for row in controls if gated_gram_triples(row)]
    assert len(surviving) == 1
    first = surviving[0]
    assert first.height == FIRST_PACKET_HEIGHT
    assert first.norm_class.norm == FIRST_PACKET_NORM
    assert first.control == FIRST_PACKET_CONTROL
    assert [
        (item.u_square, item.v_square, item.u_dot_v)
        for item in gated_gram_triples(first)
    ] == [(4, 4, 4)]


def test_all_height_zero_classes_are_empty() -> None:
    assert all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES)


def test_exact_full_shell_polar_vertex() -> None:
    assert polar_vertex_audit()
    assert t26_pair(POLAR_VERTEX, POLAR_VERTEX) == POLAR_VERTEX_SQUARE
    assert POLAR_VERTEX_SQUARE == Fraction(97, 150)
    assert POLAR_VERTEX_CONTACT_COUNT == 66
    assert POLAR_BOUNDARY_NORM == Silver(316, 0) / 25 + Silver(0, -538) / 75
    assert POLAR_BOUNDARY_NORM.sign(0) > 0


def test_strict_hyperideal_witness() -> None:
    assert polar_strict_witness_audit()
    witness = polar_strict_witness()
    assert lorentz_pair(witness, witness) == POLAR_STRICT_NORM
    assert POLAR_STRICT_NORM.sign(0) > 0
    assert POLAR_STRICT_FIRST_PACKET_MARGIN == Silver(-1, 0) / 500 + Silver(0, -3) / 1000
    assert POLAR_STRICT_FIRST_PACKET_MARGIN.sign(0) < 0


def test_frozen_shell_census_and_digest() -> None:
    data = payload()
    assert data["local_mirror_norms"]["class_count"] == 8
    assert data["height_zero"]["reflecting_wall_count"] == 0
    assert data["positive_prefix"]["control_count"] == 10
    assert data["positive_prefix"]["earlier_control_count"] == 8
    packet = data["first_packet"]
    assert packet["wall_count"] == FIRST_PACKET_WALL_COUNT == T26_KISSING_NUMBER
    shell = packet["minimal_shell_audit"]
    assert shell["signed_count"] == 117_936
    assert shell["representatives_mod_sign"] == 58_968
    assert shell["polar_absolute_maximum"] == "1"
    assert shell["polar_contact_sign_pairs"] == 66
    assert shell["representative_stream_sha256"] == (
        "69cdb97ab9d5b89e24db1b1fe94b83f5488e68f14e71e65cf04a8d991b8c4870"
    )
    assert shell["radical_three_minimum"] == 24
    assert shell["radical_three_kissing_number"] == 1_638


def test_bounded_scope_and_negative_discovery() -> None:
    data = payload()
    assert data["strict_hyperideal_witness"]["strict"] is True
    assert data["strict_hyperideal_witness"]["norm"]["signs"][0] > 0
    assert (
        data["strict_hyperideal_witness"]["maximum_first_packet_pairing"]["signs"][0]
        < 0
    )
    assert data["discovery"] is None
    assert "Later mirrors may cut it" in data["scope"]
