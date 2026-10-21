import hashlib
import json
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_coxeter_todd_squared import (
    COEFFICIENT,
    FIRST_PACKET_HEIGHT,
    FIRST_PACKET_WALL_COUNT,
    K12_GRAM,
    K12_GRAM_STREAM_SHA256,
    K12_RADICAL_THREE_MINIMUM,
    K12_SHELL_COUNTS,
    K12_SMITH_INVARIANTS,
    NORM_CLASSES,
    P17_GENERATOR,
    P7_GENERATOR,
    WITNESS_SCALE,
    coefficient_factorization_audit,
    dual_direction_audit,
    enumerate_controls,
    first_packet_chronology_audit,
    first_packet_components,
    first_packet_root,
    first_packet_root_is_reflective,
    gated_gram_triples,
    height_zero_gated_gram_triples,
    k12_basis_audit,
    lorentz_pair,
    same_height_acuteness_bounds,
    strict_hyperideal_audit,
    strict_hyperideal_witness,
    strict_witness_norm,
    strict_witness_pairing_bounds,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_coxeter_todd_squared_first_packet.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_exact_new_host_and_complete_local_alphabet() -> None:
    assert k12_basis_audit()
    assert K12_SMITH_INVARIANTS == (1,) * 6 + (3,) * 6
    assert K12_RADICAL_THREE_MINIMUM == 12
    assert coefficient_factorization_audit()
    assert COEFFICIENT == Silver(-13, 12)
    assert COEFFICIENT.norm() == -119
    assert COEFFICIENT == (Silver(0, 1) - 1) * P7_GENERATOR * P17_GENERATOR
    assert len(NORM_CLASSES) == 16
    assert len({row.norm for row in NORM_CLASSES}) == 16
    assert {
        (row.dyadic, row.ternary, row.p7, row.p17) for row in NORM_CLASSES
    } == set(
        (dyadic, ternary, p7, p17)
        for dyadic in (False, True)
        for ternary in (False, True)
        for p7 in (False, True)
        for p17 in (False, True)
    )
    stream = "\n".join(" ".join(map(str, row)) for row in K12_GRAM) + "\n"
    assert hashlib.sha256(stream.encode()).hexdigest() == K12_GRAM_STREAM_SHA256


def test_all_height_zero_and_all_earlier_controls_are_empty() -> None:
    assert all(not height_zero_gated_gram_triples(row) for row in NORM_CLASSES)
    controls = enumerate_controls()
    assert len(controls) == 50
    assert len({row.height for row in controls}) == 35
    assert sum(row.height < FIRST_PACKET_HEIGHT for row in controls) == 48
    assert all(
        not gated_gram_triples(row)
        for row in controls
        if row.height < FIRST_PACKET_HEIGHT
    )
    assert first_packet_chronology_audit()


def test_complete_two_component_packet_and_counts() -> None:
    components = first_packet_components()
    assert len(components) == 2
    assert [component.norm for component in components] == [Silver(1), Silver(2)]
    assert [component.spatial_shell_norm for component in components] == [4, 8]
    assert [component.wall_count for component in components] == [1_512, 612_360]
    assert sum(component.wall_count for component in components) == FIRST_PACKET_WALL_COUNT
    assert FIRST_PACKET_WALL_COUNT == 613_872
    assert K12_SHELL_COUNTS == {4: 756, 6: 4_032, 8: 20_412}
    assert components[0].primitive_gram_triples == ((36, 16, 24),)
    assert components[1].raw_gram_triples == (
        (64, 36, 48),
        (68, 34, 48),
        (72, 32, 48),
    )
    assert components[1].imprimitive_gram_triples == ((64, 36, 48),)
    assert components[1].impossible_gram_triples == ((68, 34, 48),)
    assert components[1].primitive_gram_triples == ((72, 32, 48),)
    assert 2 * 756 + (2 * 20_412 + 756**2) == FIRST_PACKET_WALL_COUNT


def test_explicit_roots_reflect_and_packet_is_simultaneously_acute() -> None:
    q1 = first_packet_root(4, (1,) + (0,) * 23)
    q2 = first_packet_root(8, (1, 1) + (0,) * 22)
    assert lorentz_pair(q1, q1) == 1
    assert lorentz_pair(q2, q2) == 2
    assert first_packet_root_is_reflective(q1)
    assert first_packet_root_is_reflective(q2)
    bounds = same_height_acuteness_bounds()
    assert bounds == (
        ("q1-q1", Silver(-33, -24)),
        ("q1-q2", Silver(-11, -7)),
        ("q2-q2", Silver(-32, -24)),
    )
    assert all(value.sign(0) < 0 for _name, value in bounds)


def test_strict_dual_hyperideal_witness() -> None:
    assert WITNESS_SCALE == Fraction(9, 5)
    assert dual_direction_audit()
    witness = strict_hyperideal_witness()
    assert lorentz_pair(witness, witness) == strict_witness_norm()
    assert strict_witness_norm() == Silver(Fraction(433, 25), -12)
    assert strict_witness_norm().sign(0) > 0
    bounds = strict_witness_pairing_bounds()
    assert bounds == (
        ("q1 norm-four", 2, Silver(Fraction(9, 5), Fraction(-14, 5))),
        ("q2 norm-eight", 3, Silver(Fraction(-19, 5), Fraction(9, 5))),
    )
    assert all(value.sign(0) < 0 for _name, _dot, value in bounds)
    assert strict_hyperideal_audit()


def test_frozen_artifact_scope_and_independent_shell_replay() -> None:
    data = payload()
    assert data["host"] == "<-(12*sqrt(2)-13)> + K12 + K12"
    assert data["local_mirror_norms"]["class_count"] == 16
    assert data["height_zero"]["reflecting_wall_count"] == 0
    assert data["chronology"]["control_count"] == 50
    assert data["chronology"]["earlier_surviving_control_count"] == 0
    assert data["first_packet"]["wall_count"] == FIRST_PACKET_WALL_COUNT
    shell = data["spatial_lattice"]["pari_shell_audit"]
    assert shell["ordinary_K12"]["exact_signed_shell_counts"] == {
        "4": 756,
        "6": 4_032,
        "8": 20_412,
    }
    assert shell["pairing_radical_mod_3"] == {
        "basis_index": 729,
        "minimum": 12,
        "minimal_shell_signed_count": 756,
        "representatives_mod_sign": 378,
    }
    certificate = data["strict_hyperideal_prefix_certificate"]
    assert certificate["strict"] is True
    assert certificate["finite_volume"] is False
    assert certificate["compact"] is False
    assert data["discovery"] is None
    assert "later walls may cut it" in data["scope"]
