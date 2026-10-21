import json
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_lambda22 import (
    FIRST_PACKET_HEIGHT,
    LAMBDA22_COEFFICIENT,
    LAMBDA22_MIRROR_NORMS,
    LAMBDA22_NORM10_SHELL_COUNT,
    S2_NORM20_SHELL_COUNT,
    cross_packet_pairing_from_dot,
    dyadic_packet_pairing_from_dot,
    dyadic_witness_pairing_from_dot,
    enumerate_lambda22_controls,
    first_packet_dyadic_reflection_integral,
    first_packet_pairing_from_dot,
    first_packet_reflection_integral,
    lambda22_basis_audit,
    lambda22_gram,
    lambda22_pair_silver,
    lambda22_pairing_radical_basis,
    lambda22_smith_invariants,
    necessary_components_through_first_packet,
    witness_norm,
    witness_pairing_from_dot,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_lambda22_h14.json"


def payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def test_exact_lambda22_source_and_discriminant_data() -> None:
    assert lambda22_basis_audit()
    assert lambda22_smith_invariants() == (1,) * 20 + (2, 6)
    gram = lambda22_gram()
    for modulus, index in ((2, 2**20), (3, 3**21), (6, 2**20 * 3**21)):
        basis = lambda22_pairing_radical_basis(modulus)
        assert abs(round(__import__("sympy").Matrix(basis).det())) == index
        assert all(
            sum(gram[i][j] * basis[j][column] for j in range(22)) % modulus == 0
            for i in range(22)
            for column in range(22)
        )


def test_complete_norm_alphabet_and_control_census() -> None:
    assert LAMBDA22_COEFFICIENT == Silver(-1, 1)
    assert LAMBDA22_COEFFICIENT.sign(0) > 0 > LAMBDA22_COEFFICIENT.sign(1)
    assert LAMBDA22_MIRROR_NORMS == (
        Silver(1), Silver(2), Silver(4, 2), Silver(4),
        Silver(3), Silver(6), Silver(12, 6), Silver(12),
    )
    controls = enumerate_lambda22_controls()
    assert len(controls) == 25
    assert max(row.height for row in controls) == FIRST_PACKET_HEIGHT
    assert all((2 * LAMBDA22_COEFFICIENT * row.control / row.norm).is_integral for row in controls)


def test_all_earlier_necessary_components_have_exact_resolutions() -> None:
    rows = [row for row in necessary_components_through_first_packet() if row.triples]
    assert len(rows) == 8
    assert [len(row.triples) for row in rows] == [1, 1, 5, 2, 2, 4, 1, 5]
    assert sum("realized" in row.exclusion for row in rows) == 2
    assert all(row.level.height == FIRST_PACKET_HEIGHT for row in rows[-2:])


def test_explicit_representatives_of_both_first_components_reflect() -> None:
    # q=2 component: an explicit norm-ten vector.
    w = (1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    assert lambda22_pair_silver(w, w) == 12  # calibration: not silently used below
    # Search a tiny deterministic coordinate box for a norm-ten representative.
    candidates = []
    for i in range(22):
        for j in range(i + 1, 22):
            for si in (-1, 1):
                for sj in (-1, 1):
                    v = [0] * 22
                    v[i], v[j] = si, sj
                    if lambda22_pair_silver(v, v) == 10:
                        candidates.append(tuple(v))
    assert candidates and first_packet_reflection_integral(candidates[0])

    # q=4 component: a nonzero radical class of norm twenty.
    u = [0] * 22
    for index in (1, 3, 4, 7, 10, 15):
        u[index] = 1
    assert lambda22_pair_silver(u, u) == 20
    assert first_packet_dyadic_reflection_integral(u)


def test_first_packet_is_mutually_acute_and_witness_is_strict() -> None:
    assert first_packet_pairing_from_dot(8).sign(0) < 0
    assert dyadic_packet_pairing_from_dot(14).sign(0) < 0
    assert cross_packet_pairing_from_dot(13).sign(0) < 0
    assert cross_packet_pairing_from_dot(14) == 2
    d = (1,) + (0,) * 21
    assert witness_norm(d) == Silver(13) / 9 - Silver(0, 1)
    assert witness_norm(d).sign(0) > 0
    assert witness_pairing_from_dot(5).sign(0) < 0
    assert dyadic_witness_pairing_from_dot(8).sign(0) < 0


def test_frozen_artifact_scope_and_pari_census() -> None:
    data = payload()
    assert data["control_level_count"] == 25
    assert data["discovery"] is None
    assert data["first_packet"]["dyadic_component"]["shell_count"] == S2_NORM20_SHELL_COUNT
    assert data["first_packet"]["ordinary_component"]["shell_count"] == LAMBDA22_NORM10_SHELL_COUNT
    assert data["first_packet"]["total_root_count"] == 468_965_376
    assert [
        (row["pairing_modulus"], row["signed_vector_count_at_or_below_bound"], row["maximum_norm_returned"])
        for row in data["gated_sublattices"]["pari"]["rows"]
    ] == [(1, 49896, 4), (2, 1426920, 20), (3, 1782, 24), (6, 1782, 96)]
    assert "later mirrors may cut it" in data["scope"]
