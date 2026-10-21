import hashlib
from itertools import product
import json
from pathlib import Path
import shutil
import subprocess

import pytest
from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_decomp

from coxeter_polytopes.exact.cubic7 import Cubic7
from coxeter_polytopes.search.cubic7_lambda22 import (
    GRAM_FILE_SHA256,
    GRAM_PATH,
    GRAM_STREAM_SHA256,
    LAMBDA22_UNIT_SEED,
    PARI_QFMINIM_CERTIFICATES,
    SOURCE_HTML_SHA256,
    lambda22_basis_audit,
    lambda22_control_levels,
    lambda22_gram,
    lambda22_gram_file_sha256,
    lambda22_gram_stream_sha256,
    lambda22_integer_gram_minimum_at_least,
    lambda22_level_audit,
    lambda22_norm_gate,
    lambda22_pairing_radical_basis,
    lambda22_rank3_gram_types,
    lambda22_smith_invariants,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/cubic7_lambda22_h100.json"


def _payload() -> dict[str, object]:
    return json.loads(ARTIFACT.read_text())


def _element(payload: dict[str, str]) -> Cubic7:
    return Cubic7.from_json(payload)


def _gram_tuple(value) -> tuple[int, ...]:
    return (
        value.u_square,
        value.v_square,
        value.w_square,
        value.u_dot_v,
        value.u_dot_w,
        value.v_dot_w,
        value.determinant,
    )


def _direct_box_minimum(matrix: tuple[tuple[int, ...], ...], radius: int = 5) -> int:
    values = []
    for vector in product(range(-radius, radius + 1), repeat=len(matrix)):
        if not any(vector):
            continue
        square = sum(
            matrix[i][j] * vector[i] * vector[j]
            for i in range(len(vector))
            for j in range(len(vector))
        )
        if square > 0:
            values.append(square)
    return min(values)


def test_authoritative_integral_gram_and_source_pins_are_exact() -> None:
    payload = _payload()
    lattice = payload["spatial_lattice"]
    assert lambda22_basis_audit()
    assert lambda22_gram_file_sha256() == GRAM_FILE_SHA256
    assert lambda22_gram_stream_sha256() == GRAM_STREAM_SHA256
    assert hashlib.sha256(GRAM_PATH.read_bytes()).hexdigest() == GRAM_FILE_SHA256
    assert lattice["source_html_sha256"] == SOURCE_HTML_SHA256
    assert lattice["source_pins_match"]
    assert lattice["catalogue_and_exact_matrix_invariants"] == {
        "determinant": 12,
        "even": True,
        "integral": True,
        "kissing_number": 49_896,
        "minimum": 4,
        "smith_invariants": [1] * 20 + [2, 6],
    }
    assert lambda22_smith_invariants() == (1,) * 20 + (2, 6)


def test_pairing_radical_lifts_and_exact_pari_minima_are_frozen() -> None:
    payload = _payload()
    records = {
        row["pairing_modulus"]: row
        for row in payload["pari_qfminim_replay"]["radical_lattices"]
    }
    gram = Matrix(lambda22_gram())
    expected = {
        1: (1, 4, 49_896, 12),
        2: (2**20, 12, 8_448, 36),
        3: (3**21, 24, 1_782, 72),
        6: (2**20 * 3**21, 96, 1_782, 288),
    }
    for modulus, (index, minimum, count, trace_minimum) in expected.items():
        basis = Matrix(lambda22_pairing_radical_basis(modulus))
        assert abs(int(basis.det())) == index
        assert all(int(value) % modulus == 0 for value in gram * basis)
        record = records[modulus]
        assert record["index"] == index
        assert record["ordinary"]["minimum"] == minimum
        assert record["ordinary"]["kissing_number"] == count
        assert record["trace_over_OK"]["minimum"] == trace_minimum
        certificate = PARI_QFMINIM_CERTIFICATES[modulus]
        assert record["ordinary"]["vector_stream_sha256"] == certificate[
            "ordinary_vector_stream_sha256"
        ]
        assert record["trace_over_OK"]["vector_stream_sha256"] == certificate[
            "trace_vector_stream_sha256"
        ]


def test_smith_kernel_quotient_transport_handles_dependent_gram_types() -> None:
    # Backend assertions independently check S*M*T=D, unimodularity of S,T,
    # and that the zero columns of T are an integral kernel basis.
    rank_one = ((4, 4, 0), (4, 4, 0), (0, 0, 0))
    rank_two = ((4, 0, 4), (0, 6, 6), (4, 6, 10))
    hidden_norm_two = ((4, 3, 7), (3, 4, 7), (7, 7, 14))
    for matrix in (rank_one, rank_two, hidden_norm_two):
        diagonal, left, right = smith_normal_decomp(Matrix(matrix), domain=ZZ)
        assert left * Matrix(matrix) * right == diagonal
        assert abs(int(left.det())) == abs(int(right.det())) == 1
    assert lambda22_integer_gram_minimum_at_least(rank_one, 4)
    assert lambda22_integer_gram_minimum_at_least(rank_two, 4)
    assert not lambda22_integer_gram_minimum_at_least(hidden_norm_two, 4)
    assert _direct_box_minimum(rank_one) == 4
    assert _direct_box_minimum(rank_two) == 4
    assert _direct_box_minimum(hidden_norm_two) == 2


def test_complete_six_norm_alphabet_and_height_zero_exclusion() -> None:
    payload = _payload()
    seed = LAMBDA22_UNIT_SEED
    assert seed.admissible
    assert seed.canonical_mirror_norms == tuple(
        Cubic7(value) for value in (1, 2, 4, 3, 6, 12)
    )
    assert "v_2(q) is 0,1,2" in seed.local_scope
    assert "v_3(q) is 0,1" in seed.local_scope
    expected_gates = {
        1: (1, 4, 1),
        2: (1, 4, 1),
        4: (2, 12, 2**20),
        3: (3, 24, 3**21),
        6: (3, 24, 3**21),
        12: (6, 96, 2**20 * 3**21),
    }
    for norm in seed.canonical_mirror_norms:
        gate = lambda22_norm_gate(norm)
        assert (
            gate.pairing_modulus,
            gate.coefficient_lattice_minimum,
            gate.coefficient_lattice_index,
        ) == expected_gates[int(norm.a)]
        assert not lambda22_rank3_gram_types(
            norm, gate.coefficient_lattice_minimum, gate.pairing_modulus
        )
    assert payload["height_zero"]["reflecting_wall_count"] == 0


def test_all_h100_controls_and_staged_gram_types_replay_exactly() -> None:
    payload = _payload()
    prefix = payload["positive_prefix"]
    levels = lambda22_control_levels(100)
    assert len(levels) == prefix["control_count"] == 30
    assert len({level.height for level in levels}) == prefix["distinct_height_count"] == 24
    assert [row["count"] for row in prefix["control_counts_by_norm"]] == [
        5,
        16,
        5,
        1,
        2,
        1,
    ]
    assert len(prefix["control_levels"]) == len(levels)

    ambient_count = congruence_count = gated_count = 0
    for level, row in zip(levels, prefix["control_levels"], strict=True):
        assert level.height == _element(row["height"]["exact"])
        assert level.norm == _element(row["norm"]["exact"])
        assert level.control == _element(row["control"]["exact"])
        assert level.target == _element(row["spatial_square_target"]["exact"])
        audit = lambda22_level_audit(level)
        archived_ambient = tuple(
            tuple(item["tuple_U_V_W_uv_uw_vw_det"])
            for item in row["ambient_minimum_four_gram_types"]
        )
        archived_congruence = tuple(
            tuple(item["tuple_U_V_W_uv_uw_vw_det"])
            for item in row["pairing_congruence_gram_types"]
        )
        archived_gated = tuple(
            tuple(item["tuple_U_V_W_uv_uw_vw_det"])
            for item in row["fully_gated_gram_types"]
        )
        assert tuple(map(_gram_tuple, audit.ambient_gram_types)) == archived_ambient
        assert tuple(map(_gram_tuple, audit.congruence_gram_types)) == archived_congruence
        assert tuple(map(_gram_tuple, audit.gated_gram_types)) == archived_gated
        ambient_count += len(archived_ambient)
        congruence_count += len(archived_congruence)
        gated_count += len(archived_gated)

        target = level.target
        for item in audit.ambient_gram_types:
            assert Cubic7(
                item.u_square + 2 * item.v_dot_w - item.w_square,
                2 * item.u_dot_v + 4 * item.v_dot_w - item.w_square,
                2 * item.u_dot_w
                + item.v_square
                - 2 * item.v_dot_w
                + 3 * item.w_square,
            ) == target
            assert item.determinant >= 0
            # Independent bounded scan catches every short vector in all
            # archived small rank-three types; the Smith transport test above
            # supplies the formal kernel quotient certificate.
            assert _direct_box_minimum(item.matrix) >= 4

        gate = lambda22_norm_gate(level.norm)
        for item in audit.congruence_gram_types:
            assert all(
                value % gate.pairing_modulus == 0
                for value in _gram_tuple(item)[:6]
            )
            if gate.coefficient_lattice_minimum > 4:
                # Every final rejection has a direct diagonal contradiction;
                # the zero-wall result does not rest on a heuristic embedding.
                assert any(
                    0 < square < gate.coefficient_lattice_minimum
                    for square in (
                        item.u_square,
                        item.v_square,
                        item.w_square,
                    )
                )

    assert (ambient_count, congruence_count, gated_count) == (89, 11, 0)
    assert prefix["ambient_gram_type_rank_histogram"] == {"1": 3, "2": 6, "3": 80}
    assert prefix["pairing_congruence_rank_histogram"] == {"1": 3, "2": 6, "3": 2}
    assert prefix["accepted_positive_wall_count"] == 0


def test_semidefinite_types_are_retained_and_scope_is_strictly_bounded() -> None:
    payload = _payload()
    rows = payload["positive_prefix"]["control_levels"]
    first = next(
        row
        for row in rows
        if row["norm"]["display"] == "4"
        and row["height"]["display"] == "2*theta + 1*theta^2"
    )
    assert first["ambient_minimum_four_gram_types"] == [
        {"rank": 1, "tuple_U_V_W_uv_uw_vw_det": [4, 4, 0, 4, 0, 0, 0]}
    ]
    norm12 = next(row for row in rows if row["norm"]["display"] == "12")
    assert len(norm12["ambient_minimum_four_gram_types"]) == 78
    assert norm12["pairing_congruence_gram_types"] == [
        {
            "rank": 2,
            "tuple_U_V_W_uv_uw_vw_det": [12, 144, 36, 36, 18, 72, 0],
        }
    ]
    assert not norm12["fully_gated_gram_types"]
    assert payload["frontier"]["all_controls_at_height_at_most_100_exhausted"]
    assert payload["cone_status"]["positive_wall_count_through_H100"] == 0
    assert not payload["cone_status"]["finite_volume_test_applicable"]
    assert not payload["cone_status"]["compactness_test_applicable"]
    assert "no finite-volume" in payload["scope"]
    assert "global nonreflectivity" in payload["scope"]


@pytest.mark.slow
def test_pari_independently_replays_catalogue_minimum_and_kissing_number() -> None:
    gp = shutil.which("gp")
    if gp is None:
        pytest.skip("PARI/GP is unavailable")
    matrix = "[" + ";".join(
        ",".join(str(value) for value in row) for row in lambda22_gram()
    ) + "]"
    program = (
        f"r=qfminim({matrix});"
        'print(r[1]);print(r[2]);print(matsize(r[3])[2]);print(r[3]);\n'
    )
    result = subprocess.run(
        [gp, "-q", "-s", "512000000"],
        input=program,
        text=True,
        capture_output=True,
        check=True,
        timeout=120,
    )
    lines = result.stdout.splitlines(keepends=True)
    assert [int(value) for value in lines[:3]] == [49_896, 4, 24_948]
    assert hashlib.sha256(lines[3].encode()).hexdigest() == PARI_QFMINIM_CERTIFICATES[
        1
    ]["ordinary_vector_stream_sha256"]
