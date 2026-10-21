import json
from pathlib import Path

import pytest

from coxeter_polytopes.exact.cubic7 import (
    LABEL14_PRODUCT,
    LABEL7_PRODUCT,
    THETA7,
    Cubic7,
)
from coxeter_polytopes.search.cubic7_diagonal import run_diagonal_vinberg_search
from coxeter_polytopes.verify.cubic7_alvin import (
    Cubic7AlVinError,
    alvin_lambda_integer,
    canonical_pair_product_gram,
    parse_alvin_element,
    parse_alvin_roots_text,
    parse_coxiter_text,
    parse_diagonal_form,
    verify_cubic7_alvin_candidate,
    verify_cubic7_alvin_files,
)


_PRISM_COXITER = """\
5 3
vertices labels: 1 2 3 4 5
1 2 3
2 3 4
1 4 7
3 5 1
"""

_FIXTURES = Path(__file__).with_name("fixtures") / "cubic7_alvin"


_VINNY_RC7_COXITER = """\
4 2
vertices labels: 1 2 3 4
1 3 1
2 4 1
3 4 14
"""

_VINNY_RC7_ROOTS = """\
Vectors:
    e1 = (0, -1, 0)
    e2 = (0, 0, [3,1,2])
    e3 = ([2,0,1], [1,0,0], 0)
    e4 = ([-3,-1,-3], 0, [-4,-1,-3])
Algorithm ended
"""


def _compact_prism_replay():
    search = run_diagonal_vinberg_search(dimension=3, height_bound=7)
    return verify_cubic7_alvin_candidate(
        parse_coxiter_text(_PRISM_COXITER),
        parse_diagonal_form("[-1,0,0],1,1,1"),
        search.walls,
    )


def test_alvin_trigonometric_basis_is_converted_to_power_basis_exactly():
    assert alvin_lambda_integer(1, 0, 0) == THETA7
    assert alvin_lambda_integer(0, 1, 0) == THETA7**2 - 2
    assert alvin_lambda_integer(0, 0, 1) == 1 - THETA7 - THETA7**2
    assert alvin_lambda_integer(-3, -3, -3) == 3
    assert parse_alvin_element("RCI7(1,0,0)") == THETA7
    assert parse_diagonal_form("<[-1,0,0], 1, 1, 1, 1>") == (
        -THETA7,
        Cubic7(1),
        Cubic7(1),
        Cubic7(1),
        Cubic7(1),
    )


def test_generic_pari_and_json_root_envelopes_parse():
    generic = """\
Vectors:
    e1 = (0, -1, 1)
    e2 = (RCI7(1,0,0), [0,1,0], 2)
"""
    assert parse_alvin_roots_text(generic) == (
        (Cubic7(0), Cubic7(-1), Cubic7(1)),
        (THETA7, THETA7**2 - 2, Cubic7(2)),
    )
    pari = """\
vectors = Mat( [0, -1, 1] );
vectors = concat(vectors, [RCI7(1,0,0), RCI7(0,1,0), 2] );
"""
    assert parse_alvin_roots_text(pari) == parse_alvin_roots_text(generic)
    payload = json.dumps(
        {
            "basis": "power",
            "roots": [
                [["0", "0", "0"], ["1", "0", "0"]],
                [["0", "1", "0"], ["0", "0", "1"]],
            ],
        }
    )
    assert parse_alvin_roots_text(payload) == (
        (Cubic7(0), Cubic7(1)),
        (THETA7, THETA7**2),
    )


def test_known_compact_cubic_prism_passes_every_independent_exact_stage():
    replay = _compact_prism_replay()
    assert replay.geometry_status == "compact"
    assert replay.compact is True
    assert replay.finite_volume is True
    assert len(replay.extreme_rays) == 6
    assert len(replay.finite_vertices) == 6
    assert replay.cusps == ()
    assert replay.hyperideal_rays == ()
    assert replay.ground_field == "Q(theta)"
    assert replay.pair_product_gram[0][3] == LABEL7_PRODUCT
    assert replay.canonical_fingerprint == (
        "620c8fedc0b05d7411c0470cfb035e2c5901916786b9c6d366834ff4a4018e11"
    )
    assert replay.to_json()["scope"][
        "alvin_wall_enumeration_completeness_independently_checked"
    ] is False


def test_alvin_eight_wall_diag_minus_theta_i4_control_is_exact_and_compact():
    replay = verify_cubic7_alvin_files(
        _FIXTURES / "diag_minus_theta_i4.coxiter",
        _FIXTURES / "diag_minus_theta_i4.alvin.txt",
        "[-1,0,0],1,1,1,1",
    )
    assert replay.dimension == 4
    assert len(replay.roots) == 8
    assert replay.geometry_status == "compact"
    assert replay.compact is True
    assert replay.finite_volume is True
    assert len(replay.finite_vertices) == 16
    assert replay.cusps == ()
    assert replay.hyperideal_rays == ()
    assert replay.ground_field == "Q(theta)"
    assert replay.canonical_fingerprint == (
        "468e8f2c4efe26300f8fa16289bd3b37b4e9e9d13ac28188f78a8ab226c566d2"
    )
    # The regression root fixture is exactly reconstructed; it is deliberately
    # not presented as a captured AlVin termination transcript.
    assert replay.alvin_declared_termination is None
    payload = replay.to_json()
    assert payload["vinberg_arithmeticity"]["passes"] is True
    assert payload["scope"][
        "alvin_wall_enumeration_completeness_independently_checked"
    ] is False


def test_vinny_rc7_four_wall_order_14_regression_is_exact_and_compact():
    """Replay Bottinelli's correction of the five-wall AlVin output.

    The paper's trigonometric basis is exactly AlVin's lambda basis.  Its
    fourth wall may be unit-square normalized from
    ``([-3,-1,-3],0,[-4,-1,-3])`` to
    ``([0,1,-1],0,[-1,0,-1])``; the latter has canonical norm
    ``[-6,-4,-4]`` and is the root exercising AlVin's RC7 fraction batches.
    """

    roots = parse_alvin_roots_text(_VINNY_RC7_ROOTS)
    unit = Cubic7(5, -1, -2)
    normalized_missing_wall = (
        alvin_lambda_integer(0, 1, -1),
        Cubic7(0),
        alvin_lambda_integer(-1, 0, -1),
    )
    assert unit.norm() == 1
    assert all(unit.sign(embedding) > 0 for embedding in range(3))
    assert tuple(unit * entry for entry in roots[3]) == normalized_missing_wall

    diagonal = parse_diagonal_form("[-1,0,0],1,[-2,-5,-7]")
    normalized_norm = sum(
        (
            coefficient * coordinate * coordinate
            for coefficient, coordinate in zip(diagonal, normalized_missing_wall)
        ),
        start=Cubic7(0),
    )
    assert normalized_norm == alvin_lambda_integer(-6, -4, -4)
    assert all(
        (2 * coefficient * coordinate / normalized_norm).is_integral
        for coefficient, coordinate in zip(diagonal, normalized_missing_wall)
    )
    assert normalized_missing_wall[0] ** 2 / normalized_norm == (
        Cubic7(2) - THETA7**2 / 2
    )

    replay = verify_cubic7_alvin_candidate(
        parse_coxiter_text(_VINNY_RC7_COXITER),
        diagonal,
        roots,
        transcript_text=_VINNY_RC7_ROOTS,
    )
    assert replay.geometry_status == "compact"
    assert replay.compact is True
    assert replay.finite_volume is True
    assert len(replay.finite_vertices) == 4
    assert replay.cusps == ()
    assert replay.hyperideal_rays == ()
    assert replay.pair_product_gram[2][3] == LABEL14_PRODUCT
    assert replay.canonical_fingerprint == (
        "8ec69bb39a2f6eb0ded0f6ee3332da3bd1443ad5dcce8f9ff21219a7e0c89247"
    )


def test_canonical_pair_product_gram_is_invariant_under_every_tested_permutation():
    replay = _compact_prism_replay()
    matrix = replay.pair_product_gram
    for permutation in (
        (4, 3, 2, 1, 0),
        (1, 3, 0, 4, 2),
        (2, 0, 4, 1, 3),
    ):
        permuted = tuple(
            tuple(matrix[permutation[i]][permutation[j]] for j in range(5))
            for i in range(5)
        )
        canonical, fingerprint = canonical_pair_product_gram(permuted)
        assert canonical == replay.canonical_pair_product_gram
        assert fingerprint == replay.canonical_fingerprint


def test_canonical_pair_product_gram_distinguishes_a_changed_exact_edge():
    replay = _compact_prism_replay()
    changed = [list(row) for row in replay.pair_product_gram]
    changed[0][1] = changed[1][0] = Cubic7(2)
    _, fingerprint = canonical_pair_product_gram(tuple(tuple(row) for row in changed))
    assert fingerprint != replay.canonical_fingerprint


def test_coxiter_is_checked_against_exact_cartan_products():
    search = run_diagonal_vinberg_search(dimension=3, height_bound=7)
    wrong = _PRISM_COXITER.replace("1 4 7", "1 4 3")
    with pytest.raises(Cubic7AlVinError, match="label mismatch"):
        verify_cubic7_alvin_candidate(
            parse_coxiter_text(wrong),
            parse_diagonal_form("[-1,0,0],1,1,1"),
            search.walls,
        )


def test_nonadmissible_and_noncrystallographic_inputs_fail_closed():
    search = run_diagonal_vinberg_search(dimension=3, height_bound=7)
    with pytest.raises(Cubic7AlVinError, match="not admissible"):
        verify_cubic7_alvin_candidate(
            parse_coxiter_text(_PRISM_COXITER),
            (THETA7, Cubic7(1), Cubic7(1), Cubic7(1)),
            search.walls,
        )
    malformed = list(search.walls)
    malformed[-1] = tuple(
        entry / 2 if index == 0 else entry
        for index, entry in enumerate(malformed[-1])
    )
    with pytest.raises(Cubic7AlVinError):
        verify_cubic7_alvin_candidate(
            parse_coxiter_text(_PRISM_COXITER),
            parse_diagonal_form("[-1,0,0],1,1,1"),
            malformed,
        )
