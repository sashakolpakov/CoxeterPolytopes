from __future__ import annotations

import json

import sympy as sp

from scripts.certify_roulleau_r14_c2t_pure_order4_obstruction import (
    BH_ORDER4_ROW_COUNT,
    BH_TARGETS,
    OUTPUT,
    _index_square_gate,
    build_certificate,
    smith_diagonal,
)
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


CERTIFICATE = json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_complete_pure_order4_obstruction_archive_replays() -> None:
    assert build_certificate() == CERTIFICATE
    assert CERTIFICATE["all_checks_pass"] is True
    assert CERTIFICATE["scope"]["exact_rank13_fixed_pure_order4_lift_ruled_out"]
    assert not CERTIFICATE["scope"][
        "pure_order4_lift_with_enlarged_global_fixed_lattice_ruled_out"
    ]
    assert not CERTIFICATE["scope"]["symplectic_square_order4_lift_ruled_out"]
    assert not CERTIFICATE["scope"]["higher_order_lift_ruled_out"]


def test_lefschetz_gate_and_classification_multiplicities_are_exact() -> None:
    gate = CERTIFICATE["pure_order4_lefschetz_gate"]
    assert gate["surviving_r_l_m_alpha"] == [[13, 3, 3, 2], [13, 7, 1, 1]]
    classification = CERTIFICATE["brandhorst_hofmann_classification"]
    assert classification["complete_row_count"] == BH_ORDER4_ROW_COUNT == 79
    assert [
        (row["classification_id"], row["multiplicity_in_complete_order4_file"])
        for row in classification["surviving_rows"]
    ] == [("0.4.2.9", 1), ("0.4.0.2", 1)]


def test_both_classified_fixed_lattices_fail_the_index_square_gate() -> None:
    expected = {
        (13, 3, 3): (64, [*[1] * 10, 4, 4, 4], "1/4"),
        (13, 7, 1): (256, [*[1] * 6, *[2] * 6, 4], "1/16"),
    }
    for character, (determinant, smith, ratio) in expected.items():
        gram = sp.Matrix(BH_TARGETS[character]["fixed_gram"])
        assert exact_inertia(gram) == (1, 12, 0)
        assert gram.det() == determinant
        assert smith_diagonal(gram) == smith
        gate = _index_square_gate(16, determinant)
        assert gate == {
            "det_F_over_det_S": ratio,
            "required_identity": "det(F)=[S:F]^2*det(S)",
            "is_nonnegative_integer_square": False,
            "containment_possible": False,
        }


def test_hyperbolic_handoff_remains_finite_volume_noncompact() -> None:
    handoff = CERTIFICATE["hyperbolic_polytope_handoff"]
    assert handoff == {
        "dimension": 12,
        "facets": 23,
        "classification": "finite_volume_noncompact",
        "elliptic_edges": 2308,
        "finite_vertices": 320,
        "cusps": 17,
        "bad_edges": 0,
        "bad_parabolic_neighborhoods": 0,
        "section_to_odd_unimodular_full_group_index": 8320,
    }
