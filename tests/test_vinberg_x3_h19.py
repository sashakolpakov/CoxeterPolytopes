from __future__ import annotations

import json

import sympy as sp

from scripts.certify_vinberg_x3_h19 import ARCHIVE, build_wall_gram
from scripts.verify_k3_high_rank_nef_chambers import exact_inertia


def _frozen() -> dict:
    return json.loads(ARCHIVE.read_text())


def test_exact_36_wall_gram_reconstructs_frozen_lattice() -> None:
    frozen = _frozen()
    gram = build_wall_gram()
    assert gram == sp.Matrix(frozen["wall_gram"])
    assert gram.rank() == 20
    assert exact_inertia(gram) == (1, 19, 16)
    assert frozen["lattice"]["wall_gram_smith_nonzero"] == [1] * 19 + [3]
    assert frozen["coxeter"]["pair_counts"] == {
        "2": 501,
        "3": 27,
        "6": 36,
        "infinity": 24,
        "ultraparallel": 42,
    }


def test_frozen_literal_vinberg_endpoint_certificate() -> None:
    frozen = _frozen()
    assert frozen["artifact_status"] == (
        "exact_known_finite_volume_noncompact_chamber"
    )
    assert frozen["hyperbolic_dimension"] == 19
    assert frozen["facet_count"] == 36
    assert frozen["vertices_and_cusps"]["finite_vertex_count"] == 8346
    assert frozen["vertices_and_cusps"]["ideal_vertex_count"] == 107
    assert frozen["vinberg_endpoint_audit"]["elliptic_rank18_edge_count"] == 87046
    assert frozen["vinberg_endpoint_audit"]["endpoint_type_histogram"] == {
        "finite_0_ideal_2": 487,
        "finite_1_ideal_1": 14544,
        "finite_2_ideal_0": 72015,
    }
    assert frozen["vinberg_endpoint_audit"]["extension_histogram"] == {
        "2": 87046
    }
    assert frozen["vinberg_endpoint_audit"]["bad_edge_count"] == 0
    assert frozen["volume"] == {
        "compact": False,
        "finite_volume": True,
        "reason_noncompact": "107 ideal vertices",
    }
    assert frozen["all_checks_pass"] is True
