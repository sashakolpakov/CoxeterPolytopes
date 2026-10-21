import hashlib
import json

from scripts.certify_k3_rank15_order3_fold import (
    ARCHIVE_PATH,
    ORDER_THREE_CYCLES,
    SOURCE_PATH,
    _is_order_exactly_three,
    _orbits,
    _permutation,
)
from scripts.verify_k3_high_rank_nef_chambers import build_gram, exact_inertia
from scripts.verify_k3_quotient_nef_slices import _orbit_sum_matrix


def _archive() -> dict:
    return json.loads(ARCHIVE_PATH.read_text())


def test_rank15_order_three_fold_reconstructs_quickly() -> None:
    source = json.loads(SOURCE_PATH.read_text())
    case = next(case for case in source["cases"] if case["rho"] == 15)
    wall_gram = build_gram(case)
    permutation = _permutation(wall_gram.rows, ORDER_THREE_CYCLES)
    assert _is_order_exactly_three(permutation)
    assert all(
        wall_gram[left, right]
        == wall_gram[permutation[left], permutation[right]]
        for left in range(wall_gram.rows)
        for right in range(wall_gram.rows)
    )
    orbits = _orbits(permutation)
    orbit_sums = _orbit_sum_matrix(wall_gram.rows, orbits)
    orbit_gram = orbit_sums.T * wall_gram * orbit_sums
    assert orbit_gram.rows == 17
    assert orbit_gram.rank() == 13
    assert exact_inertia(orbit_gram) == (1, 12, 4)
    assert [int(orbit_gram[index, index]) for index in range(17)] == [
        *([-2] * 13),
        -6,
        -6,
        -2,
        -2,
    ]


def test_frozen_literal_and_optimized_volume_certificates_agree() -> None:
    archived = _archive()
    assert archived["source_sha256"] == hashlib.sha256(SOURCE_PATH.read_bytes()).hexdigest()
    assert archived["all_checks_pass"] is True
    assert archived["finite_volume_certified"] is True
    assert archived["compact"] is False
    assert (
        archived["artifact_status"]
        == "exact_literal_and_optimized_finite_volume_certified"
    )
    assert (
        archived["hyperbolic_dimension"],
        archived["facet_count"],
        archived["finite_vertex_count"],
        archived["ideal_vertex_count"],
        archived["hyperideal_vertex_count"],
    ) == (12, 17, 93, 3, 0)
    assert archived["coxeter"]["edge_counts"] == {
        "2": 116,
        "3": 12,
        "4": 0,
        "6": 2,
        "infinity": 3,
        "ultraparallel": 3,
    }
    optimized = archived["optimized_vinberg"]
    literal = archived["literal_coxiter_prop_1_10"]
    assert optimized["decision"] == "finite_volume_noncompact"
    assert (
        optimized["elliptic_edge_count"],
        optimized["finite_vertex_count"],
        optimized["ideal_vertex_count"],
    ) == (623, 93, 3)
    assert optimized["finite_volume_bad_edges"] == []
    assert literal["spherical_rank_11_edge_count"] == 623
    assert literal["spherical_rank_12_vertex_count"] == 93
    assert literal["parabolic_rank_11_vertex_count"] == 3
    assert literal["edge_extension_histogram"] == {"2": 623}
    assert literal["endpoint_type_histogram"] == {
        "0,2": 3,
        "1,1": 124,
        "2,0": 496,
    }
    assert literal["bad_edges"] == []
    assert literal["finite_covolume"] is True
    assert literal["matches_optimized_endpoint_sets"] is True
    assert all(archived["checks"].values())


def test_order_three_NS_action_and_generic_realization_boundary() -> None:
    action = _archive()["neron_severi_action"]
    assert action["all_checks_pass"] is True
    assert action["determinant"] == 1
    assert action["fixed_rank"] == 13
    assert action["moving_rank"] == 2
    assert action["discriminant_action_is_identity"] is False
    assert action["discriminant_action_is_minus_identity"] is False
    assert action["discriminant_action_has_exact_order_three"] is True
    assert action["generic_K3_realization"]["status"] == (
        "not_realized_on_general_K3_special_period_unassessed"
    )
