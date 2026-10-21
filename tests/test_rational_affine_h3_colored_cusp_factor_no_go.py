"""Hostile exact tests for the affine-H3 colored cusp-factor local no-go."""

from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = (
    ROOT / "scripts/certify_rational_affine_h3_colored_cusp_factor_no_go.py"
)


def load_certificate_module():
    specification = importlib.util.spec_from_file_location(
        "rational_affine_h3_colored_cusp_factor_no_go",
        SCRIPT,
    )
    assert specification is not None and specification.loader is not None
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def test_every_even_rank_packet_has_exact_affine_null_relation():
    module = load_certificate_module()
    for rank in range(4, 52, 2):
        gram = module.affine_gram(rank)
        marks = module.affine_marks(rank)
        assert module.exact_rank(gram) == rank
        assert module.mat_vec(gram, marks) == [0] * (rank + 1)


def test_physical_divisor_is_primitive_nef_numerical_genus_three_packet():
    module = load_certificate_module()
    for rank in range(4, 52, 2):
        physical = module.physical_intersection_gram(rank)
        multiplicities = module.physical_multiplicities(rank)
        assert module.mat_vec(physical, multiplicities) == [0] * (rank + 1)
        assert min(multiplicities) == 1
        row = module.certify_rank(rank)
        assert row["canonical_mark_sum"] == 4
        assert row["positive_mark_mass"] - row["negative_mark_mass"] == 2


def test_proximity_clusters_include_the_affine_d4_special_case():
    module = load_certificate_module()
    assert module.proximity_clusters(4) == [(1, 1, 1, 1)]
    for rank in range(6, 52, 2):
        clusters = module.proximity_clusters(rank)
        assert all(sum(cluster) == 4 for cluster in clusters)
        assert sorted(len(cluster) for cluster in clusters).count(3) == 2
        assert all(len(cluster) in {2, 3} for cluster in clusters)


def test_coordinate_refinement_is_zero_colored_but_glue_parity_is_extra():
    module = load_certificate_module()
    for rank in range(4, 52, 2):
        coordinate = module.coordinate_refinement(rank)
        roots = module.finite_simple_roots(rank)
        coefficients = module.coordinate_coefficients(rank)
        assert module.dot(coordinate, coordinate) == 4
        assert all(module.dot(coordinate, root) % 2 == 0 for root in roots)
        assert (
            sum(
                coefficient * color
                for coefficient, color in zip(
                    coefficients,
                    module.canonical_colors(rank)[1:],
                )
            )
            == 0
        )

        # An ambient integral glue generator can pair oddly with a
        # primitive norm-four coordinate vector.  The factor Gram only
        # records the even pairings above and cannot exclude this value.
        outside_pairing = 1
        assert outside_pairing % 2 == 1


def test_certificate_keeps_the_surface_and_maximality_claims_out_of_scope():
    module = load_certificate_module()
    payload = module.certificate()
    assert payload["status"] == "pass"
    excluded = payload["scope"]["not_proved"]
    assert any("projective rational" in item for item in excluded)
    assert any("reflection-maximal" in item for item in excluded)
    assert any("zero-color chains" in item for item in excluded)

    report = module.REPORT.read_text()
    assert r"B_n-B_{n-1}\in2\operatorname {NS}(S)" in report
