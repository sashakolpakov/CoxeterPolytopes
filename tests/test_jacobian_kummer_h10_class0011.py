from __future__ import annotations

from fractions import Fraction
import json
from math import comb

from scripts.screen_jacobian_kummer_h10_class0011 import (
    OUTPUT,
    build_certificate,
)


def _load() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_class0011_exact_certificate_replays_from_published_stock() -> None:
    archived = _load()
    assert build_certificate() == archived


def test_class0011_is_a_certified_h10_non_coxeter_slice() -> None:
    certificate = _load()
    status = certificate["status"]

    assert certificate["group_class"]["id"] == "class_0011"
    assert certificate["group_class"]["fixed_rank"] == 11
    assert certificate["group_class"]["hyperbolic_dimension"] == 10
    assert certificate["wall_orbit_count"] == 178
    assert certificate["wall_orbit_size_histogram"] == {"1": 40, "2": 138}
    assert certificate["orbit_normal_rank"] == 11
    assert certificate["orbit_gram_inertia"] == [1, 10, 167]
    assert certificate["facet_count"] == 76
    assert certificate["redundant_orbit_count"] == 102
    assert len(certificate["facet_witnesses"]) == 76
    assert len(certificate["redundancy_certificates"]) == 102
    assert certificate["facet_normal_square_histogram"] == {
        "-1": 32,
        "-2": 16,
        "-3/4": 16,
        "-4": 12,
    }

    assert certificate["coxeter_edge_counts"] == {
        "2": 1_068,
        "3": 32,
        "4": 184,
        "infinity": 406,
        "ultraparallel": 960,
    }
    assert certificate["non_coxeter_pair_count"] == 200
    assert certificate["non_coxeter_normalized_product_histogram"] == {
        "1/3": 112,
        "1/9": 24,
        "2/3": 64,
    }
    classified_pairs = sum(certificate["coxeter_edge_counts"].values())
    assert classified_pairs + certificate["non_coxeter_pair_count"] == comb(76, 2)

    assert status["exact_fixed_slice"] is True
    assert status["exact_wall_reduction"] is True
    assert status["strict_facet_witnesses"] is True
    assert status["exact_coxeter_pair_gate"] is True
    assert status["all_facet_pairs_are_coxeter"] is False
    assert status["vinberg_volume_test_applicable"] is False
    assert status["decision"] == "not_coxeter"


def test_each_archived_forbidden_pair_replays_in_the_facet_gram() -> None:
    certificate = _load()
    gram = [
        [Fraction(entry) for entry in row]
        for row in certificate["facet_gram"]
    ]
    facet_positions = {
        orbit_index: position
        for position, orbit_index in enumerate(certificate["facet_orbit_indices"])
    }

    histogram: dict[str, int] = {}
    for pair in certificate["non_coxeter_pairs"]:
        left_orbit, right_orbit = pair["facet_orbit_indices"]
        left = facet_positions[left_orbit]
        right = facet_positions[right_orbit]
        left_square = gram[left][left]
        right_square = gram[right][right]
        product = gram[left][right]
        ratio = product * product / (left_square * right_square)

        assert str(left_square) == pair["left_square"]
        assert str(right_square) == pair["right_square"]
        assert str(product) == pair["inner_product"]
        assert str(ratio) == pair["normalized_squared_product"]
        histogram[str(ratio)] = histogram.get(str(ratio), 0) + 1

    assert histogram == certificate["non_coxeter_normalized_product_histogram"]
