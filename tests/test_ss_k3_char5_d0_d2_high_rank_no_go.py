from collections import Counter
from fractions import Fraction
import json

import numpy as np

from scripts.certify_ss_k3_char5_d0_d2_high_rank_no_go import (
    EXPECTED_HASHES,
    MODULUS,
    OUTPUT,
    _object_dot,
    parse_source,
)


def _artifact() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_source_lineage_and_lattice_are_exact() -> None:
    source = parse_source()
    artifact = _artifact()
    assert source["hashes"] == EXPECTED_HASHES
    assert artifact["lattice"] == {
        "rank": 22,
        "determinant": -25,
        "signature": [1, 21],
        "dual_form_used": "5*GramSX^-1",
        "dual_form_integral": True,
    }
    assert [len(chamber["walls"]) for chamber in source["chambers"]] == [
        16302,
        7832,
        2534,
    ]


def test_exhaustive_high_rank_counts_and_noncyclic_bound() -> None:
    artifact = _artifact()
    assert artifact["summary"] == {
        "high_rank_conjugacy_class_count": 7,
        "hyperbolic_dimensions": [13],
        "positive_examples": 0,
        "non_coxeter_casts": 7,
        "finite_volume_tests_run": 0,
        "decision": "published_AutX_D0_D2_groups_closed_negative_above_H9",
    }
    assert MODULUS == 1009
    expected = {
        "D0": (525, 1, 10),
        "D1": (315, 2, 10),
        "D2": (99, 4, 10),
    }
    for chamber in artifact["chambers"]:
        involutions, class_count, pair_max = expected[chamber["chamber"]]
        assert chamber["exhaustive_fixed_rank_screen"][
            "rank_14_involution_count"
        ] == involutions
        assert len(chamber["high_rank_conjugacy_classes"]) == class_count
        assert chamber["distinct_rank_14_pair_screen"][
            "maximum_common_fixed_rank"
        ] == pair_max
        assert sum(
            entry["conjugacy_class_size"]
            for entry in chamber["high_rank_conjugacy_classes"]
        ) == involutions


def test_every_class_has_an_exact_strict_forbidden_pair() -> None:
    artifact = _artifact()
    source = parse_source()
    scaled_dual = source["scaled_dual"]
    ratios = Counter()
    for chamber in artifact["chambers"]:
        for class_record in chamber["high_rank_conjugacy_classes"]:
            pair = class_record["strict_bad_pair"]
            left = pair["left"]
            right = pair["right"]
            left_normal = np.asarray(left["restricted_normal_covector"], dtype=object)
            right_normal = np.asarray(right["restricted_normal_covector"], dtype=object)
            assert _object_dot(left_normal, scaled_dual, left_normal) == left[
                "scaled_square"
            ]
            assert _object_dot(right_normal, scaled_dual, right_normal) == right[
                "scaled_square"
            ]
            product = _object_dot(left_normal, scaled_dual, right_normal)
            assert product == pair["scaled_product"]
            ratio = Fraction(
                product * product,
                left["scaled_square"] * right["scaled_square"],
            )
            assert str(ratio) == pair["normalized_squared_product"]
            assert ratio not in {Fraction(1, 4), Fraction(1, 2), Fraction(3, 4)}
            assert 0 < ratio < 1
            assert left["minimum_other_wall_pairing"] > 0
            assert right["minimum_other_wall_pairing"] > 0
            assert left["scaled_witness_square"] > 0
            assert right["scaled_witness_square"] > 0
            assert pair["decision"] == "not_coxeter"
            ratios[str(ratio)] += 1
    assert ratios == {"5/16": 1, "5/12": 3, "5/8": 3}


def test_certificate_scope_does_not_overclaim_full_lattice_stabilizers() -> None:
    artifact = _artifact()
    warning = artifact["scope_warning"]
    assert "not a theorem about the full lattice stabilizers Aut(D_i)" in warning
    assert "Aut_X(D_i)=Aut(D_i) intersect G_X" in warning

