from collections import Counter
import json

import sympy as sp

from scripts.search_ss_fermat_quartic_char3_folds import (
    OUTPUT,
    build_artifact,
    parse_source,
)


def test_published_frobenius_has_h11_fixed_space() -> None:
    source = parse_source()
    gram = source["gram"]
    frobenius = source["frobenius"]
    ample = source["ample"]
    assert gram.shape == (22, 22)
    assert gram.det() == -9
    assert frobenius**2 == sp.eye(22)
    assert frobenius * gram * frobenius.T == gram
    assert ample * frobenius == ample
    assert 22 - (frobenius - sp.eye(22)).rank() == 12
    assert len(source["walls"]) == 112 + 648 + 5184


def test_strict_facets_give_exact_forbidden_angle_certificate() -> None:
    artifact = build_artifact()
    gate = artifact["strict_facet_gate"]
    assert artifact["frobenius"]["hyperbolic_dimension"] == 11
    assert artifact["walls"]["frobenius_orbit_count"] == 3000
    assert artifact["walls"]["distinct_restricted_ray_count"] == 2316
    assert artifact["walls"]["restricted_gram_inertia"] == [1, 11, 2304]
    assert gate["certified_strict_facet_count"] == 384
    assert gate["bad_pair_count"] == 864
    assert gate["decision"] == "not_coxeter"
    assert Counter(
        pair["normalized_squared_product"] for pair in gate["bad_pairs"]
    ) == {"3/8": 576, "1/16": 288}

    certificates = {
        certificate["restricted_ray_index"]: certificate
        for certificate in gate["certificates"]
    }
    first = certificates[1]
    second = certificates[128]
    source = parse_source()
    inverse = source["gram"].inv()
    left = sp.Matrix(1, 22, first["primitive_normal_covector"])
    right = sp.Matrix(1, 22, second["primitive_normal_covector"])
    left_norm = sp.Rational((left * inverse * left.T)[0])
    right_norm = sp.Rational((right * inverse * right.T)[0])
    product = sp.Rational((left * inverse * right.T)[0])
    assert left_norm == -4
    assert right_norm == sp.Rational(-8, 3)
    assert product == 2
    assert sp.factor(product**2 / (left_norm * right_norm)) == sp.Rational(3, 8)


def test_archived_frobenius_census_replays() -> None:
    assert json.loads(OUTPUT.read_text(encoding="utf-8")) == build_artifact()
