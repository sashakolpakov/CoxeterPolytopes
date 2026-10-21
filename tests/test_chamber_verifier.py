from copy import deepcopy
from fractions import Fraction
from types import SimpleNamespace

import pytest

from coxeter_polytopes.exact.golden import Golden, PHI
from coxeter_polytopes.verify.chamber import CertificateError, verify_chamber_geometry
from coxeter_polytopes.verify.vinberg import verify_vinberg_prefix_stopping


def encoded(value: Golden | int | Fraction) -> dict[str, str]:
    return Golden.coerce(value).to_json()


def encoded_vector(values) -> list[dict[str, str]]:
    return [encoded(value) for value in values]


def encoded_matrix(rows) -> list[list[dict[str, str]]]:
    return [encoded_vector(row) for row in rows]


def compact_335_triangle_certificate() -> dict:
    # The arithmetic hyperbolic triangle with angles pi/3, pi/3, pi/5.
    gram = [
        [2, -1, -PHI],
        [-1, 2, -1],
        [-PHI, -1, 2],
    ]
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "compact_reflective",
        "field": "Q(sqrt(5))",
        "dimension": 2,
        "lattice_gram": encoded_matrix(gram),
        "simple_roots": encoded_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        "wall_gram": encoded_matrix(gram),
        "coxeter_relations": [[1, 3, 5], [3, 1, 3], [5, 3, 1]],
        "completeness": {
            "method": "synthetic_test_fixture",
            "evidence": {"description": "the three defining walls are the input polyhedron"},
        },
        "geometry": {
            "nonempty": True,
            "interior_point": encoded_vector([1, Fraction(1, 2), 1]),
            "vertices": [[0, 1], [0, 2], [1, 2]],
            "cusps": [],
        },
    }


def cofinite_33infinity_triangle_certificate() -> dict:
    gram = [[2, -1, -2], [-1, 2, -1], [-2, -1, 2]]
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "cofinite_noncompact",
        "field": "Q",
        "dimension": 2,
        "lattice_gram": encoded_matrix(gram),
        "simple_roots": encoded_matrix([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),
        "wall_gram": encoded_matrix(gram),
        "coxeter_relations": [
            [1, 3, "parallel"],
            [3, 1, 3],
            ["parallel", 3, 1],
        ],
        "completeness": {
            "method": "synthetic_test_fixture",
            "evidence": {"description": "the three defining walls are the input polyhedron"},
        },
        "geometry": {
            "nonempty": True,
            "interior_point": encoded_vector([1, Fraction(1, 2), 1]),
            "vertices": [[0, 1], [1, 2]],
            "cusps": [[0, 2]],
        },
    }


def isotropic_but_nonparabolic_cusp_certificate() -> dict:
    # Facets 0 and 1 form affine A1, while facet 2 is an orthogonal finite A1.
    # Their common isotropic extreme ray consequently has an unbounded
    # horospherical link and is not a finite-volume cusp.
    gram = [
        [2, -2, 0, -1],
        [-2, 2, 0, -1],
        [0, 0, 2, -1],
        [-1, -1, -1, 2],
    ]
    return {
        "schema_version": "0.1.0",
        "certificate_kind": "cofinite_noncompact",
        "field": "Q",
        "dimension": 3,
        "lattice_gram": encoded_matrix(gram),
        "simple_roots": encoded_matrix(
            [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
        ),
        "wall_gram": encoded_matrix(gram),
        "coxeter_relations": [
            [1, "parallel", 2, 3],
            ["parallel", 1, 2, 3],
            [2, 2, 1, 3],
            [3, 3, 3, 1],
        ],
        "completeness": {
            "method": "synthetic_test_fixture",
            "evidence": {"description": "negative finite-volume fixture"},
        },
        "geometry": {
            "nonempty": True,
            "interior_point": encoded_vector([Fraction(3, 2), Fraction(3, 2), 0, 1]),
            "vertices": [],
            "cusps": [[0, 1, 2]],
        },
    }


def test_exact_compact_triangle_geometry() -> None:
    result = verify_chamber_geometry(compact_335_triangle_certificate())
    assert result.dimension == 2
    assert result.facet_count == 3
    assert result.finite_vertices == ((0, 1), (0, 2), (1, 2))
    assert result.cusps == ()
    assert result.hyperideal_rays == ()


def test_exact_cofinite_triangle_with_one_cusp() -> None:
    result = verify_chamber_geometry(cofinite_33infinity_triangle_certificate())
    assert result.finite_vertices == ((0, 1), (1, 2))
    assert result.cusps == ((0, 2),)


def test_exact_cofinite_vinberg_prefix_can_trigger_stopping() -> None:
    replayed_prefix = SimpleNamespace(
        roots=(
            (Golden(1), Golden(0), Golden(0)),
            (Golden(0), Golden(1), Golden(0)),
            (Golden(0), Golden(0), Golden(1)),
        )
    )
    result = verify_vinberg_prefix_stopping(
        replayed_prefix,
        cofinite_33infinity_triangle_certificate(),
    )
    assert result.chamber.cusps == ((0, 2),)
    assert result.diagram.decision == "finite_volume_noncompact"
    assert result.diagram.elliptic_edge_count == 3


def test_cusp_cannot_be_declared_compact() -> None:
    certificate = cofinite_33infinity_triangle_certificate()
    certificate["certificate_kind"] = "compact_reflective"
    with pytest.raises(CertificateError, match="compact certificate cannot contain cusps"):
        verify_chamber_geometry(certificate)


def test_isotropic_ray_with_nonparabolic_link_is_not_a_cusp() -> None:
    with pytest.raises(CertificateError, match="componentwise affine parabolic"):
        verify_chamber_geometry(isotropic_but_nonparabolic_cusp_certificate())


def test_wrong_coxeter_label_is_rejected() -> None:
    certificate = compact_335_triangle_certificate()
    certificate["coxeter_relations"][0][2] = 4
    certificate["coxeter_relations"][2][0] = 4
    with pytest.raises(CertificateError, match="wrong Coxeter relation"):
        verify_chamber_geometry(certificate)


def test_missing_vertex_is_rejected() -> None:
    certificate = compact_335_triangle_certificate()
    certificate["geometry"]["vertices"].pop()
    with pytest.raises(CertificateError, match="declared vertices"):
        verify_chamber_geometry(certificate)


def test_nonstrict_interior_point_is_rejected() -> None:
    certificate = compact_335_triangle_certificate()
    certificate["geometry"]["interior_point"] = encoded_vector([1, 1, 1])
    with pytest.raises(CertificateError, match="every wall strictly"):
        verify_chamber_geometry(certificate)


def test_nonintegral_root_is_rejected() -> None:
    certificate = compact_335_triangle_certificate()
    certificate["simple_roots"][0] = encoded_vector([Fraction(1, 2), 0, 0])
    certificate["wall_gram"][0] = encoded_vector(
        [Fraction(1, 2), Fraction(-1, 2), -PHI / 2]
    )
    certificate["wall_gram"][1][0] = encoded(Fraction(-1, 2))
    certificate["wall_gram"][2][0] = encoded(-PHI / 2)
    with pytest.raises(CertificateError, match=r"not in Z\[phi\]"):
        verify_chamber_geometry(certificate)


def test_wrong_wall_gram_is_rejected_before_geometry() -> None:
    certificate = deepcopy(compact_335_triangle_certificate())
    certificate["wall_gram"][0][1] = encoded(-2)
    certificate["wall_gram"][1][0] = encoded(-2)
    with pytest.raises(CertificateError, match="wall_gram does not equal"):
        verify_chamber_geometry(certificate)


@pytest.mark.parametrize(
    ("mutation", "message"),
    [
        (lambda certificate: certificate.__setitem__("schema_version", "999"), "schema_version"),
        (lambda certificate: certificate.pop("completeness"), "missing"),
        (lambda certificate: certificate.__setitem__("unexpected", True), "unexpected"),
        (lambda certificate: certificate["geometry"].pop("nonempty"), "missing"),
        (
            lambda certificate: certificate["geometry"].__setitem__("nonempty", False),
            "must be true",
        ),
    ],
)
def test_closed_certificate_schema_is_enforced(mutation, message: str) -> None:
    certificate = compact_335_triangle_certificate()
    mutation(certificate)
    with pytest.raises(CertificateError, match=message):
        verify_chamber_geometry(certificate)
