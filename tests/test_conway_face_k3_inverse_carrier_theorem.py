import json
from pathlib import Path

from scripts.certify_conway_face_k3_inverse_carrier_theorem import (
    _forms_are_isometric,
    _root_gram,
    _smith_invariants,
    _standard_endpoint_picard,
    _sum_gram,
    build,
    verify,
)


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data/enumerations/conway_face_k3_inverse_carrier_theorem.json"
)


def test_rank20_endpoint_discriminant_anti_isometries() -> None:
    assert _forms_are_isometric(
        _root_gram("E6"), _root_gram("A2"), left_sign=-1
    )
    assert _forms_are_isometric(
        _root_gram("D6"),
        _sum_gram(("A1", "A1")),
        left_sign=-1,
    )


def test_rank20_endpoints_are_vinbergs_standard_picard_lattices() -> None:
    e6_endpoint = _standard_endpoint_picard("E6")
    d6_endpoint = _standard_endpoint_picard("D6")
    assert e6_endpoint.shape == (20, 20)
    assert abs(int(e6_endpoint.det())) == 3
    assert _smith_invariants(e6_endpoint) == [*([1] * 19), 3]
    assert d6_endpoint.shape == (20, 20)
    assert abs(int(d6_endpoint.det())) == 4
    assert _smith_invariants(d6_endpoint) == [*([1] * 18), 2, 2]


def test_complete_allcock_face_table_passes_embedding_gates() -> None:
    payload = build()
    assert payload["verified"] is True
    assert len(payload["face_types"]) == 21
    assert {
        row["hyperbolic_dimension"] for row in payload["face_types"]
    } == set(range(2, 20))
    assert all(row["dimension_sum"] == 19 for row in payload["face_types"])


def test_frozen_conway_face_k3_certificate() -> None:
    archived = json.loads(CERTIFICATE.read_text(encoding="utf-8"))
    summary = verify(archived)
    assert summary["face_type_count"] == 21
    assert summary["dimension_range"] == [2, 19]
    assert summary["all_dimension_sums"] == [19]
