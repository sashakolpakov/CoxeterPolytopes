import hashlib
import json
from pathlib import Path

from scripts.certify_borcherds_char2_3b_h11_22 import _order, _permutation


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data/enumerations/borcherds_char2_3b_h11_22_certificate.json"
)


def _load() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_archived_h11_22_decision_and_vinberg_counts() -> None:
    record = _load()
    assert record["all_checks_pass"] is True
    assert record["decision"] == {
        "abstract_novelty": "known: Scharlau--Kirschmer List11 entry 27",
        "algebro_geometric_novelty": (
            "priority-open explicit K3 order-3 fixed-divisor realization"
        ),
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "compact": False,
        "facet_count": 22,
        "finite_volume": True,
        "hyperbolic_dimension": 11,
    }
    assert record["vinberg"]["elliptic_edge_count"] == 14580
    assert record["vinberg"]["finite_vertex_count"] == 2304
    assert record["vinberg"]["ideal_vertex_count"] == 72
    assert record["vinberg"]["edge_endpoint_type_histogram"] == {
        "0,2": 201,
        "1,1": 3414,
        "2,0": 10965,
    }
    assert record["vinberg"]["bad_edge_count"] == 0
    assert record["vinberg"]["bad_parabolic_neighborhood_count"] == 0


def test_archived_fixed_lattice_and_facets() -> None:
    record = _load()
    fixed = record["fixed_lattice"]
    assert fixed["determinant"] == -243
    assert fixed["smith_diagonal_absolute"] == [1] * 7 + [3] * 5
    assert fixed["signature_positive_negative"] == [1, 11]
    chamber = record["restricted_chamber"]
    assert chamber["primitive_facet_type_histogram"] == {
        "orbit1_square-2_div1": 12,
        "orbit3_square-6_div3": 10,
    }
    assert chamber["coxeter_edge_counts"] == {
        "2": 186,
        "3": 35,
        "4": 0,
        "6": 10,
        "infinity": 0,
        "ultraparallel": 0,
    }
    encoded = json.dumps(
        chamber["facet_gram"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert hashlib.sha256(encoded).hexdigest() == chamber["facet_gram_sha256"]
    assert chamber["facet_gram_sha256"] == (
        "2944c4ce6ebab9a5fb1b52c9a0f7c2a105656863036c3e2badbe52b15c9763cb"
    )


def test_archived_3b_representative_has_order_three() -> None:
    permutation = _permutation()
    assert len(permutation) == 210
    assert _order(permutation) == 3

