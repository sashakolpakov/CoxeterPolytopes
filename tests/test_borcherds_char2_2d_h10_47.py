import hashlib
import json
from pathlib import Path

from scripts.certify_borcherds_char2_2d_h10_47 import _permutation
from scripts.certify_borcherds_char2_3b_h11_22 import _order


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT / "data/enumerations/borcherds_char2_2d_h10_47_certificate.json"
)


def _load() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_archived_h10_47_decision_and_vinberg_counts() -> None:
    record = _load()
    assert record["all_checks_pass"] is True
    assert record["decision"] == {
        "abstract_novelty": (
            "finite-index subgroup chamber in the classical f_2^10 "
            "commensurability class"
        ),
        "algebro_geometric_novelty": (
            "priority-open explicit K3 involution fixed-divisor realization"
        ),
        "classification": "finite_volume_noncompact_Coxeter_polytope",
        "compact": False,
        "facet_count": 47,
        "finite_volume": True,
        "hyperbolic_dimension": 10,
    }
    assert record["vinberg"]["elliptic_edge_count"] == 19910
    assert record["vinberg"]["finite_vertex_count"] == 3136
    assert record["vinberg"]["ideal_vertex_count"] == 216
    assert record["vinberg"]["edge_endpoint_type_histogram"] == {
        "0,2": 600,
        "1,1": 7260,
        "2,0": 12050,
    }
    assert record["vinberg"]["bad_edge_count"] == 0
    assert record["vinberg"]["bad_parabolic_neighborhood_count"] == 0


def test_archived_fixed_lattice_facets_and_index() -> None:
    record = _load()
    fixed = record["fixed_lattice"]
    assert fixed["determinant"] == 4096
    assert fixed["smith_diagonal_absolute"] == [2] * 10 + [4]
    assert fixed["signature_positive_negative"] == [1, 10]
    chamber = record["restricted_chamber"]
    assert chamber["primitive_facet_type_histogram"] == {
        "orbit1_square-4_div4": 6,
        "orbit2_square-2_div2": 5,
        "orbit2_square-4_div2": 16,
        "orbit2_square-4_div4": 20,
    }
    assert chamber["coxeter_edge_counts"] == {
        "2": 490,
        "3": 30,
        "4": 20,
        "6": 0,
        "infinity": 191,
        "ultraparallel": 350,
    }
    encoded = json.dumps(
        chamber["facet_gram"], sort_keys=True, separators=(",", ":")
    ).encode()
    assert hashlib.sha256(encoded).hexdigest() == chamber["facet_gram_sha256"]
    assert chamber["facet_gram_sha256"] == (
        "454e34ef2358603578a4360145ac6bdd7abfe6e0b539344323ffa5e354038881"
    )
    assert record["lineage"]["reflection_subgroup_index"] == 1566720


def test_archived_2d_representative_has_order_two() -> None:
    permutation = _permutation()
    assert len(permutation) == 210
    assert _order(permutation) == 2

