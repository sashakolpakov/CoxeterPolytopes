import json
from pathlib import Path

from scripts.certify_rational_p17_mmp import build_certificate


def test_p17_surface_mmp_has_exactly_the_diagram_orbit_of_canonical_markings():
    certificate = build_certificate()
    frozen = json.loads(
        Path("data/enumerations/rational_p17_mmp.json").read_text()
    )
    assert certificate == frozen
    assert certificate["rank_two_endpoint_types"] == {
        "F1_nonminimal": 8,
        "F2_minimal": 8,
    }
    assert certificate["rank_one_P2_endpoint_count"] == 8
    assert certificate["canonical_class_count"] == 8
    assert certificate["weighted_diagram_automorphism_order"] == 8
    assert certificate["all_diagram_automorphisms_lift_to_integral_unimodular_isometries"]
    assert certificate["p2_and_f2_reconstructions_agree"]
    assert certificate["canonical_classes_equal_standard_diagram_orbit"]
