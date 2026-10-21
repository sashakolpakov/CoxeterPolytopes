import json

from scripts.certify_coble_type_i_two_boundary_full_ns_h11 import (
    OUTPUT,
    build_certificate,
)


def test_two_boundary_type_i_coble_full_nef_h11_archive_replays() -> None:
    rebuilt = json.loads(json.dumps(build_certificate(), sort_keys=True))
    archived = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert rebuilt == archived
    assert archived["all_checks_pass"]


def test_h11_26_volume_arithmetic_and_direction_are_exact() -> None:
    certificate = json.loads(OUTPUT.read_text(encoding="utf-8"))
    assert certificate["direction"] == "SURFACE_TO_POLYTOPE_only"
    assert certificate["classification"] == (
        "finite_volume_noncompact_full_nef_H11_26_arithmetic_Coxeter_polytope"
    )
    cone = certificate["exact_double_description"]
    assert (cone["finite_vertex_count"], cone["ideal_vertex_count"], cone["hyperideal_vertex_count"]) == (508, 42, 0)
    literal = certificate["literal_vinberg"]
    assert literal["elliptic_edge_count"] == 3641
    assert literal["endpoint_type_histogram"] == {"FF": 2144, "FI": 1300, "II": 197}
    assert literal["bad_edge_count"] == literal["bad_parabolic_neighborhood_count"] == 0
    arithmetic = certificate["arithmetic"]
    assert arithmetic["all_26_reflections_integral_on_arithmetic_lattice"]
    assert arithmetic["NS_reflection_integrality"][:2] == [False, False]
    assert all(arithmetic["NS_reflection_integrality"][2:])
    assert certificate["lineage"]["worldwide_literal_priority"] == "open"
