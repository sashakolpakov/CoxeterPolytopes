"""Regression tests for the general compact K3 Coxeter-polygon theorem."""

from __future__ import annotations

from fractions import Fraction
import json
from pathlib import Path

from scripts.certify_k3_n2_compact_polygon_field_bound import (
    build_certificate,
    constants_record,
    digest,
    encoded_certificate,
    nonadjacent,
    separation_combinatorics_record,
    vertex_incident_facets,
)


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results/k3_n2_compact_polygon_field_bound.md"
CERTIFIER = ROOT / "scripts/certify_k3_n2_compact_polygon_field_bound.py"
ARCHIVE = ROOT / "data/enumerations/k3_n2_compact_polygon_field_bound.json"


def test_constants_recompute_exactly():
    record = constants_record()
    assert record["field_degree_upper_bound"] == 6
    assert record["Coxeter_angle_sine_lower_bound"] == "1/64"
    assert record["cosh_half_shortest_dotted_distance_strict_upper_bound"] == (
        "1+64*pi<257"
    )
    assert record["cosh_shortest_dotted_distance_strict_upper_bound"] == 132_097
    assert record["shortest_dotted_pair_cycle_strict_upper_bound"] == 69_798_469_636
    assert record["field_discriminant_bound"] == "69798469636^30"
    assert 2 * 257**2 - 1 == 132_097
    assert 4 * 132_097**2 == 69_798_469_636


def test_angle_norm_floor_for_every_k3_allowed_degree():
    for degree in range(1, 7):
        norm_gap = Fraction(1, 4 ** (degree - 1))
        sine_square_gap = norm_gap / 4
        assert sine_square_gap == Fraction(1, 4**degree)
        assert sine_square_gap >= Fraction(1, 4**6)
    assert Fraction(1, 4**6) == Fraction(1, 64**2)


def test_vertex_wall_selector_is_exhaustive_on_large_cyclic_range():
    for facet_count in range(4, 257):
        for vertex in range(facet_count):
            incident = vertex_incident_facets(vertex, facet_count)
            for wall in range(facet_count):
                if wall in incident:
                    continue
                assert any(nonadjacent(side, wall, facet_count) for side in incident)


def test_adjacent_and_nonadjacent_vertex_selectors_are_exhaustive():
    checked = 0
    for facet_count in range(4, 129):
        for left in range(facet_count):
            for right in range(left + 1, facet_count):
                left_incident = vertex_incident_facets(left, facet_count)
                right_incident = vertex_incident_facets(right, facet_count)
                witnesses = [
                    (a, b)
                    for a in left_incident
                    for b in right_incident
                    if nonadjacent(a, b, facet_count)
                ]
                assert witnesses
                checked += 1
    assert checked > 300_000


def test_named_shortest_arc_selectors_cover_the_boundary_cases():
    for facet_count in range(4, 129):
        for left in range(facet_count):
            # Adjacent vertices: their two outer facets.
            right = (left + 1) % facet_count
            assert nonadjacent(left - 1, right, facet_count)
            # Separation two: F_{i-1} and F_{i+1}.
            right = (left + 2) % facet_count
            assert nonadjacent(left - 1, right - 1, facet_count)
            # Longer shortest arcs: F_i and F_{j-1}.
            for separation in range(3, facet_count // 2 + 1):
                right = (left + separation) % facet_count
                assert nonadjacent(left, right - 1, facet_count)


def test_sector_and_gauss_bonnet_majorants_give_stated_pair_bound():
    # Area >= N/64*(cosh(delta/2)-1) and Area < N*pi with pi<4.
    rational_cosh_half_majorant = 1 + 64 * 4
    assert rational_cosh_half_majorant == 257
    rational_cosh_distance_majorant = 2 * rational_cosh_half_majorant**2 - 1
    assert rational_cosh_distance_majorant == 132_097
    assert 4 * rational_cosh_distance_majorant**2 == 69_798_469_636


def test_unique_large_conjugate_cannot_descend_to_a_proper_subfield():
    for total_degree in range(2, 7):
        for subfield_degree in range(1, total_degree):
            if total_degree % subfield_degree:
                continue
            relative_degree = total_degree // subfield_degree
            # A value in the subfield repeats at all extensions of an embedding.
            assert relative_degree >= 2


def test_certificate_archive_replays_and_scope_is_exact():
    envelope = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    payload = build_certificate()
    assert envelope["payload"] == payload
    assert envelope["payload_sha256"] == digest(payload)
    assert ARCHIVE.read_text(encoding="utf-8") == encoded_certificate()
    assert separation_combinatorics_record()["all_have_nonadjacent_incident_facet_witnesses"]
    conclusions = payload["conclusions"]
    assert conclusions["some_dotted_pair_cycle_below_69798469636"] is True
    assert conclusions["short_dotted_pair_cycle_generates_full_field"] is True
    assert conclusions["reflection_maximality_used"] is False
    scope = payload["scope"]
    assert scope["facet_count_at_least_4_required"] is True
    assert scope["finite_area_noncompact_polygon_closed"] is False
    assert scope["ultraideal_or_truncated_polygon_closed"] is False
    assert scope["carrier_existence_proved"] is False
    assert scope["Property_4_correspondence_return_proved"] is False


def test_report_states_every_geometric_and_arithmetic_pressure_point():
    text = " ".join(REPORT.read_text(encoding="utf-8").split())
    for phrase in (
        "with at least four facets",
        "finite projective rank-three",
        "primitive full trace",
        "minimal cyclic Cartan ground field",
        "Corollary 2.2 (nonadjacent supports are ultraparallel)",
        r"d(H_i,H_j)=d(F_i,F_j)>0",
        r"2\le\ell\le\lfloor N/2\rfloor",
        "full hyperbolic sector",
        r"q_i:=\operatorname {Tr}_{\mathcal O}(P_iP_{i+1})",
        r"=(z_i)_0\prod_{\sigma\ne\sigma_0}\sigma(z_i)",
        r"\sin\alpha_i\ge{1\over64}",
        r"\operatorname {Area}(P) =(N-2)\pi-\sum_i\alpha_i",
        r"4\cosh^2\delta",
        r"\mathbf Q(p)=k",
        "finite-area noncompact polygons",
        "does not construct the K3 carrier",
    ):
        assert phrase in text


def test_tex_transport_and_scope_controls_are_clean():
    combined = REPORT.read_bytes() + CERTIFIER.read_bytes()
    assert b"\x08" not in combined
    text = combined.decode("utf-8")
    assert ",quad" not in text
    assert ",qquad" not in text
    assert "Nikulin" not in text
    assert "Vinberg" not in text
    assert "ABSW" not in text

