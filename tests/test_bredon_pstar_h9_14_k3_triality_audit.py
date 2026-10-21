from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_bredon_pstar_h9_14_k3_triality_audit import OUTPUT, build


ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "results" / "bredon_pstar_h9_14_polytope_first_k3_triality_audit.md"


def _archived() -> dict:
    return json.loads(OUTPUT.read_text(encoding="utf-8"))


def test_polytope_first_certificate_replays_exactly() -> None:
    archived = _archived()
    assert build() == archived
    assert archived["all_checks_pass"]
    assert archived["pstar"]["dimension"] == 9
    assert archived["pstar"]["facets"] == 14
    assert archived["finite_volume"]["finite_ideal_hyperideal"] == [128, 6, 0]
    assert archived["finite_volume"]["bad_edges"] == 0
    assert archived["finite_volume"]["bad_parabolic_neighborhoods"] == 0


def test_integral_branches_and_triality_are_discriminating() -> None:
    archived = _archived()
    admissible = archived["integral_orientations"]["k3_admissible"]
    reversed_branch = archived["integral_orientations"]["reversed"]
    unfolding = archived["triality_unfolding"]

    assert admissible["rank_inertia"] == [10, [9, 1, 4]]
    assert admissible["smith_diagonal_absolute"] == [1] * 6 + [3] * 4 + [0] * 4
    assert admissible["walls_span_fixed_lattice_index"] == 1
    assert reversed_branch["rank_inertia"] == [10, [9, 1, 4]]
    assert reversed_branch["smith_diagonal_absolute"] == [1] * 4 + [3] * 6 + [0] * 4

    assert unfolding["folded_gram_equals_integral_pstar_gram"]
    assert unfolding["canonical_rank_inertia"] == [18, [17, 1, 6]]
    assert unfolding["unfolded_roots_span_ambient_lattice_index"] == 1
    assert unfolding["source_gram_equals_negative_unfolded_gram"]
    assert unfolding["triality_action_equals_archived_C07_action"]
    assert unfolding["reversed_rank_inertia"] == [28, [25, 3, 4]]


def test_exact_return_uniform_obstruction_and_inverse_score() -> None:
    archived = _archived()
    surface = archived["k3_and_surface_return"]
    uniformity = archived["adeg24_uniformity"]
    score = archived["inverse_score"]

    assert surface["exact_picard_lattice"]
    assert surface["coinvariant_minus_two_rootless"]
    assert "Nef(X) intersect Fix(g)" in surface["exact_invariant_nef_return"]
    assert uniformity["pstar"]["passes"]
    assert not uniformity["bredon_H2_triangle_with_one_m6"]["passes"]
    assert not uniformity["bredon_H2_triangle_with_three_m6"]["passes"]
    assert not uniformity["whole_24_member_uniform_simply_laced_unfolding"]
    assert score["single_Pstar"]["level"] == "L3 one-datum exact-return calibration, not L4"
    assert score["single_Pstar"]["Pushback"].startswith("FAIL")
    assert score["whole_ADEG24"]["Car"].startswith("FAIL")

    report = REPORT.read_text(encoding="utf-8")
    assert "polytope -> surface -> same polytope" in report
    assert "not an `L4` inverse theorem" in report
    assert "finite-volume and\nnoncompact" in report
