#!/usr/bin/env python3
"""Build the exact arithmetic Coxeter--dessin calibration ledger.

The producer is deliberately a small extractor, not a new geometry engine.
Every numerical entry is replayed from a pinned frozen artifact.  Missing
determinant-line, Heegner, or Borcherds computations remain ``None``/``open``.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = (
    ROOT / "data/enumerations/arithmetic_coxeter_dessin_calibration_ledger.json"
)
RSS_CEILING_BYTES = 1 << 29  # 0.5 GiB

PINNED_SHA256 = {
    "data/known_examples/bugaenko_compact_h8.json":
        "c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c",
    "data/enumerations/bugaenko_h8_determinant_line_ag.json":
        "dce7d52e301b5e215a6f4c58f55ea9fca9114b77c2b84fc929c0c36ac21e6aad",
    "data/candidates/bugaenko_h8_whittaker_gate.json":
        "1ce2771a9d2909c6a5a165b8b86341fcf049119b37676c02cd9c15cd5aaf5a06",
    "data/candidates/bugaenko_h8_whittaker_first_compensator.json":
        "3e92bfb2593abeaf1858b4f92f8877f819211f34c0bbb272d36bd783be461371",
    "data/trace_lattices/bugaenko_k3_eigen_nef.json":
        "e77bb81bcc738f4f542676748dbcf09418cab830da579a600fd1e309f03de782",
    "data/enumerations/toric_full_nef_coxeter_factory.json":
        "7bdd361d5d349b10110fa908d65e7508fc17cc653a097afe91be803da064a389",
    "data/enumerations/finite_volume_coxeter_dessin_cusp_transfer.json":
        "97259b2dd70bd48040cd998178a042d0fea675a77c93f511e8c866f78f2ed991",
    "data/enumerations/coble_one_boundary_full_ns_h10.json":
        "1e90fcd5f9b1ba8538937b91d2f26ea401f6952e8484a597624939edaf2a873e",
    "data/enumerations/coble_type_ii_b3_plane_inverse.json":
        "a8f49be994ae9bf70eae73378f4f006417be91ae827fee0a73707932f6ceac71",
    "data/enumerations/coble_type_ii_b3_plane_inverse_hostile_audit.json":
        "a2c357d26b55bcacce3578b30882a53c036cc555b3f9dcd76ac02206e03bcd04",
    "data/enumerations/polytope_first_coble_h10_hostile_audit.json":
        "12c23f9e1818f962175232a55d7497c9d18f63554db7a571539e8e8d6ce9d0f7",
    "data/enumerations/coble_type_i_two_boundary_full_ns_h11.json":
        "bcfa2de1e96bf4f41282345a64f8ee6f3a75e715287bab4962d7082584b4d985",
    "data/enumerations/coble_type_i_two_boundary_plane_inverse.json":
        "6313b11cc9ba267834231aab2469d99e9d0242f4a95f53280047351f8318da04",
    "data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json":
        "4cc26b186f43e15e8fe3ba9899a5a709bd8f9e79bc3ff86f4983c368eeea4b6f",
    "data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json":
        "8247ad85f7970d3bf1f08fe15517b70376409c17cc2784b8f14b129aaecc2ede",
    "data/enumerations/polytope_first_coble_type_vi_h14.json":
        "b1dfafffd87ff281d73c3517f18033ca7d89852e9a3e0dd39a035f890ddf26fb",
    "data/enumerations/coble_type_mii_full_ns_h17.json":
        "6ed6fba936c362e02793375fa3adfea0279d312c853a8bb4bab30c23b4b8c79c",
    "data/enumerations/coble_type_mii_h17_hostile_audit.json":
        "124a299740c6a82c7ab543ce9236c1ad38d6e988e637a0d173cbe2203756c25f",
    "data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json":
        "08d40b33637fe4d332dc2e62a368bee401c0eda385baa3f823dde3285a15dfdd",
    "data/enumerations/roulleau_rank10_two_full_nef_simplices.json":
        "11b937f0a78091452d1a50cb2b4cb554949485b1469b4ad2616ee629e24320f8",
    "data/enumerations/k3_rank11_14_literal_vinberg_volume.json":
        "f97b713d76fbbd7e18d93c018b7055164da7768f9f21355f337e8aff085a81cb",
    "results/borcherds_compact_inverse_leverage_audit.md":
        "baa3ebc265da7b4444035c8a8f659ca04b94b14426ceac2fce2280f8a46d7014",
    "results/polytope_to_automorphic_divisor_transfer_theorem.md":
        "09568e5374f915ee51350f3e8e498a925b7d56da5d71cd63a64936c6aa11a182",
    "results/determinant_line_slice_disjoint_obstruction_span.md":
        "2b2653cf40a525d05a3401186f76b080860149ef587ed7f846d6328276996b77",
}


def _rss_bytes() -> int:
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(peak if sys.platform == "darwin" else peak * 1024)


def _rss_guard() -> None:
    if _rss_bytes() > RSS_CEILING_BYTES:
        raise MemoryError("dessin calibration replay exceeded its 0.5 GiB RSS ceiling")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
            _rss_guard()
    return digest.hexdigest()


def _read(relative: str) -> bytes:
    path = ROOT / relative
    actual = _sha256(path)
    expected = PINNED_SHA256[relative]
    if actual != expected:
        raise AssertionError(f"source SHA changed for {relative}: {actual} != {expected}")
    data = path.read_bytes()
    _rss_guard()
    return data


def _load(relative: str) -> dict[str, Any]:
    return json.loads(_read(relative))


def _pin(relative: str, role: str) -> dict[str, str]:
    return {"path": relative, "sha256": PINNED_SHA256[relative], "role": role}


def _open_automorphic_status(reason: str) -> dict[str, Any]:
    return {
        "arrangement": "not_row_specifically_compiled",
        "pure_borcherds_product": None,
        "compensated_slice_exact_product": None,
        "exact_global_reflective_product": None,
        "status": "open",
        "reason": reason,
    }


def build() -> dict[str, Any]:
    # Compact benchmark.
    bug = _load("data/known_examples/bugaenko_compact_h8.json")
    bug_det = _load("data/enumerations/bugaenko_h8_determinant_line_ag.json")
    bug_whit = _load("data/candidates/bugaenko_h8_whittaker_gate.json")
    bug_first = _load("data/candidates/bugaenko_h8_whittaker_first_compensator.json")
    bug_k3 = _load("data/trace_lattices/bugaenko_k3_eigen_nef.json")
    assert bug["dimension"] == 8 and len(bug["wall_gram"]) == 11
    assert len(bug["geometry"]["vertices"]) == 41 and bug["geometry"]["cusps"] == []
    assert bug_det["all_checks_pass"] and bug_det["source"]["facet_count"] == 11
    assert bug_whit["verified"]
    assert bug_whit["verdict"]["pure_norm_two_Whittaker_principal_part"] == "obstructed"
    assert bug_first["status"] == "verified"
    assert bug_k3["eigen_nef_realization"]["restricted_chamber_equals_bugaenko_chamber"]

    # Rational cusp calibration.
    toric = _load("data/enumerations/toric_full_nef_coxeter_factory.json")
    cusp = _load("data/enumerations/finite_volume_coxeter_dessin_cusp_transfer.json")
    h6 = toric["h6_endpoint"]
    assert len(h6["wall_order"]) == 9
    assert (h6["finite_vertex_count"], h6["ideal_vertex_count"], h6["hyperideal_vertex_count"]) == (18, 3, 0)
    assert h6["decision"] == "finite_volume_noncompact"
    assert cusp["all_checks_pass"] and cusp["calibration"]["id"] == "toric-H6/9"
    assert cusp["calibration"]["finite_ideal_hyperideal"] == [18, 3, 0]

    # Horizontal Coble inverse rows.
    h10 = _load("data/enumerations/coble_one_boundary_full_ns_h10.json")
    h10_inv = _load("data/enumerations/coble_type_ii_b3_plane_inverse.json")
    h10_audit = _load("data/enumerations/coble_type_ii_b3_plane_inverse_hostile_audit.json")
    h10_census = _load("data/enumerations/polytope_first_coble_h10_hostile_audit.json")
    h10_case = next(case for case in h10["cases"].values() if len(case["wall_vectors"]) == 30)
    h10_counts = h10_case["exact_double_description"]
    assert (h10_counts["finite_vertex_count"], h10_counts["ideal_vertex_count"], h10_counts["hyperideal_vertex_count"]) == (1333, 207, 0)
    assert h10_inv["all_checks_pass"] and h10_audit["all_checks_pass"]
    assert h10_inv["same_polytope_return"]["exact_same_marked_H10_30"]
    assert h10_census["summary"]["Kondo_nonempty_complex_carrier_count"] == 1

    h11 = _load("data/enumerations/coble_type_i_two_boundary_full_ns_h11.json")
    h11_inv = _load("data/enumerations/coble_type_i_two_boundary_plane_inverse.json")
    h11_audit = _load("data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json")
    h11_counts = h11["exact_double_description"]
    assert len(h11["facets"]["labels"]) == 26
    assert (h11_counts["finite_vertex_count"], h11_counts["ideal_vertex_count"], h11_counts["hyperideal_vertex_count"]) == (508, 42, 0)
    assert h11_inv["all_checks_pass"] and h11_audit["all_checks_pass"]
    assert h11_audit["full_cone_replay"]["finite_ideal_hyperideal"] == [508, 42, 0]
    assert h11_inv["checks"]["exact_same_P_return_handoff"]

    # Vertical characteristic-three Coble rows.
    h14 = _load("data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json")
    h14_inv = _load("data/enumerations/polytope_first_coble_type_vi_h14.json")
    h17 = _load("data/enumerations/coble_type_mii_full_ns_h17.json")
    h17_audit = _load("data/enumerations/coble_type_mii_h17_hostile_audit.json")
    bridge = _load("data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json")
    h14_counts = h14["exact_double_description"]
    assert h14["literal_vinberg"]["dimension"] == 14
    assert h14["literal_vinberg"]["facet_count"] == 60
    assert (h14_counts["finite_vertex_count"], h14_counts["ideal_vertex_count"], h14_counts["hyperideal_vertex_count"]) == (27845, 1027, 0)
    assert h14_inv["all_checks_pass"] and h14_inv["exact_same_P_return"]["oriented_mirror_sets_equal"]
    assert len(h17["facets"]["labels"]) == 48
    assert h17_audit["literal_full_cone"]["finite_vertex_count"] == 190808
    assert h17_audit["literal_full_cone"]["ideal_vertex_count"] == 5244
    assert h17_audit["literal_full_cone"]["hyperideal_vertex_count"] == 0
    assert bridge["all_checks_pass"]
    assert bridge["logic_classification"]["level"] == "L3 exact return for two individual raw chambers; not L4"

    # Determinant-line positive and negative controls.
    e10_stock = _load("data/enumerations/roulleau_rank10_two_full_nef_simplices.json")
    e10 = next(case for case in e10_stock["cases"] if case["id"] == "R10-UE8-NEF-H9-10")
    ue8d4_stock = _load("data/enumerations/k3_rank11_14_literal_vinberg_volume.json")
    ue8d4 = next(case for case in ue8d4_stock["cases"] if case["key"] == "R14_UE8D4")
    auto_audit_text = _read("results/borcherds_compact_inverse_leverage_audit.md").decode()
    _read("results/polytope_to_automorphic_divisor_transfer_theorem.md")
    _read("results/determinant_line_slice_disjoint_obstruction_span.md")
    assert "weight `165`" in auto_audit_text
    assert "admits no nonconstant holomorphic modular" in auto_audit_text
    assert "form whose zero divisor is supported" in auto_audit_text
    assert (e10["hyperbolic_dimension"], e10["facet_count"]) == (9, 10)
    assert (e10["literal_vinberg"]["finite_vertex_count"], e10["literal_vinberg"]["ideal_vertex_count"]) == (9, 1)
    assert (ue8d4["dimension"], ue8d4["facet_count"]) == (13, 15)
    assert (ue8d4["finite_vertex_count"], ue8d4["ideal_vertex_count"]) == (42, 4)

    rows = [
        {
            "id": "BUGAENKO-H8-11",
            "source_pins": [
                _pin("data/known_examples/bugaenko_compact_h8.json", "complete marked Coxeter chamber"),
                _pin("data/enumerations/bugaenko_h8_determinant_line_ag.json", "determinant-line and Heegner arrangement"),
                _pin("data/candidates/bugaenko_h8_whittaker_gate.json", "pure principal-part obstruction"),
                _pin("data/candidates/bugaenko_h8_whittaker_first_compensator.json", "complete first-shell obstruction"),
                _pin("data/trace_lattices/bugaenko_k3_eigen_nef.json", "existential K3 eigen-nef carrier"),
            ],
            "polytope": {
                "dimension": 8,
                "facets": 11,
                "arithmetic_definition_field": "Q(sqrt(5))",
                "field_lane": "totally_real_compact",
                "finite_volume": True,
                "compactness": "compact",
                "vertices": {"finite": 41, "ideal_vertices_or_cusps": 0, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {
                "status": "proved_row_specific_arrangement",
                "extension": "<2*phi> + E8(-1) + <4*phi>, signature (2,8) at the identity place",
                "facet_heegner_type": "all 11 have (m,mu)=(1,0)",
                "closed_slice_shell": {"facets": 11, "lower_boundary_contacts": 309, "interior_cuts": 0},
                "associated_surface": "exists after neat cover and complete intersection; auxiliary choices are noncanonical",
            },
            "cusp_transfer": {
                "status": "not_applicable_compact",
                "ideal_vertices_or_cusps": 0,
                "surface_boundary_encoding": None,
            },
            "surface_carrier": {
                "directions": ["POLYTOPE_TO_DECORATED_AG_ARRANGEMENT", "K3_SURFACE_TO_POLYTOPE"],
                "status": "arrangement-level P-first carrier plus existential K3 eigen-nef realization",
                "ground_field_or_characteristic": "decorated carrier over a finite extension of the reflex field; K3 statement over C",
                "exact_return_level": "marked arrangement returns the same chamber; native K3 direction is a calibration, not class-wide L4",
            },
            "automorphic_status": {
                "arrangement": "proved",
                "pure_borcherds_product": "obstructed for the standard full norm-two column",
                "first_complete_compensator_shell": "obstructed",
                "compensated_slice_exact_product": "finite complex correction proved abstractly; an effective product is conditional on the stated Hilbert--Weil integral q-structure and is not explicit",
                "exact_global_reflective_product": None,
                "status": "open_beyond_first_shell",
            },
            "honest_missing_gate": "compute an explicit least safe correction and its obstruction matrix/weight/degree, or prove a Hilbert-reflective classification strong enough to return a Coxeter bound",
        },
        {
            "id": "TORIC-H6-9",
            "source_pins": [
                _pin("data/enumerations/toric_full_nef_coxeter_factory.json", "full-nef chamber and scoped inverse classifier"),
                _pin("data/enumerations/finite_volume_coxeter_dessin_cusp_transfer.json", "exact three-cusp transfer calibration"),
            ],
            "polytope": {
                "dimension": 6,
                "facets": 9,
                "arithmetic_definition_field": "Q",
                "field_lane": "rational_cuspidal",
                "finite_volume": True,
                "compactness": "noncompact",
                "vertices": {"finite": 18, "ideal_vertices_or_cusps": 3, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {
                "status": "open_not_row_specifically_compiled",
                "extension": None,
                "facet_heegner_type": None,
                "hidden_sheet_audit": None,
            },
            "cusp_transfer": {
                "status": "proved_row_specific",
                "ideal_vertices_or_cusps": 3,
                "wall_split_at_each_cusp": "7 vertical plus 2 bounded-penetration walls",
                "surface_boundary_encoding": "exceptional boundary cluster and its dual-complex marking after toroidal resolution; no single-point encoding",
            },
            "surface_carrier": {
                "directions": ["SURFACE_TO_POLYTOPE", "POLYTOPE_TO_SURFACE"],
                "status": "smooth toric weak-del-Pezzo full-nef carrier with a finite twelve-row exact recognizer",
                "ground_field_or_characteristic": "Q",
                "exact_return_level": "scoped L3 same-P classifier on the smooth-toric full-nef species, not L4",
            },
            "automorphic_status": _open_automorphic_status("the cusp theorem is proved, but no row-specific determinant support or obstruction matrix has been compiled"),
            "honest_missing_gate": "compile the determinant/Heegner support and obtain a carrier theorem on a broad P-only class rather than the toric-completable locus",
        },
        {
            "id": "COBLE-II-H10-30",
            "source_pins": [
                _pin("data/enumerations/coble_one_boundary_full_ns_h10.json", "complete H10/30 cone"),
                _pin("data/enumerations/coble_type_ii_b3_plane_inverse.json", "strict raw-P plane inverse"),
                _pin("data/enumerations/coble_type_ii_b3_plane_inverse_hostile_audit.json", "independent hostile replay"),
                _pin("data/enumerations/polytope_first_coble_h10_hostile_audit.json", "six-chamber frozen denominator"),
            ],
            "polytope": {
                "dimension": 10,
                "facets": 30,
                "arithmetic_definition_field": "Q",
                "field_lane": "rational_cuspidal",
                "finite_volume": True,
                "compactness": "noncompact",
                "vertices": {"finite": 1333, "ideal_vertices_or_cusps": 207, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {"status": "open_not_row_specifically_compiled", "extension": None, "facet_heegner_type": None, "hidden_sheet_audit": None},
            "cusp_transfer": {
                "status": "eligible_but_not_instantiated",
                "ideal_vertices_or_cusps": 207,
                "surface_boundary_encoding": "exceptional boundary clusters with dual-complex and cusp-wall markings",
                "missing_input": "fixed determinant extension, finite Heegner support, and literal support-wise parabolic replay",
            },
            "surface_carrier": {
                "directions": ["SURFACE_TO_POLYTOPE", "POLYTOPE_TO_SURFACE_TO_SAME_POLYTOPE"],
                "status": "raw chamber forces two B3 plane branches and a Q-defined strong Type-II Coble carrier",
                "ground_field_or_characteristic": "equations over Q; classification and nef return over C",
                "exact_return_level": "L3 strict single raw-P inverse realization, not variable-class L4",
            },
            "automorphic_status": _open_automorphic_status("the native Coble inverse is exact, but its orthogonal determinant dessin has not been compiled"),
            "honest_missing_gate": "prove the contraction/incidence/Hilbert grammar exhaustive on a variable P-only class and then classify its realization schemes",
        },
        {
            "id": "COBLE-I-H11-26",
            "source_pins": [
                _pin("data/enumerations/coble_type_i_two_boundary_full_ns_h11.json", "complete H11/26 cone"),
                _pin("data/enumerations/coble_type_i_two_boundary_plane_inverse.json", "repaired line-quintic inverse"),
                _pin("data/enumerations/coble_type_i_two_boundary_plane_inverse_raw_hostile_audit.json", "independent missing-curve repair and return"),
            ],
            "polytope": {
                "dimension": 11,
                "facets": 26,
                "arithmetic_definition_field": "Q",
                "field_lane": "rational_cuspidal",
                "finite_volume": True,
                "compactness": "noncompact",
                "vertices": {"finite": 508, "ideal_vertices_or_cusps": 42, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {"status": "open_not_row_specifically_compiled", "extension": None, "facet_heegner_type": None, "hidden_sheet_audit": None},
            "cusp_transfer": {
                "status": "eligible_but_not_instantiated",
                "ideal_vertices_or_cusps": 42,
                "surface_boundary_encoding": "exceptional boundary clusters with dual-complex and cusp-wall markings",
                "missing_input": "fixed determinant extension, finite Heegner support, and literal support-wise parabolic replay",
            },
            "surface_carrier": {
                "directions": ["SURFACE_TO_POLYTOPE", "POLYTOPE_TO_SURFACE_TO_SAME_POLYTOPE"],
                "status": "repaired raw-P line-plus-quintic construction realizes all 26 curves and excludes hidden negative curves",
                "ground_field_or_characteristic": "Q-model, interpreted as a complex Type-I Coble surface",
                "exact_return_level": "proved one-input exact return; not a variable-class L4 theorem",
            },
            "automorphic_status": _open_automorphic_status("no row-specific determinant/Heegner or vector-valued obstruction computation exists"),
            "honest_missing_gate": "promote the line-quintic determinant grammar from one chamber to an exhaustive variable family with bounded Hilbert data",
        },
        {
            "id": "COBLE-VI3-H14-60",
            "source_pins": [
                _pin("data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json", "complete H14/60 cone and lineage"),
                _pin("data/enumerations/polytope_first_coble_type_vi_h14.json", "raw boundary/CM/cusp recognizer"),
                _pin("data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json", "P-derived plane incidence and characteristic-three existence"),
            ],
            "polytope": {
                "dimension": 14,
                "facets": 60,
                "arithmetic_definition_field": "Q",
                "field_lane": "rational_cuspidal_with_vertical_char3_native_carrier",
                "finite_volume": True,
                "compactness": "noncompact",
                "vertices": {"finite": 27845, "ideal_vertices_or_cusps": 1027, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {"status": "open_not_row_specifically_compiled", "extension": None, "facet_heegner_type": None, "hidden_sheet_audit": None},
            "cusp_transfer": {
                "status": "eligible_but_not_instantiated",
                "ideal_vertices_or_cusps": 1027,
                "surface_boundary_encoding": "exceptional boundary clusters with dual-complex and cusp-wall markings",
                "missing_input": "fixed determinant extension, finite Heegner support, and literal support-wise parabolic replay",
            },
            "surface_carrier": {
                "directions": ["SURFACE_TO_POLYTOPE", "POLYTOPE_TO_SURFACE_TO_SAME_POLYTOPE"],
                "status": "raw P uniquely recovers the five-boundary package; five plane blowdowns force the conic-plus-four-tangents incidence point in characteristic 3",
                "ground_field_or_characteristic": "algebraically closed field of characteristic 3",
                "exact_return_level": "L3 exact return for an individual raw chamber, not L4",
            },
            "automorphic_status": _open_automorphic_status("the rational Shimura lane and the characteristic-three native carrier are distinct; no bridge between their divisor theories is compiled"),
            "honest_missing_gate": "replace the singled-out template and Kondo no-hidden-curve input by a P-only exhaustive carrier class and classify its characteristic strata",
        },
        {
            "id": "COBLE-MII3-H17-48",
            "source_pins": [
                _pin("data/enumerations/coble_type_mii_full_ns_h17.json", "complete H17/48 marked chamber"),
                _pin("data/enumerations/coble_type_mii_h17_hostile_audit.json", "full cone and curve exhaustion"),
                _pin("data/enumerations/coble_polytope_first_blowdown_incidence_bridge.json", "raw P1xP1 grid compiler and characteristic-three existence"),
            ],
            "polytope": {
                "dimension": 17,
                "facets": 48,
                "arithmetic_definition_field": "Q",
                "field_lane": "rational_cuspidal_with_vertical_char3_native_carrier",
                "finite_volume": True,
                "compactness": "noncompact",
                "vertices": {"finite": 190808, "ideal_vertices_or_cusps": 5244, "hyperideal": 0},
            },
            "determinant_line_and_heegner": {"status": "open_not_row_specifically_compiled", "extension": None, "facet_heegner_type": None, "hidden_sheet_audit": None},
            "cusp_transfer": {
                "status": "eligible_but_not_instantiated",
                "ideal_vertices_or_cusps": 5244,
                "surface_boundary_encoding": "exceptional boundary clusters with dual-complex and cusp-wall markings",
                "missing_input": "fixed determinant extension, finite Heegner support, and literal support-wise parabolic replay",
            },
            "surface_carrier": {
                "directions": ["SURFACE_TO_POLYTOPE", "POLYTOPE_TO_SURFACE_TO_SAME_POLYTOPE"],
                "status": "raw P uniquely recovers the eight-boundary frame, P1xP1 grid, all 24 graph curves, and the characteristic-three incidence point",
                "ground_field_or_characteristic": "algebraically closed field of characteristic 3",
                "exact_return_level": "L3 exact return for an individual raw chamber, not L4",
            },
            "automorphic_status": _open_automorphic_status("no row-specific determinant/Heegner support or obstruction calculation exists"),
            "honest_missing_gate": "turn the unique-frame and projectivity equations into a non-singleton P-only classification with an intrinsic no-hidden-wall theorem",
        },
    ]

    determinant_controls = [
        {
            "id": "E10-H9-10-POSITIVE",
            "source_pins": [
                _pin("data/enumerations/roulleau_rank10_two_full_nef_simplices.json", "complete H9/10 chamber"),
                _pin("results/borcherds_compact_inverse_leverage_audit.md", "determinant-line 2-reflective control theorem"),
            ],
            "polytope": {
                "dimension": 9,
                "facets": 10,
                "arithmetic_definition_field": "Q",
                "compactness": "noncompact",
                "vertices": {"finite": 9, "ideal_vertices_or_cusps": 1, "hyperideal": 0},
            },
            "hyperbolic_lattice": "II_(1,9)=U+E8(-1)",
            "determinant_extension": "II_(1,9)+<2> isometric to 2U+E7(-1)",
            "automorphic_status": {
                "arrangement": "proved",
                "pure_borcherds_product": "proved",
                "compensated_slice_exact_product": "not_needed_for_exact_complete_arrangement",
                "exact_global_reflective_product": "proved complete 2-reflective, multiplicity one, weight 165",
            },
            "surface_carrier_direction": "SURFACE_TO_POLYTOPE full-nef K3 control; no native P-first classifier claimed here",
        },
        {
            "id": "UE8D4-H13-15-NEGATIVE",
            "source_pins": [
                _pin("data/enumerations/k3_rank11_14_literal_vinberg_volume.json", "complete H13/15 chamber"),
                _pin("results/borcherds_compact_inverse_leverage_audit.md", "determinant-line non-2-reflectivity theorem"),
            ],
            "polytope": {
                "dimension": 13,
                "facets": 15,
                "arithmetic_definition_field": "Q",
                "compactness": "noncompact",
                "vertices": {"finite": 42, "ideal_vertices_or_cusps": 4, "hyperideal": 0},
            },
            "hyperbolic_lattice": "U+E8(-1)+D4(-1)",
            "determinant_extension": "U+E8(-1)+D4(-1)+<8>, signature (2,13)",
            "automorphic_status": {
                "arrangement": "available by determinant-line transfer",
                "pure_borcherds_product": None,
                "compensated_slice_exact_product": "general rational interior compensation is available, but no row-specific product is archived",
                "exact_global_reflective_product": "impossible: Wang's complete classification excludes every 2-reflective form in signature (2,13)",
            },
            "surface_carrier_direction": "SURFACE_TO_POLYTOPE full-nef K3 control; determinant exact-reflectivity is negative",
        },
        {
            "id": "BUGAENKO-H8-11-OBSTRUCTED",
            "source_pins": [
                _pin("data/candidates/bugaenko_h8_whittaker_gate.json", "pure obstruction"),
                _pin("data/candidates/bugaenko_h8_whittaker_first_compensator.json", "first-shell obstruction"),
            ],
            "automorphic_status": {
                "arrangement": "proved",
                "pure_borcherds_product": "obstructed on the standard full-group column",
                "first_complete_compensator_shell": "obstructed",
                "compensated_slice_exact_product": "abstract complex correction exists; effective Hilbert product is conditional and not explicit",
                "exact_global_reflective_product": None,
            },
        },
    ]

    synthesis = {
        "general_compiler_support": {
            "direct_determinant_and_heegner_calibration": ["BUGAENKO-H8-11"],
            "direct_cusp_calibration": ["TORIC-H6-9"],
            "native_exact_P_to_surface_experiments": [
                "TORIC-H6-9",
                "COBLE-II-H10-30",
                "COBLE-I-H11-26",
                "COBLE-VI3-H14-60",
                "COBLE-MII3-H17-48",
            ],
            "surface_to_polytope_only_controls": ["E10-H9-10-POSITIVE", "UE8D4-H13-15-NEGATIVE"],
            "interpretation": "The determinant/Heegner and cusp rows test the general decorated compiler.  The Coble and toric inverse rows test stronger native divisor-surface realization, but only in scoped or individual carrier lanes.  They do not prove that a general arithmetic chamber has a native nef-cone surface.",
        },
        "single_missing_classification_invariant": {
            "name": "minimal faithful Coxeter-dessin Hilbert degree T_dess(P)",
            "definition": "the least marked projective degree/stratum-degree package of a faithful decorated surface carrier, including all ordinary face strata, the exceptional boundary clusters and their dual complexes at every cusp, the sheet-separating level, and the marked Heegner curves",
            "why_it_matters": "a P-computable bound in terms of n, [k:Q], discriminant, lattice level, and facet count would place every carrier in finitely many explicit Hilbert/Hurwitz stacks, making surface classification capable of pushing a theorem back to P",
            "current_status": "open; no calibration row computes a class-uniform bound, and the native inverse rows rely on singled-out contraction templates or surface classification",
        },
        "boundary_terminology": "For finite volume, boundary-at-infinity data are ideal vertices/cusps.  A faithful surface encoding uses resolved exceptional boundary clusters and their dual-complex incidence, not one point carrying all cusp-wall curves.",
    }

    artifact = {
        "schema": "coxeter.arithmetic_dessin_calibration_ledger.v1",
        "status_date": "2026-08-09",
        "all_checks_pass": True,
        "scope": {
            "claim": "exact calibration ledger and theorem-facing gap analysis",
            "not_claimed": [
                "a new abstract Coxeter polytope",
                "a general native Neron--Severi or nef-cone carrier theorem",
                "a completed Borcherds product for rows recorded as open",
                "a class-wide L4 classification pushback",
            ],
        },
        "rows": rows,
        "determinant_line_controls": determinant_controls,
        "synthesis": synthesis,
        "resource_contract": {
            "in_process_rss_ceiling_bytes": RSS_CEILING_BYTES,
            "outer_watchdog_command": "python3 scripts/run_with_rss_limit.py --memory-gib 0.5 -- python3 scripts/certify_arithmetic_coxeter_dessin_calibration_ledger.py --check-output data/enumerations/arithmetic_coxeter_dessin_calibration_ledger.json",
        },
    }
    encoded = json.dumps(artifact, sort_keys=True)
    forbidden_boundary_phrase = "ideal " + "faces"
    assert forbidden_boundary_phrase not in encoded.lower()
    assert len(rows) == 6 and len({row["id"] for row in rows}) == 6
    assert all(row["polytope"]["vertices"]["hyperideal"] == 0 for row in rows)
    _rss_guard()
    return artifact


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--check-output", type=Path)
    args = parser.parse_args()
    artifact = build()
    payload = json.dumps(artifact, indent=2, sort_keys=True) + "\n"
    if args.check_output is not None:
        expected = args.check_output
        if not expected.is_absolute():
            expected = ROOT / expected
        if expected.read_text() != payload:
            raise AssertionError(f"archive mismatch: {expected}")
    output = args.output
    if output is not None:
        if not output.is_absolute():
            output = ROOT / output
        output.write_text(payload)
    elif args.check_output is None:
        DEFAULT_OUTPUT.write_text(payload)
    print(
        "arithmetic Coxeter-dessin calibration PASS: "
        f"rows={len(artifact['rows'])}, controls={len(artifact['determinant_line_controls'])}, "
        f"peak_rss_gib={_rss_bytes() / (1 << 30):.3f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
