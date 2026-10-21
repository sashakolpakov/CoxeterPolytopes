import hashlib
import json
from pathlib import Path

import sympy as sp

from coxeter_polytopes.known_examples.bugaenko import bugaenko_simple_roots
from scripts.bugaenko_k3_eigen_nef import multiplication_by_phi
from scripts.bugaenko_trace_lattice import (
    basis_change_certificate,
    expected_standard_gram,
)
from scripts.certify_rational_surface_eigencarrier_universality import (
    _trace_coordinate_vector,
)


ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "data/enumerations/six_carrier_integral_packet_rank_calibration.json"


def load(path):
    return json.loads((ROOT / path).read_text(encoding="utf-8"))


def max_abs(values):
    if isinstance(values, dict):
        values = list(values.values())
    flattened = []

    def visit(value):
        if isinstance(value, list):
            for item in value:
                visit(item)
        else:
            flattened.append(abs(int(value)))

    visit(values)
    return max(flattened)


def pairing_degrees(gram, ample, walls):
    if isinstance(walls, dict):
        walls = list(walls.values())
    return [
        sum(
            ample[i] * gram[i][j] * wall[j]
            for i in range(len(ample))
            for j in range(len(ample))
        )
        for wall in walls
    ]


def test_source_hashes_are_pinned():
    ledger = load("data/enumerations/six_carrier_integral_packet_rank_calibration.json")
    for relative, expected in ledger["source_pins"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == expected, relative


def test_trace_rank_and_surface_rank_table():
    rows = load("data/enumerations/six_carrier_integral_packet_rank_calibration.json")[
        "rows"
    ]
    assert len(rows) == 6
    assert [row["trace"]["rank_d_times_n_plus_1"] for row in rows] == [
        row["d"] * (row["n"] + 1) for row in rows
    ]
    assert [row["carrier"]["rho"] for row in rows] == [20, 7, 11, 12, 15, 18]
    assert [row["carrier"]["b2"] for row in rows] == [20, 7, 11, 12, 15, 18]
    assert all(not row["trace"]["external_only"] for row in rows)


def test_source_wall_spans_and_dimensions():
    toric = load("data/enumerations/toric_full_nef_coxeter_factory.json")
    assert toric["h6_endpoint"]["wall_gram_rank"] == 7
    assert toric["h6_endpoint"]["wall_span_smith_nonzero"] == [1] * 7

    h10 = load("data/enumerations/coble_one_boundary_full_ns_h10.json")["cases"][
        "II"
    ]
    assert h10["arithmetic"]["wall_rank"] == 11
    assert h10["arithmetic"]["roots_span_whole_lattice"]

    h11 = load("data/enumerations/coble_type_i_two_boundary_full_ns_h11.json")
    assert sp.Matrix(h11["facets"]["NS_wall_vectors"]).rank() == 12
    assert h11["arithmetic"]["walls_span_whole_arithmetic_lattice"]

    h14 = load("data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json")
    assert h14["full_NS_lattice"]["rank"] == 15
    assert sp.Matrix(h14["facets"]["NS_wall_vectors"]).rank() == 15

    h17 = load("data/enumerations/coble_type_mii_full_ns_h17.json")
    assert h17["full_NS"]["visible_curves_span_NS_with_index"] == 1
    assert sp.Matrix(h17["facets"]["physical_NS_wall_vectors"]).rank() == 18


def test_two_certified_total_ample_facet_degrees():
    toric = load("data/enumerations/toric_full_nef_coxeter_factory.json")[
        "h6_endpoint"
    ]
    toric_degrees = pairing_degrees(
        toric["intersection_matrix"], toric["strict_interior"], toric["walls"]
    )
    assert sorted(toric_degrees) == [1] * 6 + [2] * 3
    assert sum(toric_degrees) == 12

    h17 = load("data/enumerations/coble_type_mii_full_ns_h17.json")
    h17_degrees = pairing_degrees(
        h17["full_NS"]["gram"],
        h17["facets"]["strict_interior"],
        h17["facets"]["physical_NS_wall_vectors"],
    )
    assert sorted(h17_degrees) == [1] * 24 + [6] * 24
    assert sum(h17_degrees) == 168


def test_six_exact_packet_coordinate_ceilings():
    change = sp.Matrix(basis_change_certificate())
    standard_action = change.inv() * sp.Matrix(multiplication_by_phi()) * change
    standard_walls = [
        change.inv() * _trace_coordinate_vector(root)
        for root in bugaenko_simple_roots()
    ]
    bugaenko_entries = list(sp.Matrix(expected_standard_gram())) + list(
        2 * standard_action
    )
    for wall in standard_walls:
        bugaenko_entries.extend(list(wall))
    assert all(entry.is_Integer for entry in bugaenko_entries)
    assert max(abs(int(entry)) for entry in bugaenko_entries) == 48

    toric = load("data/enumerations/toric_full_nef_coxeter_factory.json")
    h10 = load("data/enumerations/coble_one_boundary_full_ns_h10.json")["cases"][
        "II"
    ]
    h11 = load("data/enumerations/coble_type_i_two_boundary_full_ns_h11.json")
    h14 = load("data/enumerations/coble_type_vi_char3_five_boundary_full_ns_h14.json")
    h17 = load("data/enumerations/coble_type_mii_full_ns_h17.json")
    ceilings = [
        48,
        max(
            max_abs(toric["h6_endpoint"]["intersection_matrix"]),
            max_abs(toric["h6_endpoint"]["walls"]),
        ),
        max(max_abs(h10["lattice_gram"]), max_abs(h10["wall_vectors"])),
        max(
            max_abs(h11["full_NS_lattice"]["gram"]),
            max_abs(h11["facets"]["NS_wall_vectors"]),
        ),
        max(
            max_abs(h14["full_NS_lattice"]["gram"]),
            max_abs(h14["facets"]["NS_wall_vectors"]),
        ),
        max(
            max_abs(h17["full_NS"]["gram"]),
            max_abs(h17["facets"]["physical_NS_wall_vectors"]),
        ),
    ]
    assert ceilings == [48, 1, 17, 21, 8, 1]
    ledger_ceilings = [
        row["complexity"]["packet_coordinate_ceiling_C"]
        for row in load(
            "data/enumerations/six_carrier_integral_packet_rank_calibration.json"
        )["rows"]
    ]
    assert ledger_ceilings == ceilings


def test_interval_generators_and_no_manufactured_monodromy():
    rows = load("data/enumerations/six_carrier_integral_packet_rank_calibration.json")[
        "rows"
    ]
    sqrt5 = sp.sqrt(5)
    assert sp.simplify((1 + sqrt5) ** 2 - 2 * (1 + sqrt5) - 4) == 0
    assert (1 + sqrt5 >= 2) is sp.true
    assert (1 + sqrt5 <= 14) is sp.true
    assert (1 - sqrt5 >= -2) is sp.true
    assert (1 - sqrt5 <= 1) is sp.true
    assert [row["field_gate"]["t"] for row in rows] == [1, 0, 0, 0, 0, 0]
    assert all(row["complexity"]["monodromy_rank"] is None for row in rows)
    assert [
        row["complexity"]["total_ample_facet_degree_D"] for row in rows
    ] == [None, 12, None, None, None, 168]
