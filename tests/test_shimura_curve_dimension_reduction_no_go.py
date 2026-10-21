from __future__ import annotations

import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/certify_shimura_curve_dimension_reduction_no_go.py"
REPORT = ROOT / "results/shimura_curve_dimension_reduction_no_go.md"


def load_module():
    spec = importlib.util.spec_from_file_location("shimura_curve_dimension_no_go", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_all_exact_checks_pass() -> None:
    cert = load_module().build_certificate()
    assert cert["all_checks_pass"] is True
    assert all(cert["checks"].values())


def test_absolute_maximal_tetrahedron_kills_three_natural_extractions() -> None:
    cert = load_module().build_certificate()
    h6 = cert["maximal_arithmetic_tetrahedron_control"]
    assert h6["three_node_principal_determinants"] == {
        "012": "3 - sqrt(5)",
        "013": "3 - sqrt(5)",
        "023": "4",
        "123": "4",
    }
    assert all(not facet["is_coxeter_polygon"] for facet in h6["facets"].values())
    topology = h6["reflection_centralizer_topology"]
    assert topology["odd_graph_first_betti_number"] == 1
    assert topology["nonreflection_group"] == "Z"
    assert topology["oriented_double_compactified_genus"] == 1


def test_same_curve_occurs_in_dimensions_four_and_five() -> None:
    cert = load_module().build_certificate()
    shared = cert["same_curve_different_dimension_control"]
    assert shared["coxeter_triangle"] == [2, 3, 10]
    assert shared["ground_field"] == "Q(sqrt(5))"
    assert shared["real_DM_orbifold_curve"] == "P1_(2,3,10) with chosen real half"
    assert "five-prism" in shared["same_face_in_H5_control"]


def test_report_states_the_exact_scope_and_sources() -> None:
    raw = REPORT.read_bytes()
    assert all(byte not in raw for byte in (0, 8, 11, 13))
    text = raw.decode("utf-8")
    assert "conditional face-curve theorem" in text
    assert "Felikson" in text
    assert "Long--Maclachlan--Reid" in text
    assert "genus one" in text
    assert "does not extend" in text
    assert "certify_shimura_curve_dimension_reduction_no_go.py" in text
    assert r"\(H_6\)" in text
    assert r"\mathbb H^n" in text
    assert r"\Gamma_\Omega\simeq\pi_1(C_4)" in text
