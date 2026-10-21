from __future__ import annotations

import json
from pathlib import Path

from scripts.certify_a20_c2a_c2b_k3_automorphism_lifts import build_certificate


ROOT = Path(__file__).resolve().parents[1]
CERTIFICATE = (
    ROOT
    / "data"
    / "enumerations"
    / "a20_c2a_c2b_k3_automorphism_lifts.json"
)
REPORT = ROOT / "results" / "a20_c2a_c2b_k3_automorphism_lifts.md"
LEDGER = ROOT / "results" / "higher_dimensional_polytope_lineage_ledger.md"


def _archived() -> dict:
    return json.loads(CERTIFICATE.read_text(encoding="utf-8"))


def test_certificate_replays_exactly() -> None:
    archived = _archived()
    assert build_certificate() == archived
    assert archived["all_checks_pass"]
    assert archived["certificate_payload_sha256"] == (
        "be9baceabfeb2c12364fc546815d62b09f714580e998685fac6214a8286388c1"
    )


def test_common_k3_lattice_period_and_discriminant_glue() -> None:
    archived = _archived()
    assert archived["ambient_picard_lattice"] == {
        "determinant": -4,
        "gram_sha256": "9f734cb8358ff8733df90957b71a4cec12c7c08c74dc08ad0b38b1ac47032c0c",
        "invariants": [20, 2, 1],
        "rank": 20,
        "signature": [1, 19],
        "source_wall_count": 50,
    }
    glue = archived["discriminant_gluing"]
    assert glue["transcendental_lattice"] == "<2>+<2>"
    assert glue["both_actions_equivariant"]
    assert len(glue["quadratic_form_table"]) == 4
    assert all(
        (row["4q_N_mod_8"] + row["4q_T_mod_8"]) % 8 == 0
        for row in glue["quadratic_form_table"]
    )
    k3 = archived["common_marked_k3_lattice"]
    assert k3["signature"] == [3, 19]
    assert k3["unimodular_determinant"] == -1
    assert k3["common_period"] == ["1", "i"]
    assert k3["resulting_picard_lattice_is_exactly_N"]


def test_two_automorphisms_and_exact_fixed_ns() -> None:
    archived = _archived()
    cases = {row["id"]: row for row in archived["cases"]}
    c01 = cases["A20_2_1_C01"]
    c02 = cases["A20_2_1_C02"]

    assert c01["k3_action_order"] == 4
    assert c01["period_eigenvalue"] == "i"
    assert c01["k3_action_characteristic_polynomial"] == (
        "(x-1)^13*(x+1)^7*(x^2+1)"
    )
    assert c01["fixed_neron_severi_rank"] == 13
    assert c01["fixed_neron_severi_is_archived_fixed_lattice"]
    assert c01["primitive_invariant_polarization_square"] == "32068"

    assert c02["k3_action_order"] == 2
    assert c02["period_eigenvalue"] == "-1"
    assert c02["k3_action_characteristic_polynomial"] == (
        "(x-1)^12*(x+1)^10"
    )
    assert c02["fixed_neron_severi_rank"] == 12
    assert c02["fixed_neron_severi_is_archived_fixed_lattice"]
    assert c02["primitive_invariant_polarization_square"] == "5386"


def test_literal_vinberg_volume_is_complete() -> None:
    cases = {row["id"]: row for row in _archived()["cases"]}
    assert cases["A20_2_1_C01"]["hyperbolic_fixed_section"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "compact": False,
        "dimension": 12,
        "elliptic_edge_count": 3846,
        "facet_count": 19,
        "finite_vertex_count": 527,
        "finite_volume": True,
        "ideal_vertex_count": 25,
    }
    assert cases["A20_2_1_C02"]["hyperbolic_fixed_section"] == {
        "bad_edge_count": 0,
        "bad_parabolic_neighborhood_count": 0,
        "compact": False,
        "dimension": 11,
        "elliptic_edge_count": 2229,
        "facet_count": 20,
        "finite_vertex_count": 311,
        "finite_volume": True,
        "ideal_vertex_count": 29,
    }


def test_report_and_ledger_preserve_scope() -> None:
    report = " ".join(REPORT.read_text(encoding="utf-8").lower().split())
    ledger = LEDGER.read_text(encoding="utf-8")
    assert "non-symplectic order four" in report
    assert "non-symplectic involution" in report
    assert "not claim" in report and "abstract" in report
    assert "50-wall chamber itself is the literal nef cone" in report
    assert "A20-C2a-H12-19" in ledger
    assert "A20-C2b-H11-20" in ledger
    assert "genuine" in ledger
