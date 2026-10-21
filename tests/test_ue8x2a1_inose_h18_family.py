import json
from pathlib import Path

from scripts.certify_ue8x2a1_inose_h18_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ue8x2a1_inose_h18_family.json"


def load() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == load()


def test_equation_open_and_fibers() -> None:
    payload = load()
    family = payload["family"]
    assert family["affine_equation"] == (
        "y^2=x^3-3*s^2*t^4*x+t^5*(t^2-2*(1-s^3)*t+1)"
    )
    assert family["parameter_open"] == "s^3*(s^3-1) != 0"
    assert [fiber["kodaira"] for fiber in family["fibers"]] == [
        "II*",
        "II*",
        "I2",
        "I1",
    ]
    assert family["fiber_euler_sum"] == 24
    assert family["open_nonempty_witness"] == {
        "residual_discriminant": 896,
        "residual_quadratic": "t**2 + 30*t + 1",
        "s": 2,
    }


def test_exact_very_general_ns() -> None:
    ns = load()["neron_severi"]
    assert ns["very_general_exact_lattice"] == "U+E8(-1)^2+A1(-1)"
    assert ns["rank_signature_determinant"] == [19, [1, 18, 0], 2]
    assert ns["parameter_recovery"] == "s^3=([t]Q+2)/4"
    assert "very generally rho=19" in ns["generic_exactness"]
    assert "Mordell--Weil rank and torsion zero" in ns[
        "saturation_and_mordell_weil"
    ]


def test_complete_24_wall_handoff() -> None:
    payload = load()
    chamber = payload["complete_nef_handoff"]
    assert payload["positive_object_id"] == "R19-UE8x2A1-NEF-H18-24"
    assert chamber["wall_count"] == 24
    assert chamber["wall_gram_rank_inertia"] == [19, [1, 18, 5]]
    assert chamber["full_rank_basis_determinant"] == 2
    assert chamber["wall_analysis"]["all_checks_pass"] is True


def test_literal_vinberg_gate() -> None:
    payload = load()
    vinberg = payload["vinberg_handoff"]
    assert payload["all_checks_pass"] is True
    assert vinberg["classification"] == "finite_volume_noncompact"
    assert (vinberg["hyperbolic_dimension"], vinberg["facets"]) == (18, 24)
    assert (vinberg["finite_vertices"], vinberg["ideal_vertices"]) == (525, 10)
    assert vinberg["elliptic_edges"] == 5160
    assert vinberg["bad_edges"] == 0
    assert vinberg["connected_affine_neighborhoods"] == 22
    assert vinberg["bad_full_parabolic_neighborhoods"] == 0
    assert vinberg["simple_roots_complete"] is True
    assert vinberg["canonical_full_audit_sha256"] == (
        "dfead85c337bc4e6e58a2925a4f905949c73f0e271c559333fd48adf6d0715be"
    )

