import json
from pathlib import Path

from scripts.certify_ue8x2_inose_h17_family import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/ue8x2_inose_h17_family.json"


def load() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == load()


def test_inose_equation_open_and_fibers() -> None:
    payload = load()
    family = payload["family"]
    assert family["affine_equation"] == (
        "y^2 = x^3 - 3*A*t^4*x + t^5*(t^2 - 2*B*t + 1)"
    )
    assert family["parameter_open"] == (
        "A*(A^3-(B+1)^2)*(A^3-(B-1)^2) != 0"
    )
    assert [fiber["kodaira"] for fiber in family["fibers"]] == [
        "II*",
        "I1",
        "II*",
    ]
    assert family["fiber_euler_sum"] == 24
    assert family["open_nonempty_witness"]["residual_discriminant"] == 184320


def test_exact_very_general_neron_severi_lattice() -> None:
    payload = load()
    ns = payload["neron_severi"]
    assert ns["very_general_exact_lattice"] == "U+E8(-1)+E8(-1)"
    assert ns["rank_signature_determinant"] == [18, [1, 17, 0], -1]
    assert ns["generic_exactness"]["parameter_moduli_dimension"] == 2
    assert ns["generic_exactness"]["coefficient_recovery"] == {
        "A^3": "(4*B^2+2-[t^2]R)/4",
        "B": "-[t^3]R/4",
    }
    assert "very general member has rho=18" in ns["generic_exactness"][
        "period_argument"
    ]
    assert "Mordell--Weil rank and torsion zero" in ns["unimodular_closure"]


def test_complete_nef_wall_handoff() -> None:
    payload = load()
    chamber = payload["complete_nef_handoff"]
    assert payload["positive_object_id"] == "R18-UE8x2-NEF-H17-19"
    assert chamber["wall_count"] == 19
    assert chamber["wall_gram_rank_inertia"] == [18, [1, 17, 1]]
    assert chamber["full_rank_basis_determinant"] == -1
    assert chamber["wall_analysis"]["all_checks_pass"] is True


def test_literal_vinberg_finite_volume_gate() -> None:
    payload = load()
    vinberg = payload["vinberg_handoff"]
    assert payload["all_checks_pass"] is True
    assert vinberg["classification"] == "finite_volume_noncompact"
    assert (vinberg["hyperbolic_dimension"], vinberg["facets"]) == (17, 19)
    assert (vinberg["finite_vertices"], vinberg["ideal_vertices"]) == (80, 2)
    assert vinberg["elliptic_edges"] == 729
    assert vinberg["bad_edges"] == 0
    assert vinberg["connected_affine_neighborhoods"] == 3
    assert vinberg["bad_full_parabolic_neighborhoods"] == 0
    assert vinberg["simple_roots_complete"] is True
    assert vinberg["canonical_full_audit_sha256"] == (
        "6e77968ffeb645ad008cc52bcc4995c9a983d928a88d6a4a34e5dc35ad4d0139"
    )

