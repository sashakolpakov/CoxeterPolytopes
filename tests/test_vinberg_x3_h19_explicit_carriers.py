import json
from pathlib import Path

from scripts.certify_vinberg_x3_h19_explicit_carriers import build_certificate


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = (
    ROOT / "data/enumerations/vinberg_x3_h19_explicit_carriers.json"
)


def load() -> dict:
    return json.loads(ARTIFACT.read_text(encoding="utf-8"))


def test_archive_rebuilds_exactly() -> None:
    assert build_certificate() == load()


def test_two_explicit_algebraic_carriers() -> None:
    payload = load()
    theorem = payload["positive_theorem"]
    assert theorem["explicit_carrier_count_certified_here"] == 2
    assert theorem["neron_severi"] == "U+2E8(-1)+A2(-1)"
    assert payload["inose_carrier"]["affine_equation"] == (
        "y^2=x^3+t^5*(t-1)^2"
    )
    assert payload["abelian_quotient_carrier"][
        "weierstrass_equation_after_constant_complex_rescaling"
    ] == "y^2=x^3+(t^3-1)^4"


def test_inose_and_quotient_lattice_exactness() -> None:
    payload = load()
    inose = payload["inose_carrier"]
    quotient = payload["abelian_quotient_carrier"]
    assert [fiber["kodaira"] for fiber in inose["fibers"]] == [
        "II*",
        "IV",
        "II*",
    ]
    assert inose["trivial_lattice"] == "U+2E8(-1)+A2(-1)"
    assert inose["mordell_weil"] == "0"
    assert quotient["quotient_singularities"] == (
        "9 rational double points of type A2"
    )
    assert quotient["primitive_closure_index"] == 3
    assert quotient["identified_ns"] == "U+2E8(-1)+A2(-1)"
    assert payload["lattice_identification"][
        "rank_signature_determinant"
    ] == [20, [1, 19, 0], -3]


def test_mixed_divisor_wall_interpretation_is_scoped() -> None:
    geometry = load()["algebro_geometric_interpretation"]
    assert "36 primitive wall normals" in geometry["divisor_carrier"]
    assert "24 square-minus-two" in geometry["long_walls"]
    assert "12 square-minus-six" in geometry["short_walls"]
    assert "not asserted to be irreducible curves" in geometry["short_walls"]
    assert "contained in the nef cone" in geometry["nef_cone_relation"]


def test_complete_h19_vinberg_handoff() -> None:
    payload = load()
    hyperbolic = payload["hyperbolic_handoff"]
    assert payload["all_checks_pass"] is True
    assert (hyperbolic["dimension"], hyperbolic["facets"]) == (19, 36)
    assert hyperbolic["root_norm_counts"] == {
        "-2_divisibility_1": 24,
        "-6_divisibility_3": 12,
    }
    assert (hyperbolic["finite_vertices"], hyperbolic["cusps"]) == (8346, 107)
    assert hyperbolic["elliptic_rank18_edges"] == 87046
    assert hyperbolic["bad_edges"] == 0
    assert hyperbolic["finite_volume"] is True
    assert hyperbolic["compact"] is False

