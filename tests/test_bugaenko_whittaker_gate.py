import json
from pathlib import Path

from scripts.certify_bugaenko_whittaker_gate import (
    archived_chamber_data,
    certify,
    discriminant_residue_counts,
    norm_four_gauss_histograms,
)


def test_bugaenko_determinant_discriminant_module() -> None:
    certificate = certify()
    chamber = archived_chamber_data()
    module = certificate["discriminant_weil_module"]
    trace = certificate["determinant_line_extension"]["trace_certificate"]

    assert trace["standard_gram_verified"] is True
    assert trace["determinant"] == "64"
    assert trace["inertia"] == [18, 2, 0]
    assert module["field_module"] == "O/(2) + O/(4)"
    assert module["order"] == 64
    assert module["isotropic_elements"] == 1
    assert module["even_overlattice_status"] == "none_nontrivial"
    assert sum(discriminant_residue_counts().values()) == 64
    assert chamber["facet_count"] == chamber["facet_sheets"] == 11
    assert chamber["all_facet_norms_are_two"] is True
    assert chamber["lower_face_contact_sheets"] == 309
    assert chamber["closed_chamber_contact_sheets"] == 320
    assert chamber["interior_cut_sheets"] == 0


def test_bugaenko_norm_four_gauss_row_vanishes() -> None:
    expected = {0: 16, 2: 16, 4: 16, 6: 16}
    histograms = norm_four_gauss_histograms()

    assert len(histograms) == 3
    assert all(histogram == expected for histogram in histograms.values())


def test_bugaenko_pure_whittaker_gate_is_obstructed() -> None:
    certificate = certify()

    assert certificate["verified"] is True
    assert certificate["weights"] == {
        "Whittaker_input": [-3, 5],
        "cusp_obstruction": [5, 5],
        "Bessel_order": 4,
    }
    assert certificate["principal_part"]["independent_columns_for_eleven_facets"] == 1
    assert certificate["principal_part"]["facet_sheets"] == 11
    assert certificate["principal_part"]["lower_face_contact_sheets"] == 309
    assert certificate["principal_part"]["interior_cut_sheets"] == 0
    assert certificate["verdict"]["Poincare_series_P_m0_0_nonzero"] is True
    assert certificate["verdict"]["pure_norm_two_Whittaker_principal_part"] == "obstructed"
    assert certificate["verdict"]["larger_compensating_principal_parts"] == "not_decided"

    frozen = json.loads(
        (
            Path(__file__).resolve().parents[1]
            / "data/candidates/bugaenko_h8_whittaker_gate.json"
        ).read_text(encoding="utf-8")
    )
    assert frozen["verified"] is True
    assert frozen["identity_term_certificate"]["strictly_below_one"] is True
    assert frozen["verdict"]["pure_norm_two_Whittaker_principal_part"] == (
        "obstructed"
    )
