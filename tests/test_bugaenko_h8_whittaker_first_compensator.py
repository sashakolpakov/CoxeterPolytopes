import hashlib
import json
from fractions import Fraction
from pathlib import Path

import pytest


pytest.importorskip("flint")

from scripts.certify_bugaenko_h8_whittaker_first_compensator import (  # noqa: E402
    MINIMAL_COSET,
    M_MINIMAL,
    W_MINIMAL,
    _minimal_shell_data,
    certify,
)


def test_complete_globally_minimal_shell() -> None:
    result = certify()
    shell = _minimal_shell_data()

    assert result["verified"] is True
    assert all(result["checks"].values())
    assert sum(len(values) for values in shell["residues"].values()) == 64
    assert len(shell["residue_keys"]) == 40
    assert shell["norm_residues"] == {0, 1, 3, 4}
    assert W_MINIMAL.norm() == -4
    assert M_MINIMAL.norm() == Fraction(1, 80)
    assert shell["minimal_cosets"] == shell["coset_orbit"]
    assert MINIMAL_COSET in shell["coset_orbit"]
    assert len(shell["coset_orbit"]) == 3


def test_minimal_shell_is_nonempty_and_slice_disjoint() -> None:
    result = certify()
    geometry = result["real_slice_classification"]
    shell = result["complete_minimal_shell"]

    assert shell["represented_Hhat_value"] == "1/2"
    assert shell["component_orbit_count_modulo_units_and_sign"] == 1
    assert shell["components_in_full_orbit"] == 3
    assert shell["every_component_is_self_opposite"] is True
    assert geometry == {
        "facet_sheets": 0,
        "lower_face_contact_sheets": 0,
        "interior_cut_sheets": 0,
        "disposition": "entire nonempty Heegner divisor is disjoint",
        "proof": (
            "every lambda in mu+L has auxiliary coefficient y/4+a; "
            "it cannot vanish when y is nonzero modulo 4"
        ),
    }


def test_exact_two_column_diagonal_dominance() -> None:
    result = certify()
    bounds = result["scaled_fourier_matrix_bounds"]

    assert bounds["matrix_size"] == 2
    assert Fraction(bounds["target_diagonal_error_upper"]) == Fraction(3, 4)
    assert Fraction(bounds["minimal_diagonal_error_upper"]) == Fraction(1, 500)
    assert Fraction(bounds["each_cross_error_upper"]) == Fraction(1, 6)
    assert Fraction(bounds["target_row_dominance_margin_lower"]) == Fraction(
        1, 12
    )
    assert Fraction(bounds["minimal_row_dominance_margin_lower"]) == Fraction(
        1247, 1500
    )
    assert "no weakly holomorphic principal part" in result[
        "bounded_conclusion"
    ]


def test_archived_certificate_digests() -> None:
    root = Path(__file__).resolve().parents[1]
    artifact = json.loads(
        (
            root
            / "data/candidates/bugaenko_h8_whittaker_first_compensator.json"
        ).read_text(encoding="utf-8")
    )
    digests = artifact["digests"]

    for key, relative_path in {
        "executable_sha256": (
            "scripts/certify_bugaenko_h8_whittaker_first_compensator.py"
        ),
        "report_sha256": (
            "results/bugaenko_h8_whittaker_first_compensator.md"
        ),
        "test_sha256": (
            "tests/test_bugaenko_h8_whittaker_first_compensator.py"
        ),
    }.items():
        assert hashlib.sha256((root / relative_path).read_bytes()).hexdigest() == (
            digests[key]
        )

    stdout = json.dumps(certify(), indent=2, sort_keys=True) + "\n"
    assert hashlib.sha256(stdout.encode("utf-8")).hexdigest() == digests[
        "certificate_stdout_sha256"
    ]
