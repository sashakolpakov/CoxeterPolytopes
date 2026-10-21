import hashlib
import json
from fractions import Fraction
from pathlib import Path

import pytest


pytest.importorskip("flint")

from scripts.certify_og10_whittaker_first_compensator import (  # noqa: E402
    MINIMAL_COSET_A,
    MINIMAL_COSET_B,
    M_MINIMAL,
    W_MINIMAL,
    _minimal_shell_data,
    certify,
)


def test_complete_minimal_off_slice_shell():
    result = certify()

    assert result["verified"]
    assert all(result["checks"].values())
    assert result["discriminant_module"]["order"] == 176
    assert result["minimal_shell"]["index_norm"] == "9/880"
    assert result["minimal_shell"]["unit_orbit_count_modulo_mu_sign"] == 2
    assert M_MINIMAL.norm() == Fraction(9, 880)
    assert W_MINIMAL.norm() == -99


def test_finite_discriminant_enumeration_and_two_orbits():
    shell = _minimal_shell_data()

    assert sum(len(cosets) for cosets in shell["allowed"].values()) == 132
    assert len(shell["allowed_keys"]) == 54
    assert shell["norm_residues"] == {0, 11, 33}
    assert shell["smaller_absolute_norms"] == (11, 33, 44, 55, 77, 88)
    assert all(
        not intersection
        for intersections in shell["possible_intersections"].values()
        for intersection in intersections
    )
    assert len(shell["minimal_intersections"][0]) == 15
    assert not shell["minimal_intersections"][1]
    assert MINIMAL_COSET_A in shell["orbit_a"]
    assert MINIMAL_COSET_B in shell["orbit_b"]
    assert not (shell["orbit_a"] & shell["orbit_b"])
    assert shell["minimal_allowed_cosets"] == (
        shell["orbit_a"] | shell["orbit_b"]
    )


def test_exact_scaled_coefficient_matrix_diagonal_dominance():
    result = certify()
    bounds = result["scaled_fourier_matrix_bounds"]

    target_error = Fraction(bounds["target_diagonal_error_upper"])
    minimal_error = Fraction(bounds["minimal_diagonal_error_upper"])
    cross = Fraction(bounds["target_minimal_cross_upper"])
    target_margin = Fraction(bounds["target_row_dominance_margin_lower"])
    minimal_margin = Fraction(bounds["minimal_row_dominance_margin_lower"])

    assert bounds["matrix_size"] == 3
    assert target_margin == 1 - target_error - 2 * cross > Fraction(2, 5)
    assert minimal_margin == 1 - 2 * minimal_error - cross > Fraction(99, 100)
    assert "complete globally minimal norm-9/880 off-slice shell" in result[
        "bounded_conclusion"
    ]


def test_archived_certificate_digests():
    root = Path(__file__).resolve().parents[1]
    artifact = json.loads(
        (root / "data/candidates/og10_whittaker_first_compensator.json").read_text()
    )
    digests = artifact["digests"]

    for key, relative_path in {
        "executable_sha256": "scripts/certify_og10_whittaker_first_compensator.py",
        "report_sha256": "results/og10_whittaker_first_compensator.md",
        "test_sha256": "tests/test_og10_whittaker_first_compensator.py",
    }.items():
        assert hashlib.sha256((root / relative_path).read_bytes()).hexdigest() == digests[
            key
        ]

    stdout = json.dumps(certify(), indent=2, sort_keys=True) + "\n"
    assert hashlib.sha256(stdout.encode()).hexdigest() == digests[
        "certificate_stdout_sha256"
    ]
