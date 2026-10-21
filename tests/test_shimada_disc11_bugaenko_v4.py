import os
from pathlib import Path

import pytest

from scripts.certify_shimada_disc11_bugaenko_v4 import (
    EXPECTED_FIELD_MATRIX_SHA256,
    EXPECTED_FOLD_WORD_SHA256,
    _bugaenko_reversal_isometry,
    _diagram_automorphisms,
    certify,
)
from coxeter_polytopes.known_examples.bugaenko import bugaenko_coxeter_relations


DISC11 = Path(
    os.environ.get(
        "SHIMADA_DISC11", "/tmp/shimada_threediscs/threediscs/disc11.txt"
    )
)


def test_bugaenko_diagram_has_only_the_reversal() -> None:
    automorphisms = _diagram_automorphisms(bugaenko_coxeter_relations())
    assert len(automorphisms) == 2
    assert tuple(range(11)) in automorphisms
    assert (8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9) in automorphisms
    reversal_bug, reversal_std = _bugaenko_reversal_isometry()
    assert reversal_bug
    assert reversal_std


@pytest.mark.skipif(
    not DISC11.is_file(), reason="Shimada's 45 MB disc11.txt supplement is not vendored"
)
def test_exact_bugaenko_fixed_slice_v4_certificate() -> None:
    result = certify(DISC11)
    assert result["verified"] is True
    assert result["A_fixed_lattice"]["nikulin_invariants_r_a_delta"] == [18, 2, 1]
    assert result["A_fixed_lattice"]["fixed_locus"] == {
        "unique_positive_genus_component_genus": 1,
        "rational_curve_components": 8,
    }
    assert result["folds"]["exact_words"]["plus_root_count"] == 88
    assert result["folds"]["exact_words"]["minus_root_count"] == 88
    assert result["folds"]["exact_words"]["plus_exact_word_sha256"] == (
        EXPECTED_FOLD_WORD_SHA256["B_plus"]
    )
    assert result["folds"]["exact_words"]["minus_exact_word_sha256"] == (
        EXPECTED_FOLD_WORD_SHA256["B_minus"]
    )
    residuals = {item["name"]: item for item in result["residuals"]}
    assert {
        name: item["field_row_action_sha256"] for name, item in residuals.items()
    } == EXPECTED_FIELD_MATRIX_SHA256
    assert result["multiplication_table_checks"]["plus_minus_equals_A"] is True
    assert result["multiplication_table_checks"]["minus_plus_equals_A"] is True
