import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_dual_extremal_28_3 import (
    FIRST_CUTTING_HEIGHT,
    FIRST_CUTTING_TRIPLE,
    EARLIER_ACCEPTANCE_SPECS,
    REPLACEMENT_WITNESS_SCALE,
    controls,
    earlier_acceptance_bound,
    first_cutting_old_witness_pairing,
    first_cutting_root_audit,
    necessary_triples,
    norm_eight_pairing_maximum,
    replacement_margin,
    replacement_witness_norm,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_dual_extremal_28_3_h362_replacement.json"


def test_complete_corrected_chronology_through_h362() -> None:
    rows = controls(FIRST_CUTTING_HEIGHT)
    assert len(rows) == 103
    assert len({row.height for row in rows}) == 70
    assert tuple(
        (row.height, row.norm_class.norm, len(necessary_triples(row)))
        for row in rows if necessary_triples(row)
    ) == (
        (Silver(93, 66), Silver(3), 1),
        (Silver(93, 66), Silver(6), 5),
        (Silver(124, 88), Silver(2), 1),
        (Silver(181, 128), Silver(2), 6),
        (Silver(186, 132), Silver(6), 4),
        (Silver(362, 256), Silver(2), 8),
    )


def test_first_threat_is_realized_primitive_and_cuts_old_witness() -> None:
    assert first_cutting_root_audit()
    maximum = norm_eight_pairing_maximum(FIRST_CUTTING_TRIPLE)
    assert maximum.maximum == Silver(53, 37)
    assert maximum.maximizing_pair == (53, 37)
    assert maximum.augmented_determinant == 0
    assert first_cutting_old_witness_pairing() == Silver(1, 1) / 2


def test_replacement_is_strict_hyperideal_against_complete_superset() -> None:
    assert REPLACEMENT_WITNESS_SCALE == Fraction(12, 25)
    assert replacement_witness_norm() == Silver(1777) / 625 - Silver(0, 2)
    assert replacement_witness_norm().sign(0) > 0
    margins = tuple(
        replacement_margin(row, triple)
        for row in controls(FIRST_CUTTING_HEIGHT)
        for triple in necessary_triples(row)
    )
    assert all(value.sign(0) > 0 for value in margins)
    assert min(margins) == Silver(14, 6) / 25


def test_cross_gram_acceptance_and_two_short_vector_exclusions() -> None:
    bounds = tuple(earlier_acceptance_bound(spec) for spec in EARLIER_ACCEPTANCE_SPECS)
    assert len(bounds) == 17
    assert all(bound.admitted_maximum.sign(0) <= 0 for bound in bounds)
    assert [bound.relaxed_maximum for bound in bounds[7:13]] == [
        Silver(0), Silver(-1), Silver(1), Silver(0), Silver(1), Silver(0)
    ]
    forbidden = tuple(row for bound in bounds for row in bound.forbidden_positive_extensions)
    assert len(forbidden) == 2
    assert [(row[0], row[1], row[3]) for row in forbidden] == [
        ((9, -7, 0, 15), Silver(1), 6),
        ((6, -6, 9, 13), Silver(1), 2),
    ]


def test_frozen_scope_and_replay() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["chronology"]["control_count"] == 103
    assert [row["short_square"] for row in payload["chronology"]["binary_minimum_gate"]["removed_spurious_types"]] == [6, 4, 6, 2, 6, 2]
    assert payload["first_cutting_wall"]["accepted"] is True
    assert payload["replacement"]["strict_against_every_necessary_type_through_bound"] is True
    completed = subprocess.run(
        [sys.executable, str(ROOT / "scripts/audit_silver_dual_extremal_28_3_h362_replacement.py"), "--verify"],
        cwd=ROOT, text=True, capture_output=True, check=True,
    )
    assert "PASS" in completed.stdout
