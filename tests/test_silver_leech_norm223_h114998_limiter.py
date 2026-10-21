import json
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_leech_norm223_limiter import (
    IMPOSSIBLE_TRIPLE,
    LIMITING_D_SCALED,
    LIMITING_E_SCALED,
    LIMITING_MAXIMUM,
    REALIZED_TRIPLE,
    norm223_limiting_root_audit,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_leech_norm223_h114998_limiter.json"


def test_complete_maximizing_locus_and_exact_obstruction() -> None:
    audit = norm223_limiting_root_audit()
    assert audit.surviving_triple_count == 216
    assert tuple(row.triple for row in audit.maximizing_types) == (
        IMPOSSIBLE_TRIPLE,
        REALIZED_TRIPLE,
    )
    assert tuple(
        (row.d_dot_u, row.d_dot_v, row.augmented_determinant)
        for row in audit.maximizing_types
    ) == ((679, 479, 14), (679, 479, 0))
    assert audit.maximizing_types[0].shifted_difference_square == 2
    assert not audit.maximizing_types[0].realized
    assert audit.maximizing_types[1].shifted_difference_square == 0
    assert audit.maximizing_types[1].realized


def test_explicit_root_is_primitive_reflective_and_realizes_maximum() -> None:
    audit = norm223_limiting_root_audit()
    assert audit.explicit_d_scaled == LIMITING_D_SCALED
    assert audit.explicit_e_scaled == LIMITING_E_SCALED
    assert audit.d_e_gram == ((4, -1), (-1, 6))
    assert audit.realized_u_v_gram == REALIZED_TRIPLE
    assert audit.direction_pairing == LIMITING_MAXIMUM
    assert audit.root_square == 2
    assert audit.primitive
    assert audit.reflection_integral
    assert audit.height_zero_empty


def test_complete_earlier_schur_and_same_height_acceptance() -> None:
    audit = norm223_limiting_root_audit()
    assert audit.earlier_component_count == 34
    assert all(
        row.perpendicular_residual.sign(0) > 0
        and row.threshold_gap.sign(0) > 0
        and row.monotonicity_discriminant.sign(0) > 0
        and row.rejection_discriminant.sign(0) > 0
        for row in audit.earlier_schur_certificates
    )
    tight = min(
        audit.earlier_schur_certificates,
        key=lambda row: row.rejection_discriminant,
    )
    assert tight.control_ordinal == audit.tight_earlier_control_ordinal == 565
    assert tight.threshold_gap == Silver(Fraction(63, 4), Fraction(45, 4))
    assert tight.perpendicular_residual == Silver(Fraction(55, 4), 9)
    assert tight.rejection_discriminant == Silver(57, 41)
    assert audit.same_height_q1_empty
    assert set(audit.positive_distinct_same_height_pairings) == {
        Silver(1),
        Silver(0, 1),
    }
    assert set(audit.same_height_difference_targets) == {
        Silver(2),
        Silver(4, -2),
    }
    assert audit.same_height_acute
    assert audit.accepted


def test_frozen_limiter_scope_and_replay() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["audit"]["accepted"] is True
    assert payload["conclusion"]["chronologically_vinberg_accepted"] is True
    assert payload["conclusion"]["finite_volume_through_prefix"] is False
    assert payload["conclusion"]["vinberg_termination"] is False
    completed = subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_silver_leech_norm223_h114998_limiter.py"),
            "--verify",
        ],
        cwd=ROOT,
        text=True,
        capture_output=True,
        check=True,
    )
    assert "PASS" in completed.stdout
