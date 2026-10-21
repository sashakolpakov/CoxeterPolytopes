from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from coxeter_polytopes.exact.silver import SQRT2, Silver
from coxeter_polytopes.search.silver_dual_extremal_28_3 import (
    DETERMINANT,
    FIRST_NECESSARY_HEIGHT,
    RADICAL_THREE_INDEX,
    RADICAL_THREE_KISSING_NUMBER,
    RADICAL_THREE_MINIMUM,
    catalogue_audit,
    first_prefix_audit,
)


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_dual_extremal_28_3_first_prefix.json"


def test_catalogue_and_radical_invariants() -> None:
    assert catalogue_audit()
    payload = json.loads(ARTIFACT.read_text())
    spatial = payload["spatial_lattice"]
    assert spatial["determinant"] == str(DETERMINANT)
    pari = spatial["pari"]
    assert pari["pairing_radical_preimage_index"] == RADICAL_THREE_INDEX
    assert pari["pairing_radical_preimage_minimum"] == RADICAL_THREE_MINIMUM
    assert (
        pari["pairing_radical_preimage_minimum_shell_signed_count"]
        == RADICAL_THREE_KISSING_NUMBER
    )


def test_complete_first_necessary_prefix_has_strict_hyperideal_witness() -> None:
    audit = first_prefix_audit()
    assert audit.control_count == 51
    assert audit.height_count == 35
    assert audit.earlier_nonempty_controls == 0
    assert len(audit.endpoint_components) == 2
    assert sum(len(row[2]) for row in audit.endpoint_components) == 6
    assert audit.witness_scale == Silver(1) / 2
    assert audit.witness_excess == 3 - 2 * SQRT2
    assert audit.witness_excess.sign(0) > 0
    assert all(margin.sign(0) > 0 for margin in audit.margins)
    assert tuple(maximum.maximum for maximum in audit.maxima) == (
        24 + 18 * SQRT2,
        45 + 30 * SQRT2,
        44 + 31 * SQRT2,
        44 + 32 * SQRT2,
        40 + 31 * SQRT2,
        33 + 27 * SQRT2,
    )
    assert FIRST_NECESSARY_HEIGHT == 93 + 66 * SQRT2


def test_frozen_certificate_regenerates_exactly() -> None:
    subprocess.run(
        [
            sys.executable,
            str(ROOT / "scripts/audit_silver_dual_extremal_28_3_first_prefix.py"),
            "--verify",
        ],
        cwd=ROOT,
        check=True,
    )


def test_scope_is_bounded_and_does_not_claim_a_compact_example() -> None:
    payload = json.loads(ARTIFACT.read_text())
    assert payload["verdict"] == {
        "bounded_prefix_has_finite_volume": False,
        "bounded_prefix_has_open_hyperideal_patch": True,
        "compact_example_found": False,
        "eventual_chamber": "OPEN",
    }
    assert "bounded prefix theorem" in payload["scope"]
