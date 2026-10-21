from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from coxeter_polytopes.exact.silver import Silver
from coxeter_polytopes.search.silver_modular_rank24_sequence import (
    N4_HEIGHT_BOUND,
    N4_NORM_CLASSES,
    family_matrix_identity,
    level,
    minimum,
    n4_first_packet_audit,
    n4_radical_basis,
    n4_scaled_radical_gram,
)
from coxeter_polytopes.search.silver_modular_rank24_sequence import _rref_kernel_mod_prime, sequence_gram


ROOT = Path(__file__).resolve().parents[1]
ARTIFACT = ROOT / "data/enumerations/silver_modular_rank24_sequence_n4.json"


def test_exact_family_and_n4_local_structure() -> None:
    assert family_matrix_identity()
    assert level(4) == 41 and minimum(4) == 24
    kernel, pivots = _rref_kernel_mod_prime(sequence_gram(4), 41)
    assert len(kernel) == len(pivots) == 12
    assert len(n4_radical_basis()) == len(n4_scaled_radical_gram()) == 24
    assert tuple(row.norm for row in N4_NORM_CLASSES) == (
        Silver(1), Silver(7, -2), Silver(7, 2), Silver(41),
        Silver(2), Silver(14, -4), Silver(14, 4), Silver(82),
    )


def test_complete_n4_prefix_has_strict_hyperideal_witness() -> None:
    audit = n4_first_packet_audit()
    assert N4_HEIGHT_BOUND == 5000
    assert audit.control_count == 196
    assert audit.distinct_height_count == 132
    assert audit.first_control.height == Silver(806, 570)
    assert [(row.u_square, row.v_square, row.u_dot_v) for row in audit.triples] == [
        (798, 408, 570), (794, 410, 570),
    ]
    assert [row.maximum for row in audit.maxima] == [Silver(17, 12), Silver(14, 10)]
    assert audit.witness_norm == Silver(339, 0) / 41 - Silver(0, 5)
    assert audit.margins == (Silver(7, 4), Silver(16, 10))
    assert audit.strict_hyperideal


def test_frozen_certificate_and_scope() -> None:
    subprocess.run(
        [sys.executable, str(ROOT / "scripts/audit_silver_modular_rank24_sequence_n4.py"), "--verify"],
        cwd=ROOT,
        check=True,
    )
    payload = json.loads(ARTIFACT.read_text())
    assert payload["verdict"] == {
        "bounded_prefix_has_finite_volume": False,
        "bounded_prefix_has_open_hyperideal_patch": True,
        "compact_example_found": False,
        "eventual_chamber": "OPEN",
    }
    assert "bounded" in payload["scope"].lower()

