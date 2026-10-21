from __future__ import annotations

from hashlib import sha256

import pytest

from scripts.certify_h9_14_record30_a0_dual_exhaustion import (
    verify_frozen_replay_tools,
)


def test_deep_replay_accepts_only_declared_binary_digests(tmp_path) -> None:
    solver = tmp_path / "solver"
    checker = tmp_path / "checker"
    solver.write_bytes(b"frozen solver")
    checker.write_bytes(b"frozen checker")
    document = {
        "proof": {
            "solver_binary_sha256": sha256(solver.read_bytes()).hexdigest(),
            "checker_binary_sha256": sha256(checker.read_bytes()).hexdigest(),
        }
    }
    verify_frozen_replay_tools(document, solver=solver, checker=checker)

    checker.write_bytes(b"substituted checker")
    with pytest.raises(ValueError, match="DRAT checker binary digest mismatch"):
        verify_frozen_replay_tools(document, solver=solver, checker=checker)


def test_deep_replay_rejects_missing_binary(tmp_path) -> None:
    missing = tmp_path / "missing"
    document = {
        "proof": {
            "solver_binary_sha256": "0" * 64,
            "checker_binary_sha256": "0" * 64,
        }
    }
    with pytest.raises(ValueError, match="solver binary is not a file"):
        verify_frozen_replay_tools(document, solver=missing, checker=missing)
