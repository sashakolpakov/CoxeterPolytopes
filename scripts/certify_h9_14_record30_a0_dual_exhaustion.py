#!/usr/bin/env python3
"""Verify the frozen dual-CNF/DRAT theorem for Record-30 ``a=0``."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Mapping


ROOT = Path(__file__).resolve().parents[1]
ARCHIVE = ROOT / "data/enumerations/h9_14_record30_a0_dual_exhaustion.json"
CNF_ZST = ROOT / "data/enumerations/h9_14_record30_a0_dual_exhaustion.cnf.zst"
PROOF_ZST = ROOT / "data/enumerations/h9_14_record30_a0_dual_exhaustion.drat.zst"
TRANSCRIPT = ROOT / "data/enumerations/h9_14_record30_a0_dual_exhaustion_drat_trim.txt"
PRODUCTION = ROOT / "data/enumerations/h9_14_record30_a0_connected_ridge_census.json"
LABELS = ROOT / "data/enumerations/h9_14_record30_a0_all_target_labels.json"
GENERATOR = ROOT / "scripts/independent_h9_14_record30_a0_dual_cnf.py"
HOSTILE = ROOT / "scripts/verify_h9_14_record30_a0_dual_cnf.py"
TEST_SOURCE = ROOT / "tests/test_h9_14_record30_a0_dual_cnf.py"

EXPECTED_CERTIFICATE_SHA256 = "436b9bc834529f7d2d94754e43c3f785bb937ec0f7157570ad9a06d56858ab2d"
EXPECTED_GENERATOR_SHA256 = "af5d5b1d3429647a6aacef11c0ff4781f821b90384305bed73171ea686a6c1f0"
EXPECTED_HOSTILE_SHA256 = "c8461a1cd1f715f9845ed08af38e0cd192470ef691060863760ab81fea6d53d9"
EXPECTED_TEST_SHA256 = "7acd3937b6803cd45f5c5a35eafa64085537bb719021a1700c943a53b34dafd9"


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as stream:
        while block := stream.read(1 << 20):
            value.update(block)
    return value.hexdigest()


def internal_digest(document: Mapping[str, object]) -> str:
    payload = dict(document)
    payload.pop("certificate_sha256", None)
    return sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify(document: Mapping[str, object]) -> dict[str, object]:
    if document.get("certificate_sha256") != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("Record-30 a0 theorem certificate digest changed")
    if internal_digest(document) != EXPECTED_CERTIFICATE_SHA256:
        raise ValueError("Record-30 a0 theorem internal digest is invalid")
    expected_sources = {
        GENERATOR: EXPECTED_GENERATOR_SHA256,
        HOSTILE: EXPECTED_HOSTILE_SHA256,
        TEST_SOURCE: EXPECTED_TEST_SHA256,
    }
    for path, expected in expected_sources.items():
        if digest(path) != expected:
            raise ValueError(f"frozen source changed: {path.name}")

    dual = document.get("dual_cnf")
    proof = document.get("proof")
    label = document.get("label_gate")
    scope = document.get("scope")
    hostile = document.get("hostile_encoder_audit")
    if not all(isinstance(item, Mapping) for item in (dual, proof, label, scope, hostile)):
        raise ValueError("Record-30 a0 theorem has malformed sections")
    assert isinstance(dual, Mapping) and isinstance(proof, Mapping)
    assert isinstance(label, Mapping) and isinstance(scope, Mapping)
    assert isinstance(hostile, Mapping)

    archived_files = (
        (CNF_ZST, dual["cnf_zstd_sha256"], dual["cnf_zstd_size_bytes"]),
        (PROOF_ZST, proof["proof_zstd_sha256"], proof["proof_zstd_size_bytes"]),
        (TRANSCRIPT, proof["checker_transcript_sha256"], TRANSCRIPT.stat().st_size),
        (PRODUCTION, document["production_crosscheck"]["connected_ridge_census_file_sha256"], PRODUCTION.stat().st_size),
        (LABELS, label["all_target_label_replay_file_sha256"], LABELS.stat().st_size),
    )
    for path, expected_hash, expected_size in archived_files:
        if digest(path) != expected_hash or path.stat().st_size != expected_size:
            raise ValueError(f"frozen artifact changed: {path.name}")
    if "s VERIFIED" not in TRANSCRIPT.read_text(encoding="utf-8"):
        raise ValueError("normalized DRAT transcript lost VERIFIED")

    production = json.loads(PRODUCTION.read_text(encoding="utf-8"))
    labels = json.loads(LABELS.read_text(encoding="utf-8"))
    if production.get("target_orbit_count") != 8 or len(
        production.get("target_orbit_representatives", [])
    ) != 8:
        raise ValueError("production a0 target census changed")
    if (
        labels.get("target_orbit_count") != 8
        or labels.get("sparse_label_survivor_count") != 0
        or labels.get("summary_sha256")
        != label["all_target_label_replay_summary_sha256"]
    ):
        raise ValueError("all-target a0 label no-go changed")
    if (
        dual.get("variable_count") != 226358
        or dual.get("clause_count") != 984935
        or dual.get("blocked_seed_form_count") != 63960
        or dual.get("seed_stabilizer_lex_encoded") is not True
        or proof.get("checker_status") != "VERIFIED"
        or proof.get("solver_status") != "UNSATISFIABLE"
        or hostile.get("positive_target_witness_count") != 8
        or hostile.get("all_blocking_clauses_exactly_reconstructed") is not True
        or hostile.get("all_positive_witnesses_literally_satisfy_unblocked_cnf") is not True
        or scope.get("whole_record30_a0_schedule_excluded") is not True
        or scope.get("finite_volume_noncompact_statement") is not False
    ):
        raise ValueError("Record-30 a0 theorem gate or scope changed")
    return {
        "target_orbit_count": 8,
        "blocked_seed_form_count": 63960,
        "cnf_variable_count": 226358,
        "cnf_clause_count": 984935,
        "drat_verified": True,
        "positive_witness_count": 8,
        "compact_label_survivor_count": 0,
        "finite_volume_noncompact_statement": False,
    }


def decompress(zstd: Path, source: Path, target: Path) -> None:
    with target.open("wb") as output:
        process = subprocess.run(
            [str(zstd), "-d", "-c", str(source)], stdout=output, check=False
        )
    if process.returncode != 0:
        raise ValueError(f"failed to decompress {source.name}")


def verify_frozen_replay_tools(
    document: Mapping[str, object], *, solver: Path, checker: Path
) -> None:
    """Require the exact solver/checker binaries named by the certificate.

    The archived DRAT proof is independently checkable, so accepting a newer
    trusted checker would not invalidate the theorem.  ``--deep`` is defined
    more narrowly, however: it reproduces the *frozen* run and therefore must
    not silently substitute different executables.
    """
    proof = document.get("proof")
    if not isinstance(proof, Mapping):
        raise ValueError("Record-30 a0 theorem has malformed proof provenance")
    expected = (
        (solver, proof.get("solver_binary_sha256"), "solver"),
        (checker, proof.get("checker_binary_sha256"), "DRAT checker"),
    )
    for path, expected_hash, label in expected:
        if not isinstance(expected_hash, str):
            raise ValueError(f"frozen {label} digest is missing")
        if not path.is_file():
            raise ValueError(f"frozen {label} binary is not a file: {path}")
        actual_hash = digest(path)
        if actual_hash != expected_hash:
            raise ValueError(
                f"{label} binary digest mismatch: expected {expected_hash}, "
                f"got {actual_hash}"
            )


def deep_replay(
    document: Mapping[str, object], *, zstd: Path, solver: Path, checker: Path
) -> None:
    verify_frozen_replay_tools(document, solver=solver, checker=checker)
    with tempfile.TemporaryDirectory(prefix="record30_a0_dual_deep_") as raw:
        temporary = Path(raw)
        cnf = temporary / "formula.cnf"
        proof = temporary / "proof.drat"
        decompress(zstd, CNF_ZST, cnf)
        decompress(zstd, PROOF_ZST, proof)
        if digest(cnf) != document["dual_cnf"]["cnf_raw_sha256"]:
            raise ValueError("decompressed CNF digest changed")
        if digest(proof) != document["proof"]["proof_raw_sha256"]:
            raise ValueError("decompressed proof digest changed")
        subprocess.run(
            [
                "python3",
                str(HOSTILE),
                "--cnf",
                str(cnf),
                "--proof",
                str(proof),
                "--solver",
                str(solver),
                "--drat-trim",
                str(checker),
            ],
            cwd=ROOT,
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deep", action="store_true")
    parser.add_argument("--zstd", type=Path, default=Path("zstd"))
    parser.add_argument("--solver", type=Path)
    parser.add_argument("--drat-trim", type=Path)
    args = parser.parse_args()
    document = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    summary = verify(document)
    if args.deep:
        if args.solver is None or args.drat_trim is None:
            raise ValueError("--deep requires --solver and --drat-trim")
        deep_replay(
            document,
            zstd=args.zstd,
            solver=args.solver,
            checker=args.drat_trim,
        )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
