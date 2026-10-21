#!/usr/bin/env python3
"""Run the deliberately bounded RC7 p7 first-wall regression.

Unlike ``run_cubic7_alvin_live.py``, this checker intentionally passes
``-maxv 6``.  Its only claim is that AlVin reaches the exact sixth wall which
was lost by the historical nonintegral-R1 coefficient bound.  It makes no
termination or finite-volume claim.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from typing import Any


REPOSITORY = Path(__file__).resolve().parents[1]
EXPECTED_COMMIT = "d92006a842cb4668e7c89e8c8e2d30337e4b0f26"
QF = "[-1,0,0],1,1,1,[-3,-2,-2]"
EXPECTED_ROOTS = (
    "(0, -1, 0, 0, 0)",
    "(0, 0, -1, 1, 0)",
    "(0, 0, 0, -1, 1)",
    "(0, 0, 0, 0, -1)",
    "(1, 0, [0,0,-1], 0, 0)",
    "([0,1,-1], [0,-1,-1], 0, 0, 0)",
)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _run_text(command: list[str]) -> str:
    return subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    ).stdout.strip()


def _image_metadata(image: str) -> dict[str, Any]:
    inspected = json.loads(_run_text(["docker", "image", "inspect", image]))[0]
    commit = _run_text(
        [
            "docker",
            "run",
            "--rm",
            "--entrypoint",
            "git",
            image,
            "-C",
            "/opt/AlVin",
            "rev-parse",
            "HEAD",
        ]
    )
    return {
        "id": inspected["Id"],
        "created": inspected["Created"],
        "architecture": inspected["Architecture"],
        "os": inspected["Os"],
        "alvin_commit": commit,
    }


def parse_roots(stdout: str) -> tuple[str, ...]:
    rows: dict[int, str] = {}
    for match in re.finditer(r"^\s*e(\d+)\s*=\s*(.+?)\s*$", stdout, re.MULTILINE):
        rows[int(match.group(1))] = match.group(2)
    return tuple(rows[index] for index in sorted(rows))


def assess_prefix(stdout: str, stderr: str, return_code: int) -> dict[str, bool]:
    roots = parse_roots(stdout)
    combined = stdout + "\n" + stderr
    return {
        "return_code_zero": return_code == 0,
        "exact_first_six_roots_in_order": roots == EXPECTED_ROOTS,
        "repaired_wall_is_e6": len(roots) == 6 and roots[-1] == EXPECTED_ROOTS[-1],
        "repaired_norm_class_listed": "[-6,-4,-4]" in stdout,
        "no_unknown_weight": "Unknown weight" not in combined,
        "bounded_nontermination_marker": (
            "The algorithm did not terminate; the polyhedron may be of infinite volume"
            in stdout
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--image", default="coxeter-alvin:local")
    arguments = parser.parse_args(argv)

    case_directory = arguments.output_root.resolve() / "p7-spatial-r1-prefix"
    if case_directory.exists():
        parser.error(f"refusing to overwrite existing case directory: {case_directory}")
    case_directory.mkdir(parents=True)

    image = _image_metadata(arguments.image)
    if image["alvin_commit"] != EXPECTED_COMMIT:
        parser.error(
            f"image contains AlVin {image['alvin_commit']}, expected {EXPECTED_COMMIT}"
        )
    command = [
        "docker",
        "run",
        "--rm",
        "--network=none",
        arguments.image,
        "-k",
        "rc7",
        "-qf",
        QF,
        "-noimage",
        "-maxv",
        "6",
        "-debug",
    ]
    request = {
        "schema": "coxeter.k7_alvin.p7_prefix_request.v1",
        "claim": "bounded exact first-six-wall regression only",
        "qf": QF,
        "max_vectors": 6,
        "contains_maxv": True,
        "command_argv": command,
        "image": image,
        "expected_alvin_commit": EXPECTED_COMMIT,
        "runner_sha256": _sha256(Path(__file__)),
        "started_at_utc": _utc_now(),
    }
    (case_directory / "request.json").write_text(
        json.dumps(request, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    completed = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    stdout_path = case_directory / "stdout.txt"
    stderr_path = case_directory / "stderr.txt"
    stdout_path.write_text(completed.stdout, encoding="utf-8")
    stderr_path.write_text(completed.stderr, encoding="utf-8")
    roots = parse_roots(completed.stdout)
    checks = assess_prefix(completed.stdout, completed.stderr, completed.returncode)
    result = {
        "schema": "coxeter.k7_alvin.p7_prefix_result.v1",
        "claim": "bounded exact first-six-wall regression only",
        "complete_termination_claimed": False,
        "finite_volume_claimed": False,
        "return_code": completed.returncode,
        "roots": roots,
        "checks": checks,
        "acceptance_passes": all(checks.values()),
        "stdout_sha256": _sha256(stdout_path),
        "stderr_sha256": _sha256(stderr_path),
        "finished_at_utc": _utc_now(),
    }
    (case_directory / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["acceptance_passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
