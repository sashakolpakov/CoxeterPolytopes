#!/usr/bin/env python3
"""Run one pinned AlVin RC7 case and archive an auditable transcript.

This wrapper deliberately does not offer AlVin's ``-maxv`` option.  A seed is
reported as terminating only when the container exits, AlVin prints
``Algorithm ended``, and a CoxIter graph is present.  Killing this wrapper or
using its external timeout produces an explicitly inconclusive record.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading
import time
from typing import Any, BinaryIO


REPOSITORY = Path(__file__).resolve().parents[1]
ALVIN_CONTEXT = REPOSITORY / "tools" / "alvin"
EXPECTED_COMMIT = "d92006a842cb4668e7c89e8c8e2d30337e4b0f26"
CONTEXT_FILES = (
    "Dockerfile",
    "vendor.sha256",
    "alvin-norm-powers.patch",
    "alvin-aarch64.patch",
    "alvin-rc7-exact-bounds.patch",
    "alvin-rc7-angle14.patch",
    "alvin-rc7-fraction-conjugation.patch",
    "alvin-rc7-spatial-r1.patch",
    "gaol-aarch64.patch",
)


HOSTILE_MATRIX = (
    (1, 2, 1, 2),
    (2, 1, 2, 1),
    (1, 2, 1, 14),
    (2, 1, 14, 1),
)


def _high_dimensional_qf(hyperbolic_dimension: int, last_coefficient: str) -> str:
    """Build ``diag(-theta, 1^(d-1), t)`` as an exact AlVin argument."""

    if hyperbolic_dimension < 2:
        raise ValueError("a hyperbolic diagonal form must have dimension at least two")
    return ",".join(
        (
            "[-1,0,0]",
            *("1" for _ in range(hyperbolic_dimension - 1)),
            last_coefficient,
        )
    )


CASES: dict[str, dict[str, Any]] = {
    "control": {
        "qf": "[-1,0,0],1,1,1,1",
        "kind": "regression",
        "expected_facets": 8,
    },
    "hostile-m14": {
        "qf": "[-1,0,0],1,[-2,-5,-7]",
        "kind": "regression",
        "expected_facets": 4,
        "expected_matrix": HOSTILE_MATRIX,
    },
    "ramified-7-rank5": {
        "qf": "[-1,0,0],1,1,1,[-3,-2,-2]",
        "kind": "seed",
        "priority": 1,
        "support": "p7",
        "radical_norm": 7,
    },
    "inert-2-rank5": {
        "qf": "[-1,0,0],1,1,1,[-2,-2,-2]",
        "kind": "seed",
        "priority": 1,
        "support": "(2)",
        "radical_norm": 8,
    },
    "ramified-7-rank6": {
        "qf": "[-1,0,0],1,1,1,1,[-3,-2,-2]",
        "kind": "seed",
        "priority": "rank6-stretch",
        "support": "p7",
        "radical_norm": 7,
        "lane": "H5-stretch",
    },
    "inert-2-rank6": {
        "qf": "[-1,0,0],1,1,1,1,[-2,-2,-2]",
        "kind": "seed",
        "priority": "rank6-stretch",
        "support": "(2)",
        "radical_norm": 8,
        "lane": "H5-stretch",
    },
    "ramified-7-h22": {
        "qf": _high_dimensional_qf(22, "[-3,-2,-2]"),
        "kind": "seed",
        "priority": "headline-H22",
        "support": "p7",
        "radical_norm": 7,
        "lane": "high-dimensional-BnA1",
        "hyperbolic_dimension": 22,
        "ambient_rank": 23,
    },
    "inert-2-h22": {
        "qf": _high_dimensional_qf(22, "[-2,-2,-2]"),
        "kind": "seed",
        "priority": "headline-H22",
        "support": "(2)",
        "radical_norm": 8,
        "lane": "high-dimensional-BnA1",
        "hyperbolic_dimension": 22,
        "ambient_rank": 23,
    },
    "ramified-7-h29": {
        "qf": _high_dimensional_qf(29, "[-3,-2,-2]"),
        "kind": "seed",
        "priority": "headline-H29",
        "support": "p7",
        "radical_norm": 7,
        "lane": "high-dimensional-BnA1",
        "hyperbolic_dimension": 29,
        "ambient_rank": 30,
    },
    "inert-2-h29": {
        "qf": _high_dimensional_qf(29, "[-2,-2,-2]"),
        "kind": "seed",
        "priority": "headline-H29",
        "support": "(2)",
        "radical_norm": 8,
        "lane": "high-dimensional-BnA1",
        "hyperbolic_dimension": 29,
        "ambient_rank": 30,
    },
    "timelike-minus-2theta-rank5": {
        "qf": "[-2,0,0],1,1,1,1",
        "kind": "seed",
        "priority": "1b",
        "lane": "timelike-scaled",
        "description": "diag(-2theta,1^4)",
    },
    "timelike-minus-theta-times-2-minus-theta-rank5": {
        "qf": "[-4,-1,-2],1,1,1,1",
        "kind": "seed",
        "priority": "1b",
        "lane": "timelike-scaled",
        "description": "diag(-theta(2-theta),1^4)",
    },
    "split-13-r7-rank5": {
        "qf": "[-1,0,0],1,1,1,[-5,-2,-4]",
        "kind": "seed",
        "priority": 2,
        "support": "p13(theta=7)",
        "radical_norm": 13,
    },
    "split-13-r10-rank5": {
        "qf": "[-1,0,0],1,1,1,[-4,-2,-3]",
        "kind": "seed",
        "priority": 2,
        "support": "p13(theta=10)",
        "radical_norm": 13,
    },
    "split-13-r8-rank5": {
        "qf": "[-1,0,0],1,1,1,[-3,-2,-3]",
        "kind": "seed",
        "priority": 2,
        "support": "p13(theta=8)",
        "radical_norm": 13,
    },
    "p7-p13-r7-rank5": {
        "qf": "[-1,0,0],1,1,1,[-5,-4,-5]",
        "kind": "seed",
        "priority": 3,
        "support": "p7*p13(theta=7)",
        "radical_norm": 91,
        "lane": "odd-squarefree-product",
    },
    "p7-p13-r10-rank5": {
        "qf": "[-1,0,0],1,1,1,[-9,-4,-8]",
        "kind": "seed",
        "priority": 3,
        "support": "p7*p13(theta=10)",
        "radical_norm": 91,
        "lane": "odd-squarefree-product",
    },
    "p7-p13-r8-rank5": {
        "qf": "[-1,0,0],1,1,1,[-10,-4,-7]",
        "kind": "seed",
        "priority": 3,
        "support": "p7*p13(theta=8)",
        "radical_norm": 91,
        "lane": "odd-squarefree-product",
    },
}


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
        "repo_digests": inspected.get("RepoDigests", []),
        "alvin_commit": commit,
    }


def _split_qf(text: str) -> tuple[str, ...]:
    parts: list[str] = []
    start = 0
    depth = 0
    for index, character in enumerate(text):
        if character == "[":
            depth += 1
        elif character == "]":
            depth -= 1
        elif character == "," and depth == 0:
            parts.append(re.sub(r"\s+", "", text[start:index]))
            start = index + 1
    parts.append(re.sub(r"\s+", "", text[start:]))
    return tuple(parts)


def _canonical_qf_element(text: str) -> tuple[int, int, int]:
    """Convert an AlVin RC7 diagonal entry to the integral power basis.

    AlVin may print a rational integer ``n`` either as ``n`` or in its
    trigonometric basis as ``[-n,-n,-n]``.  Comparing the source spelling
    would therefore report a harmless canonicalization as a form change.
    For ``theta = 2*cos(2*pi/7)``, AlVin's ``[u,v,w]`` equals
    ``(w-2v) + (u-w) theta + (v-w) theta^2``.
    """

    value = re.sub(r"\s+", "", text)
    if re.fullmatch(r"[+-]?\d+", value):
        return (int(value), 0, 0)
    match = re.fullmatch(r"\[([+-]?\d+),([+-]?\d+),([+-]?\d+)\]", value)
    if match is None:
        raise ValueError(f"unsupported integral RC7 coefficient: {text!r}")
    u, v, w = (int(coefficient) for coefficient in match.groups())
    return (w - 2 * v, u - w, v - w)


def _canonical_qf(text: str) -> tuple[tuple[int, int, int], ...]:
    return tuple(_canonical_qf_element(entry) for entry in _split_qf(text))


def _printed_qf(stdout: str) -> str | None:
    match = re.search(r"^Quadratic form \([^\n]+\):\s*(.+)$", stdout, re.MULTILINE)
    return match.group(1).strip() if match else None


def _parse_coxiter(path: Path) -> tuple[int, int, tuple[tuple[int, ...], ...]]:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    facets, dimension = (int(value) for value in lines[0].split())
    labels = [[2 for _ in range(facets)] for _ in range(facets)]
    for index in range(facets):
        labels[index][index] = 1
    cursor = 1
    names = tuple(str(index + 1) for index in range(facets))
    if lines[cursor].lower().startswith("vertices labels:"):
        names = tuple(lines[cursor].split(":", 1)[1].split())
        cursor += 1
    indices = {name: index for index, name in enumerate(names)}
    for line in lines[cursor:]:
        left, right, label = line.split()[:3]
        i, j = indices[left], indices[right]
        labels[i][j] = labels[j][i] = int(label)
    return facets, dimension, tuple(tuple(row) for row in labels)


def _pump(source: BinaryIO, destination: BinaryIO, terminal: BinaryIO | None) -> None:
    # ``Popen(..., stdout=PIPE)`` exposes a ``BufferedReader``.  Its ordinary
    # ``read(n)`` is allowed to wait for all ``n`` bytes, which hid hours of
    # live AlVin output until either 64 KiB accumulated or the child exited.
    # ``read1`` returns the bytes currently available from the underlying
    # pipe, while retaining the same byte-for-byte archive semantics.  Keep
    # the fallback for small BinaryIO test doubles without ``read1``.
    read_available = getattr(source, "read1", source.read)
    for block in iter(lambda: read_available(64 * 1024), b""):
        destination.write(block)
        destination.flush()
        if terminal is not None:
            terminal.write(block)
            terminal.flush()


def _stop_owned_container(name: str) -> None:
    subprocess.run(
        ["docker", "stop", "--time", "10", name],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", required=True, choices=tuple(CASES))
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--image", default="coxeter-alvin:local")
    parser.add_argument(
        "--cpus",
        type=float,
        help="optional Docker CPU limit for this case (for parallel sweeps)",
    )
    parser.add_argument(
        "--external-timeout-seconds",
        type=float,
        help="optional wall-clock stop; any such stop is recorded as inconclusive",
    )
    parser.add_argument("--echo-program-output", action="store_true")
    arguments = parser.parse_args(argv)
    if arguments.cpus is not None and arguments.cpus <= 0:
        parser.error("--cpus must be positive")

    case = dict(CASES[arguments.case])
    case_directory = arguments.output_root.resolve() / arguments.case
    if case_directory.exists():
        parser.error(f"refusing to overwrite existing case directory: {case_directory}")
    output_directory = case_directory / "output"
    output_directory.mkdir(parents=True)

    image = _image_metadata(arguments.image)
    if image["alvin_commit"] != EXPECTED_COMMIT:
        parser.error(
            f"image contains AlVin {image['alvin_commit']}, expected {EXPECTED_COMMIT}"
        )

    context = {
        name: _sha256(ALVIN_CONTEXT / name)
        for name in CONTEXT_FILES
    }
    container_name = f"k7-alvin-{arguments.case}-{time.time_ns()}"
    command = [
        "docker",
        "run",
        "--rm",
        "--name",
        container_name,
        "--network=none",
        "--volume",
        f"{output_directory}:/opt/AlVin/build/output",
    ]
    if arguments.cpus is not None:
        command.extend(("--cpus", str(arguments.cpus)))
    command.extend(
        [
            arguments.image,
            "-k",
            "rc7",
            "-qf",
            case["qf"],
            "-ip",
            "-noimage",
        ]
    )
    if "-maxv" in command:
        raise AssertionError("the live sweep must never use -maxv")

    request = {
        "schema": "coxeter.k7_alvin.live_request.v1",
        "case": arguments.case,
        **case,
        "command_argv": command,
        "contains_maxv": False,
        "image": image,
        "expected_alvin_commit": EXPECTED_COMMIT,
        "runner_sha256": _sha256(Path(__file__)),
        "context_sha256": context,
        "started_at_utc": _utc_now(),
        "external_timeout_seconds": arguments.external_timeout_seconds,
        "docker_cpu_limit": arguments.cpus,
    }
    (case_directory / "request.json").write_text(
        json.dumps(request, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    started = time.monotonic()
    externally_stopped: str | None = None
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert process.stdout is not None and process.stderr is not None
    with (case_directory / "stdout.txt").open("wb") as stdout_file, (
        case_directory / "stderr.txt"
    ).open("wb") as stderr_file:
        stdout_thread = threading.Thread(
            target=_pump,
            args=(
                process.stdout,
                stdout_file,
                sys.stdout.buffer if arguments.echo_program_output else None,
            ),
        )
        stderr_thread = threading.Thread(
            target=_pump,
            args=(
                process.stderr,
                stderr_file,
                sys.stderr.buffer if arguments.echo_program_output else None,
            ),
        )
        stdout_thread.start()
        stderr_thread.start()
        try:
            while True:
                try:
                    return_code = process.wait(timeout=30)
                    break
                except subprocess.TimeoutExpired:
                    elapsed = time.monotonic() - started
                    print(
                        f"[{_utc_now()}] {arguments.case} still running; "
                        f"elapsed={elapsed:.1f}s stdout_bytes={stdout_file.tell()}",
                        flush=True,
                    )
                    if (
                        arguments.external_timeout_seconds is not None
                        and elapsed >= arguments.external_timeout_seconds
                    ):
                        externally_stopped = "external_timeout"
                        _stop_owned_container(container_name)
        except KeyboardInterrupt:
            externally_stopped = "runner_interrupted"
            _stop_owned_container(container_name)
            return_code = process.wait()
        stdout_thread.join()
        stderr_thread.join()

    elapsed = time.monotonic() - started
    stdout_path = case_directory / "stdout.txt"
    stderr_path = case_directory / "stderr.txt"
    stdout = stdout_path.read_text(encoding="utf-8", errors="replace")
    printed_qf = _printed_qf(stdout)
    coxiter_files = sorted(output_directory.glob("*.coxiter"))
    graph: dict[str, Any] | None = None
    if len(coxiter_files) == 1:
        facets, dimension, matrix = _parse_coxiter(coxiter_files[0])
        graph = {
            "coxiter_file": coxiter_files[0].name,
            "facets": facets,
            "dimension": dimension,
            "coxeter_matrix": matrix,
        }

    algorithm_ended = "Algorithm ended" in stdout
    vertices_at_infinity_match = re.search(
        r"^Number of vertices at infinity:\s*(\d+)\s*$", stdout, re.MULTILINE
    )
    vertices_at_infinity = (
        int(vertices_at_infinity_match.group(1))
        if vertices_at_infinity_match is not None
        else None
    )
    finite_volume_reported = (
        algorithm_ended and "Information about the polyhedron:" in stdout
    )
    compact_reported = finite_volume_reported and vertices_at_infinity == 0
    complete_termination = (
        externally_stopped is None
        and return_code == 0
        and algorithm_ended
        and graph is not None
    )
    qf_exact_order = printed_qf is not None and _canonical_qf(printed_qf) == _canonical_qf(
        case["qf"]
    )
    qf_same_multiset = printed_qf is not None and sorted(_canonical_qf(printed_qf)) == sorted(
        _canonical_qf(case["qf"])
    )
    checks: dict[str, bool] = {
        "complete_termination": complete_termination,
        "printed_qf_no_squarefactor_change": qf_same_multiset,
    }
    if "expected_facets" in case:
        checks["expected_facets"] = graph is not None and graph["facets"] == case["expected_facets"]
    if "expected_matrix" in case:
        checks["expected_matrix"] = graph is not None and tuple(
            tuple(row) for row in graph["coxeter_matrix"]
        ) == case["expected_matrix"]

    verifier: dict[str, Any] | None = None
    if complete_termination and printed_qf is not None and len(coxiter_files) == 1:
        verifier_path = case_directory / "independent_verifier.json"
        labelg = shutil.which("labelg")
        verifier_command = [
            sys.executable,
            str(REPOSITORY / "tools" / "verify_cubic7_alvin.py"),
            "--coxiter",
            str(coxiter_files[0]),
            "--roots",
            str(stdout_path),
            "--qf",
            printed_qf,
        ]
        if labelg is not None:
            verifier_command.extend(("--labelg", labelg))
        verifier_command.extend(("--output", str(verifier_path)))
        verified = subprocess.run(
            verifier_command,
            cwd=REPOSITORY,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False,
        )
        verifier_payload = (
            json.loads(verifier_path.read_text(encoding="utf-8"))
            if verifier_path.exists()
            else {"passes": False, "error": verified.stderr.strip()}
        )
        verifier = {
            "command_argv": verifier_command,
            "return_code": verified.returncode,
            "passes": verified.returncode == 0 and verifier_payload.get("passes") is True,
            "report": verifier_path.name,
            "geometry": verifier_payload.get("geometry"),
            "scope": verifier_payload.get("scope"),
        }
        checks["independent_exact_replay"] = verifier["passes"]

    result = {
        "schema": "coxeter.k7_alvin.live_result.v1",
        "case": arguments.case,
        "kind": case["kind"],
        "qf_requested": case["qf"],
        "qf_printed": printed_qf,
        "printed_qf_exact_input_order": qf_exact_order,
        "printed_qf_same_diagonal_multiset": qf_same_multiset,
        "return_code": return_code,
        "externally_stopped": externally_stopped,
        "elapsed_seconds": elapsed,
        "finished_at_utc": _utc_now(),
        "algorithm_ended_marker": algorithm_ended,
        "complete_termination": complete_termination,
        "finite_volume_reported_by_alvin": finite_volume_reported,
        "vertices_at_infinity_reported_by_alvin": vertices_at_infinity,
        "compact_reported_by_alvin": compact_reported,
        "graph": graph,
        "independent_verifier": verifier,
        "checks": checks,
        "acceptance_passes": all(checks.values()),
        "stdout_sha256": _sha256(stdout_path),
        "stderr_sha256": _sha256(stderr_path),
        "output_artifacts": [
            {
                "path": str(path.relative_to(case_directory)),
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in sorted(output_directory.iterdir())
            if path.is_file()
        ],
    }
    (case_directory / "result.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2, sort_keys=True), flush=True)
    return 0 if result["acceptance_passes"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
