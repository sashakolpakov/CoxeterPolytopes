"""Build or independently verify the H9/14 missing-edge graph manifest."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil
import subprocess
from typing import Sequence

from coxeter_polytopes.combinatorics.missing_edges import (
    build_manifest,
    verify_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = ROOT / "data" / "enumerations" / "h9_14_missing_edges.json"
GENG_COMMAND = ["geng", "-q", "14", "2:10"]


def _file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _run_geng(executable: str) -> tuple[list[str], dict[str, object]]:
    resolved_text = shutil.which(executable)
    if resolved_text is None:
        raise RuntimeError(f"cannot find nauty generator {executable!r}")
    resolved = Path(resolved_text).resolve()

    command = [str(resolved), *GENG_COMMAND[1:]]
    completed = subprocess.run(
        command,
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    version_run = subprocess.run(
        [str(resolved), "--version"],
        check=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    version = (version_run.stdout + version_run.stderr).decode("utf-8").strip()
    graph6_values = completed.stdout.decode("ascii").splitlines()
    provenance: dict[str, object] = {
        "name": "nauty geng",
        "version": version,
        "command": GENG_COMMAND,
        "resolved_command": command,
        "executable": str(resolved),
        "executable_sha256": _file_sha256(resolved),
        "stdout_sha256": sha256(completed.stdout).hexdigest(),
        "stderr": completed.stderr.decode("utf-8").strip(),
    }
    return graph6_values, provenance


def generate(output: Path, geng: str = "geng") -> dict[str, object]:
    graph6_values, provenance = _run_geng(geng)
    manifest = build_manifest(graph6_values, provenance)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def verify(path: Path) -> dict[str, object]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    return verify_manifest(manifest)


def main(argv: Sequence[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--geng", default="geng")
    parser.add_argument(
        "--verify-only",
        action="store_true",
        help="verify the archived output without invoking nauty",
    )
    args = parser.parse_args(argv)

    if args.verify_only:
        summary = verify(args.output)
        mode = "verified"
    else:
        manifest = generate(args.output, args.geng)
        summary = verify_manifest(manifest)
        mode = "generated and verified"
    print(
        json.dumps(
            {"mode": mode, "output": str(args.output), "counts": summary},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
