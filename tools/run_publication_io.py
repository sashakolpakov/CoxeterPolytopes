#!/usr/bin/env python3
"""Monitor small-memory publication I/O; never a proof/Lean runner.

macOS does not reliably support lowering RLIMIT_DATA. This parent instead
measures the entire owned process group, stops it above 256 MiB RSS, and
keeps the existing mathematical proof admission policy untouched. Only
publication packing, restoration, its synthetic tests and Sphinx are allowed.
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import re
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
GIB = 1024**3


def resources() -> dict:
    free = shutil.disk_usage(ROOT).free
    if free < 20 * GIB:
        raise RuntimeError("Publication I/O stopped: disk free below 20 GiB")
    ram = None
    if shutil.which("vm_stat"):
        output = subprocess.run(["vm_stat"], capture_output=True, text=True,
                                check=True, timeout=10).stdout
        page = int(re.search(r"page size of (\d+) bytes", output).group(1))
        ram = page * sum(int(re.search(re.escape(k) + r":\s*(\d+)", output).group(1))
                         for k in ("Pages free", "Pages inactive", "Pages speculative"))
        if ram < 3 * GIB:
            raise RuntimeError("Publication I/O stopped: available RAM below 3 GiB")
    return {"disk_free_bytes": free, "ram_available_bytes": ram}


def rss(group: int) -> int:
    output = subprocess.run(["ps", "-axo", "pgid=,rss="], check=True,
                            capture_output=True, text=True, timeout=10).stdout
    return sum(int(fields[1]) * 1024 for line in output.splitlines()
               if len(fields := line.split()) == 2 and int(fields[0]) == group)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seconds", type=int, default=1800)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    allowed = {"tools/prepare_public_data.py", "scripts/fetch_public_data.py",
               "tools/audit_public_tree.py",
               "tools/upload_public_data.py",
               "tools/audit_public_data_contents.py",
               "tools/revise_public_data_reference_exclusion.py",
               "tests/test_public_data_io_20261021.py"}
    suffix = command[1:]
    if suffix[:1] == ["-B"]:
        suffix = suffix[1:]
    if not command or Path(command[0]).name not in {"python3", "python", "python3.13", "python3.14"} or not suffix or not (suffix[0] in allowed or suffix[:2] == ["-m", "sphinx"]):
        parser.error("Only the listed publication I/O tools and Sphinx are allowed")
    if args.seconds <= 0:
        parser.error("--seconds must be positive")
    # Establish monitoring access before starting any child.
    start_resources = resources()
    rss(os.getpgrp())
    environment = dict(os.environ, COXETER_PUBLICATION_IO_GUARD="1",
                       PYTHONDONTWRITEBYTECODE="1", OMP_NUM_THREADS="1",
                       OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    started = time.monotonic()
    child = subprocess.Popen(command, cwd=ROOT, env=environment,
                             start_new_session=True)
    peak, reason = 0, None
    try:
        while child.poll() is None:
            resources()
            peak = max(peak, rss(child.pid))
            if peak > 256 * 1024**2:
                raise RuntimeError("Publication I/O exceeded 256 MiB process-group RSS")
            if time.monotonic() - started > args.seconds:
                raise RuntimeError("Publication I/O exceeded its time limit")
            time.sleep(1)
    except BaseException as error:
        reason = str(error) or type(error).__name__
        if child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
    print(json.dumps({"publication_io_only": True, "command": command,
                      "exit_code": child.returncode, "stop_reason": reason,
                      "observed_peak_group_rss_bytes": peak,
                      "rss_sampling_seconds": 1, "start_resources": start_resources,
                      "seconds": round(time.monotonic() - started, 3)}), flush=True)
    return 1 if reason else child.returncode


if __name__ == "__main__":
    raise SystemExit(main())
