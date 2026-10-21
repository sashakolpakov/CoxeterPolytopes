#!/usr/bin/env python3
"""Run one command with sampled process-tree RSS and wall-time ceilings.

macOS may expose ``RLIMIT_AS`` while rejecting attempts to lower it.  This
wrapper therefore watches the resident memory of the child and, when the OS
permits process enumeration, all of its descendants from the parent.  In a
managed macOS sandbox it reports ``scope=direct`` and the command must itself
be the memory-owning process.  If either ceiling is crossed, only the process
group created by this wrapper is killed.

Exit status 124 means timeout; 125 means the RSS ceiling was crossed.  On a
normal exit the child's status is returned unchanged.
"""

from __future__ import annotations

import argparse
import os
import signal
import subprocess
import sys
import time

import psutil


def _resident_tree_bytes(pid: int) -> tuple[int, str]:
    """Return sampled RSS and whether descendant enumeration was available."""

    try:
        process = psutil.Process(pid)
        total = process.memory_info().rss
        try:
            descendants = process.children(recursive=True)
        except (PermissionError, psutil.AccessDenied):
            # The managed macOS sandbox blocks the sysctl-based all-PID scan
            # used by psutil.children().  The spawned command is itself the
            # expensive process in the intended uses, so direct RSS remains a
            # real ceiling rather than silently disabling the watchdog.
            return total, "direct"
        for child in descendants:
            try:
                total += child.memory_info().rss
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                pass
        return total, "tree"
    except (psutil.NoSuchProcess, psutil.AccessDenied):
        return 0, "exited"


def _kill_group(pid: int) -> None:
    """Kill exactly the isolated process group started by this wrapper."""

    try:
        os.killpg(pid, signal.SIGKILL)
    except ProcessLookupError:
        pass


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--memory-gib", type=float, required=True)
    parser.add_argument("--timeout", type=float)
    parser.add_argument("--poll-seconds", type=float, default=0.05)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = list(args.command)
    if command and command[0] == "--":
        command.pop(0)
    if not command:
        parser.error("a command is required after --")
    if args.memory_gib <= 0:
        parser.error("--memory-gib must be positive")
    if args.timeout is not None and args.timeout <= 0:
        parser.error("--timeout must be positive")
    if args.poll_seconds <= 0:
        parser.error("--poll-seconds must be positive")

    limit = int(args.memory_gib * (1 << 30))
    started = time.monotonic()
    peak = 0
    process = subprocess.Popen(command, start_new_session=True)
    disposition = "child_exit"
    measurement_scope = "direct"
    while process.poll() is None:
        resident, sampled_scope = _resident_tree_bytes(process.pid)
        if sampled_scope != "exited":
            measurement_scope = sampled_scope
        peak = max(peak, resident)
        if peak > limit:
            disposition = "rss_limit"
            _kill_group(process.pid)
            break
        if args.timeout is not None and time.monotonic() - started > args.timeout:
            disposition = "timeout"
            _kill_group(process.pid)
            break
        time.sleep(args.poll_seconds)
    returncode = process.wait()
    elapsed = time.monotonic() - started
    print(
        "rss-watchdog: "
        f"disposition={disposition} "
        f"scope={measurement_scope} "
        f"peak_gib={peak / (1 << 30):.3f} "
        f"elapsed_seconds={elapsed:.3f}",
        file=sys.stderr,
    )
    if disposition == "timeout":
        return 124
    if disposition == "rss_limit":
        return 125
    return returncode


if __name__ == "__main__":
    raise SystemExit(main())
