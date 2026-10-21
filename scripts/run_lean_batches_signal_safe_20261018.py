#!/usr/bin/env python3
"""Capped Lean batches with cleanup that needs no signal-zero group probe.

Keep the earlier runner byte-identical for frozen execution provenance.
All its resource, cache, lock and compilation checks remain in use.
"""
import os
import signal
import subprocess
import run_lean_batches_20261018 as base

resources = base.resources
STATE = base.STATE


def stop_group(child, grace_seconds=5):
    """Reap the leader and kill remaining owned group members after TERM."""
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        child.wait(timeout=grace_seconds)
    except subprocess.TimeoutExpired:
        pass
    # A leader may exit while a descendant ignores TERM. Always send the
    # final signal to its owned group; a signal-zero probe is unnecessary.
    try:
        os.killpg(child.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    child.wait()
    if child.stdout is not None:
        child.stdout.close()


def guarded_run(*args, **kwargs):
    base.stop_group = stop_group
    return base.guarded_run(*args, **kwargs)


def main():
    base.stop_group = stop_group
    base.main()


if __name__ == '__main__':
    main()
