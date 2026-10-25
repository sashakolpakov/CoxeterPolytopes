#!/usr/bin/env python3
"""Focused disappearance-race regressions; no mathematical worker is imported."""
import os
import sys

if not __debug__ or sys.flags.optimize or os.environ.get('PYTHONOPTIMIZE', '') not in ('', '0'):
    raise SystemExit('ASSERTIONS_REQUIRED')

import ctypes
import errno
import json
import subprocess
from unittest.mock import patch

import monitor_owned_proof_tree_20261024 as monitor


def main():
    if sys.argv[1:] != ['--parent-lock']:
        raise SystemExit('Run through the bounded public entry with --parent-lock')
    reader = monitor.ProcessReader.__new__(monitor.ProcessReader)
    pid = 987654321
    cases = 0

    def run_case(error, length, answer, expected):
        nonlocal cases

        def proc(*args):
            ctypes.set_errno(error)
            return length

        reader.proc = proc
        with patch.object(monitor.sys, 'platform', 'darwin'), patch.object(
                monitor.subprocess, 'run', return_value=answer) as probe:
            try:
                value = reader.process(pid)
            except RuntimeError as exc:
                assert expected == 'REFUSED', (expected, str(exc))
                assert 'OWNED_TREE_IDENTITY_DENIED_OR_ABI' in str(exc)
            else:
                assert expected == 'ABSENT' and value is None
            if error == errno.EPERM and length == 0:
                probe.assert_called_once_with(
                    ['ps', '-p', str(pid), '-o', 'pid='],
                    capture_output=True, text=True, timeout=3)
            else:
                probe.assert_not_called()
        cases += 1

    def answer(code, stdout='', stderr=''):
        return subprocess.CompletedProcess([], code, stdout, stderr)

    run_case(errno.EPERM, 0, answer(1), 'ABSENT')
    run_case(errno.EPERM, 0, answer(0, str(pid)+'\n'), 'REFUSED')
    run_case(errno.EPERM, 0, answer(1, str(pid)+'\n'), 'REFUSED')
    run_case(errno.EPERM, 0, answer(1, '', 'permission denied'), 'REFUSED')
    run_case(errno.EPERM, 0, answer(2), 'REFUSED')
    run_case(errno.EPERM, 0, answer(0), 'REFUSED')
    run_case(errno.ESRCH, 0, None, 'ABSENT')
    run_case(errno.ENOENT, 0, None, 'ABSENT')
    run_case(0, 135, None, 'REFUSED')
    print(json.dumps({'status': 'PASS_DISAPPEARED_EPERM_REGRESSIONS',
                      'cases': cases, 'scope': 'mocked PID identity disappearance only',
                      'proof_resource_gates_changed': False,
                      'live_permission_denial_still_refused': True,
                      'mathematical_verification': False}))


if __name__ == '__main__':
    main()
