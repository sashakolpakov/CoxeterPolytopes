#!/usr/bin/env python3
"""Run one repository checker with a live59s/1GiB watchdog.

Example: python3 -B scripts/run_bounded_research_check_20261018.py --
         scripts/check_compact_arithmetic20_composition_20261018.py --verify
Use --parent-lock for older checkers which do not acquire the shared lock.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import selectors
import signal
import subprocess
import sys
import time
import build_illustrated_ra_manuscript_20261018 as monitor
import run_lean_batches_signal_safe_20261018 as safe

ROOT = Path(__file__).resolve().parents[1]

# The child retains an inherited lock and stops its own process group if its
# owner disappears. This closes the abrupt-parent-exit hole for legacy scripts.
BOOTSTRAP = r'''
import os, runpy, signal, sys, threading, time
owner = int(os.environ['COXETER_CHECK_OWNER'])
def watch_owner():
    while True:
        time.sleep(.2)
        if os.getppid() != owner:
            os.killpg(os.getpgrp(), signal.SIGKILL)
threading.Thread(target=watch_owner, daemon=True).start()
signal.alarm(59)
sys.argv = sys.argv[1:]
sys.path.insert(0, os.path.dirname(sys.argv[0]))
runpy.run_path(sys.argv[0], run_name='__main__')
'''


def main():
    owner = os.getppid()
    runner_original = monitor.digest(Path(__file__))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger', help='New compact execution ledger in results/')
    parser.add_argument('--parent-lock', action='store_true',
                        help='Acquire and retain the lock for a non-locking checker')
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    assert command
    source = (ROOT/command[0]).resolve()
    assert source.parent == ROOT/'scripts' and source.suffix == '.py' and source.is_file()
    ledger = (ROOT/args.ledger).resolve() if args.ledger else None
    if ledger:
        assert ledger.parent == ROOT/'results' and not ledger.exists()
    if (ROOT/'repro/RESEARCH_PAUSED').exists():
        raise RuntimeError('Research is paused; no checker was started')
    before = safe.resources()
    monitor.rss()  # Check watchdog access before spawning a child.
    original = monitor.digest(source)
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, monitor.interrupt)
    started, peak = time.monotonic(), 0
    output, child, lock = bytearray(), None, None
    signal.alarm(59)
    try:
        if args.parent_lock:
            lock = (ROOT/'formal/.proof-batches/exclusive.lock').open('a')
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        env = dict(os.environ)
        env['COXETER_CHECK_OWNER'] = str(os.getpid())
        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS',
                     'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
            env[name] = '1'
        assert os.getppid() == owner, 'Owner exited before checker launch'
        child = subprocess.Popen([sys.executable, '-B', '-c', BOOTSTRAP,
                                  str(source), *command[1:]],
            cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            start_new_session=True, pass_fds=(lock.fileno(),) if lock else ())
        with selectors.DefaultSelector() as poll:
            poll.register(child.stdout, selectors.EVENT_READ)
            while poll.get_map() or child.poll() is None:
                for key, _ in poll.select(.1):
                    data = os.read(key.fd, 65536)
                    if data:
                        output.extend(data)
                        assert len(output) < 524288, 'Checker output limit'
                    else:
                        poll.unregister(key.fileobj)
                peak = max(peak, monitor.rss())
                assert peak < 1024**3, 'Process-tree RSS limit'
                assert time.monotonic()-started < 59, 'Wall-clock limit'
                assert os.getppid() == owner, 'Owner exited'
                safe.resources(start=False)
        assert child.wait() == 0, output[-7000:].decode(errors='replace')
        safe.stop_group(child, .3)
        child = None
        assert monitor.digest(source) == original, 'Checker changed during execution'
        assert monitor.digest(Path(__file__)) == runner_original, 'Runner changed during execution'
        text = output.decode()
        result = {'status': 'PASS_BOUNDED_RESEARCH_CHECK', 'command': command,
            'elapsed_seconds': round(time.monotonic()-started, 3),
            'peak_process_tree_rss_bytes': peak, 'initial_resources': before,
            'wall_limit_seconds': 59, 'rss_limit_bytes': 1024**3, 'threads': 1,
            'stdout_bytes': len(output), 'stdout': text,
            'runner_sha256': runner_original,
            'checker_sha256': original,
            'lock_owner': 'wrapper' if args.parent_lock else 'checker'}
        if ledger:
            with ledger.open('x') as stream:
                json.dump(result, stream, indent=2)
                stream.write('\n')
        print(text, end='' if text.endswith('\n') else '\n')
        print(json.dumps({k: v for k, v in result.items() if k != 'stdout'}))
    finally:
        signal.alarm(0)
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(signum, signal.SIG_IGN)
        try:
            if child is not None:
                safe.stop_group(child, .3)
        finally:
            if lock is not None:
                lock.close()


if __name__ == '__main__':
    main()
