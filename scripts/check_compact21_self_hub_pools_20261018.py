#!/usr/bin/env python3
"""Fresh exact replay of the self-hub pool diagnostic, without optimization."""
import fcntl
from hashlib import file_digest
import json
import os
import resource
import signal
import sys
import threading
import time
import fresh_compact21_self_hub_pool_dual_probe_20261018 as M


def main():
    started = time.monotonic()
    scale = 1 if sys.platform == 'darwin' else 1024
    owner, pid, stop = os.getppid(), os.getpid(), threading.Event()
    def monitor():
        next_resources = 0
        while not stop.wait(.1):
            try:
                assert resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*scale < 1024**3
                assert os.getppid() == owner
                if time.monotonic() >= next_resources:
                    M.C.H.S.G.resources(start=False); next_resources = time.monotonic()+1
            except BaseException:
                os.kill(pid, signal.SIGTERM); return
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, M.C.H.S.G.interrupt)
    signal.alarm(59)
    with (M.ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = M.C.H.S.G.resources()
        watcher = threading.Thread(target=monitor, daemon=True); watcher.start()
        try:
            target = M.ROOT/f'results/{M.TAG}.json'
            assert target.stat().st_size < 1024**2
            with target.open() as stream:
                saved = json.load(stream)
            for path, digest in (saved['sources'] | saved['supplemental_sources']).items():
                with (M.ROOT/path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
            fresh = M.run()
            for key, value in fresh.items():
                assert saved[key] == value, key
            M.C.H.S.G.resources(start=False)
        finally:
            stop.set(); watcher.join(timeout=2)
    signal.alarm(0)
    print(json.dumps({'status': 'PASS_FRESH_EXACT_SELF_HUB_POOL_REPLAY',
                      'strictly_breaking_columns': fresh['strictly_breaking_columns'],
                      'elapsed_seconds': round(time.monotonic()-started, 3),
                      'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*scale,
                      'initial_resources': before, 'numeric_optimizer_run': False}))


if __name__ == '__main__':
    main()
