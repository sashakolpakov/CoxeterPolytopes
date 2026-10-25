#!/usr/bin/env python3
"""Storage-only transparent compression of ONE fixed historical JSON file.

No JSON parsing, scientific worker, publication upload or proof-gate change.
The original path is replaced only after the copied data stream matches it.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import signal
import stat
import subprocess
import sys
import time
import types

ROOT = Path(__file__).resolve().parents[1]
RELATIVE = 'results/fresh_sqrt5_rank13_all_selected_nonade_separation_20261017.json'
SIZE = 1283771506
GIB = 1024**3

# Reuse the reviewed birth-identity monitor, without invoking its proof main.
MONITOR = ROOT/'scripts/monitor_owned_proof_tree_20261024.py'
require_monitor_sha = 'b07f84fbcc5b35e484696f0515b8e233ea975f13c2da6ff75c0d3b54b99ed5db'
with MONITOR.open('rb') as monitor_stream:
    monitor_bytes = monitor_stream.read(65537)
if len(monitor_bytes) > 65536 or hashlib.sha256(monitor_bytes).hexdigest() != require_monitor_sha:
    raise SystemExit('Reviewed owned-tree monitor identity required')
owned = types.ModuleType('_coxeter_storage_owned_monitor')
owned.__file__ = str(MONITOR)
sys.modules[owned.__name__] = owned
exec(compile(monitor_bytes, str(MONITOR), 'exec'), owned.__dict__)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def main():
    started = time.monotonic()
    parent = os.getppid()
    group = os.getpgrp()
    child = None
    tree = None
    cleanup = None
    peak = 0
    last = 0

    def guard(force=False):
        nonlocal peak, last
        now = time.monotonic()
        require(now-started < 59 and os.getppid() == parent,
                'Storage deadline or owner interruption')
        if not force and now-last < .2:
            return
        last = now
        free = shutil.disk_usage(ROOT).free
        require(free >= 20*GIB, 'Storage operation requires20GiB free')
        vm = subprocess.run(['vm_stat'], capture_output=True, text=True,
                            check=True, timeout=3).stdout
        page = int(re.search(r'page size of (\d+) bytes', vm).group(1))
        ram = page*sum(int(re.search(re.escape(k)+r':\s*(\d+)', vm).group(1))
                       for k in ('Pages free', 'Pages inactive', 'Pages speculative'))
        require(ram >= 3*GIB, 'Storage operation requires3GiB RAM')
        table = subprocess.run(['ps', '-axo', 'pgid=,rss='], capture_output=True,
                               text=True, check=True, timeout=3).stdout
        groups = {group}
        if tree is not None and cleanup is None:
            groups.update(p.pgid for p in tree.sample().values())
        rss = sum(int(row[1])*1024 for line in table.splitlines()
                  if len(row := line.split()) == 2 and int(row[0]) in groups)
        peak = max(peak, rss)
        require(rss < 256*1024**2, 'Storage process-group RSS limit')

    def digest(path):
        h = hashlib.sha256()
        size = 0
        with path.open('rb') as inp:
            for block in iter(lambda: inp.read(1024**2), b''):
                guard()
                h.update(block)
                size += len(block)
        return size, h.hexdigest()

    def interrupted(signum, frame):
        raise RuntimeError('Storage operation interrupted by signal'+str(signum))

    source = ROOT/RELATIVE
    target = source.with_name(source.name+'.owned-apfs-copy-20261025')
    receipt = ROOT/'results/preserved_json_apfs_compression_20261025.json'
    guard(True)  # Required monitoring access before any data mutation.
    reader = owned.ProcessReader()
    reader.probe()
    require(not target.exists() and not receipt.exists(), 'New owned output paths required')
    require(source.is_file() and not source.is_symlink() and source.stat().st_size == SIZE,
            'Exact fixed historical file required')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = identity(source)
        info = source.stat()
        original = digest(source)
        require(original[0] == SIZE and identity(source) == before, 'Source changed during hashing')
        for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
            signal.signal(signum, interrupted)
        signal.alarm(59)
        try:
            env = dict(os.environ, DITTOABORT='1', OMP_NUM_THREADS='1',
                       OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
            env.pop('DITTO_TEST_OPTIONS', None)
            child = subprocess.Popen(['/usr/bin/ditto', '--nocache', '--noclone',
                '--hfsCompression', str(source), str(target)], env=env,
                start_new_session=True)
            tree = owned.OwnedTree(child, reader)
            while child.poll() is None:
                guard(True)
                time.sleep(.1)
            require(child.wait() == 0, 'ditto failed; original preserved')
            cleanup = tree.stop()
            require(not cleanup['signals'], 'Unexpected copy descendants; original preserved')
            child = None
            require(target.is_file() and not target.is_symlink(), 'Regular copied file required')
            copied_before = identity(target)
            require(digest(target) == original, 'Copied byte stream differs; original preserved')
            require(identity(target) == copied_before, 'Copy changed during verification')
            require(identity(source) == before, 'Original changed during copy; preserved')
            copied = target.stat()
            require((stat.S_IMODE(copied.st_mode), copied.st_uid, copied.st_gid) ==
                    (stat.S_IMODE(info.st_mode), info.st_uid, info.st_gid),
                    'Copied mode/owner mismatch; original preserved')
            require(copied.st_blocks < info.st_blocks, 'No allocation saving; original preserved')
            guard(True)
            require(MONITOR.read_bytes() == monitor_bytes, 'Monitor source changed during storage operation')
            require(identity(source) == before and identity(target) == copied_before,
                    'Original or copy changed before replacement; original preserved')
            os.replace(target, source)
            result = {'status': 'PRESERVED_JSON_TRANSPARENT_FILESYSTEM_COMPRESSION',
                      'source': RELATIVE, 'bytes': SIZE, 'sha256': original[1],
                      'allocated_bytes_before': info.st_blocks*512,
                      'allocated_bytes_after': source.stat().st_blocks*512,
                      'copied_byte_stream_verified': True,
                      'path_and_byte_contents_preserved': True,
                      'peak_group_rss_bytes': peak,
                      'owned_copy_process_cleanup': cleanup,
                      'seconds': round(time.monotonic()-started, 3),
                      'scientific_or_publication_execution': False,
                      'proof_admission_policy_changed': False,
                      'failure_after_replace': 'An interruption or receipt-write failure after atomic replacement can leave no receipt; the verified byte-identical source remains in place. No power-loss durability or hostile concurrent-filesystem guarantee is claimed.',
                      'monitor': 'Sampled256MiB group RSS and20GiB/3GiB storage floors; OS probes3s each; streaming1MiB blocks.'}
            with receipt.open('x') as out:
                json.dump(result, out, indent=2)
                out.write('\n')
            print(json.dumps(result))
        finally:
            signal.alarm(0)
            for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
                signal.signal(signum, signal.SIG_IGN)
            if tree is not None and cleanup is None:
                tree.stop()
            elif child is not None and child.poll() is None:
                child.terminate()
                try:
                    child.wait(timeout=3)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
            # A failed partial copy is an owned, unused storage product.
            if target.exists():
                target.unlink()


if __name__ == '__main__':
    main()
