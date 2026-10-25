#!/usr/bin/env python3
"""Byte-preserving storage maintenance for TWO unused generated Lean objects.

This is neither publication packing nor a scientific/Lean/TeX launcher. It
does not change the40GiB/6GiB proof admission gate. It retains the originals
until a streamed gzip round trip matches their fixed content identities.
"""
import fcntl
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
GIB = 1024**3
FILES = {
    'Ideal59': (107735496, '4f941c500c527933c6de005af23c373718f2f830dcf991e5c24fed8e72ff6588'),
    'Finite59': (63746968, 'ed42354fb5332a5fb470fa5a21529971311edbc7fed1bfb06b66f9627153b500'),
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def main():
    started = time.monotonic()
    peak = 0
    group = os.getpgrp()

    def guard():
        nonlocal peak
        require(time.monotonic()-started < 59, 'Storage-maintenance59s limit')
        disk = shutil.disk_usage(ROOT).free
        require(disk >= 20*GIB, 'Storage maintenance requires20GiB free')
        vm = subprocess.run(['vm_stat'], capture_output=True, text=True,
                            check=True, timeout=3).stdout
        page = int(re.search(r'page size of (\d+) bytes', vm).group(1))
        ram = page*sum(int(re.search(re.escape(k)+r':\s*(\d+)', vm).group(1))
                       for k in ('Pages free', 'Pages inactive', 'Pages speculative'))
        require(ram >= 3*GIB, 'Storage maintenance requires3GiB RAM')
        table = subprocess.run(['ps', '-axo', 'pgid=,rss='], capture_output=True,
                               text=True, check=True, timeout=3).stdout
        rss = sum(int(row[1])*1024 for line in table.splitlines()
                  if len(row := line.split()) == 2 and int(row[0]) == group)
        peak = max(peak, rss)
        require(rss < 256*1024**2, 'Storage-maintenance group RSS limit')
        return {'free_disk_bytes': disk, 'available_ram_bytes': ram}

    initial = guard()  # Monitoring access is mandatory before touching files.
    archive = ROOT/'formal/.proof-batches/cache_archive_20261025'
    receipt = ROOT/'results/lean_generated_cache_byte_preserving_cleanup_20261025.json'
    require(not archive.exists() and not receipt.exists(), 'New maintenance targets required')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        sources = []
        for name, (size, digest) in FILES.items():
            p = ROOT/f'formal/.lake/build/lib/lean/CoxeterProofs/Generated/{name}.olean'
            require(p.is_file() and not p.is_symlink() and p.stat().st_size == size,
                    'Exact unused generated-object path/size required')
            git = subprocess.run(['git', '-c', 'core.packedGitWindowSize=8m',
                '-c', 'core.packedGitLimit=16m', '-c', 'core.deltaBaseCacheLimit=4m',
                'ls-files', '--error-unmatch', '--', str(p.relative_to(ROOT))],
                cwd=ROOT, capture_output=True, text=True, timeout=3)
            require(git.returncode == 1 and not git.stdout,
                    'Only untracked generated cache objects may be archived')
            sources.append((name, p, size, digest, identity(p)))
        archive.mkdir()
        records = []
        for name, source, size, expected, before in sources:
            target = archive/(name+'.olean.gz')
            original_hash = hashlib.sha256()
            count = 0
            with source.open('rb') as inp, target.open('xb') as out:
                with gzip.GzipFile(filename='', mode='wb', fileobj=out,
                                   compresslevel=1, mtime=0) as packed:
                    for block in iter(lambda: inp.read(1024**2), b''):
                        guard()
                        original_hash.update(block)
                        count += len(block)
                        packed.write(block)
            require(count == size and original_hash.hexdigest() == expected,
                    'Original object content identity mismatch; original retained')
            restored_hash = hashlib.sha256()
            count = 0
            with gzip.open(target, 'rb') as inp:
                for block in iter(lambda: inp.read(1024**2), b''):
                    guard()
                    restored_hash.update(block)
                    count += len(block)
            require(count == size and restored_hash.hexdigest() == expected,
                    'Gzip round trip failed; original retained')
            require(identity(source) == before, 'Original object changed; retained')
            record = {'original': str(source.relative_to(ROOT)),
                      'archive': str(target.relative_to(ROOT)), 'bytes': size,
                      'sha256': expected, 'archive_bytes': target.stat().st_size,
                      'roundtrip_verified': True}
            records.append(record)
        # Verify BOTH archive streams before removing EITHER unused object.
        for name, source, size, expected, before in sources:
            guard()
            require(identity(source) == before, 'Object changed before removal')
        for name, source, size, expected, before in sources:
            source.unlink()
        final = guard()
        result = {'status': 'BYTE_PRESERVING_UNUSED_LEAN_CACHE_ARCHIVE',
                  'records': records, 'initial': initial, 'final': final,
                  'peak_group_rss_bytes': peak,
                  'seconds': round(time.monotonic()-started, 3),
                  'source_files_or_proof_evidence_deleted': False,
                  'scientific_or_lean_execution': False,
                  'proof_admission_policy_changed': False,
                  'restore': 'Stream gzip decompression to each original cache path; verify its recorded SHA256 before use.'}
        with receipt.open('x') as out:
            json.dump(result, out, indent=2)
            out.write('\n')
        print(json.dumps(result))


if __name__ == '__main__':
    main()
