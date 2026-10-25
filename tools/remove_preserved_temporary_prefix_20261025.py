#!/usr/bin/env python3
"""Remove one stale temporary prefix only after its complete byte preservation.

The full historical gzip is independently checked in both compressed and
restored form. Every byte of the temporary file must equal its prefix.
This is bounded storage I/O, not a mathematical check or a proof launcher.
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
TEMP = 'results/.BC.T_YvOsCu'
TEMP_BYTES = 175767552
ARCHIVE = ('formal/.proof-batches/historical_json_archive_20261025/'
           'fresh_sqrt5_rank13_all_selected_nonade_separation_20261017.json.gz')
ARCHIVE_BYTES = 89001268
ARCHIVE_SHA = 'd969c8721beb8417f2b4283dc17c439bc6d2371158e740222ac0e4d9bf1e8094'
RESTORED_BYTES = 1283771506
RESTORED_SHA = '4328a9600faec034bd59069e609a9d9323aceea4b9cc00383472c85f883c9437'
VERIFIED = 'results/temporary_prefix_preservation_verified_20261025.json'
REMOVED = 'results/temporary_prefix_removal_20261025.json'
GIB = 1024**3


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def identity(path):
    s = path.stat()
    return s.st_dev, s.st_ino, s.st_size, s.st_mtime_ns, s.st_ctime_ns


def main():
    started = time.monotonic()
    owner, group = os.getppid(), os.getpgrp()
    last, peak = 0, 0

    def guard(force=False):
        nonlocal last, peak
        now = time.monotonic()
        require(os.getppid() == owner, 'Storage owner exited')
        require(now-started < 59, 'Storage-maintenance59s limit')
        if not force and last and now-last < .2:
            return
        last = now
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

    initial = guard(True)
    temporary, archive = ROOT/TEMP, ROOT/ARCHIVE
    for name in (TEMP, ARCHIVE):
        path = ROOT/name
        require(not any(p.is_symlink() for p in [path, *path.parents]
                        if p.is_relative_to(ROOT)), 'No symlink paths')
        require(path.resolve() == path and path.is_file(), 'Exact existing path')
    require(temporary.stat().st_size == TEMP_BYTES, 'Fixed temporary size')
    require(archive.stat().st_size == ARCHIVE_BYTES, 'Fixed preserved archive size')
    require(not (ROOT/VERIFIED).exists() and not (ROOT/REMOVED).exists(),
            'Refuse overwriting preservation evidence')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        original_identity, archive_identity = identity(temporary), identity(archive)
        git = subprocess.run(['git', '-c', 'core.packedGitWindowSize=8m',
            '-c', 'core.packedGitLimit=16m', '-c', 'core.deltaBaseCacheLimit=4m',
            'ls-files', '--error-unmatch', '--', TEMP], cwd=ROOT,
            capture_output=True, text=True, timeout=3)
        require(git.returncode == 1 and not git.stdout, 'Only an untracked temporary file')
        digest = hashlib.sha256()
        with archive.open('rb') as inp:
            for block in iter(lambda: inp.read(1024**2), b''):
                guard()
                digest.update(block)
        require(digest.hexdigest() == ARCHIVE_SHA, 'Preserved archive digest mismatch')
        restored_digest, prefix_digest = hashlib.sha256(), hashlib.sha256()
        total, compared = 0, 0
        with gzip.open(archive, 'rb') as preserved, temporary.open('rb') as candidate:
            for block in iter(lambda: preserved.read(1024**2), b''):
                guard()
                restored_digest.update(block)
                total += len(block)
                length = min(len(block), TEMP_BYTES-compared)
                if length:
                    prefix = block[:length]
                    require(candidate.read(length) == prefix,
                            'Temporary contains bytes absent from preserved prefix')
                    prefix_digest.update(prefix)
                    compared += length
            require(candidate.read(1) == b'', 'Uncompared temporary tail')
        require(total == RESTORED_BYTES and restored_digest.hexdigest() == RESTORED_SHA,
                'Complete historical archive round trip failed')
        require(compared == TEMP_BYTES, 'Every temporary byte must be compared')
        require(identity(temporary) == original_identity and identity(archive) == archive_identity,
                'Storage inputs changed during verification')
        verification = {
            'status': 'EVERY_TEMPORARY_BYTE_PRESERVED_IN_VERIFIED_FULL_ARCHIVE',
            'temporary': TEMP, 'temporary_bytes': TEMP_BYTES,
            'temporary_sha256': prefix_digest.hexdigest(),
            'archive': ARCHIVE, 'archive_bytes': ARCHIVE_BYTES,
            'archive_sha256': ARCHIVE_SHA, 'full_restored_bytes': RESTORED_BYTES,
            'full_restored_sha256': RESTORED_SHA,
            'all_temporary_bytes_equal_initial_full_archive_prefix': True,
            'mathematical_verification_claimed': False,
            'restore': 'Stream exactly the first175767552 decompressed archive bytes to the temporary path; verify temporary_sha256.'}
        with (ROOT/VERIFIED).open('x') as out:
            json.dump(verification, out, indent=2)
            out.write('\n')
        guard(True)
        require(identity(temporary) == original_identity and identity(archive) == archive_identity,
                'Storage input changed before duplicate removal')
        temporary.unlink()
        final = guard(True)
        result = {'status': 'REMOVED_BYTE_VERIFIED_STALE_TEMPORARY_PREFIX',
                  'preservation_receipt': VERIFIED, 'removed_path': TEMP,
                  'logical_bytes_reclaimed': TEMP_BYTES, 'initial': initial, 'final': final,
                  'peak_group_rss_bytes': peak, 'seconds': round(time.monotonic()-started, 3),
                  'scientific_bytes_preserved': True, 'proof_resource_policy_changed': False,
                  'scientific_lean_or_tex_worker_started': False,
                  'source_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
        with (ROOT/REMOVED).open('x') as out:
            json.dump(result, out, indent=2)
            out.write('\n')
        print(json.dumps(result))


if __name__ == '__main__':
    main()
