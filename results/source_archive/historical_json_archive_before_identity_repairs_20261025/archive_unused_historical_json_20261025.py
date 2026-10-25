#!/usr/bin/env python3
"""Byte-preserving streamed archive of ONE unused historical JSON result.

This is neither publication packing nor a scientific/Lean/TeX launcher. It
does not change the40GiB/6GiB proof admission gate. It retains the original
until a streamed gzip round trip matches its fixed content identity.
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
    'fresh_sqrt5_rank13_all_selected_nonade_separation_20261017':
        (1283771506, '4328a9600faec034bd59069e609a9d9323aceea4b9cc00383472c85f883c9437'),
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
    last = 0
    owner = os.getppid()
    group = os.getpgrp()

    def guard(force=False):
        nonlocal peak, last
        require(os.getppid() == owner, 'Storage owner exited')
        require(time.monotonic()-started < 59, 'Storage-maintenance59s limit')
        now = time.monotonic()
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

    initial = guard(True)  # Monitoring access is mandatory before touching files.
    with (ROOT/'results/finite57_terminal_adapter_inputs_small_affine_20261025.json').open('rb') as inp:
        manifest = inp.read(65537)
    require(hashlib.sha256(manifest).hexdigest() ==
            '8a0973d7b7ef08231e199bee9f04f45c6a8dc01ba66ca8de8acb2d8e513ee7b7',
            'Current terminal manifest changed')
    require(not any(x['path'] == 'results/'+name+'.json'
                    for x in json.loads(manifest)['inputs'] for name in FILES),
            'A current terminal input must not be archived')
    archive = ROOT/'formal/.proof-batches/historical_json_archive_20261025'
    receipt = ROOT/'results/historical_json_streamed_archive_20261025.json'
    require(not archive.exists() and not receipt.exists(), 'New maintenance targets required')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        sources = []
        for name, (size, digest) in FILES.items():
            p = ROOT/f'results/{name}.json'
            require(p.is_file() and not p.is_symlink() and p.stat().st_size == size,
                    'Exact unused historical-result path/size required')
            git = subprocess.run(['git', '-c', 'core.packedGitWindowSize=8m',
                '-c', 'core.packedGitLimit=16m', '-c', 'core.deltaBaseCacheLimit=4m',
                'ls-files', '--error-unmatch', '--', str(p.relative_to(ROOT))],
                cwd=ROOT, capture_output=True, text=True, timeout=3)
            require(git.returncode == 1 and not git.stdout,
                    'Only this untracked historical result may be archived')
            sources.append((name, p, size, digest, identity(p)))
        archive.mkdir()
        records = []
        for name, source, size, expected, before in sources:
            target = archive/(name+'.json.gz')
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
        # Verify all archive streams before removing any original.
        # This immutable byte manifest also survives interruption during removal.
        (archive/'manifest.json').write_text(json.dumps(records, indent=2)+'\n')
        for name, source, size, expected, before in sources:
            guard(True)
            require(identity(source) == before, 'Object changed before removal')
        for name, source, size, expected, before in sources:
            source.unlink()
        final = guard(True)
        result = {'status': 'BYTE_PRESERVING_UNUSED_HISTORICAL_JSON_ARCHIVE',
                  'records': records, 'initial': initial, 'final': final,
                  'peak_group_rss_bytes': peak,
                  'seconds': round(time.monotonic()-started, 3),
                  'original_raw_path_archived': True,
                  'scientific_byte_contents_preserved': True,
                  'current_terminal_manifest_110_inputs_unchanged': True,
                  'scientific_or_lean_execution': False,
                  'proof_admission_policy_changed': False,
                  'restore': 'Stream gzip decompression to its original results path; verify the recorded SHA256 before use. The current110-input terminal manifest does not use this historical result.',
                  'monitor_precision': 'Sampled256MiB group RSS and20GiB/3GiB storage floors every0.2s; each1MiB block checks owner/deadline. No mathematical admission policy changed.'}
        with receipt.open('x') as out:
            json.dump(result, out, indent=2)
            out.write('\n')
        print(json.dumps(result))


if __name__ == '__main__':
    main()
