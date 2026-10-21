#!/usr/bin/env python3
"""Build the included rank20 appendix with every declared source pinned.

The existing bounded builder owns the shared lock and resource gates.
This wrapper adds before/after binding of BOTH TeX files and the48-file
internal supplement. It starts no parallel job and imports no science.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT/'manuscripts/supplements/compact_arithmetic_eight_field19_sources_final_20261019.json'
MANIFEST_SHA = '4e3a560d4cebded5725c528652d68b12a385b59458b5a1a8be970d70deb4ae6c'
BUILDER = ROOT/'scripts/build_illustrated_ra_manuscript_20261018.py'
BUILDER_SHA = '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15'


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inspect_inputs():
    assert digest(MANIFEST) == MANIFEST_SHA
    assert digest(BUILDER) == BUILDER_SHA
    manifest = json.loads(MANIFEST.read_text())
    files = manifest['files']
    assert len(files) == 48
    for name, info in files.items():
        path = (ROOT/name).resolve()
        assert path.is_relative_to(ROOT) and path.is_file()
        assert path.stat().st_size == info['bytes'] and digest(path) == info['sha256'], name
    for name in ('manuscripts/compact_arithmetic20.tex',
                 'manuscripts/appendices/compact_arithmetic_eight_field19.tex'):
        assert name in files
    return {name: info['sha256'] for name, info in files.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--replace-unfrozen', action='store_true')
    parser.add_argument('--ledger', type=Path)
    args = parser.parse_args()
    ledger = (ROOT/args.ledger).resolve() if args.ledger else None
    if ledger:
        assert ledger.parent == ROOT/'results' and not ledger.exists()
    original = inspect_inputs()
    command = [sys.executable, '-B', str(BUILDER), 'manuscripts/compact_arithmetic20.tex']
    if args.check_only:
        command.append('--check-only')
    if args.replace_unfrozen:
        command.append('--replace-unfrozen')
    completed = subprocess.run(command, cwd=ROOT, text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=62)
    assert len(completed.stdout.encode()) < 524288
    if completed.returncode:
        print(completed.stdout[-5000:], end='')
        raise SystemExit(completed.returncode)
    build = json.loads(completed.stdout)
    assert build['status'] == 'PASS_TEX_BUILD'
    assert inspect_inputs() == original, 'An included proof source changed during build'
    record = {'status': 'PASS_INCLUDED_TEX_AND_48_SOURCE_BINDING',
              'source_pins': original, 'manifest_sha256': MANIFEST_SHA,
              'builder_sha256': BUILDER_SHA, 'wrapper_sha256': digest(Path(__file__)),
              'build': build, 'mathematical_verification_claimed_by_build': False}
    if ledger:
        with ledger.open('x') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')
    print(json.dumps({'status': record['status'], 'pages': build['pages'],
                      'warnings': build['warnings'], 'pdf_sha256': build['pdf_sha256'],
                      'bound_sources': len(original)}))


if __name__ == '__main__':
    main()
