#!/usr/bin/env python3
"""Build the nine-field appendix package with before/after source binding.

The pinned generic builder alone owns the global lock, resource gates,
59-second/1-GiB TeX batch and intermediate cleanup. This wrapper starts
one such builder, bounds its output, and binds the full named supplement.
It imports no scientific code and never launches a Lean job.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import selectors
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT/'manuscripts/supplements/compact_arithmetic_nine_field19_sources_20261020.json'
MANIFEST_SHA = '5cf6afc51440ae1510e6b6e54ec96308749dd11713b29a866c5d9edb96065dae'
PREDECESSOR = ROOT/'manuscripts/supplements/compact_arithmetic_eight_field19_sources_final_20261019.json'
PREDECESSOR_SHA = '4e3a560d4cebded5725c528652d68b12a385b59458b5a1a8be970d70deb4ae6c'
BUILDER = ROOT/'scripts/build_illustrated_ra_manuscript_20261018.py'
BUILDER_SHA = '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15'
PARENT = 'manuscripts/compact_arithmetic20.tex'
REQUIRED_TEX = {PARENT,
    'manuscripts/appendices/compact_arithmetic_eight_field19.tex',
    'manuscripts/appendices/compact_arithmetic_k7_field19.tex'}
MAX_OUTPUT = 512*1024


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def inspect_inputs():
    assert isinstance(MANIFEST_SHA, str) and len(MANIFEST_SHA) == 64, \
        'Nine-field manifest is not sealed; no builder may start'
    assert digest(MANIFEST) == MANIFEST_SHA, 'Nine-field manifest changed'
    assert digest(PREDECESSOR) == PREDECESSOR_SHA, 'Eight-field manifest changed'
    assert digest(BUILDER) == BUILDER_SHA, 'Generic bounded builder changed'
    manifest = json.loads(MANIFEST.read_text())
    assert manifest['status'] == 'FINAL_SOURCE_BOUND_NINE_FIELD_PACKAGE'
    assert manifest['preserved_eight_field_manifest_sha256'] == PREDECESSOR_SHA
    files = manifest['files']
    assert len(files) == manifest['file_count'] and len(files) >= 48
    assert REQUIRED_TEX <= files.keys(), 'An included TeX source is unbound'
    predecessor = json.loads(PREDECESSOR.read_text())
    assert len(predecessor['files']) == 48
    for name, old in predecessor['files'].items():
        assert name in files, ('Missing old dependency', name)
        if name != PARENT:
            assert files[name]['bytes'] == old['bytes'], name
            assert files[name]['sha256'] == old['sha256'], name
    for name, info in files.items():
        assert not Path(name).is_absolute()
        path = (ROOT/name).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), name
        assert path.stat().st_size == info['bytes'], name
        assert digest(path) == info['sha256'], name
    return {name: info['sha256'] for name, info in files.items()}


def interrupted(signum, _frame):
    raise RuntimeError(f'Interrupted included-source build: signal {signum}')


def stop_builder(child):
    if child.poll() is None:
        # The generic builder handles SIGTERM and cleans its separately
        # grouped TeX child in a finally block. Allow that cleanup to finish.
        child.terminate()
        try:
            child.wait(timeout=3)
        except subprocess.TimeoutExpired:
            child.kill()
            child.wait(timeout=3)


def run_builder(command):
    output = bytearray()
    started = time.monotonic()
    child = subprocess.Popen(command, cwd=ROOT, stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT, start_new_session=True)
    try:
        with selectors.DefaultSelector() as poll:
            poll.register(child.stdout, selectors.EVENT_READ)
            while poll.get_map() or child.poll() is None:
                for key, _ in poll.select(.1):
                    data = os.read(key.fd, 65536)
                    if data:
                        assert len(output)+len(data) <= MAX_OUTPUT, \
                            'Bounded builder output exceeded512KiB'
                        output.extend(data)
                    else:
                        poll.unregister(key.fileobj)
                assert time.monotonic()-started < 63, 'Builder wrapper timeout'
        return child.wait(), bytes(output)
    finally:
        stop_builder(child)
        child.stdout.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--replace-unfrozen', action='store_true')
    parser.add_argument('--ledger', type=Path)
    args = parser.parse_args()
    ledger = (ROOT/args.ledger).resolve() if args.ledger else None
    if ledger:
        assert ledger.parent == ROOT/'results' and not ledger.exists()
    for sig in (signal.SIGTERM, signal.SIGINT):
        signal.signal(sig, interrupted)
    original = inspect_inputs()
    wrapper_original = digest(Path(__file__))
    command = [sys.executable, '-B', str(BUILDER), PARENT]
    if args.check_only:
        command.append('--check-only')
    if args.replace_unfrozen:
        command.append('--replace-unfrozen')
    try:
        returncode, output = run_builder(command)
    finally:
        assert inspect_inputs() == original, 'An included source changed during build'
        assert digest(Path(__file__)) == wrapper_original, 'Build wrapper changed'
    if returncode:
        print(output[-5000:].decode(errors='replace'), end='')
        raise SystemExit(returncode)
    build = json.loads(output)
    assert build['status'] == 'PASS_TEX_BUILD'
    record = {'status': 'PASS_INCLUDED_TEX_AND_NINE_FIELD_SOURCE_BINDING',
        'source_pins': original, 'manifest_sha256': MANIFEST_SHA,
        'preserved_eight_field_manifest_sha256': PREDECESSOR_SHA,
        'builder_sha256': BUILDER_SHA, 'wrapper_sha256': wrapper_original,
        'bound_sources': len(original), 'build': build,
        'mathematical_verification_claimed_by_build': False}
    if ledger:
        with ledger.open('x') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')
    print(json.dumps({'status': record['status'], 'pages': build['pages'],
        'warnings': build['warnings'], 'pdf_sha256': build['pdf_sha256'],
        'bound_sources': len(original)}))


if __name__ == '__main__':
    main()
