#!/usr/bin/env python3
"""Build one local manuscript in a59s/1GiB batch; discard TeX intermediates."""
import argparse
import fcntl
from hashlib import file_digest
import json
import os
from pathlib import Path
import re
import selectors
import shutil
import signal
import subprocess
import tempfile
import time
import run_lean_batches_signal_safe_20261018 as safe

ROOT = Path(__file__).resolve().parents[1]
LIMIT = 1024**3


def digest(path):
    with path.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def rss():
    output = subprocess.check_output(
        ['ps', '-axo', 'pid=,ppid=,rss='], text=True, timeout=3)
    rows = [tuple(map(int, line.split())) for line in output.splitlines()
            if line.strip()]
    owned = {os.getpid()}
    while True:
        more = {p for p, parent, _ in rows if parent in owned}-owned
        if not more:
            break
        owned.update(more)
    return 1024*sum(size for p, _, size in rows if p in owned)


def frozen(path):
    name = str(path.relative_to(ROOT))
    for manifest in (ROOT/'results').glob('*freeze.json'):
        with manifest.open() as stream:
            if any(name in line for line in stream):
                return True
    return False


def interrupt(signum, frame):
    raise RuntimeError(f'Interrupted bounded TeX batch: signal{signum}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--replace-unfrozen', action='store_true')
    args = parser.parse_args()
    source = (ROOT/args.source).resolve()
    assert source.parent == ROOT/'manuscripts' and source.suffix == '.tex'
    assert source.is_file()
    target = ROOT/'manuscripts/build'/f'{source.stem}.pdf'
    ledger = ROOT/'results'/f'{source.stem}_build_20261018.json'
    if not args.check_only:
        for path in (target, ledger):
            if path.exists():
                assert args.replace_unfrozen and not frozen(path), path
    engine = shutil.which('pdflatex')
    assert engine, 'pdflatex is required'
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, interrupt)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = safe.resources()
        rss()  # Verify watchdog availability before spawning.
        original = digest(source)
        owner, started, peak, passes = os.getppid(), time.monotonic(), 0, []
        signal.alarm(59)
        try:
            with tempfile.TemporaryDirectory(prefix='coxeter-ra-tex-',
                                             dir='/private/tmp') as temporary:
                for iteration in range(3):
                    child, output = None, bytearray()
                    try:
                        env = dict(os.environ)
                        for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS',
                                     'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
                            env[name] = '1'
                        command = [engine, '-interaction=nonstopmode',
                                   '-halt-on-error', '-file-line-error',
                                   '-output-directory='+temporary, str(source)]
                        child = subprocess.Popen(command, cwd=ROOT, env=env,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            start_new_session=True)
                        with selectors.DefaultSelector() as poll:
                            poll.register(child.stdout, selectors.EVENT_READ)
                            while poll.get_map() or child.poll() is None:
                                for key, _ in poll.select(.1):
                                    data = os.read(key.fd, 65536)
                                    if data:
                                        output.extend(data)
                                        assert len(output) < 524288
                                    else:
                                        poll.unregister(key.fileobj)
                                peak = max(peak, rss())
                                assert peak < LIMIT
                                assert time.monotonic()-started < 59
                                assert os.getppid() == owner
                                safe.resources(start=False)
                        assert child.wait() == 0, output[-5000:].decode(errors='replace')
                        log = (Path(temporary)/(source.stem+'.log')).read_text(
                            errors='replace')
                        passes.append({'pass': iteration+1,
                                       'stdout_bytes': len(output)})
                        if iteration >= 1 and 'Rerun to get' not in log:
                            break
                    finally:
                        if child is not None:
                            safe.stop_group(child, .3)
                assert digest(source) == original, 'Source changed during build'
                pdf = Path(temporary)/(source.stem+'.pdf')
                assert pdf.is_file()
                warnings = [line.strip() for line in log.splitlines() if re.search(
                    r'LaTeX Warning|Overfull|Underfull|Package .* Warning', line)]
                pages = re.search(r'Output written on .*?\((\d+) pages?', log, re.S)
                result = {'status': 'PASS_TEX_BUILD',
                    'source': str(source.relative_to(ROOT)),
                    'source_sha256': original,
                    'pdf': str(target.relative_to(ROOT)),
                    'pdf_sha256': digest(pdf), 'pdf_bytes': pdf.stat().st_size,
                    'pages': int(pages.group(1)) if pages else None,
                    'warnings': warnings, 'passes': passes,
                    'elapsed_seconds': round(time.monotonic()-started, 3),
                    'peak_rss_bytes': peak, 'initial_resources': before,
                    'temporary_build_files_removed': True,
                    'mathematical_verification_claimed': False,
                    'check_only': args.check_only,
                    'builder_sha256': digest(Path(__file__))}
                if not args.check_only:
                    shutil.copyfile(pdf, target)
                    ledger.write_text(json.dumps(result, indent=2)+'\n')
        finally:
            signal.alarm(0)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
