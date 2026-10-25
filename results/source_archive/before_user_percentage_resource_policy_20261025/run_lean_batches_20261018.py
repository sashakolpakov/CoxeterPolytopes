#!/usr/bin/env python3
"""One module per batch, bounded Lean memory, resource gates, automatic cleanup.

Legacy objects may be reused, but are never labelled newly capped verification.
This guard is not a mathematical proof and never starts a parallel Lake build.
"""
import argparse
from collections import deque
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import selectors
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
FORMAL = ROOT / 'formal'
STATE = FORMAL / '.proof-batches'
GIB = 1024**3
MEMORY_MB = 2048
OUTPUT_LIMIT = 524288
TIMEOUT_SECONDS = 900
_COMPILER_HASHES = {}


def available_memory():
    output = subprocess.check_output(['vm_stat'], text=True, timeout=5)
    page = int(re.search(r'page size of (\d+) bytes', output).group(1))
    counts = {k: int(v) for k, v in re.findall(r'^([^:\n]+):\s+(\d+)', output, re.M)}
    return page * sum(counts.get(k, 0) for k in
                      ('Pages free', 'Pages inactive', 'Pages speculative'))


def resources(start=True, mark_cleanup=True):
    disk = shutil.disk_usage(ROOT).free
    memory = available_memory()
    if disk < 20*GIB and mark_cleanup:
        STATE.mkdir(exist_ok=True)
        marker = STATE / 'cleanup_needed.json'
        if not marker.exists():
            marker.write_text(json.dumps(dict(status='CLEANUP_AGENT_REQUEST',
                free_disk_gib=round(disk/GIB, 2), threshold_gib=20,
                action='Coordinator must start a cleanup agent before further jobs'))+'\n')
            print(f'CLEANUP_AGENT_REQUEST: {marker}', file=sys.stderr, flush=True)
    disk_floor, memory_floor = (40*GIB, 6*GIB) if start else (8*GIB, 3*GIB)
    if disk < disk_floor or memory < memory_floor:
        raise RuntimeError(f'Resource gate: disk={disk/GIB:.1f}GiB '
                           f'memory={memory/GIB:.1f}GiB; '
                           f'require {disk_floor/GIB:.0f}/{memory_floor/GIB:.0f}GiB')
    return {'free_disk_gib': round(disk/GIB, 2),
            'available_memory_gib': round(memory/GIB, 2)}


def source(module):
    if not re.fullmatch(r'CoxeterProofs(?:\.[A-Za-z0-9_]+)*', module):
        raise ValueError(f'Invalid local module: {module}')
    path = FORMAL / (module.replace('.', '/')+'.lean')
    if not path.is_file():
        raise ValueError(f'Missing source: {module}')
    return path


def local_imports(module):
    # Generated certificates can be large: inspect their import lines as a stream.
    with source(module).open() as stream:
        for line in stream:
            match = re.fullmatch(r'\s*import\s+(CoxeterProofs(?:\.[A-Za-z0-9_]+)*)\s*', line)
            if match:
                yield match.group(1)


def dependency_order(target):
    seen, active, ordered = set(), set(), []
    def visit(module):
        if module in seen:
            return
        if module in active:
            raise RuntimeError(f'Import cycle: {module}')
        active.add(module)
        for name in local_imports(module):
            visit(name)
        active.remove(module)
        seen.add(module)
        ordered.append(module)
    visit(target)
    return ordered


def object_path(module):
    return FORMAL / '.lake/build/lib/lean' / (module.replace('.', '/')+'.olean')


def ledger_path(module, audit=False):
    return STATE / (module.replace('.', '_')+('_axiom_audit' if audit else '')+'.json')


def file_identity(path):
    stat = path.stat()
    return (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)


def file_hash(path, cache=None):
    before = file_identity(path)
    if cache is not None and cache.get(str(path), (None,))[0] == before:
        return cache[str(path)][1]
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            digest.update(chunk)
    if file_identity(path) != before:
        raise RuntimeError(f'Input changed while hashing: {path}')
    result = digest.hexdigest()
    if cache is not None:
        cache[str(path)] = (before, result)
    return result


def compiler_identity():
    configs = ('lean-toolchain', 'lakefile.toml', 'lakefile.lean', 'lake-manifest.json')
    result = {name: file_hash(FORMAL/name) for name in configs if (FORMAL/name).exists()}
    toolchain = os.environ.get('ELAN_TOOLCHAIN') or (FORMAL/'lean-toolchain').read_text().strip()
    elan = Path(os.environ.get('ELAN_HOME', str(Path.home()/'.elan')))
    install = Path(toolchain) if Path(toolchain).is_absolute() else elan/'toolchains'/toolchain.replace('/', '--').replace(':', '---')
    compiler = install/'bin/lean'
    if not compiler.is_file():
        raise RuntimeError(f'Cannot fingerprint pinned Lean compiler: {compiler}')
    binaries = [compiler, *sorted((install/'lib/lean').glob('libleanshared*'))]
    result['compiler_files'] = {str(path): file_hash(path, _COMPILER_HASHES)
                                for path in binaries if path.is_file()}
    result['environment'] = {name: os.environ.get(name, '') for name in
                             ('ELAN_TOOLCHAIN', 'LEAN_PATH', 'LEAN_SRC_PATH')}
    return result


def input_snapshot(module):
    return dict(source_sha256=file_hash(source(module)),
                dependency_olean_sha256={name: file_hash(object_path(name))
                    for name in dependency_order(module) if name != module},
                compiler=compiler_identity(), memory_cap_mib=MEMORY_MB, threads=1)


def clean_unused():
    removed = 0
    targets = ((FORMAL/'.lake/build', ('*.ilean', '*.ilean.hash', '*.c', '*.c.o',
                '*.bc', '*.ir', '*.setup.json', '*.log')),
               (ROOT/'manuscripts/build', ('*.aux', '*.fls', '*.fdb_latexmk',
                '*.log', '*.out', '*.toc', '*.synctex.gz', '*.nav', '*.snm', '*.vrb')))
    for directory, patterns in targets:
        for pattern in patterns:
            for path in directory.rglob(pattern):
                if path.is_file():
                    removed += path.stat().st_size
                    path.unlink()
    work = STATE / 'work'
    if work.exists():
        removed += sum(p.stat().st_size for p in work.rglob('*') if p.is_file())
        shutil.rmtree(work)
    return removed


def cached(module):
    """Return truthful provenance, or None when the object must be rebuilt."""
    path, olean = source(module), object_path(module)
    if not olean.exists() or olean.stat().st_mtime_ns < path.stat().st_mtime_ns:
        return None
    for name in local_imports(module):
        dep = object_path(name)
        if not dep.exists() or dep.stat().st_mtime_ns > olean.stat().st_mtime_ns:
            return None
    ledger = ledger_path(module)
    if not ledger.exists():
        return 'LEGACY_CACHED_UNCAPPED'
    try:
        with ledger.open() as stream:
            record = json.load(stream)
        if (record.get('status') == 'KERNEL_COMPILED'
                and record.get('inputs') == input_snapshot(module)
                and record.get('olean_sha256') == file_hash(olean)):
            return 'KERNEL_COMPILED'
    except (OSError, ValueError, RuntimeError):
        pass
    return None


def interrupt(signum, frame):
    raise KeyboardInterrupt(f'signal {signum}')


def stop_group(child, grace_seconds=5):
    try:
        os.killpg(child.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    deadline = time.monotonic()+grace_seconds
    while time.monotonic() < deadline:
        child.poll()
        try:
            os.killpg(child.pid, 0)
        except ProcessLookupError:
            break
        time.sleep(0.1)
    else:
        try:
            os.killpg(child.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
    child.wait()
    if child.stdout is not None:
        child.stdout.close()


def guarded_run(command, cwd=FORMAL, timeout_seconds=None, termination_grace=5):
    """Monitor until exit even if stdout closes; retain at most 512KiB."""
    before, parent = resources(), os.getppid()
    timeout = TIMEOUT_SECONDS if timeout_seconds is None else timeout_seconds
    chunks, digest = deque(), hashlib.sha256()
    started, size, retained, child = time.monotonic(), 0, 0, None
    try:
        child = subprocess.Popen(command, cwd=cwd, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        with selectors.DefaultSelector() as poll:
            poll.register(child.stdout, selectors.EVENT_READ)
            while poll.get_map() or child.poll() is None:
                if not poll.get_map():
                    time.sleep(0.25)
                for key, _ in poll.select(timeout=2):
                    data = os.read(key.fd, 65536)
                    if data:
                        size += len(data)
                        digest.update(data)
                        chunks.append(data)
                        retained += len(data)
                        while retained > OUTPUT_LIMIT:
                            excess = retained-OUTPUT_LIMIT
                            first = chunks.popleft()
                            if len(first) > excess:
                                chunks.appendleft(first[excess:])
                                retained -= excess
                            else:
                                retained -= len(first)
                    else:
                        poll.unregister(key.fileobj)
                resources(start=False)
                if os.getppid() != parent:
                    raise RuntimeError('Batch owner exited; stopping its child')
                if timeout and time.monotonic()-started > timeout:
                    raise RuntimeError('Batch exceeded 15 minutes; split this module further')
        if child.wait():
            raise RuntimeError(f'Child exited {child.returncode}; see bounded last_failure.txt')
        return dict(output=b''.join(chunks), stdout_sha256=digest.hexdigest(),
                    output_bytes=size, elapsed_seconds=round(time.monotonic()-started, 2),
                    resources=before)
    except BaseException:
        if STATE.exists():
            (STATE/'last_failure.txt').write_bytes(b''.join(chunks)[-OUTPUT_LIMIT:])
        raise
    finally:
        if child is not None:
            stop_group(child, termination_grace)


def batch(module, audit=False):
    path = source(module)
    inputs = input_snapshot(module)  # Capture before launch, then verify before publish.
    work = STATE/'work'
    work.mkdir(parents=True, exist_ok=True)
    temporary = work/'batch.olean'
    command = ['lake', 'env', 'lean', '-M', str(MEMORY_MB), '-j', '1',
               str(path.relative_to(FORMAL))]
    if not audit:
        command += ['-o', str(temporary)]
    try:
        result = guarded_run(command)
        if input_snapshot(module) != inputs:
            raise RuntimeError('Source, dependencies or toolchain changed; object not published')
        record = dict(status='KERNEL_AXIOM_AUDIT' if audit else 'KERNEL_COMPILED',
                      module=module, inputs=inputs,
                      **{key: value for key, value in result.items() if key != 'output'})
        if audit:
            if result['output_bytes'] > len(result['output']):
                raise RuntimeError('Axiom audit exceeded bounded output; audit incomplete')
            sys.stdout.write(result['output'].decode(errors='replace'))
        else:
            if not temporary.is_file():
                raise RuntimeError('Lean did not produce the requested object')
            record['olean_sha256'] = file_hash(temporary)
            output = object_path(module)
            output.parent.mkdir(parents=True, exist_ok=True)
            os.replace(temporary, output)
        ledger = ledger_path(module, audit)
        pending = ledger.with_suffix('.json.tmp')
        pending.write_text(json.dumps(record)+'\n')
        os.replace(pending, ledger)
        (STATE/'last_failure.txt').unlink(missing_ok=True)
        print(json.dumps(record), flush=True)
    finally:
        clean_unused()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('module', nargs='?', default='CoxeterProofs')
    parser.add_argument('--check-only', action='store_true')
    parser.add_argument('--audit-axioms', action='store_true')
    args = parser.parse_args()
    order = dependency_order(args.module)
    if args.check_only:
        states = {module: cached(module) or 'NEEDS_BATCH' for module in order}
        print(json.dumps(dict(modules=len(order), parallel_jobs=1, memory_cap_mib=MEMORY_MB,
            threads=1, resources=resources(mark_cleanup=False), cache_provenance=states)))
        return
    if (ROOT/'repro/RESEARCH_PAUSED').exists():
        raise RuntimeError('Research is paused; no proof batch was started')
    STATE.mkdir(exist_ok=True)
    with (STATE/'exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
            signal.signal(sig, interrupt)
        legacy, compiled = [], []
        try:
            clean_unused()  # Only after obtaining the exclusive lock.
            resources(start=False)
            for module in order:
                state = cached(module)
                if state == 'LEGACY_CACHED_UNCAPPED':
                    legacy.append(module)
                elif state is None:
                    batch(module)
                    compiled.append(module)
            if args.audit_axioms:
                batch(args.module, audit=True)
            print(json.dumps(dict(status='PASS_SEQUENTIAL_MODULE_CHECK',
                newly_capped_modules=compiled, legacy_cached_uncapped_modules=legacy,
                all_dependencies_have_capped_ledger=not legacy)), flush=True)
        finally:
            clean_unused()  # Includes an entirely cached run and every failure path.


if __name__ == '__main__':
    try:
        main()
    except (RuntimeError, ValueError, OSError, KeyboardInterrupt) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
