#!/usr/bin/env python3
"""Resume CURRENT-source distance9 receipts; no mathematical worker is imported.

  --ledger results/my_replay.json --plan [--width isolation=32]
  --ledger results/my_replay.json --max-jobs 2 [--width overlaps=8]
  --self-test

Plan/self-test execute no mathematical checks. Execution uses unchanged59s/1GiB
wrapper gates and lock, stops on failure/defer, and credits only parsed receipts.
Unsigned local receipts are not authenticated attestations. Widths need timing
pilots; no adaptive retry, timeout splitting, or deletion of old rows occurs.
"""
import argparse
import fcntl
from hashlib import sha256
import importlib.util
import json
import os
from pathlib import Path
import signal
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
COVERAGE = 'scripts/check_partitioned_replay_coverage_20261020.py'
DRIVER = 'scripts/verify_manuscript.py'
PINS = {COVERAGE: '456abc2bfe92262412106ebd6f94b07e48d01757cd1c76d8cef952c8bf3d3b0e',
        DRIVER: 'f6a0b39f981c41f3c7d75b184ee3d274e309616bcfa72b693a3b018bfc7394a1'}
DEFAULTS = dict(zip(('metadata', 'domains', 'overlap-domain', 'disjoint',
                     'inertia', 'overlaps', 'isolation', 'profiles'), (1, 1, 8, 8, 8, 4, 8, 8)))


def load(relative):
    spec = importlib.util.spec_from_file_location(Path(relative).stem, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path):
    value = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            value.update(block)
    return value.hexdigest()


def select(missing, widths):
    for phase in DEFAULTS:
        if missing.get(phase):
            lo, hi = missing[phase][0]
            return phase, lo, min(hi, lo + widths[phase])
    return None


def frozen(path):
    name = str(path.relative_to(ROOT))
    for manifest in (ROOT/'results').glob('*freeze.json'):
        with manifest.open() as stream:
            if any(name in line for line in stream):
                return True
    return False


def self_test(c):
    missing = {p: c.gaps([False]*n) for p, n in c.DOMAINS.items()}
    assert select(missing, DEFAULTS) == ('metadata', 0, 1)
    del missing['metadata'], missing['domains']
    assert select(missing, DEFAULTS) == ('overlap-domain', 0, 8)
    assert select({'overlaps': [[2, 5], [8, 9]]}, DEFAULTS) == ('overlaps', 2, 5)
    assert select({'profiles': [[7, 8]]}, DEFAULTS) == ('profiles', 7, 8)
    for phase, size in c.DOMAINS.items():
        assert select({phase: [[size-1, size]]}, DEFAULTS) == (phase, size-1, size)
    assert select({p: c.gaps([True]*n) for p, n in c.DOMAINS.items()}, DEFAULTS) is None
    assert sum(c.DOMAINS.values()) == 2196
    print(json.dumps(dict(status='PASS_INTERVAL_SELECTION_SELF_TEST_ONLY',
                         mathematical_execution_or_actual_coverage=False)))


def main():
    if not __debug__:
        raise RuntimeError('disabled assertions are not permitted')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ledger')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--max-jobs', type=int, default=0)
    parser.add_argument('--width', action='append', default=[], metavar='PHASE=N')
    args = parser.parse_args()
    for relative, expected in PINS.items():
        if digest(ROOT/relative) != expected:
            raise ValueError('coordinator input changed: '+relative)
    c, v = load(COVERAGE), load(DRIVER)  # Both definition-only, stdlib modules.
    c.need(tuple(c.DOMAINS) == tuple(DEFAULTS), 'phase contract changed')
    widths = dict(DEFAULTS)
    for item in args.width:
        phase, separator, raw = item.partition('=')
        c.need(separator and phase in widths and raw.isdecimal(), 'invalid width: '+item)
        widths[phase] = int(raw)
        c.need(0 < widths[phase] <= c.DOMAINS[phase], 'width outside phase domain')
    if args.self_test:
        c.need(not args.ledger and not args.plan and not args.max_jobs and not args.width,
               'self-test is separate from actual ledgers')
        self_test(c)
        return 0
    c.need(args.ledger and (args.plan or 1 <= args.max_jobs <= 2196),
           'supply --ledger and --plan or --max-jobs 1..2196')
    path = v.local(args.ledger, 'results', '.json')
    own, owner = digest(Path(__file__).resolve()), os.getppid()
    fixed = dict(c.FIXED, **PINS)

    def pins():
        c.need(os.getppid() == owner and owner > 1, 'coordinator parent lost')
        c.need(digest(Path(__file__).resolve()) == own, 'coordinator source changed')
        for relative, expected in fixed.items():
            c.need(c.digest(ROOT/relative) == expected, 'current pin mismatch: '+relative)

    def read():
        if not path.exists():
            return dict(jobs=[]), None
        c.need(path.stat().st_size <= c.MAX_LEDGER, 'ledger exceeds coverage-reader cap')
        with path.open('rb') as stream:
            raw = stream.read(c.MAX_LEDGER + 1)
        c.need(len(raw) <= c.MAX_LEDGER, 'ledger grew beyond cap')
        data = c.decode(raw)
        c.need(isinstance(data, dict) and isinstance(data.get('jobs'), list), 'expected jobs ledger')
        c.need(len(data['jobs']) < c.MAX_RECORDS, 'ledger record cap reached')
        stamp = sha256(raw).hexdigest()
        c.need(digest(path) == stamp, 'ledger changed during read')
        return data, stamp

    def summary(data):
        return c.analyze([('distance9', data)], {'distance9'})

    def report(data, state):
        answer = summary(data)
        answer.update(status=state, next_interval=select(answer['distance9']['missing'], widths),
                      widths=widths, ledger=str(path.relative_to(ROOT)), coordinator_sha256=own,
                      scope='Recorded current-source coverage only; not an additional math proof.')
        print(json.dumps(answer, sort_keys=True))
        return answer

    pins()
    if args.plan:
        data, stamp = read()
        report(data, 'PLAN_ONLY_NOT_EXECUTED')
        c.need((digest(path) if path.exists() else None) == stamp, 'ledger changed during plan')
        pins()
        return 0
    c.need(not frozen(path), 'refusing to modify a frozen ledger')
    lockpath = ROOT/'formal/.proof-batches/distance9_resume.lock'
    with lockpath.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        data, stamp = read()

        def save(state):
            nonlocal stamp
            c.need(not frozen(path), 'refusing to modify a frozen ledger')
            c.need((digest(path) if path.exists() else None) == stamp, 'concurrent ledger writer')
            data['distance9_resume'] = dict(status=state, coordinator_sha256=own,
                                           coverage_verifier_sha256=PINS[COVERAGE], widths=widths)
            raw = (json.dumps(data, ensure_ascii=False, allow_nan=False, indent=2)+'\n').encode()
            c.need(len(raw) <= c.MAX_LEDGER, 'ledger exceeds coverage-reader cap')
            temporary = None
            try:
                with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.distance9_resume_',
                                                 suffix='.tmp', delete=False) as stream:
                    temporary = Path(stream.name)
                    stream.write(raw); stream.flush(); os.fsync(stream.fileno())
                os.replace(temporary, path)
                stamp = digest(path)
            finally:
                if temporary is not None:
                    temporary.unlink(missing_ok=True)

        def interrupt(signum, frame):
            signal.setitimer(signal.ITIMER_REAL, 0)
            for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
                signal.signal(sig, signal.SIG_IGN)
            raise KeyboardInterrupt('parent lost' if signum == signal.SIGALRM else 'signal '+str(signum))

        def parent_poll(signum, frame):
            if os.getppid() != owner:
                interrupt(signum, frame)

        for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
            signal.signal(sig, interrupt)
        signal.signal(signal.SIGALRM, parent_poll)
        signal.setitimer(signal.ITIMER_REAL, .25, .25)
        state, interval = 'JOB_LIMIT_REACHED_INCOMPLETE', None
        try:
            for _ in range(args.max_jobs):
                pins()
                interval = select(summary(data)['distance9']['missing'], widths)
                if interval is None:
                    state = 'COMPLETE_RECORDED_CURRENT_COVERAGE'
                    break
                # Reserve space for a worst-case escaped run_job output before launch.
                c.need(len(json.dumps(data, ensure_ascii=False, indent=2).encode()) +
                       6*v.CAP + 65536 <= c.MAX_LEDGER,
                       'ledger needs a larger interval/new combined-ledger strategy; no job launched')
                c.need(len(data['jobs']) < c.MAX_RECORDS - 1, 'ledger record cap reached')
                phase, lo, hi = interval
                command = [sys.executable, '-B', str(ROOT/c.RUNNER), '--parent-lock', '--',
                           c.D9, '--verify', '--phase', phase, '--start', str(lo), '--stop', str(hi)]
                job = dict(id=f'distance9-{phase}-{lo}-{hi}-{len(data["jobs"])}', kind='script',
                           command=command, source=c.D9, coverage=f'{phase} [{lo},{hi}) only')
                result = v.run_job(job)
                data['jobs'].append(result)  # Preserve failures/stale history verbatim.
                state = result['status']
                if state == 'PASS_DECLARED_SCOPE':
                    c.need(c.parse_record(result, 'distance9') == ('distance9', interval), 'invalid new receipt')
                    state = 'JOB_LIMIT_REACHED_INCOMPLETE'
                save(state)
                print(json.dumps(dict(phase=phase, start=lo, stop=hi, status=result['status'],
                                      elapsed_seconds=result['elapsed_seconds'])), flush=True)
                pins()
                if result['status'] != 'PASS_DECLARED_SCOPE':
                    break
            if summary(data)['complete']:
                state = 'COMPLETE_RECORDED_CURRENT_COVERAGE'
        except (Exception, KeyboardInterrupt) as error:
            state = 'INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'COORDINATOR_ERROR'
            data.setdefault('coordinator_events', []).append(
                dict(status=state, interval=interval, error=str(error)[:400]))
        finally:
            signal.setitimer(signal.ITIMER_REAL, 0)
            for sig in (signal.SIGINT, signal.SIGTERM, signal.SIGHUP):
                signal.signal(sig, signal.SIG_IGN)
            save(state)
        report(data, state)
        return 0 if state in ('COMPLETE_RECORDED_CURRENT_COVERAGE', 'JOB_LIMIT_REACHED_INCOMPLETE') else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, TypeError, KeyError, OverflowError, RuntimeError) as error:
        print(json.dumps(dict(status='COORDINATOR_REFUSED', error=str(error)[:400])))
        sys.exit(2)
