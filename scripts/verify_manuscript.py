#!/usr/bin/env python3
"""Replay a manuscript's declared checks sequentially, with explicit coverage.

--list and --plan perform no mathematical computation. A successful run verifies
only the scopes in the catalogue; assumed geometry is never promoted to a Lean
theorem by this orchestration script.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import selectors
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CATALOGUE = ROOT/'manuscripts/catalogue.json'
CAP = 524288


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def local(value, parent=None, suffix=None):
    path = (ROOT/value).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError('Path leaves repository: '+str(value))
    if parent is not None and path.parent != ROOT/parent:
        raise ValueError('Unexpected parent: '+str(value))
    if suffix is not None and path.suffix != suffix:
        raise ValueError('Unexpected suffix: '+str(value))
    return path


def registry(value):
    return value if isinstance(value, dict) else {item['id']: item for item in value}


def plan_for(manuscript, checks, args):
    jobs = []
    if args.mode in ('scripts', 'core', 'full'):
        for identity in manuscript.get('verification', []):
            check = checks[identity]
            if args.mode != 'full' and check.get('tier', 'core') != 'core':
                continue
            argv = check['argv']
            if not isinstance(argv, list) or not argv or not all(
                    isinstance(item, str) for item in argv):
                raise ValueError('Check argv must be a nonempty string array')
            # Catalogue commands are script-relative; tolerate an explicit
            # Python prefix for convenient human-readable inventories.
            if Path(argv[0]).name in ('python', 'python3'):
                argv = argv[1:]
                if argv[:1] == ['-B']:
                    argv = argv[1:]
            source = local(argv[0], 'scripts', '.py')
            if not source.is_file():
                raise ValueError('Missing checker: '+str(source))
            command = [sys.executable, '-B',
                       str(ROOT/'scripts/run_bounded_research_check_20261018.py')]
            if not check.get('owns_lock', False):
                command.append('--parent-lock')
            command += ['--', *argv]
            jobs.append(dict(id=identity, kind='script', command=command,
                             coverage=check.get('coverage', ''),
                             source=str(source.relative_to(ROOT))))
    if args.mode in ('logic', 'core', 'full'):
        for module in manuscript.get('lean', {}).get('modules', []):
            if not re.fullmatch(r'CoxeterProofs(?:\.[A-Za-z0-9_]+)*', module):
                raise ValueError('Invalid Lean module: '+module)
            module_source = 'formal/'+module.replace('.', '/')+'.lean'
            audit = False
            path = local(module_source)
            if path.is_file():
                with path.open() as stream:
                    audit = any(line.lstrip().startswith('#print axioms ') for line in stream)
            command = [sys.executable, '-B',
                str(ROOT/'scripts/run_lean_batches_signal_safe_20261018.py'), module]
            if audit:
                command.append('--audit-axioms')
            jobs.append(dict(id=module, kind='conditional_lean', command=command,
                printed_axiom_audit_requested=audit,
                coverage=manuscript.get('lean', {}).get('coverage', ''),
                source=module_source))
    if args.pdf:
        override = manuscript.get('pdf_builder_argv')
        if override is not None:
            if not isinstance(override, list) or not override or not all(
                    isinstance(item, str) for item in override):
                raise ValueError('PDF builder argv must be a nonempty string array')
            builder = local(override[0], 'scripts', '.py')
            if not builder.is_file() or '--check-only' not in override[1:]:
                raise ValueError('Custom PDF builder must exist and request check-only')
            command = [sys.executable, '-B', str(builder), *override[1:]]
        else:
            command = [sys.executable, '-B', str(ROOT/
                       'scripts/build_illustrated_ra_manuscript_20261018.py'),
                       manuscript['tex'], '--check-only']
        jobs.append(dict(id=manuscript['id']+'-pdf', kind='typesetting',
                         command=command,
                         coverage='Fresh PDF compilation; preserves published PDF and ledger.',
                         source=manuscript['tex']))
    return jobs


def run_job(job):
    output, child = bytearray(), None
    source = local(job['source'])
    original = digest(source)
    execution_runner = Path(job['command'][2])
    runner_original = digest(execution_runner)
    started = time.monotonic()
    try:
        child = subprocess.Popen(job['command'], cwd=ROOT, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT, start_new_session=True)
        with selectors.DefaultSelector() as poll:
            poll.register(child.stdout, selectors.EVENT_READ)
            while poll.get_map() or child.poll() is None:
                for key, _ in poll.select(.2):
                    data = os.read(key.fd, 65536)
                    if data:
                        output.extend(data)
                        if len(output) > CAP:
                            raise RuntimeError('Batch output exceeded 512KiB')
                    else:
                        poll.unregister(key.fileobj)
        code = child.wait()
        if digest(source) != original:
            raise RuntimeError('Job source changed during verification: '+job['source'])
        if digest(execution_runner) != runner_original:
            raise RuntimeError('Execution runner changed during verification')
        decoded = output.decode(errors='replace')
        state = 'PASS_DECLARED_SCOPE' if code == 0 else (
            'DEFERRED_RESOURCE_GATE' if 'Resource gate:' in decoded else 'FAILED')
        record = dict(id=job['id'], kind=job['kind'], status=state,
                      coverage=job['coverage'], exit_code=code,
                      elapsed_seconds=round(time.monotonic()-started, 3),
                      stdout_sha256=hashlib.sha256(output).hexdigest(),
                      stdout_bytes=len(output),
                      source_sha256=original, execution_runner_sha256=runner_original)
        if job['kind'] == 'conditional_lean':
            record['printed_axiom_audit_requested'] = job['printed_axiom_audit_requested']
        # Retain compilation/cache provenance. A printed axiom audit occurs
        # only when the target contains actual #print axioms commands.
        if code == 0:
            record['output'] = decoded
        else:
            record['failure_tail'] = decoded[-7000:]
        return record
    finally:
        if child is not None:
            if child.poll() is None:
                child.terminate()
                try:
                    child.wait(timeout=8)
                except subprocess.TimeoutExpired:
                    child.kill()
                    child.wait()
            if child.stdout:
                child.stdout.close()


def interrupt(signum, frame):
    raise KeyboardInterrupt('Manuscript verification interrupted')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('manuscript', nargs='?')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--mode', choices=['scripts', 'logic', 'core', 'full'],
                        default='core')
    parser.add_argument('--pdf', action='store_true')
    parser.add_argument('--ledger', help='New JSON file directly under results/')
    args = parser.parse_args()
    catalogue_hash = digest(CATALOGUE)
    catalogue = json.loads(CATALOGUE.read_text())
    manuscripts = registry(catalogue['manuscripts'])
    checks = registry(catalogue.get('checks', []))
    if args.list:
        print(json.dumps(dict(status='INVENTORY_ONLY', manuscripts=[
            dict(id=item['id'], tex=item['tex'], role=item.get('role'),
                 bundle=item.get('bundle')) for item in manuscripts.values()]), indent=2))
        return 0
    if not args.manuscript or args.manuscript not in manuscripts:
        parser.error('Choose a catalogue manuscript ID; use --list')
    manuscript = manuscripts[args.manuscript]
    jobs = plan_for(manuscript, checks, args)
    gaps = manuscript.get('lean', {}).get('gaps', [])
    if args.plan:
        print(json.dumps(dict(status='PLAN_ONLY_NOT_VERIFIED',
            manuscript=args.manuscript, jobs=jobs, remaining_formal_gaps=gaps), indent=2))
        return 0
    if not jobs:
        raise ValueError('No checks declared for this selection; no pass is claimed')
    source = local(manuscript['tex'], suffix='.tex')
    if not source.is_relative_to(ROOT/'manuscripts'):
        raise ValueError('External reference is an inventory entry, not a local manuscript')
    if not source.is_file():
        raise ValueError('Manuscript source is not yet available')
    ledger = local(args.ledger, 'results', '.json') if args.ledger else None
    if ledger and ledger.exists():
        raise ValueError('Refusing to replace an existing verification ledger')
    source_hash = digest(source)
    assumptions = ROOT/'docs/manuscript_assumptions.md'
    if not assumptions.is_file():
        raise ValueError('Assumption register is not yet available; no pass is claimed')
    assumptions_hash = digest(assumptions)
    record = dict(status='INCOMPLETE', manuscript=args.manuscript, mode=args.mode,
        created_utc=datetime.now(timezone.utc).isoformat(),
        catalogue_sha256=catalogue_hash, manuscript_sha256=source_hash,
        runner_sha256=digest(Path(__file__)), jobs=[], planned_jobs=len(jobs),
        remaining_formal_gaps=gaps,
        assumption_document='docs/manuscript_assumptions.md',
        assumption_document_sha256=assumptions_hash,
        end_to_end_geometric_lean_proof_claimed=False)
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(signum, interrupt)
    try:
        for job in jobs:
            if (digest(CATALOGUE) != catalogue_hash or digest(source) != source_hash
                    or digest(assumptions) != assumptions_hash
                    or digest(Path(__file__)) != record['runner_sha256']):
                raise RuntimeError('Catalogue, manuscript, assumptions or runner changed')
            print('Checking '+job['id'], file=sys.stderr, flush=True)
            result = run_job(job)
            record['jobs'].append(result)
            if result['status'] != 'PASS_DECLARED_SCOPE':
                record['status'] = result['status']
                break
        else:
            if (digest(CATALOGUE) != catalogue_hash or digest(source) != source_hash
                    or digest(assumptions) != assumptions_hash
                    or digest(Path(__file__)) != record['runner_sha256']):
                raise RuntimeError('Catalogue, manuscript, assumptions or runner changed')
            record['status'] = 'PASS_DECLARED_CHECKS_CONDITIONAL_ON_LISTED_ASSUMPTIONS'
    except BaseException as error:
        record['status'] = 'INTERRUPTED' if isinstance(error, KeyboardInterrupt) else 'FAILED'
        record['error'] = str(error)
    if ledger:
        with ledger.open('x') as stream:
            json.dump(record, stream, indent=2)
            stream.write('\n')
    print(json.dumps(record, indent=2))
    return 0 if record['status'].startswith('PASS_') else 2


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, OSError, KeyError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(2)
