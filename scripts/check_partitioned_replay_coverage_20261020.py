#!/usr/bin/env python3
"""Check recorded D13/name and distance9/interval replay coverage, not math.

Examples (repeat either ledger option to combine independent receipts):
  --verify --d13-ledger results/finite57_d13_batched_replay_runtime_20261020.json
  --verify --distance9-ledger results/my_distance9_replay.json
  --self-test

Only valid current-source successful records contribute. Failed, deferred,
stale or malformed records contribute nothing; duplicate successes are
allowed. Top-level status/plans and nested historical attempts never prove
coverage. A complete current conjunction may coexist with historical failed
records. Empty or incomplete ledgers exit1, malformed input/current-pin
failures exit2. Self-test success is NOT actual replay coverage.

No worker is imported or executed. Local unsigned ledgers are trusted as
execution records, not authenticated attestations: this cannot detect a
fabricated record, prove the correctness of a watchdog, rerun input-domain
computations, or establish geometric hypotheses.
"""
import argparse
from collections import Counter
from hashlib import sha256
import json
import math
from pathlib import Path
import signal
import sys

ROOT = Path(__file__).resolve().parents[1]
D13 = 'scripts/check_d13_transfer_batches_20261020.py'
D9 = 'scripts/check_selective_distance9_batches_20261020.py'
PARSER = 'scripts/check_d13_json_projection_contract_20261020.py'
RUNNER = 'scripts/run_bounded_research_check_20261018.py'
FIXED = {
    D13: 'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    D9: 'f5f0b47084135f8bb1886a37e0478d1f6b8c9c4327499f986a8ec9e173ecd843',
    PARSER: '2f6b63daa564b1689215d105d215bd2fb707f5873b95604cf19b8be90509de94',
    RUNNER: '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/verify_manuscript.py': 'f6a0b39f981c41f3c7d75b184ee3d274e309616bcfa72b693a3b018bfc7394a1',
    'scripts/fresh_d13_three_rank15_context_transfer_20261018.py':
        'a95ea679854ba907ead8640964a258c3baeb3e998052a3e2e3ff20f230cce486',
    'scripts/fresh_selective_distance9_dimension108_20261017.py':
        'f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a',
    'results/fresh_d13_three_rank15_context_transfer_20261018.json':
        '2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
    'results/fresh_selective_distance9_dimension108_20261017.json':
        'e7c6d1a1e56d8788e7d7a7d02bd5dd71837e3a3aef4ae768c1efc0a187bcacd9',
}
CERT13 = FIXED['results/fresh_d13_three_rank15_context_transfer_20261018.json']
CERT9 = FIXED['results/fresh_selective_distance9_dimension108_20261017.json']
BATCHES = (
    'pins-0', 'pins-1', 'pins-2', 'pins-3', 'metadata', 'support', 'rank15',
    'one-sided', 'quiet', 'literal', 'later-0', 'later-1', 'later-2', 'later-3',
    'newer', 'variants', 'd12', 'universal', 'finite',
)
DOMAINS = {'metadata': 1, 'domains': 1, 'overlap-domain': 384,
           'disjoint': 384, 'inertia': 85, 'overlaps': 384,
           'isolation': 949, 'profiles': 8}
MAX_FILES, MAX_RECORDS = 4096, 8192
MAX_LEDGER, MAX_TOTAL, MAX_OUTPUT = 8*1024**2, 64*1024**2, 512*1024
RSS_CAP, TIME_CAP = 1024**3, 59


def need(condition, message):
    if not condition:
        raise ValueError(message)


def pairs(items):
    result = {}
    for key, value in items:
        need(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def decode(raw):
    def invalid(_):
        raise ValueError('nonfinite JSON constant')
    def finite_float(token):
        value = float(token)
        need(math.isfinite(value), 'nonfinite JSON number')
        return value
    return json.loads(raw, object_pairs_hook=pairs, parse_constant=invalid, parse_float=finite_float)


def number(value):
    return type(value) in (int, float) and math.isfinite(value)


def digest(path, cap=40*1024**2):
    need(path.is_file() and path.stat().st_size <= cap, 'file missing/too large: '+str(path))
    h = sha256()
    size = 0
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            size += len(block)
            need(size <= cap, 'file grew beyond cap')
            h.update(block)
    return h.hexdigest()


def wrapper(record):
    """Accept direct wrapper JSON or actual verify_manuscript.run_job JSON."""
    need(isinstance(record, dict), 'record is not an object')
    if record.get('status') == 'PASS_BOUNDED_RESEARCH_CHECK':
        meta, text = record, record.get('stdout')
    else:
        need(record.get('status') == 'PASS_DECLARED_SCOPE' and
             type(record.get('exit_code')) is int and record['exit_code'] == 0 and
             record.get('kind') == 'script', 'not a successful script record')
        output = record.get('output')
        need(isinstance(output, str), 'missing complete outer output')
        raw = output.encode('utf-8')
        need(len(raw) <= MAX_OUTPUT and record.get('stdout_bytes') == len(raw) and
             record.get('stdout_sha256') == sha256(raw).hexdigest(), 'outer output byte/hash mismatch')
        lines = output.splitlines(keepends=True)
        need(len(lines) >= 2 and lines[-1].strip(), 'missing wrapper metadata line')
        meta, text = decode(lines[-1]), ''.join(lines[:-1])
        need(isinstance(meta, dict) and 'stdout' not in meta, 'not printed wrapper metadata')
        need(record.get('source_sha256') == meta.get('checker_sha256') and
             record.get('execution_runner_sha256') == FIXED[RUNNER], 'outer source/runner mismatch')
        need(number(record.get('elapsed_seconds')) and record['elapsed_seconds'] >= 0,
             'invalid outer elapsed time')
    need(isinstance(text, str), 'missing worker output')
    need(isinstance(meta, dict) and meta.get('status') == 'PASS_BOUNDED_RESEARCH_CHECK',
         'no actual wrapper PASS')
    need(meta.get('runner_sha256') == FIXED[RUNNER] and meta.get('lock_owner') == 'wrapper',
         'wrapper pin/lock mismatch')
    need(meta.get('wall_limit_seconds') == TIME_CAP and
         meta.get('rss_limit_bytes') == RSS_CAP and
         type(meta.get('threads')) is int and meta['threads'] == 1,
         'resource cap/thread declaration mismatch')
    need(number(meta.get('elapsed_seconds')) and 0 <= meta['elapsed_seconds'] <= TIME_CAP,
         'elapsed time outside cap')
    if 'output' in record:
        need(record['elapsed_seconds'] + 0.01 >= meta['elapsed_seconds'],
             'outer/wrapper elapsed time inconsistent')
    need(type(meta.get('peak_process_tree_rss_bytes')) is int and
         0 < meta['peak_process_tree_rss_bytes'] < RSS_CAP, 'RSS outside cap')
    resources = meta.get('initial_resources')
    need(isinstance(resources, dict) and number(resources.get('free_disk_gib')) and
         resources['free_disk_gib'] >= 40 and number(resources.get('available_memory_gib')) and
         resources['available_memory_gib'] >= 6, 'launch resource gate not recorded')
    raw = text.encode('utf-8')
    need(len(raw) <= MAX_OUTPUT and meta.get('stdout_bytes') == len(raw), 'worker output bytes mismatch')
    return meta, text


def parse_record(record, family):
    meta, text = wrapper(record)
    argv = meta.get('command')
    need(isinstance(argv, list) and all(isinstance(s, str) for s in argv), 'invalid argv')
    if family == 'd13' and argv == [PARSER, '--verify']:
        need(meta.get('checker_sha256') == FIXED[PARSER], 'stale parser source')
        expected = ('PASS_D13_JSON_PROJECTION_CONTRACT valid_projections 19 required_rejections 15 '
                    'actual_certificate_projections 1 scientific_replay False temporary_fixtures_removed True\n')
        need(text == expected, 'parser marker mismatch')
        return 'parser', None
    if family == 'd13':
        need(len(argv) == 4 and argv[:3] == [D13, '--verify', '--batch'] and
             argv[3] in BATCHES, 'D13 argv/domain mismatch')
        name = argv[3]
        need(meta.get('checker_sha256') == FIXED[D13], 'stale D13 source')
        expected = (f'PASS_D13_FROZEN_TRANSFER_REPLAY_BATCH {name} outer_pins 484 '
                    f'saved_certificate {CERT13} complete_replay_requires_all 19\n')
        need(text == expected, 'D13 source/certificate/batch marker mismatch')
        return 'd13', name
    need(family == 'distance9' and argv and argv[0] == D9, 'distance9 worker mismatch')
    need(meta.get('checker_sha256') == FIXED[D9], 'stale distance9 source')
    flags, index = {}, 1
    while index < len(argv):
        flag = argv[index]
        need(flag in ('--verify', '--phase', '--start', '--stop') and flag not in flags,
             'unknown/duplicate distance9 argument')
        if flag == '--verify':
            flags[flag], index = True, index + 1
        else:
            need(index + 1 < len(argv), 'missing argument value')
            flags[flag], index = argv[index + 1], index + 2
    need(flags.get('--verify') is True and flags.get('--phase') in DOMAINS, 'missing phase/verify')
    phase = flags['--phase']
    lo = int(flags.get('--start', '0'))
    hi = int(flags.get('--stop', str(lo + 1)))
    need(str(lo) == flags.get('--start', '0') and str(hi) == flags.get('--stop', str(lo + 1)) and
         0 <= lo < hi <= DOMAINS[phase], 'invalid phase interval')
    proof = decode(text)
    need(isinstance(proof, dict) and set(proof) == {
        'status', 'phase', 'start', 'stop', 'detail', 'complete_phase_domains', 'source', 'certificate'},
        'distance9 marker schema mismatch')
    need(proof['status'] == 'PASS_DISTANCE9_REPLAY_PARTITION' and proof['phase'] == phase and
         type(proof['start']) is int and type(proof['stop']) is int and
         (proof['start'], proof['stop']) == (lo, hi), 'argv/printed interval mismatch')
    need(proof['complete_phase_domains'] == DOMAINS and
         all(type(n) is int for n in proof['complete_phase_domains'].values()) and
         isinstance(proof['detail'], dict), 'complete domain map/detail mismatch')
    need(proof['source'] == FIXED[D9] and proof['certificate'] == CERT9, 'distance9 source/CERT marker mismatch')
    return 'distance9', (phase, lo, hi)


def records(document):
    need(isinstance(document, dict), 'ledger root is not an object')
    if 'jobs' in document:
        need(isinstance(document['jobs'], list), 'jobs is not a list')
        return document['jobs']
    if document.get('status') in ('PASS_BOUNDED_RESEARCH_CHECK', 'PASS_DECLARED_SCOPE'):
        return [document]
    # No implicit traversal of previous_attempt, operational failures, plans,
    # summary coverage claims, or other historical nested records.
    return []


def gaps(bits):
    result, index = [], 0
    while index < len(bits):
        if bits[index]:
            index += 1
            continue
        end = index + 1
        while end < len(bits) and not bits[end]:
            end += 1
        result.append([index, end])
        index = end
    return result


def analyze(inputs, requested):
    covered13, covered9 = set(), {p: [False]*n for p, n in DOMAINS.items()}
    accepted, ignored, parser_passes, seen = Counter(), Counter(), 0, 0
    for family, document in inputs:
        for record in records(document):
            seen += 1
            need(seen <= MAX_RECORDS, 'too many records')
            try:
                kind, value = parse_record(record, family)
            except (ValueError, TypeError, KeyError, OverflowError, UnicodeError) as error:
                reason = str(error)[:160]
                if reason not in ignored and len(ignored) >= 40:
                    reason = 'other invalid records (reason list capped)'
                ignored[reason] += 1
                continue
            accepted[family] += 1
            if kind == 'parser':
                parser_passes += 1
            elif kind == 'd13':
                covered13.add(value)
            else:
                phase, lo, hi = value
                covered9[phase][lo:hi] = [True]*(hi-lo)
    answer = dict(records_examined=seen, accepted_records=dict(accepted),
                  ignored_records=sum(ignored.values()), ignored_reasons=dict(ignored),
                  current_parser_contract_pass_records=parser_passes)
    complete = True
    if 'd13' in requested:
        missing = [name for name in BATCHES if name not in covered13]
        answer['d13'] = dict(required_batches=len(BATCHES), covered_batches=len(covered13), missing=missing)
        complete &= not missing
    if 'distance9' in requested:
        missing = {p: gaps(bits) for p, bits in covered9.items() if not all(bits)}
        answer['distance9'] = dict(required_domains=DOMAINS,
                                  covered_units={p: sum(v) for p, v in covered9.items()}, missing=missing)
        complete &= not missing
    answer['complete'] = bool(complete)
    return answer


def oracle_self_test():
    """Synthetic receipt tests only; fabricated fixtures never count as runs."""
    def receipt(command, text):
        return dict(status='PASS_BOUNDED_RESEARCH_CHECK', command=command,
                    elapsed_seconds=0.1, peak_process_tree_rss_bytes=1000000,
                    initial_resources=dict(free_disk_gib=50, available_memory_gib=7),
                    wall_limit_seconds=59, rss_limit_bytes=1024**3, threads=1,
                    stdout_bytes=len(text.encode()), stdout=text,
                    runner_sha256=FIXED[RUNNER], checker_sha256=FIXED[command[0]], lock_owner='wrapper')

    def outer(value):
        text = value['stdout'] + json.dumps({k: v for k, v in value.items() if k != 'stdout'}) + '\n'
        return dict(id='synthetic', kind='script', status='PASS_DECLARED_SCOPE', exit_code=0,
                    elapsed_seconds=0.2, stdout_sha256=sha256(text.encode()).hexdigest(),
                    stdout_bytes=len(text.encode()), source_sha256=value['checker_sha256'],
                    execution_runner_sha256=FIXED[RUNNER], output=text)

    def batch(name):
        text = (f'PASS_D13_FROZEN_TRANSFER_REPLAY_BATCH {name} outer_pins 484 '
                f'saved_certificate {CERT13} complete_replay_requires_all 19\n')
        return receipt([D13, '--verify', '--batch', name], text)

    def interval(phase, lo, hi):
        proof = dict(status='PASS_DISTANCE9_REPLAY_PARTITION', phase=phase, start=lo, stop=hi,
                     detail={}, complete_phase_domains=DOMAINS, source=FIXED[D9], certificate=CERT9)
        return receipt([D9, '--verify', '--phase', phase, '--start', str(lo), '--stop', str(hi)],
                       json.dumps(proof, sort_keys=True)+'\n')

    def check(items, requested, expected):
        need(analyze(items, requested)['complete'] is expected, 'coverage oracle failed')

    complete13 = [outer(batch(n)) if i % 2 else batch(n) for i, n in enumerate(BATCHES)]
    complete9 = [outer(interval(p, 0, n)) for p, n in DOMAINS.items()]
    check([('d13', {'status': 'PASS', 'jobs': []})], {'d13'}, False)
    check([('d13', {'jobs': [], 'previous_attempt': {'jobs': complete13}})], {'d13'}, False)
    check([('d13', {'jobs': complete13})], {'d13'}, True)
    check([('d13', {'jobs': complete13[:-1]})], {'d13'}, False)
    check([('d13', {'jobs': complete13 + [complete13[0], {'status': 'FAILED'}]})], {'d13'}, True)
    check([('d13', {'status': 'FAILED', 'jobs': complete13})], {'d13'}, True)
    check([('distance9', {'jobs': complete9})], {'distance9'}, True)
    check([('distance9', {'jobs': complete9[:-1]})], {'distance9'}, False)
    split = [interval('isolation', 0, 400), interval('isolation', 400, 949)]
    rest = [interval(p, 0, n) for p, n in DOMAINS.items() if p != 'isolation']
    check([('distance9', {'jobs': rest + split + [interval('isolation', 200, 500)]})], {'distance9'}, True)
    check([('distance9', {'jobs': rest + [interval('isolation', 0, 948)]})], {'distance9'}, False)
    mutations = (
        ('status', 'FAILED'), ('checker_sha256', '0'*64), ('runner_sha256', '0'*64),
        ('wall_limit_seconds', 60), ('rss_limit_bytes', 2*1024**3), ('threads', 2),
        ('threads', True),
        ('elapsed_seconds', 60), ('peak_process_tree_rss_bytes', 1024**3),
        ('stdout_bytes', 1), ('lock_owner', 'checker'),
        ('initial_resources', dict(free_disk_gib=39, available_memory_gib=7)),
        ('command', [D13, '--verify', '--batch', 'unknown']),
        ('stdout', batch('finite')['stdout'].replace(CERT13, '0'*64)),
    )
    for key, value in mutations:
        bad = dict(batch('finite')); bad[key] = value
        check([('d13', {'jobs': complete13[:-1] + [bad]})], {'d13'}, False)
    bad = outer(batch('finite')); bad['stdout_sha256'] = '0'*64
    check([('d13', {'jobs': complete13[:-1] + [bad]})], {'d13'}, False)
    bad = interval('profiles', 0, 8); bad['command'][-1] = '7'
    check([('distance9', {'jobs': complete9[:-1] + [bad]})], {'distance9'}, False)
    bad = interval('profiles', 0, 8)
    proof = decode(bad['stdout']); proof['source'] = '0'*64
    bad['stdout'] = json.dumps(proof)+'\n'; bad['stdout_bytes'] = len(bad['stdout'].encode())
    check([('distance9', {'jobs': complete9[:-1] + [bad]})], {'distance9'}, False)
    check([('d13', {'jobs': complete13}), ('distance9', {'jobs': complete9})], {'d13', 'distance9'}, True)
    check([('d13', {'jobs': complete13}), ('distance9', {'status': 'PASS', 'jobs': []})],
          {'d13', 'distance9'}, False)
    try:
        decode('{"jobs":[],"jobs":[]}')
    except ValueError:
        pass
    else:
        raise ValueError('duplicate-key fixture not rejected')
    for malformed in ('{"x":NaN}', '{"x":1e9999}'):
        try:
            decode(malformed)
        except ValueError:
            pass
        else:
            raise ValueError('nonfinite-number fixture not rejected')
    return dict(status='PASS_COVERAGE_ORACLE_SELF_TEST_ONLY', synthetic=True,
                mathematical_execution_or_actual_coverage=False, adversarial_mutations=len(mutations))


def main():
    if not __debug__:
        raise RuntimeError('disabled assertions are not permitted')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--self-test', action='store_true')
    parser.add_argument('--d13-ledger', action='append', default=[], metavar='PATH')
    parser.add_argument('--distance9-ledger', action='append', default=[], metavar='PATH')
    args = parser.parse_args()
    signal.alarm(45)
    if args.self_test:
        need(not args.verify and not args.d13_ledger and not args.distance9_ledger,
             'self-test is separate from actual ledger verification')
        print(json.dumps(oracle_self_test(), sort_keys=True))
        return 0
    need(args.verify and (args.d13_ledger or args.distance9_ledger), 'supply --verify and at least one ledger')
    own = Path(__file__).resolve()
    own_hash = digest(own)
    for relative, expected in FIXED.items():
        need(digest(ROOT/relative) == expected, 'current pin mismatch: '+relative)
    selected = [('d13', p) for p in args.d13_ledger] + [('distance9', p) for p in args.distance9_ledger]
    need(len(selected) <= MAX_FILES, 'too many ledger files')
    snapshots, total = {}, 0

    def inputs():
        nonlocal total
        for family, name in selected:
            path = Path(name).resolve()
            need(path.is_file() and path.stat().st_size <= MAX_LEDGER, 'ledger missing/too large')
            with path.open('rb') as stream:
                raw = stream.read(MAX_LEDGER + 1)
            total += len(raw)
            need(len(raw) <= MAX_LEDGER and total <= MAX_TOTAL, 'ledger byte cap exceeded')
            stamp = sha256(raw).hexdigest()
            need(path not in snapshots or snapshots[path] == stamp, 'ledger changed between reads')
            snapshots[path] = stamp
            yield family, decode(raw)

    answer = analyze(inputs(), {f for f, _ in selected})
    for path, expected in snapshots.items():
        need(digest(path, MAX_LEDGER) == expected, 'ledger changed during check')
    for relative, expected in FIXED.items():
        need(digest(ROOT/relative) == expected, 'current pin changed: '+relative)
    need(digest(own) == own_hash, 'coverage verifier changed')
    answer.update(status=('PASS_RECORDED_PARTITION_COVERAGE' if answer['complete'] else
                          'INCOMPLETE_RECORDED_PARTITION_COVERAGE'),
                  scope='Recorded current-source execution/name/interval coverage only; no mathematical replay.',
                  ledger_files=len(snapshots), verifier_sha256=own_hash,
                  provenance_limit='Unsigned local ledgers are not authenticated execution attestations.')
    print(json.dumps(answer, sort_keys=True))
    return 0 if answer['complete'] else 1


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (ValueError, TypeError, KeyError, OSError, OverflowError, RecursionError) as error:
        print(json.dumps(dict(status='INVALID_RECORDED_PARTITION_COVERAGE', error=str(error)[:300])))
        sys.exit(2)
