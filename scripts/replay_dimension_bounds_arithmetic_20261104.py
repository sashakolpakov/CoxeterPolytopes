#!/usr/bin/env python3
"""Replay one frozen arithmetic combiner in an isolated, exact input snapshot.

--plan reads metadata only. Execution requires the CURRENT bounded parent
runner and its shared lock. Archived policy code is staged as inert hash input;
it is never imported or executed. No frozen source or live policy is modified.
"""
import argparse
import ast
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / 'repro/dimension_bounds/manifest.json'
WORKERS = {
    'A1': 'scripts/check_a17_two_isolated_quadratic_endpoints_20261022.py',
    'A2': 'scripts/check_a17_finite_cubic_endpoint_20261022.py',
    'A3': 'scripts/check_d17_one_isolated_quadratic_fork_20261022.py',
    'A4': 'scripts/check_arithmetic16_b16_isolated_quadratic_head_short_pairs_20261023.py',
}
IMPORTS = {'argparse', 'dataclasses', 'fractions', 'hashlib', 'json', 'math', 'pathlib'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def relative(value):
    path = Path(value)
    require(not path.is_absolute() and '..' not in path.parts and path.parts,
            'Unsafe relative path')
    return path


def source_path(value):
    path = (ROOT / relative(value)).resolve()
    require(ROOT in path.parents and path.is_file(), 'Missing local input: ' + value)
    require(path.stat().st_size <= 2 * 1024 * 1024, 'Snapshot file exceeds 2 MiB')
    return path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def inspect(item):
    data = json.loads(MANIFEST.read_text())
    entry = data['arithmetic_replays'][item]
    require(entry['worker'] == WORKERS[item], 'Unexpected scientific worker')
    worker = source_path(entry['worker'])
    require(sha(worker) == entry['worker_sha256'], 'Frozen worker hash mismatch')
    tree = ast.parse(worker.read_text())
    pins = None
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            require(all(x.name in IMPORTS for x in node.names), 'Non-stdlib worker import')
        elif isinstance(node, ast.ImportFrom):
            require(node.level == 0 and node.module in IMPORTS, 'Non-stdlib worker import')
        elif isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
            require(node.func.id not in {'exec', 'eval', '__import__', 'compile'},
                    'Dynamic worker execution refused')
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(
                isinstance(t, ast.Name) and t.id == 'PINS' for t in node.targets):
            require(pins is None, 'Duplicate PINS')
            pins = ast.literal_eval(node.value)
    require(isinstance(pins, dict), 'Missing literal PINS')
    inputs = entry['snapshot_inputs']
    require(set(pins) <= set(inputs), 'Incomplete pin closure')
    size = worker.stat().st_size
    for logical, binding in inputs.items():
        relative(logical)
        path = source_path(binding['source'])
        require(sha(path) == binding['sha256'], 'Snapshot input hash mismatch: ' + logical)
        if logical in pins:
            require(pins[logical] == binding['sha256'], 'Historical pin mismatch')
        size += path.stat().st_size
    require(size <= 4 * 1024 * 1024, 'Snapshot exceeds 4 MiB')
    command = entry['historical_command']
    require(command[0] == entry['worker'] and command.count('--output') == 1,
            'Unexpected command shape')
    require(command[1:3] == ['--expected-source-sha256', entry['worker_sha256']]
            and command[3] == 'combine', 'Only recorded combine commands are supported')
    part_count = 2 if item == 'A2' else 3
    require(len(command) == 6 + 4 * part_count and command[-2] == '--output',
            'Unexpected combine arguments')
    for i in range(4, len(command) - 2, 4):
        require(command[i] == '--part' and command[i + 2] == '--part-sha256',
                'Unexpected certificate argument order')
        require(command[i + 1] in inputs, 'Unbound certificate part')
        require(inputs[command[i + 1]]['sha256'] == command[i + 3],
                'Certificate argument hash differs from snapshot binding')
    require(len({command[i + 1] for i in range(4, len(command) - 2, 4)}) == part_count,
            'Duplicate certificate part')
    return entry, {'item': item, 'worker': entry['worker'],
                   'worker_sha256': entry['worker_sha256'],
                   'snapshot_files': len(inputs) + 1, 'snapshot_bytes': size,
                   'replay_scope': entry['replay_scope'],
                   'historical_runner_is_hash_input_only': True}


def new_output(value):
    path = (ROOT / relative(value)).resolve()
    require(path.parent == ROOT / 'results' and not path.exists(),
            'Choose a NEW file directly under results/')
    return path


def require_current_parent():
    require(os.environ.get('COXETER_CHECK_OWNER') == str(os.getppid()),
            'Use scripts/run_proof_batch_20261024.py check --ledger NEW -- ...')
    lock = (ROOT / 'formal/.proof-batches/exclusive.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
    else:
        fcntl.flock(lock, fcntl.LOCK_UN)
        lock.close()
        raise RuntimeError('Current bounded parent must hold exclusive.lock')
    require(not (ROOT / 'repro/RESEARCH_PAUSED').exists(), 'Research is paused')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('item', choices=tuple(WORKERS))
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--output', help='NEW replay certificate under results/')
    parser.add_argument('--receipt', help='NEW source/snapshot receipt under results/')
    args = parser.parse_args()
    def interrupted(signum, _frame):
        raise InterruptedError('Snapshot interrupted by signal ' + str(signum))
    for signum in (signal.SIGTERM, signal.SIGINT, signal.SIGHUP, signal.SIGALRM):
        signal.signal(signum, interrupted)
    entry, report = inspect(args.item)
    report['manifest_sha256'] = sha(MANIFEST)
    report['launcher_sha256'] = sha(Path(__file__))
    if args.plan:
        require(not args.output and not args.receipt, 'Plan does not create outputs')
        report.update(status='SNAPSHOT_METADATA_READY_NO_SCIENTIFIC_EXECUTION',
                      scientific_workers_executed=0)
        print(json.dumps(report, indent=2))
        return
    require(args.output and args.receipt, 'Execution needs --output and --receipt')
    require_current_parent()
    output, receipt = new_output(args.output), new_output(args.receipt)
    require(output != receipt, 'Output and receipt must differ')
    with tempfile.TemporaryDirectory(prefix='dimension-bounds-',
                                     dir=ROOT / 'formal/.proof-batches') as name:
        snapshot = Path(name)
        for logical, binding in entry['snapshot_inputs'].items():
            target = snapshot / relative(logical)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(source_path(binding['source']).read_bytes())
        worker = snapshot / entry['worker']
        worker.parent.mkdir(parents=True, exist_ok=True)
        worker.write_bytes(source_path(entry['worker']).read_bytes())
        (snapshot / 'results').mkdir(exist_ok=True)
        command = list(entry['historical_command'])
        command[command.index('--output') + 1] = 'results/local_replay.jsonl'
        require(not (snapshot / command[-1]).exists(), 'Output would overwrite an input')
        # The exact scientific source alone runs. Pinned historical runners are
        # read by its digest function, never selected as subprocesses/imports.
        subprocess.run([sys.executable, '-I', '-B', str(worker), *command[1:]],
                       cwd=snapshot, check=True, timeout=50)
        result = snapshot / 'results/local_replay.jsonl'
        require(result.is_file() and result.stat().st_size <= 2 * 1024 * 1024,
                'Missing or oversized replay certificate')
        with output.open('xb') as stream:
            stream.write(result.read_bytes())
        report.update(status='PASS_SCOPED_FROZEN_ARITHMETIC_COMBINE',
                      scientific_workers_executed=1, command=command,
                      output={'path': args.output, 'sha256': sha(output),
                              'bytes': output.stat().st_size},
                      full_geometric_formalization=False,
                      independent_implementation=False)
        with receipt.open('x') as stream:
            json.dump(report, stream, indent=2)
            stream.write('\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
