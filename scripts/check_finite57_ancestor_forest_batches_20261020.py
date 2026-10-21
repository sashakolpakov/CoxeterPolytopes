#!/usr/bin/env python3
"""Replay one complete inherited finite57 forest phase, without eager corpora.

Use --verify --phase coarse|increasing through the unchanged bounded wrapper
with --parent-lock.  The two phases check the frozen C.finite and
N.conditional_upper computations; they do not replay either producer's local
geometry, actual-type refinement or other output fields.  No generation mode.
"""
import argparse
import ast
import fcntl
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import signal


ROOT = Path(__file__).resolve().parents[1]
C = 'fresh_a16_d16_witness_context_ambient58_frontier_20261018'
N = 'fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018'
PROJECTOR = 'check_d13_transfer_batches_20261020'
FIXED = {
    f'scripts/{C}.py':
        '69dfbe890dd3de9c7bd54a5a79d50a9cbf38a515741e3ea464d3ef90c4090c9b',
    f'results/{C}.json':
        '9e5838e5eb9afec9c5b050d4d753e2ec28dc962c788af0fc7e7d5cc144ccc70e',
    f'scripts/{N}.py':
        'cf43d27af47f9f23216b89f6ef76fe35ccb5973e33a1f4216ff12fe60293cdac',
    f'results/{N}.json':
        'f80edecac4af2021e39806a194befebcba0ab83c748138a6e40556a7dedc117c',
    f'scripts/{PROJECTOR}.py':
        'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}
# SHA-256 of the canonical JSON map {relative .py path: source SHA-256}.
# This is the complete static local import closure of C and N, including
# imports inside function bodies.  It was computed by metadata-only AST
# inspection.  Neither native closure contains dynamic Python imports.
SOURCE_CLOSURE_COUNT = 120
SOURCE_CLOSURE_SHA256 = (
    '97a5f76d8017a6e9b3f5210815bc3b7de89869e793c3b2a85219ac2ba492e505')
COARSE_FLAGS = {
    'actual_maximum_claimed', 'actual_attainment_of_global_upper_claimed'}
COARSE_APPENDED = {'exact_actual_maximum', 'actual_attaining_forest'}
COARSE_EXCLUDED = COARSE_FLAGS | COARSE_APPENDED
STATUSES = {
    'coarse': 'PASS_CONDITIONAL_A16_D16_WITNESS_CONTEXT_AMBIENT58_FRONTIER',
    'increasing': 'PASS_D17_A9_A2_WITNESS_CONTEXT_ZERO_AMBIENT58_FRONTIER',
}


def digest(relative):
    path = (ROOT / relative).resolve()
    assert path.is_relative_to(ROOT) and path.is_file(), relative
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def canonical_hash(value):
    return sha256(json.dumps(value, sort_keys=True,
                             separators=(',', ':')).encode()).hexdigest()


def parent_lock():
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 1 and os.getppid() == owner, 'bounded parent wrapper required'
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('r+') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise AssertionError('--parent-lock required')
    return owner


def native_source_pins():
    """Identify/hash local definitions without importing scientific code."""
    seen = set()
    pending = [C, N]
    pins = {}
    while pending:
        name = pending.pop()
        if name in seen:
            continue
        seen.add(name)
        relative = f'scripts/{name}.py'
        path = ROOT / relative
        pins[relative] = digest(relative)
        tree = ast.parse(path.read_text(), filename=relative)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert node.level == 0, ('unexpected relative import', relative)
                names = [node.module]
            else:
                continue
            for imported in names:
                if (ROOT / 'scripts' / f'{imported}.py').is_file():
                    pending.append(imported)
    assert len(pins) == SOURCE_CLOSURE_COUNT, 'local source-domain changed'
    assert canonical_hash(pins) == SOURCE_CLOSURE_SHA256, 'local source pins'
    return pins


def equal(native, expected, encode):
    # Same enc -> JSON normalization used by the unchanged saved producers.
    actual = json.loads(json.dumps(encode(native)))
    assert actual == expected, 'complete saved phase output differs'


def replay_coarse(native, expected):
    actual = native.finite(expected)
    # Original run changes these two flags only after a separate actual-type
    # refinement and then appends two actual-refinement fields.  This phase
    # excludes exactly those four fields; every coarse field is compared.
    assert COARSE_EXCLUDED <= set(expected), 'coarse output schema'
    assert COARSE_FLAGS <= set(actual), 'native coarse output schema'
    assert set(expected) - set(actual) == COARSE_APPENDED
    assert not set(actual) - set(expected)
    equal({k: v for k, v in actual.items() if k not in COARSE_EXCLUDED},
          {k: v for k, v in expected.items() if k not in COARSE_EXCLUDED},
          native.enc)
    assert actual['fixed_large_vectors'] == 1295
    assert actual['fixed_large_and_A5_vectors'] == 3298
    assert len(actual['tables']) == 6596
    assert actual['Bellman_upper_inequalities'] == 8287874
    assert actual['attaining_predecessors'] == 382568
    assert actual['context_options'] == 10677
    return dict(fixed_vectors=1295, expanded_choices=3298, tables=6596,
                transition_inequalities=8287874, tight_predecessors=382568,
                context_options=10677, coarse_upper=str(actual['certified_upper']),
                excluded_actual_refinement_fields=sorted(COARSE_EXCLUDED))


def replay_increasing(native, expected):
    # conditional_upper itself compares every normalized reconstructed row
    # with the corresponding saved row, after validating its full key domain.
    actual = native.conditional_upper(expected)
    equal(actual, expected, native.enc)
    assert actual['complete_domain_rows'] == len(actual['rows']) == 759
    assert sum(len(row['tables']) for row in actual['rows']) == 1518
    assert actual['Bellman_inequalities'] == 1907367
    assert actual['attaining_predecessors'] == 88044
    assert actual['item_identities'] == 37191
    assert actual['context_options'] == 4453
    return dict(rows=759, tables=1518, transition_inequalities=1907367,
                tight_predecessors=88044, item_identities=37191,
                context_options=4453, upper=str(actual['maximum']))


def main():
    if not __debug__:
        raise RuntimeError('assertions must remain enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--phase', choices=('coarse', 'increasing'), required=True)
    args = parser.parse_args()
    owner = parent_lock()
    signal.alarm(45)  # Outer wrapper retains its unchanged59s/1GiB watchdog.
    try:
        self_path = Path(__file__).relative_to(ROOT).as_posix()
        self_hash = digest(self_path)
        for path, expected in FIXED.items():
            assert digest(path) == expected, ('fixed pin', path)
        source_pins = native_source_pins()
        projector = importlib.import_module(PROJECTOR)
        assert Path(projector.__file__).resolve() == ROOT / 'scripts' / f'{PROJECTOR}.py'
        stem = C if args.phase == 'coarse' else N
        field = 'finite' if args.phase == 'coarse' else 'increasing_family_upper'
        certificate = f'results/{stem}.json'
        # The selected objects are20,639,009 and5,815,948 serialized bytes,
        # below JsonProjection's unchanged24MiB selected-value cap.  All other
        # certificate fields are structurally skipped, never decoded.
        saved = projector.JsonProjection(certificate).load(
            (('status',), ('sources',), (field,)))
        assert saved['status'] == STATUSES[args.phase]
        assert isinstance(saved['sources'], dict) and saved['sources']
        expected_pins = dict(FIXED)
        for path, expected in saved['sources'].items():
            assert isinstance(expected, str) and len(expected) == 64
            assert path not in expected_pins or expected_pins[path] == expected
            expected_pins[path] = expected
        for path, expected in source_pins.items():
            assert path not in expected_pins or expected_pins[path] == expected
            expected_pins[path] = expected
        before = {path: digest(path) for path in expected_pins}
        assert before == expected_pins, 'saved input/source pins'
        # Definition imports only.  Neither unchanged original run() is called;
        # in particular N's149-document eager loader is never entered.
        native = importlib.import_module(stem)
        assert Path(native.__file__).resolve() == ROOT / 'scripts' / f'{stem}.py'
        result = (replay_coarse(native, saved[field]) if args.phase == 'coarse'
                  else replay_increasing(native, saved[field]))
        assert os.getppid() == owner, 'parent exited'
        assert {path: digest(path) for path in before} == before, 'input/source changed'
        assert native_source_pins() == source_pins, 'local import domain changed'
        assert digest(self_path) == self_hash, 'worker changed'
        print('PASS_FINITE57_ANCESTOR_FOREST_PHASE', args.phase,
              json.dumps(dict(result, certificate_sha256=FIXED[certificate],
                              source_pins=len(before),
                              native_source_closure_sha256=SOURCE_CLOSURE_SHA256,
                              full_producer_replay_claimed=False,
                              mathematical_jobs_run=1), sort_keys=True), flush=True)
    finally:
        signal.alarm(0)


if __name__ == '__main__':
    main()
