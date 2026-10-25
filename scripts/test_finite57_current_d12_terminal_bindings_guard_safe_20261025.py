#!/usr/bin/env python3
"""Structural regressions for the guard-safe selected-D12 payment checker.

Run only through the admitted serial proof entry. These small synthetic
finite forests test the validator; they are not geometric receipt evidence.
No mathematical producer or repository worker is imported.
"""
import argparse
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import tempfile
import types

ROOT = Path(__file__).resolve().parents[1]
CHECKER = 'scripts/check_finite57_current_d12_terminal_bindings_guard_safe_20261025.py'
LIBRARY = 'scripts/finite57_literal_weights_audit_20261024.py'
PINS = {
    CHECKER: 'cb24da35f459dc3b54240ddbf84048c0cfcbeec21f58f83152c8edf8ea9a5b79',
    LIBRARY: 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def captured(relative):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(128 * 1024 + 1)
    require(0 < len(data) <= 128 * 1024, 'source read cap')
    require(sha256(data).hexdigest() == PINS[relative], 'source pin: ' + relative)
    return data


def load(name, relative, data):
    module = types.ModuleType(name)
    module.__file__ = str(ROOT / relative)
    exec(compile(data, module.__file__, 'exec'), module.__dict__)
    return module


def gram(kind, rank, weights, isolated=0):
    graph = weights.classical_graph(kind, rank)
    result = [[2 * int(i == j) for j in range(rank + isolated)]
              for i in range(rank + isolated)]
    for i, row in enumerate(graph):
        for j, label in row.items():
            require(label == 3, 'synthetic test uses only simple bonds')
            result[i][j] = -1
    return result


def row(kind, rank, marks, images, payment):
    return dict(marks=marks, payment=payment,
                roles=[dict(nodes=list(range(rank)), type=kind + str(rank))],
                images=images)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixture-parent', required=True)
    args = parser.parse_args()
    parent = Path(args.fixture_parent).resolve()
    require(parent.parent == Path('/private/tmp'), 'owned fixture must be in /private/tmp')
    require(parent.name.startswith('coxeter-finite57-d12-structural-'), 'owned fixture prefix')
    require(not parent.exists(), 'fixture parent already exists')
    checker_bytes, library_bytes = captured(CHECKER), captured(LIBRARY)
    checker = load('finite57_d12_structural_subject', CHECKER, checker_bytes)
    weights = load('finite57_d12_structural_literal', LIBRARY, library_bytes)
    # The real batch is still protected by the wrapper's actual shared lock.
    checker.lock()
    cases = []

    def positive(name, action):
        action()
        cases.append(dict(name=name, expected='success', status='PASS'))

    def rejects(name, clause, action):
        try:
            action()
        except ValueError as error:
            require(clause in str(error), 'wrong rejection for ' + name + ': ' + str(error))
        else:
            raise ValueError('invalid input accepted: ' + name)
        cases.append(dict(name=name, expected=clause, status='PASS'))

    small = row('A', 3, [0, 1], [[0, 1], [1, 2]], '1')
    small_gram = gram('A', 3, weights)
    positive('all A3 native marked images',
             lambda: checker.bind(small, small_gram, [2], weights, [0]))

    missing = deepcopy(small)
    missing['images'] = [[0, 1]]
    rejects('missing marked image', 'same-component pair binding',
            lambda: checker.bind(missing, small_gram, [2], weights, [0]))
    repeated = deepcopy(small)
    repeated['images'] = [[0, 1], [0, 1], [1, 2]]
    rejects('duplicate marked image', 'sorted and distinct',
            lambda: checker.bind(repeated, small_gram, [2], weights, [0]))
    expensive = deepcopy(small)
    expensive['payment'] = '2'
    rejects('chosen lower above current value', 'below chosen terminal payment',
            lambda: checker.bind(expensive, small_gram, [2], weights, [0]))

    guarded = row('A', 13, [0, 12], [[0, 12]], '1')
    no_witness = gram('A', 13, weights)
    guarded_white = list(range(1, 12))
    rejects('unearned A13 endpoint guard', 'unconditional near-end pairs',
            lambda: checker.bind(guarded, no_witness, guarded_white, weights, [0]))
    with_witness = gram('A', 13, weights, isolated=1)
    rejects('displayed A1 earns no G0 guard credit', 'unconditional near-end pairs',
            lambda: checker.bind(guarded, with_witness,
                                 guarded_white + [13], weights, [0]))
    near_end = row('A', 13, [0, 11], [[0, 11], [1, 12]], '1')
    positive('unconditional A13 near-end unit',
             lambda: checker.bind(near_end, no_witness,
                                  [i for i in range(13) if i not in (0, 11)], weights, [0]))
    false_witness = deepcopy(with_witness)
    false_witness[12][13] = false_witness[13][12] = -1
    rejects('attached node is not a whole witness', 'entire displayed component',
            lambda: checker.bind(guarded, false_witness,
                                 guarded_white + [13], weights, [0]))

    negative_guard = row('D', 13, [0, 11], [[0, 11], [0, 12]], '1')
    d13 = gram('D', 13, weights)
    d13_white = [i for i in range(13) if i not in (0, 11)]
    rejects('changed D13 pair refuses the G0 domain', 'outside this actual332 domain',
            lambda: checker.bind(negative_guard, d13, d13_white, weights, [0]))
    conservative = deepcopy(negative_guard)
    conservative['payment'] = '1/2'
    rejects('even a half does not bypass the guarded-role domain', 'outside this actual332 domain',
            lambda: checker.bind(conservative, d13, d13_white, weights, [0]))
    partial_a2 = gram('D', 13, weights, isolated=2)
    partial_a2[13][14] = partial_a2[14][13] = -1
    rejects('unsealed A2 cannot dismiss a possible rank15 witness', 'outside this actual332 domain',
            lambda: checker.bind(negative_guard, partial_a2,
                                 d13_white + [13, 14], weights, [0]))
    center_unit = row('D', 13, [0, 10], [[0, 10]], '1')
    positive('D13 long-end versus fork-center ordinary unit',
             lambda: checker.bind(center_unit, d13,
                                  [i for i in range(13) if i not in (0, 10)], weights, [0]))

    extended_row = row('A', 4, [0, 1], [[0, 1], [2, 3]], '1')
    extended_gram = gram('A', 4, weights)
    document = dict(gram=small_gram, white=[2],
                    extended_gram=extended_gram, extended_white=[2, 3],
                    before=small, after=extended_row)

    def correct_extended_context():
        actual = {path: (gg, ww, gp, wp)
                  for path, _, gg, ww, gp, wp in checker.terminals(document, fuel=[0])}
        require(set(actual) == {'/before', '/after'}, 'exact synthetic terminal paths')
        require(actual['/before'] == (small_gram, [2], '/gram', '/white'), 'old context')
        require(actual['/after'] == (extended_gram, [2, 3],
                                     '/extended_gram', '/extended_white'), 'extended context')
        for gg, ww, gp, wp in actual.values():
            require(checker.pointer(document, gp) == gg, 'Gram pointer')
            require(checker.pointer(document, wp) == ww, 'white pointer')

    positive('before and after retain distinct physical provenance', correct_extended_context)
    rejects('using old context for an extended receipt', 'entire displayed component',
            lambda: checker.bind(extended_row, small_gram, [2], weights, [0]))
    rejects('duplicate JSON keys', 'duplicate JSON key',
            lambda: checker.decode(b'{"a":1,"a":2}'))
    rejects('nonfinite JSON number', 'nonfinite JSON number',
            lambda: checker.decode(b'{"a":NaN}'))
    rejects('Boolean is not a rational', 'exact rational type',
            lambda: checker.fraction(True))
    rejects('excessive rational length', 'bounded exact rational',
            lambda: checker.fraction('1' * 33))
    rejects('fuel refuses instead of truncating an image set', 'fuel exhausted',
            lambda: checker.images(small_gram, [0, 1, 2], [0, 1], ('A', 3),
                                   weights, [checker.MAX_VISITS]))

    parent.mkdir()
    try:
        with tempfile.TemporaryDirectory(prefix='owned-', dir=parent) as temporary:
            fixture = Path(temporary)
            (fixture / 'formal/.proof-batches').mkdir(parents=True)
            (fixture / 'formal/.proof-batches/exclusive.lock').touch()
            (fixture / 'oversized').write_bytes(b'123456789')
            previous_root = checker.ROOT
            checker.ROOT = fixture
            try:
                rejects('missing parent lock', 'must already be held', checker.lock)
                rejects('actual input read cap', 'bounded input',
                        lambda: checker.raw('oversized', maximum=8))
            finally:
                checker.ROOT = previous_root
    finally:
        parent.rmdir()

    require(captured(CHECKER) == checker_bytes and captured(LIBRARY) == library_bytes,
            'source changed during structural checks')
    print(json.dumps(dict(status='PASS_GUARD_SAFE_D12_TERMINAL_BINDER_STRUCTURAL_REGRESSIONS',
        source_pins=PINS, positive_cases=sum(c['expected'] == 'success' for c in cases),
        negative_cases=sum(c['expected'] != 'success' for c in cases), cases=cases,
        scope='Synthetic validator behavior only. No fixed selected-D12 receipt, geometric whole-component proof or full finite57 claim is certified.'),
        sort_keys=True))


if __name__ == '__main__':
    main()
