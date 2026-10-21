#!/usr/bin/env python3
"""Fixed parser-contract checks for the source-only D13 replay replacement.

Compare exact projections against independent stdlib json.loads on valid
adversarial fixtures, verify required rejections, and compare a nested
projection of the unchanged frozen D13 certificate.  No scientific module
is imported, no certificate is generated, and temporary fixtures are owned
by one TemporaryDirectory under /private/tmp.  Invoke with the shared
bounded wrapper --parent-lock; root owns execution.
"""
import argparse
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import signal
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
HELPER = 'scripts/check_d13_transfer_batches_20261020.py'
CERT = 'results/fresh_d13_three_rank15_context_transfer_20261018.json'
FIXED = {
    HELPER: 'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    CERT: '2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}


def digest(path):
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def baseline(raw, paths, optional=()):
    """Independent oracle: decode the complete JSON, then ordinary lookups.

    No scanner, parser tokenization, or projection implementation is reused.
    These fixtures use nonoverlapping object paths by construction.
    """
    decoded = json.loads(raw)
    if paths == ((),):
        return decoded
    result = {}
    for path in paths:
        source = decoded
        try:
            for key in path:
                source = source[key]
        except KeyError:
            assert len(path) == 1 and path[0] in optional
            continue
        target = result
        for key in path[:-1]:
            target = target.setdefault(key, {})
        target[path[-1]] = source
    return result


def valid_fixtures():
    # Literal fixtures deliberately exercise wire spellings, not merely
    # json.dumps roundtrips: exponent/negative-zero, escapes, and Unicode.
    literals = (
        (b'{"object":{"empty":{},"nested":{"v":7}},"array":[1,{"k":[true,false,null]}],"zero":-0.0,"exp":1.25e+2}',
         ((('object', 'nested', 'v'), ('array',), ('zero',), ('exp',)), ((),))),
        (b'{"skip":["} ] { [ , :",{"x":"\\\"quoted\\\" and \\\\ slash"}],"wanted":{"text":"\\\"[{}]\\\"\\\\","n":42},"after":null}',
         ((('wanted', 'text'), ('wanted', 'n'), ('after',)), (('skip',),))),
        ('{"κλειδί":{"雪":{"😀":"a\\tb\\nc\\rd\\bf\\f"}},"пусто":[],"é":"\\u005b\\u007b\\u005d\\u007d"}'.encode('utf-8'),
         ((('κλειδί', '雪', '😀'), ('пусто',), ('é',)), ((),))),
        (b' \n\t { "a\\\"[\\\\]" : {"b}":["\\uD83D\\uDE00",0]}, "keep":false, "unused": "\\\\\\\"{]" } \r\n',
         ((('a"[\\]', 'b}'), ('keep',)), (('keep',),))),
        (b'{"emptyObject":{},"emptyArray":[],"string":"","bool":true,"integer":123456789012345678901234567890,"null":null}',
         ((('emptyObject',), ('emptyArray',), ('string',), ('bool',), ('integer',), ('null',)),)),
        (b'{"optionalPresent":{"v":[1,2]},"required":"ok"}',
         ((('optionalPresent',), ('required',)),)),
        (b'{"required":{"v":"ok"},"skip":{"deep":[{"more":[[],{},"[\\\"{}\\\"]"]}]}}',
         ((('required', 'v'),),)),
    )
    fixtures = []
    for number, (raw, selections) in enumerate(literals):
        for selected in selections:
            fixtures.append((f'literal-{number}', raw, selected, ()))
    # Values in different branches and large skipped arrays must not alter
    # the projection.  Allocation remains small even for full stdlib decode.
    array = [{'i': i, 'text': '[}\\"' * (i % 5), 'nested': [i, None, False]}
             for i in range(400)]
    generated = json.dumps({'before': array, 'selected': {'x': ['雪', 3]},
                            'after': list(reversed(array))}, ensure_ascii=False).encode('utf-8')
    fixtures.append(('skipped-arrays', generated, (('selected', 'x'),), ()))
    fixtures.append(('optional-absent', b'{"required":7}',
                     (('required',), ('absent',)), ('absent',)))
    fixtures.append(('optional-present', b'{"required":7,"opt":null}',
                     (('required',), ('opt',)), ('opt',)))
    fixtures.append(('optional-only-empty', b'{}', (('absent',),), ('absent',)))
    for label, raw in (
        ('root-array', b'[1,{"q":"[}"},null]'),
        ('root-string', b'"scalar"'),
        ('root-number', b'-12.5e-2'),
        ('root-null', b'null'),
    ):
        fixtures.append((label, raw, ((),), ()))
    return fixtures


def rejection_fixtures():
    return (
        ('duplicate-root', b'{"a":1,"a":2}', (('a',),), ()),
        ('duplicate-ancestor', b'{"a":{"b":1,"b":2}}', (('a', 'b'),), ()),
        ('duplicate-unused-at-visited-root', b'{"a":1,"skip":[],"skip":{}}', (('a',),), ()),
        ('overlap-parent-first', b'{"a":{"b":1}}', (('a',), ('a', 'b')), ()),
        ('overlap-child-first', b'{"a":{"b":1}}', (('a', 'b'), ('a',)), ()),
        ('duplicate-selector', b'{"a":1}', (('a',), ('a',)), ()),
        ('missing-mandatory-root', b'{"other":1}', (('a',),), ()),
        ('missing-mandatory-nested', b'{"a":{"other":1}}', (('a', 'b'),), ()),
        ('optional-does-not-hide-mandatory', b'{"opt":1}', (('required',), ('opt',)), ('opt',)),
        ('trailing-data', b'{"a":1} {"b":2}', (('a',),), ()),
        ('trailing-scalar', b'{"a":1} false', (('a',),), ()),
        ('nonobject-ancestor', b'{"a":[1,2]}', (('a', 'b'),), ()),
        ('mismatched-skipped-delimiters', b'{"skip":[1},"a":7}', (('a',),), ()),
        ('invalid-selected-scalar', b'{"a":tru}', (('a',),), ()),
        ('truncated-string', b'{"skip":"unterminated,"a":7}', (('a',),), ()),
    )


def main():
    if not __debug__:
        raise RuntimeError('assertions must remain enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.parse_args()
    assert sys.dont_write_bytecode, '-B required'
    own = Path(__file__)
    own_hash = digest(own)
    for relative, expected in FIXED.items():
        assert digest(ROOT / relative) == expected, ('fixed input', relative)
    # Import definitions only: this module has no native scientific imports,
    # and its guarded main()/replay() are not called by this test.
    helper = importlib.import_module('check_d13_transfer_batches_20261020')
    owner = helper.parent_lock()
    def interrupt(signum, _frame):
        raise TimeoutError(('parser-contract interrupted', signum))
    for signum in (signal.SIGALRM, signal.SIGTERM, signal.SIGINT, signal.SIGHUP):
        signal.signal(signum, interrupt)
    signal.alarm(55)
    valid, rejected = [], []
    with TemporaryDirectory(prefix='coxeter_d13_json_contract_', dir='/private/tmp') as directory:
        parent = Path(directory)
        for index, (label, raw, paths, optional) in enumerate(valid_fixtures()):
            fixture = parent / f'valid_{index}.json'
            fixture.write_bytes(raw)
            expected = baseline(raw, paths, optional)
            actual = helper.JsonProjection(fixture).load(paths, optional)
            assert actual == expected, ('stdlib comparison', label, paths)
            assert fixture.read_bytes() == raw, ('fixture changed', label)
            valid.append(label)
        for index, (label, raw, paths, optional) in enumerate(rejection_fixtures()):
            fixture = parent / f'invalid_{index}.json'
            fixture.write_bytes(raw)
            try:
                helper.JsonProjection(fixture).load(paths, optional)
            except (AssertionError, json.JSONDecodeError):
                rejected.append(label)
            else:
                raise AssertionError(('required rejection did not occur', label))
        # Actual immutable file, independently decoded in full (1.1MiB).
        # Its saved scientific facts are NOT replayed by this parser test.
        certificate = ROOT / CERT
        assert certificate.stat().st_size < 2 * 1024**2
        raw = certificate.read_bytes()
        paths = (('status',), ('selected_support', 'rows'),
                 ('current_interfaces', 'later_corners'),
                 ('finite', 'complete_increasing_families'))
        assert helper.JsonProjection(certificate).load(paths) == baseline(raw, paths)
        assert sha256(raw).hexdigest() == FIXED[CERT]
    assert not parent.exists(), 'temporary fixtures not removed'
    assert os.getppid() == owner, 'parent exited'
    for relative, expected in FIXED.items():
        assert digest(ROOT / relative) == expected, ('fixed input changed', relative)
    assert digest(own) == own_hash, 'contract checker changed'
    print('PASS_D13_JSON_PROJECTION_CONTRACT', 'valid_projections', len(valid),
          'required_rejections', len(rejected), 'actual_certificate_projections', 1,
          'scientific_replay', False, 'temporary_fixtures_removed', True, flush=True)


if __name__ == '__main__':
    main()
