#!/usr/bin/env python3
"""Independent stdlib fixtures for the bounded distance9 JSON row scanner.

No mathematical producer or replay is invoked. Use the existing parent-lock
wrapper. Temporary module rebinding is process-local and always restored;
only files in this checker's owned TemporaryDirectory are written.
"""
import argparse
import ast
from contextlib import contextmanager
from copy import deepcopy
import fcntl
from hashlib import sha256
import importlib
import json
import os
from pathlib import Path
import signal
import tempfile

ROOT = Path(__file__).resolve().parents[1]
WORKER = 'scripts/check_selective_distance9_batches_20261020.py'
HELPER = 'scripts/check_d13_transfer_batches_20261020.py'
CERT = 'results/fresh_selective_distance9_dimension108_20261017.json'
PINS = {
    WORKER: 'f5f0b47084135f8bb1886a37e0478d1f6b8c9c4327499f986a8ec9e173ecd843',
    HELPER: 'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    CERT: 'e7c6d1a1e56d8788e7d7a7d02bd5dd71837e3a3aef4ae768c1efc0a187bcacd9',
}
MAX_FIXTURE = 32 * 1024
MAX_REAL_PREFIX = 2 * 1024**2
FAILURES = (AssertionError, ValueError, IndexError, KeyError, TypeError)


def digest(path):
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def parent_lock():
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 1 and os.getppid() == owner, 'bounded parent wrapper required'
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise AssertionError('--parent-lock required')
    return owner


def must_fail(call, label):
    try:
        call()
    except FAILURES:
        return
    raise AssertionError(('unexpected acceptance', label))


@contextmanager
def fixture(worker, directory, payload, **limits):
    assert isinstance(payload, bytes) and len(payload) <= MAX_FIXTURE
    names = ('ROOT', 'CERT', 'MAX_FILE', 'MAX_ROW')
    before = {name: getattr(worker, name) for name in names}
    path = directory / 'fixture.json'
    path.write_bytes(payload)
    worker.ROOT, worker.CERT = directory, path.name
    for name, value in limits.items():
        assert name in ('MAX_FILE', 'MAX_ROW') and 0 < value <= MAX_FIXTURE
        setattr(worker, name, value)
    try:
        yield
    finally:
        for name, value in before.items():
            setattr(worker, name, value)
        path.unlink()
        assert {name: getattr(worker, name) for name in names} == before


def oracle_path(document, path):
    value = document
    for key in path:
        assert isinstance(value, dict) and key in value
        value = value[key]
    assert isinstance(value, list)
    return value


def valid_case(worker, directory, payload, path, label):
    # A complete, independent stdlib parse supplies the fixture oracle.
    oracle = oracle_path(json.loads(payload), path)
    with fixture(worker, directory, payload):
        visited = []
        total = worker.Arrays.walk(path, lambda i, row: visited.append((i, row)))
        assert total == len(oracle), ('row count', label)
        assert visited == list(enumerate(oracle)), ('walk order/value', label)
        predicates = (
            ('all', lambda i, row: True),
            ('none', lambda i, row: False),
            ('odd indices', lambda i, row: i % 2 == 1),
            ('bounded interval', lambda i, row: 1 <= i < 3),
            ('row predicate', lambda i, row: isinstance(row, dict)),
        )
        for name, predicate in predicates:
            rows, count = worker.Arrays.select(path, predicate)
            expected = [row for i, row in enumerate(oracle) if predicate(i, row)]
            assert count == len(oracle) and rows == expected, ('selection', label, name)


def positive_fixtures(worker, directory):
    odd = 'quotes " \\ slash / brackets ] } [ { newline\n tab\t 雪 é 😀'
    examples = [
        ('empty nested array', b'{"scope":{"rows":[]}}', ('scope', 'rows')),
        ('empty root array', b'[]', ()),
        ('empty root array trailing whitespace', b' \t[]\n ', ()),
        ('root trailing whitespace regression', b' \n[1,{"two":[2,3]},null]\t \n', ()),
        ('mixed rows', b'{"rows":[null,true,false,0,-12,1.25,"x",[],{"a":[1]}]}', ('rows',)),
        ('escaped unicode', json.dumps({'scope': {'rows': [odd, {'v': odd}, [odd]]}},
                                      ensure_ascii=True).encode(), ('scope', 'rows')),
        ('literal unicode', json.dumps({'雪': {'é😀': [odd, {'nested': ['é', '\u2028']}, '雪']}},
                                       ensure_ascii=False).encode(), ('雪', 'é😀')),
        ('exact path versus decoy', json.dumps({
            'skip': {'scope': {'rows': ['wrong']}},
            'scope': {'text': '"rows":[null] } ]', 'rows': [1, 2, 3, 4]},
            'after': [{'scope': {'rows': ['also wrong']}}],
        }).encode(), ('scope', 'rows')),
        ('punctuation in path keys', json.dumps({'a.b': {'[0]': {'rows': [{'x': '}"['}, 7]}}}).encode(),
         ('a.b', '[0]', 'rows')),
        ('object root surrounding whitespace', b' \t{"scope":{"rows":[0,1,2]}} \r\n', ('scope', 'rows')),
        # Duplicate keys inside a selected item follow json.loads semantics.
        # Only object-path containers, not recursively decoded item objects,
        # are advertised as rejecting duplicates by the production reader.
        ('selected row duplicates match stdlib', b'{"scope":{"rows":[{"x":1,"x":2}]}}', ('scope', 'rows')),
        ('unvisited duplicate subtree is not validated',
         b'{"skip":{"x":1,"x":2},"scope":{"rows":[3]}}', ('scope', 'rows')),
    ]
    for label, payload, path in examples:
        valid_case(worker, directory, payload, path, label)
    payload = b'{"scope":{"rows":[],"other":0},"outside":false}'
    oracle = json.loads(payload)
    with fixture(worker, directory, payload):
        assert worker.Arrays.keys() == set(oracle)
        assert worker.Arrays.keys(('scope',)) == set(oracle['scope'])
    return len(examples) + 1


def negative_fixtures(worker, directory):
    count = 0
    # These malformed bytes affect the selected array, a visited container,
    # or the root's structural/trailing boundary; skipped subtree grammar is
    # deliberately not asserted to be fully validated.
    malformed = [
        b'', b'{"scope":{"rows":[1,]}}', b'{"scope":{"rows":[1 2]}}',
        b'{"scope":{"rows":[01]}}', b'{"scope":{"rows":[{"a":}]}}',
        b'{"scope":{"rows":["bad\\q"]}}',
        b'{"scope":{"rows":["unterminated]}}',
        b'{"scope":{"rows":[1,2]false}}',
        b'{"scope":{"rows":[1,2]}} trailing',
        b'{"scope":{"rows":[1,2]}} {}',
    ]
    full = b'{"scope":{"rows":[1,{"x":"y"},[2,3]]}}'
    malformed.extend(full[:-cut] for cut in (1, 2, 5, 11, 20))
    for payload in malformed:
        must_fail(lambda: json.loads(payload), 'stdlib malformed oracle')
        with fixture(worker, directory, payload):
            must_fail(lambda: worker.Arrays.select(('scope', 'rows')), 'malformed selected array')
        count += 1
    for payload in (b'[1,]', b'[1 2]', b'[1]true', b'[1]\n{}', b'[1,2'):
        must_fail(lambda: json.loads(payload), 'stdlib malformed root array')
        with fixture(worker, directory, payload):
            must_fail(lambda: worker.Arrays.select(()), 'malformed root array')
        count += 1
    # Visited path objects reject duplicates even though json.loads itself
    # accepts them. The oracle value is not used to turn duplicates into a
    # valid selected domain.
    duplicates = (
        b'{"scope":{"rows":[0]},"scope":{"rows":[1]}}',
        b'{"scope":{"rows":[0],"rows":[1]}}',
        b'{"before":0,"before":1,"scope":{"rows":[2]}}',
    )
    for payload in duplicates:
        with fixture(worker, directory, payload):
            must_fail(lambda: worker.Arrays.select(('scope', 'rows')), 'visited duplicate object key')
        count += 1
    payload = b'{"scope":{"rows":[1],"scalar":2},"nothing":null}'
    for path in (('missing',), ('scope', 'absent'), ('scope',),
                 ('nothing',), ('scope', 'scalar'), ('scope', 'rows', '0')):
        with fixture(worker, directory, payload):
            must_fail(lambda: oracle_path(json.loads(payload), path), 'oracle wrong/missing path')
            must_fail(lambda: worker.Arrays.select(path), 'scanner wrong/missing path')
        count += 1
    with fixture(worker, directory, payload):
        must_fail(lambda: worker.Arrays.keys(('scope', 'rows')), 'keys requires object')
    count += 1
    with fixture(worker, directory, b'{"rows":["0123456789"]}', MAX_ROW=8):
        must_fail(lambda: worker.Arrays.select(('rows',)), 'item allocation cap')
    count += 1
    with fixture(worker, directory, b'{"rows":[1]}', MAX_FILE=8):
        must_fail(lambda: worker.Arrays.select(('rows',)), 'file allocation cap')
    return count + 1


def optimization_guard(relative):
    # Inspect and execute ONLY the initial guard, not the worker main body.
    # optimize=1 makes __debug__ False without starting a second process or
    # importing any mathematical producer. An assert would disappear here.
    path = ROOT / relative
    data = path.read_bytes()
    assert len(data) < 256 * 1024
    tree = ast.parse(data, filename=relative)
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == 'main')
    guard = main.body[0]
    assert isinstance(guard, ast.If) and isinstance(guard.test, ast.UnaryOp)
    assert isinstance(guard.test.op, ast.Not) and isinstance(guard.test.operand, ast.Name)
    assert guard.test.operand.id == '__debug__' and not guard.orelse
    assert len(guard.body) == 1 and isinstance(guard.body[0], ast.Raise)
    error = guard.body[0].exc
    assert isinstance(error, ast.Call) and isinstance(error.func, ast.Name) and error.func.id == 'RuntimeError'
    fragment = ast.Module(body=[deepcopy(guard)], type_ignores=[])
    code = compile(ast.fix_missing_locations(fragment), relative + ':guard-only', 'exec', optimize=1)
    try:
        exec(code, {'RuntimeError': RuntimeError})
    except RuntimeError as exc:
        assert str(exc) == 'assertions must remain enabled'
    else:
        raise AssertionError(('optimization guard disappeared', relative))


def real_prefix_oracle(path, object_path):
    """Independent stdlib decoder; read at most a2MiB certificate prefix.

    Descend JSON objects by decoding their keys and complete preceding values
    with JSONDecoder.raw_decode. Stop once the selected array ends, then use
    json.loads on that exact substring. No production scanner/regex is used.
    This is not a parse of the entire21MiB scientific certificate.
    """
    with path.open('rb') as stream:
        prefix = stream.read(MAX_REAL_PREFIX)
    text = prefix.decode('utf-8')
    decoder = json.JSONDecoder()

    def whitespace(i):
        while i < len(text) and text[i] in ' \r\n\t':
            i += 1
        return i

    def descend(i, keys):
        i = whitespace(i)
        assert text[i] == '{'
        i = whitespace(i + 1)
        seen = set()
        while True:
            key, i = decoder.raw_decode(text, i)
            assert isinstance(key, str) and key not in seen
            seen.add(key)
            i = whitespace(i)
            assert text[i] == ':'
            i = whitespace(i + 1)
            if key == keys[0]:
                if len(keys) > 1:
                    return descend(i, keys[1:])
                _, end = decoder.raw_decode(text, i)
                assert text[whitespace(end)] in ',}'
                result = json.loads(text[i:end])
                assert isinstance(result, list)
                return result
            _, i = decoder.raw_decode(text, i)
            i = whitespace(i)
            assert text[i] == ',', ('requested prefix path not found', object_path)
            i = whitespace(i + 1)
    assert object_path
    return descend(0, object_path)


def main():
    if not __debug__:
        raise RuntimeError('assertions must remain enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.parse_args()
    owner = parent_lock()
    signal.alarm(45)
    own_path = Path(__file__)
    own_hash = digest(own_path)
    for relative, expected in PINS.items():
        assert digest(ROOT / relative) == expected, ('fixture pin', relative)
    worker = importlib.import_module('check_selective_distance9_batches_20261020')
    original = {name: getattr(worker, name) for name in ('ROOT', 'CERT', 'MAX_FILE', 'MAX_ROW')}
    assert worker.ROOT == ROOT and worker.CERT == CERT
    optimization_guard(WORKER)
    optimization_guard(HELPER)
    with tempfile.TemporaryDirectory(prefix='distance9-scanner-fixtures-') as temporary:
        directory = Path(temporary)
        positive = positive_fixtures(worker, directory)
        negative = negative_fixtures(worker, directory)
        assert not any(directory.iterdir()), 'owned fixture files leaked'
    assert not directory.exists(), 'owned temp directory leaked'
    assert {name: getattr(worker, name) for name in original} == original
    oracle = real_prefix_oracle(ROOT / CERT, ('census', 'ports'))
    visited = []
    total = worker.Arrays.walk(('census', 'ports'), lambda i, row: visited.append((i, row)))
    assert total == len(oracle) == 1262 and visited == list(enumerate(oracle))
    selected, count = worker.Arrays.select(('census', 'ports'),
                                           lambda i, row: row['allowed'] and i % 3 == 1)
    assert count == 1262
    assert selected == [row for i, row in enumerate(oracle) if row['allowed'] and i % 3 == 1]
    assert os.getppid() == owner, 'parent exited'
    for relative, expected in PINS.items():
        assert digest(ROOT / relative) == expected, ('fixture input changed', relative)
    assert digest(own_path) == own_hash, 'fixture checker changed'
    print(json.dumps(dict(status='PASS_DISTANCE9_SCANNER_FIXTURES',
                          positive_cases=positive, negative_cases=negative,
                          real_array_rows=1262, real_prefix_cap=MAX_REAL_PREFIX,
                          worker_sha256=PINS[WORKER], checker_sha256=own_hash,
                          scope='Object-path and row projection tests, not a full JSON validator of skipped subtrees.'),
                     sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
