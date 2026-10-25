#!/usr/bin/env python3
"""Prepared source ONLY: bind the332 selected-D12 finite outputs to the literal rule.

No repository mathematical worker is imported. Root must review this source and
launch one named batch through the unchanged guarded bounded runner. This is
neither a geometric whole-component proof nor the complete terminal inventory.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization or inherited PYTHONOPTIMIZE is forbidden')

import argparse
import fcntl
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import re
import types

ROOT = Path(__file__).resolve().parents[1]
LIB = 'scripts/finite57_literal_weights_audit_20261024.py'
TEX = 'manuscripts/finite_volume57.tex'
REFERENCE = 'results/fresh_a1_d5_witness_context_retention_20261018.json'
CASES = {
    'd12-first': ('results/fresh_d12_radical_saturated_half_replacements_20261018.json',
                  '2d39a231e9d65d295017695913c32951df792edcf1cfc876ab0c272f199b379c', 280),
    'd12-rest': ('results/fresh_d12_remaining_selected_half_replacements_20261018.json',
                 '931c5abc3ff4807a8bc31afb6a8c54ad8f8aeabf4c4b7cfe24288c715fc755a3', 52),
}
PINS = {LIB: 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0',
        TEX: 'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0',
        REFERENCE: 'f5bfafea703eb81e1f57a458534b7ca138a1e993e1db56cc243f9f823d39e907'}
MAX_FILE = 4 * 1024**2
MAX_NODES = 28
MAX_VISITS = 200000
TYPE = re.compile(r'([ABD])(\d{1,2})\Z|E6\Z')
NUMBER = re.compile(r'-?\d{1,32}(?:/\d{1,32})?\Z')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def raw(relative, maximum=MAX_FILE):
    path = (ROOT / relative).resolve()
    require(path.is_relative_to(ROOT) and path.is_file(), 'input path')
    with path.open('rb') as stream:
        data = stream.read(maximum + 1)
    require(0 < len(data) <= maximum, 'bounded input: ' + relative)
    return data


def unique(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def fraction(value):
    require(type(value) in (str, int), 'exact rational type')
    require(NUMBER.fullmatch(str(value)) is not None, 'bounded exact rational')
    return F(value)


def indices(value, upper, length=None):
    require(type(value) is list and (length is None or len(value) == length), 'index list')
    require(all(type(i) is int and 0 <= i < upper for i in value), 'index type/range')
    require(len(value) == len(set(value)), 'duplicate node')
    return value


def matrix(value):
    require(type(value) is list and 0 < len(value) <= MAX_NODES, 'matrix order')
    n = len(value)
    require(all(type(row) is list and len(row) == n for row in value), 'square matrix')
    g = [[fraction(x) for x in row] for row in value]
    require(all(g[i][i] > 0 for i in range(n)), 'positive root norms')
    require(all(g[i][j] == g[j][i] <= 0 for i in range(n) for j in range(i)), 'acute symmetric Gram')
    return g


def components(g, ids):
    remaining = set(ids)
    result = []
    while remaining:
        seen = {min(remaining)}
        todo = list(seen)
        while todo:
            i = todo.pop()
            for j in remaining - seen:
                if g[i][j]:
                    seen.add(j)
                    todo.append(j)
        remaining -= seen
        result.append(sorted(seen))
    return result


def native(g, nodes):
    labels = [[0] * len(nodes) for _ in nodes]
    for a, i in enumerate(nodes):
        for b, j in enumerate(nodes[:a]):
            if g[i][j]:
                square = 4 * g[i][j]**2 / (g[i][i] * g[j][j])
                require(square in (1, 2), 'unsupported native bond in a paying component')
                labels[a][b] = labels[b][a] = 3 if square == 1 else 4
    return labels


def typename(name):
    require(type(name) is str and TYPE.fullmatch(name) is not None, 'paying finite type')
    return ('E', 6) if name == 'E6' else (name[0], int(name[1:]))


def charge(fuel):
    fuel[0] += 1
    require(fuel[0] <= MAX_VISITS, 'global traversal/isomorphism fuel exhausted')


def images(g, nodes, marks, typ, weights, fuel):
    """All full graph isomorphisms, including every edge and nonedge."""
    k, n = typ
    require(n == len(nodes) and 1 <= n <= MAX_NODES, 'component size/type')
    left = native(g, nodes)
    if k == 'E':
        right = [dict() for _ in range(6)]
        for i, j in ((0, 1), (0, 2), (2, 3), (0, 4), (4, 5)):
            right[i][j] = right[j][i] = 3
    else:
        right = weights.classical_graph(k, n)
    degree = [sum(bool(x) for x in row) for row in left]
    order = sorted(range(n), key=lambda i: (-degree[i], i))
    assignment = {}
    used = set()
    out = set()

    def visit(depth):
        charge(fuel)
        if depth == n:
            out.add(tuple(sorted(assignment[nodes.index(m)] for m in marks)))
            return
        i = order[depth]
        for j in range(n):
            if j in used or degree[i] != len(right[j]):
                continue
            if any(left[i][a] != right[j].get(b, 0) for a, b in assignment.items()):
                continue
            assignment[i] = j
            used.add(j)
            visit(depth + 1)
            used.remove(j)
            del assignment[i]

    visit(0)
    require(bool(out), 'claimed finite component has no full native isomorphism')
    return sorted(out)


def witness_type(g, nodes):
    """Only positive-guard witnesses are retained; unrecognized types are omitted."""
    if len(nodes) not in (1, 2, 3, 9, 10, 11, 15):
        return None
    try:
        labels = native(g, nodes)
    except ValueError:
        return None
    if any(x not in (0, 3) for row in labels for x in row):
        return None
    degree = [sum(bool(x) for x in row) for row in labels]
    if sum(degree) != 2 * (len(nodes) - 1):
        return None
    if max(degree) <= 2:
        return ('A', len(nodes))
    centers = [i for i, d in enumerate(degree) if d == 3]
    if len(nodes) == 15 and len(centers) == 1 and max(degree) == 3:
        if sum(degree[j] == 1 for j, label in enumerate(labels[centers[0]]) if label) == 2:
            return ('D', 15)
    return None


def supplied_images(value, n, arity):
    require(type(value) is list and bool(value), 'nonempty image list')
    rows = [tuple(indices(row, n, arity)) for row in value]
    require(all(tuple(sorted(row)) == row for row in rows), 'unordered image convention')
    require(rows == sorted(set(rows)), 'complete image list must be sorted and distinct')
    return rows


def bind(row, gg, ww, weights, fuel):
    g = matrix(gg)
    white = indices(ww, len(g))
    marks = indices(row['marks'], len(g), 2)
    require(not set(white) & set(marks), 'marks are outside white support')
    ids = sorted(set(white) | set(marks))
    parts = components(g, ids)
    require(type(row['roles']) is list and len(row['roles']) in (1, 2), 'marked roles schema')
    selected = [next(part for part in parts if mark in part) for mark in marks]
    shared = selected[0] == selected[1]
    require(len(row['roles']) == (1 if shared else 2), 'roles versus actual connectivity')
    roles = []
    for part in [selected[0]] if shared else selected:
        hits = [r for r in row['roles'] if r.get('nodes') == part]
        require(len(hits) == 1, 'role must be an entire displayed component')
        role = hits[0]
        typ = typename(role['type'])
        local_marks = marks if shared else [next(m for m in marks if m in part)]
        actual = images(g, part, local_marks, typ, weights, fuel)
        if shared:
            require(actual == supplied_images(row['images'], typ[1], 2), 'same-component pair binding')
        elif typ == ('E', 6):
            require(type(role.get('long_endpoint')) is bool, 'E6 role Boolean')
            require(role['long_endpoint'] == all(x[0] in weights.E6_LONG_ENDS for x in actual), 'E6 long-end role')
        else:
            require(actual == supplied_images(role['images'], typ[1], 1), 'cross-component node binding')
        roles.append((typ, actual))
    spectators = [witness_type(g, p) for p in parts if p not in selected]
    witnesses = [t for t in spectators if t is not None]
    forest = [r[0] for r in roles] + witnesses
    lower = []
    guard_possible = False
    for choice in product(*(r[1] for r in roles)):
        first = (0, choice[0][0])
        second = (0, choice[0][1]) if shared else (1, choice[1][0])
        value = weights.finite_corner_weight(forest, first, second)
        # Positive extra witnesses are ignored. The sole final negative guard
        # is relaxed to any rank-feasible three A15/D15 spectators; no padding.
        if shared and roles[0][0] == ('D', 13) and choice[0] in ((0, 11), (0, 12)):
            count = sum(t in (('A', 15), ('D', 15)) for t in witnesses)
            guard_possible = count + (58 - len(ids)) // 15 >= 3
            if guard_possible:
                value = min(value, F(1, 2))
        lower.append(value)
    require(len(ids) <= 58, 'terminal cannot lie in an ambient58 finite corner')
    require(len(set(lower)) == 1, 'literal payment is not invariant under component maps')
    chosen = fraction(row['payment'])
    require(chosen > 0 and lower[0] >= chosen, 'current literal bound below chosen terminal payment')
    return dict(marks=marks, chosen=str(chosen), current_lower=str(lower[0]),
                same_component=shared, types=[list(t) for t, _ in roles],
                known_witnesses=[list(t) for t in witnesses],
                possible_negative_D13_guard=guard_possible)


def terminals(value, path='', gram=None, white=None, gram_path=None,
              white_path=None, depth=0, fuel=None):
    require(fuel is not None and depth <= 64, 'bounded traversal context')
    charge(fuel)
    if isinstance(value, dict):
        if 'gram' in value:
            gram_path = path + '/gram'
        if 'white' in value:
            white_path = path + '/white'
        gram = value.get('gram', gram)
        white = value.get('white', white)
        if {'marks', 'payment', 'roles'} <= value.keys() and fraction(value['payment']) > 0:
            require(gram is not None and white is not None and gram_path and white_path,
                    'missing inherited physical context/provenance')
            yield path, value, gram, white, gram_path, white_path
        for key, child in value.items():
            use_extended = key in ('after', 'extended_corners')
            gg = value.get('extended_gram', gram) if use_extended else gram
            ww = value.get('extended_white', white) if use_extended else white
            gp = path + '/extended_gram' if use_extended and 'extended_gram' in value else gram_path
            wp = path + '/extended_white' if use_extended and 'extended_white' in value else white_path
            yield from terminals(child, path + '/' + key, gg, ww, gp, wp, depth + 1, fuel)
    elif isinstance(value, list):
        for i, child in enumerate(value):
            yield from terminals(child, path + '/' + str(i), gram, white,
                                 gram_path, white_path, depth + 1, fuel)


def pointer(doc, path):
    require(type(path) is str and path.startswith('/'), 'physical JSON pointer')
    value = doc
    for part in path[1:].split('/'):
        require(part and '~' not in part, 'unambiguous fixed-schema JSON pointer')
        value = value[int(part)] if isinstance(value, list) else value[part]
    return value


def decode(encoded):
    return json.loads(encoded, object_pairs_hook=unique,
                      parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON number')))


def lock():
    path = ROOT / 'formal/.proof-batches/exclusive.lock'
    require(path.is_file(), 'shared proof lock is missing')
    with path.open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            return
        fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
    raise ValueError('shared parent proof lock must already be held')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(CASES), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    lock()
    relative, expected, count = CASES[args.batch]
    inputs = dict(PINS, **{relative: expected})
    before = {}
    source = raw(str(Path(__file__).resolve().relative_to(ROOT)))
    for name, pin in inputs.items():
        data = raw(name)
        before[name] = sha256(data).hexdigest()
        require(before[name] == pin, 'frozen input identity: ' + name)
    library = raw(LIB)
    require(sha256(library).hexdigest() == PINS[LIB], 'literal library changed before loading')
    weights = types.ModuleType('finite57_literal_terminal_weights')
    exec(compile(library, str(ROOT / LIB), 'exec'), weights.__dict__)
    encoded = raw(relative)
    require(sha256(encoded).hexdigest() == expected, 'terminal input changed before parsing')
    doc = decode(encoded)
    reference_bytes = raw(REFERENCE)
    require(sha256(reference_bytes).hexdigest() == PINS[REFERENCE], 'reference changed before parsing')
    reference = decode(reference_bytes)
    packets = reference['interfaces']['D12_positive_outputs']
    require(type(packets) is list and len(packets) == 2, 'complete reference source domain')
    require({p['source'] for p in packets} == {Path(v[0]).stem for v in CASES.values()},
            'reference source keys')
    packet = next(p for p in packets if p['source'] == Path(relative).stem)
    require(type(packet['positive_outputs']) is list and len(packet['positive_outputs']) == count,
            'bounded expected positive path list')
    expected_paths = {}
    for row in packet['positive_outputs']:
        key = row['path']
        require(type(key) is str and key not in expected_paths, 'duplicate expected positive path')
        expected_paths[key] = row
    require(len(expected_paths) == count, 'expected complete positive path domain')
    result = []
    seen = set()
    fuel = [0]
    for path, row, gram, white, gram_path, white_path in terminals(doc, fuel=fuel):
        require(path not in seen, 'duplicate terminal path')
        require(path in expected_paths, 'unexpected positive terminal path')
        seen.add(path)
        expected_row = expected_paths[path]
        require(row['payment'] == expected_row['payment'] and row['reason'] == expected_row['reason']
                and [r['type'] for r in row['roles']] == expected_row['marked_types'],
                'positive-output provenance payload')
        require(pointer(doc, gram_path) == gram and pointer(doc, white_path) == white,
                'physical context source binding')
        result.append(dict(path=path, gram_path=gram_path, white_path=white_path,
                           **bind(row, gram, white, weights, fuel)))
    require(seen == set(expected_paths), 'complete exact positive-output path set')
    for name, pin in inputs.items():
        require(sha256(raw(name)).hexdigest() == pin, 'input changed during batch: ' + name)
    require(raw(str(Path(__file__).resolve().relative_to(ROOT))) == source, 'checker source changed during batch')
    print(json.dumps(dict(status='PASS_CURRENT_D12_TERMINAL_LITERAL_BINDING',
        batch=args.batch, positive_output_paths=count, inputs=before,
        source_sha256=sha256(source).hexdigest(), global_fuel_used=fuel[0], rows=result,
        scope='Current literal lower dominates each selected-D12 chosen finite payment. Whole-component seals, complete selected-domain maps, inertia and geometric funding remain explicit premises; no other terminal family is certified.'), sort_keys=True))


if __name__ == '__main__':
    main()
