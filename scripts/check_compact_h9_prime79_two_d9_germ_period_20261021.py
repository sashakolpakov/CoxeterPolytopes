#!/usr/bin/env python3
"""Prepared exact constant-size two-D9 period check; no root enumeration.

Run only through the unchanged bounded runner with --parent-lock. Actual
full-germ and arithmetic-covolume interpretations remain mathematical inputs.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROOF = 'results/compact_h9_prime79_two_actual_d9_germ_period_candidate_20261021.md'
CERT = 'results/prime79_actual_a9_end0_interval_20261021.jsonl'
LEDGER = 'results/prime79_actual_a9_end0_interval_execution_20261021.json'
PINS = {
    'results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md':
        '33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3',
    'results/compact_h9_prime79_actual_a9_end_edge_upper_and_complete_contract_candidate_20261021.md':
        '615d2915b1f95ffe5dc4d0ec634840b67991a7afb7b620aaca63d498bb10c2d3',
    'results/compact_h9_prime79_actual_a9_end_interval_checker_spec_20261021.md':
        'b21fe0b47693d0767b9b27283bdc7ba78e5e141857577b7ebd153aad7df94d4d',
    CERT: '7f35075ed06766a4d78bc328ac9f8b8a819e9a860ec8a32b92d13d2b123319f0',
    LEDGER: 'ff8d43f1dd228712db84147e2ce23e8813d286a23f4c4dbc00fd39ca016146a5',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def q(a=0, b=0):
    return F(a), F(b)


ZERO, ONE, PHI = q(), q(1), q(0, 1)


def add(x, y):
    return x[0]+y[0], x[1]+y[1]


def neg(x):
    return -x[0], -x[1]


def sub(x, y):
    return add(x, neg(y))


def mul(x, y):
    return x[0]*y[0]+x[1]*y[1], x[0]*y[1]+x[1]*y[0]+x[1]*y[1]


def inv(x):
    norm = x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
    require(norm != 0, 'zero field inverse')
    return (x[0]+x[1])/norm, -x[1]/norm


def sign(x):
    a, b = 2*x[0]+x[1], x[1]
    if not b:
        return (a > 0)-(a < 0)
    if a >= 0 and b > 0:
        return 1
    if a <= 0 and b < 0:
        return -1
    d = a*a-5*b*b
    return ((d > 0)-(d < 0)) * (1 if a > 0 else -1)


def total(items):
    value = ZERO
    for item in items:
        value = add(value, item)
    return value


def vadd(x, y):
    return [add(a, b) for a, b in zip(x, y)]


def scale(a, x):
    return [mul(a, b) for b in x]


def vsum(items):
    value = [ZERO]*10
    for item in items:
        value = vadd(value, item)
    return value


def transpose(a):
    return [list(row) for row in zip(*a)]


def matmul(a, b):
    bt = transpose(b)
    return [[total(mul(x, y) for x, y in zip(row, col)) for col in bt] for row in a]


def apply(a, v):
    return [total(mul(x, y) for x, y in zip(row, v)) for row in a]


def identity(n):
    return [[ONE if i == j else ZERO for j in range(n)] for i in range(n)]


def inverse(a):
    n = len(a)
    b = [list(row)+tail for row, tail in zip(a, identity(n))]
    for j in range(n):
        pivot = next((i for i in range(j, n) if b[i][j] != ZERO), None)
        require(pivot is not None, 'singular finite matrix')
        b[j], b[pivot] = b[pivot], b[j]
        b[j] = [mul(x, inv(b[j][j])) for x in b[j]]
        for i in range(n):
            if i != j:
                factor = b[i][j]
                b[i] = [sub(x, mul(factor, y)) for x, y in zip(b[i], b[j])]
    return [row[n:] for row in b]


def decode(v):
    require(len(v) == 10 and all(len(x) == 2 for x in v), 'bad original vector')
    return [q(*x) for x in v]


def integral(a):
    return all(x.denominator == 1 for row in a for value in row for x in value)


def encode(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, (tuple, list)):
        return [encode(x) for x in value]
    if isinstance(value, dict):
        return {key: encode(x) for key, x in value.items()}
    return value


def check_pins(pins):
    for name, expected in sorted(pins.items()):
        path = (ROOT/name).resolve()
        require(path.is_relative_to(ROOT), 'input outside repository')
        require(digest(path) == expected, 'changed input: '+name)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--expected-source-sha', required=True)
    parser.add_argument('--proof-sha', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    require(owner > 1 and os.getppid() == owner, 'requires live bounded owner')
    target = (ROOT/args.out).resolve()
    require(target.parent == ROOT/'results' and not target.exists(), 'new results output only')
    pins = dict(PINS)
    pins[str(Path(__file__).resolve().relative_to(ROOT))] = args.expected_source_sha
    pins[PROOF] = args.proof_sha
    check_pins(pins)
    require((ROOT/CERT).stat().st_size == 16733, 'changed END0 size')
    with (ROOT/CERT).open() as stream:
        records = [json.loads(line) for line in stream]
    require(len(records) == 9, 'changed finite record count')
    header, setup, result = records[0], records[1], records[-1]
    require(header['record'] == 'header' and setup['record'] == 'setup'
            and result['record'] == 'result', 'bad finite record structure')
    for name, value in header['pins'].items():
        require(name not in pins or pins[name] == value, 'conflicting nested input')
        pins[name] = value
    check_pins(pins)
    require(result['complete'] and result['all_six_profiles_complete']
            and result['unique_simple_germ'] and result['end'] == 0
            and result['profile_selection'] == 'all'
            and result['simple_mirror_ids'] == [4]
            and [row['id'] for row in result['minimum_ties']] == [4], 'not complete unique END0')
    require(sum(row['raw_trials'] for row in result['counts'].values()) == 810436,
            'incomplete inherited finite interval')
    ledger = json.loads((ROOT/LEDGER).read_text())
    stdout = json.loads(ledger['stdout'])
    require(ledger['status'] == 'PASS_BOUNDED_RESEARCH_CHECK'
            and ledger['wall_limit_seconds'] == 59 and ledger['rss_limit_bytes'] == 1073741824
            and ledger['threads'] == 1 and ledger['lock_owner'] == 'wrapper'
            and stdout['output_sha256'] == pins[CERT]
            and stdout['output_bytes'] == 16733 and stdout['raw_scalar_trials'] == 810436,
            'bad inherited execution boundary')

    basis = identity(10)
    bx, by = basis[:2]
    e = basis[2:]
    d, a = q(18, 31), q(3, 3)
    gram = [[ZERO]*10 for _ in range(10)]
    gram[0][0], gram[0][1], gram[1][0], gram[1][1] = q(0, -2), a, a, q(2)
    for i in range(2, 10):
        gram[i][i] = q(2)
    for i, j in [(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7)]:
        gram[i+2][j+2] = gram[j+2][i+2] = q(-1)

    def dot(x, y):
        return total(mul(u, v) for u, v in zip(x, apply(gram, y)))

    c = scale(q(-1), by)
    p1 = vadd(scale(q(2), bx), scale(neg(a), by))
    dual0 = [ZERO, ZERO]+[q(x) for x in (4, 5, 7, 10, 8, 6, 4, 2)]
    dual7 = [ZERO, ZERO]+[q(x) for x in (2, 3, 4, 6, 5, 4, 3, 2)]
    n = vsum([scale(q(0, F(1, 2)), p1), scale(q(F(-1, 2)), c),
              dual7, scale(q(1, 2), dual0)])
    known = {row['name']: decode(row['normal']) for row in setup['known_walls']}
    require(known['n'] == n and known['c'] == c, 'changed OLD actual roots')
    s, h = known['observed_s'], decode(next(row for row in records
                                         if row.get('record') == 'mirror' and row['id'] == 4)['normal'])
    require(h == decode(result['minimum_ties'][0]['normal']), 'changed minimum normal')
    point, omega = decode(setup['corner']), decode(setup['dual'])
    require(dot(point, point) == mul(q(-10), d) and dot(omega, omega) == q(F(9, 10))
            and dot(point, omega) == ZERO
            and sub(mul(q(11, 9), q(11, 9)), mul(q(9), d)) == q(40),
            'wrong NEW corner/dual normalizations')
    u = vadd(scale(q(3, 4), p1), scale(mul(q(F(1, 2)), mul(d, PHI)), dual0))
    v = scale(q(F(1, 20)), vadd(scale(q(11, 9), point), scale(mul(q(-10), d), omega)))
    old = [c, n]+[scale(q(-1), e[j]) for j in (7, 6, 5, 4, 3, 1, 2)]
    new = [scale(q(-1), e[j]) for j in (2, 3, 4, 5, 6, 7)]+[s, c, h]
    require(dot(u, u) == dot(v, v) == neg(d) and sign(dot(u, p1)) < 0
            and sign(dot(v, point)) < 0, 'wrong equal future normalization')
    for roots, corner in [(old, u), (new, v)]:
        require(all(dot(root, corner) == ZERO for root in roots), 'corner not perpendicular')
    old_gram = [[dot(x, y) for y in old] for x in old]
    d9 = [[q(2) if i == j else ZERO for j in range(9)] for i in range(9)]
    for i, j in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7), (6, 8)]:
        d9[i][j] = d9[j][i] = q(-1)
    require(old_gram == d9 == [[dot(x, y) for y in new] for x in new],
            'unequal ordered full D9 Grams')
    old_matrix = transpose(old+[u])
    old_inverse = inverse(old_matrix)
    cases = []
    for swap in (False, True):
        roots = list(new)
        if swap:
            roots[-2], roots[-1] = roots[-1], roots[-2]
        transform = matmul(transpose(roots+[v]), old_inverse)
        undo = inverse(transform)
        require(matmul(transpose(transform), matmul(gram, transform)) == gram,
                'not an ambient isometry')
        require(matmul(transform, undo) == matmul(undo, transform) == identity(10), 'bad inverse')
        require(apply(transform, n) == scale(q(-1), e[3])
                and apply(transform, apply(transform, n)) == s, 'wrong two-step root orbit')
        require(dot(n, n) == q(2) and dot(n, apply(transform, n)) == ZERO
                and dot(n, s) == q(-80, -132) and sign(add(dot(n, s), q(6))) < 0,
                'wrong finite-order obstruction')
        images = transpose(transform)
        residual = vsum([scale(q(2), v), scale(neg(mul(q(2), u[1])), images[1])]
                        +[scale(neg(mul(q(2), u[j+2])), images[j+2]) for j in range(1, 8)])
        require(all(x.denominator == 1 for value in residual for x in value), 'nonintegral residual')
        residues = [[int(x) % 4 for x in value] for value in residual]
        is_integral = integral(transform) and integral(undo)
        require(is_integral == (not swap), 'manual leaf-assignment theorem failed')
        require(all(value == [0, 0] for value in residues) if not swap
                else residues[2] == [0, 2], 'manual complete mod4 residue theorem failed')
        cases.append({'swap_target_leaves': swap, 'original_integral_and_inverse': is_integral,
                      'matrix': transform, 'inverse': undo, 'residual_mod4': residues,
                      'trace': total(transform[j][j] for j in range(10))})
    check_pins(pins)
    require(os.getppid() == owner, 'bounded owner exited')
    output = {'status': 'PASS_EXACT_TWO_D9_ORIGINAL_MODULE_PERIOD_ALGEBRA', 'pins': pins,
              'cases': cases, 'old_corner': u, 'new_corner': v,
              'new_root_enumeration': False, 'n_dot_T2n': q(-80, -132),
              'actual_germs_and_finite_covolume_are_external_inputs': True,
              'literal_nonreflectivity_requires_accepted_geometric_boundary': True,
              'no_new_example_dimension_or_form_class_exclusion': True}
    with target.open('x') as stream:
        json.dump(encode(output), stream, indent=2, sort_keys=True)
        stream.write('\n')
        stream.flush()
        os.fsync(stream.fileno())
    print(json.dumps({'status': output['status'], 'output': str(target.relative_to(ROOT)),
                      'output_sha256': digest(target), 'output_bytes': target.stat().st_size,
                      'integral_leaf_assignments': 1}, sort_keys=True))


if __name__ == '__main__':
    main()
