#!/usr/bin/env python3
"""Prepared exact necessary census: ONE isolated quadratic D17 fork.

UNEXECUTED source. Root-owned guarded execution only; no repository imports.
The pinned contract proves the complete36-field domain, full integer bases,
both signed-norm square branches, actual raw cycles and strict scope.
Three disjoint twelve-field batches and a separate streaming combine.
At most8 head and17 fork candidates per field,4896 pair tests in total.
No output is a polytope, entire-D17 exclusion or dimension theorem.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__).resolve()
PINS = {
    'results/arithmetic17_whole_d17_dotted_head_degree_three_and_nonisolated_quadratic_candidate_20261022.md':
        '4baa131a7609f1af7254e154fae99834f0ab38f9e7a31bc2a49a186b0916ceba',
    'results/arithmetic17_whole_d17_fresh_long_head_and_all_finite_head_exclusion_candidate_20261022.md':
        '9e2cf7d75f12c73b69e0bc9c39044acb56c37be36f94fe707efb808ee5ae9eb6',
    'results/arithmetic_uniform_classical_weak_helper_obstruction_20261021.md':
        '535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905',
    'results/arithmetic17_d17_one_isolated_quadratic_fork_exact_census_contract_20261022.md':
        '991c97fab25b01a5d7aadc5a2bb220e469525647fac331cf74c02665af01df21',
}
EXPECTED_DISCRIMINANTS = (
    5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40,
    41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77,
    85, 88, 89, 92, 93, 97, 101, 104, 105, 109, 113, 120,
)
FAILURES = (
    'time_not_field_square', 'g_nonintegral', 'L_nonintegral',
    'nonphysical_cycle_outside_open_interval', 'parallel_or_nonnative_product',
)
SCOPE = 'necessary_one_isolated_quadratic_dotted_fork_only'
MAX_JSONL_LINE = 65536


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def validate_pins(source_sha):
    require(len(source_sha) == 64 and all(c in '0123456789abcdef' for c in source_sha),
            'invalid expected source SHA256')
    require(digest(SOURCE) == source_sha, 'source SHA256 mismatch')
    for rel, expected in PINS.items():
        require(digest(ROOT / rel) == expected, 'input SHA256 mismatch: ' + rel)


def integer(value):
    return type(value) is int


@dataclass(frozen=True)
class Q:
    d: int
    a: F
    b: F = F(0)

    def __post_init__(self):
        object.__setattr__(self, 'a', F(self.a))
        object.__setattr__(self, 'b', F(self.b))

    def lift(self, other):
        if isinstance(other, Q):
            require(self.d == other.d, 'different quadratic fields')
            return other
        require(isinstance(other, (int, F)), 'non-exact scalar')
        return Q(self.d, F(other))

    def __add__(self, other):
        other = self.lift(other)
        return Q(self.d, self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(self.d, -self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        return Q(self.d, self.a * other.a + self.d * self.b * other.b,
                 self.a * other.b + self.b * other.a)

    __rmul__ = __mul__

    def conjugate(self):
        return Q(self.d, self.a, -self.b)

    def norm(self):
        return self.a * self.a - self.d * self.b * self.b

    def __truediv__(self, other):
        other = self.lift(other)
        n = other.norm()
        require(n != 0, 'division by zero')
        p = self * other.conjugate()
        return Q(self.d, p.a / n, p.b / n)

    def sign(self):
        if self.b == 0:
            return (self.a > 0) - (self.a < 0)
        if self.a == 0:
            return (self.b > 0) - (self.b < 0)
        if self.a * self.b > 0:
            return (self.a > 0) - (self.a < 0)
        n = self.norm()
        require(n != 0, 'nonzero quadratic pair cannot have zero norm')
        return ((self.a > 0) - (self.a < 0)) * ((n > 0) - (n < 0))

    def integral(self):
        if self.d % 4 == 1:
            return (2 * self.b).denominator == 1 and (self.a - self.b).denominator == 1
        return self.a.denominator == 1 and self.b.denominator == 1

    def encode(self):
        return {'d': self.d,
                'a': [self.a.numerator, self.a.denominator],
                'b': [self.b.numerator, self.b.denominator]}


def rational_sqrt(value):
    value = F(value)
    if value < 0:
        return None
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    return F(n, d) if n * n == value.numerator and d * d == value.denominator else None


def field_sqrt(value):
    """Return the physical NONNEGATIVE K-root, or None; all branches exact."""
    if value.b == 0:
        r = rational_sqrt(value.a)
        if r is not None:
            return Q(value.d, r)
        r = rational_sqrt(value.a / value.d)
        if r is not None:
            root = Q(value.d, 0, r)
            require(root * root == value, 'pure-radical square identity')
            return root
        return None
    n = rational_sqrt(value.norm())
    if n is None:
        return None
    for signed_norm in (n, -n):
        x = rational_sqrt((value.a + signed_norm) / 2)
        if x is not None and x != 0:
            root = Q(value.d, x, value.b / (2 * x))
            if root * root == value:
                if root.sign() < 0:
                    root = -root
                require(root.sign() > 0, 'mixed square root sign')
                return root
    return None


def arithmetic_fixtures():
    """Prepared exact fixtures; execution belongs to each eventual guarded unit."""
    require(Q(2, 2, -1).sign() > 0 and Q(2, 1, -1).sign() < 0, 'sign cancellation')
    require(Q(2, -2, 1).sign() < 0 and Q(2, -1, 1).sign() > 0, 'opposite sign cases')
    require(field_sqrt(Q(2, 0)) == Q(2, 0), 'zero square')
    require(field_sqrt(Q(2, 4)) == Q(2, 2), 'rational square')
    require(field_sqrt(Q(2, 2)) == Q(2, 0, 1), 'pure radical square')
    require(field_sqrt(Q(2, F(4, 3))) is None, 'rational nonsquare')
    require(field_sqrt(Q(2, -1)) is None, 'negative square')
    # The second root has NEGATIVE norm; testing only +sqrt(norm(value)) fails.
    for root in (Q(2, 2, 1), Q(2, 1, 1), Q(2, -1, 1)):
        found = field_sqrt(root * root)
        require(found is not None and found * found == root * root and found.sign() > 0,
                'both signed-norm mixed branches')
    require(Q(5, F(1, 2), F(1, 2)).integral(), 'full half-integral basis')
    require(not Q(5, F(1, 2), 0).integral(), 'nonintegral half scalar')
    require(not Q(2, F(3, 2), 1).integral(), 'sqrt basis nonintegrality')
    return ['exact_signs', 'zero', 'rational_square', 'pure_radical',
            'rational_nonsquare', 'negative', 'both_mixed_norm_signs', 'full_integer_bases']


def squarefree(d):
    return all(d % (p * p) != 0 for p in range(2, isqrt(d) + 1))


def fields():
    out = []
    for d in range(2, 121):
        if squarefree(d):
            D = d if d % 4 == 1 else 4 * d
            if D < 121:
                out.append((D, d))
    out.sort()
    require(tuple(D for D, _ in out) == EXPECTED_DISCRIMINANTS, 'complete field domain')
    return out


def candidates(D, d, which):
    require(D == (d if d % 4 == 1 else 4 * d), 'fundamental discriminant')
    half = d % 4 == 1
    omega = Q(d, F(1, 2), F(1, 2)) if half else Q(d, 0, 1)
    answer = []
    for b in range(1, 5 if which == 'U' else 18):
        if which == 'U' and b * b * D >= 121:
            continue
        if which == 'W' and 289 * b * b * D >= 676 * 676:
            continue
        floor = (isqrt(b * b * d) - b) // 2 if half else isqrt(b * b * d)
        for a in (floor + 1, floor + 2) if which == 'U' else (floor + 1,):
            value = a + b * omega
            conjugate = value.conjugate()
            low, high, conjugate_high = ((4, 11, 2) if which == 'U'
                                         else (2, F(676, 17), F(4, 17)))
            if ((value - low).sign() > 0 and (high - value).sign() > 0
                    and conjugate.sign() > 0 and (conjugate_high - conjugate).sign() > 0):
                require(value.integral() and value.b > 0, 'full-ring generator')
                answer.append((a, b, value))
    require(len(answer) <= (8 if which == 'U' else 17), 'static candidate bound')
    require(len(set(answer)) == len(answer), 'duplicate coefficient candidate')
    return answer


def product_class(L):
    if (L - 4).sign() > 0:
        return 'dotted'
    if (L - 4).sign() == 0:
        return None
    alphabet = [(Q(L.d, 1), 'native3'), (Q(L.d, 2), 'native4'), (Q(L.d, 3), 'native6')]
    if L.d == 5:
        alphabet.extend([(Q(5, F(3, 2), F(1, 2)), 'native5'),
                         (Q(5, F(5, 2), F(1, 2)), 'native10')])
    if L.d == 2:
        alphabet.append((Q(2, 2, 1), 'native8'))
    if L.d == 3:
        alphabet.append((Q(3, 2, 1), 'native12'))
    return next((label for value, label in alphabet if L == value), None)


def survivor_record(batch_index, D, iu, iw, uitem, witem, root):
    ua, ub, U = uitem
    wa, wb, W = witem
    H0 = U * (U - 2)
    Hf = W * (17 * W - 4)
    R = H0 * Hf
    require(root * root == R and root.sign() > 0, 'positive exact time witness')
    g = U * W - root
    L = g * g / (2 * U * W)
    label = product_class(L)
    require(g.sign() < 0 and g.integral() and L.integral(), 'actual raw cycle witnesses')
    require(L.sign() > 0 and L.conjugate().sign() > 0
            and (4 - L.conjugate()).sign() > 0 and label is not None, 'actual product bounds')
    return {'kind': 'survivor', 'batch_index': batch_index, 'D': D,
            'index_u': iu, 'index_w': iw, 'U_coefficients': [ua, ub], 'W_coefficients': [wa, wb],
            'U': U.encode(), 'W': W.encode(), 'head_time': H0.encode(), 'fork_time': Hf.encode(),
            'R': R.encode(), 'positive_sqrt_R': root.encode(), 'g': g.encode(), 'L': L.encode(),
            'classification': label}


def pair_test(batch_index, D, iu, iw, uitem, witem):
    U, W = uitem[2], witem[2]
    R = U * (U - 2) * W * (17 * W - 4)
    require(R.sign() > 0 and R.conjugate().sign() > 0, 'positive time radicand at both places')
    root = field_sqrt(R)
    if root is None:
        return 'time_not_field_square', None
    g = U * W - root
    require(g.sign() < 0, 'domain forces strict head/fork acuteness')
    if not g.integral():
        return 'g_nonintegral', None
    L = g * g / (2 * U * W)
    if not L.integral():
        return 'L_nonintegral', None
    if L.conjugate().sign() <= 0 or (4 - L.conjugate()).sign() <= 0:
        return 'nonphysical_cycle_outside_open_interval', None
    if product_class(L) is None:
        return 'parallel_or_nonnative_product', None
    return None, survivor_record(batch_index, D, iu, iw, uitem, witem, root)


def emit(stream, value):
    stream.write(json.dumps(value, sort_keys=True, separators=(',', ':')) + '\n')
    stream.flush()


def output_file(argument):
    path = (ROOT / argument).resolve()
    require((ROOT / 'results').resolve() in path.parents, 'output must be under results/')
    require(path.parent.is_dir(), 'output parent must already exist')
    return path.open('x', encoding='utf-8')


def blank_totals():
    return {'fields': 0, 'U_total': 0, 'W_total': 0, 'pair_tests': 0, 'survivors': 0,
            'rejects': {key: 0 for key in FAILURES}}


def add_field(totals, record):
    for target, source in (('fields', None), ('U_total', 'U_count'), ('W_total', 'W_count'),
                           ('pair_tests', 'pair_tests'), ('survivors', 'survivors')):
        totals[target] += 1 if source is None else record[source]
    for key in FAILURES:
        totals['rejects'][key] += record['rejects'][key]


def batch_header(index, source_sha):
    selected = fields()[12 * index:12 * (index + 1)]
    return {'kind': 'header', 'schema': 1, 'unit': 'batch', 'batch_index': index,
            'discriminants': [D for D, _ in selected], 'source_sha256': source_sha,
            'fixed_pins': PINS, 'scope': SCOPE}


def batch(args, fixture_names):
    selected = fields()[12 * args.index:12 * (args.index + 1)]
    totals = blank_totals()
    with output_file(args.output) as stream:
        emit(stream, batch_header(args.index, args.expected_source_sha256))
        for D, d in selected:
            us, ws = candidates(D, d, 'U'), candidates(D, d, 'W')
            rejects = {key: 0 for key in FAILURES}
            surviving = 0
            for iu, uitem in enumerate(us):
                for iw, witem in enumerate(ws):
                    failure, record = pair_test(args.index, D, iu, iw, uitem, witem)
                    if failure is None:
                        emit(stream, record)
                        surviving += 1
                    else:
                        rejects[failure] += 1
            record = {'kind': 'field', 'batch_index': args.index, 'D': D, 'd': d,
                      'U_count': len(us), 'W_count': len(ws), 'pair_tests': len(us) * len(ws),
                      'survivors': surviving, 'rejects': rejects}
            require(sum(rejects.values()) + surviving == record['pair_tests'], 'pair partition')
            emit(stream, record)
            add_field(totals, record)
        validate_pins(args.expected_source_sha256)
        summary = {'kind': 'summary', 'unit': 'batch', 'batch_index': args.index,
                   'status': 'PASS', 'scope': SCOPE, 'arithmetic_fixtures': fixture_names, **totals}
        emit(stream, summary)
    print(json.dumps(summary, sort_keys=True))


def decode_q(value, d):
    require(type(value) is dict and set(value) == {'d', 'a', 'b'} and value['d'] == d,
            'quadratic witness schema')
    parts = []
    for key in ('a', 'b'):
        pair = value[key]
        require(type(pair) is list and len(pair) == 2 and all(integer(x) for x in pair)
                and pair[1] > 0, 'rational witness schema')
        q = F(pair[0], pair[1])
        require([q.numerator, q.denominator] == pair, 'noncanonical rational witness')
        parts.append(q)
    return Q(d, *parts)


def read_record(stream):
    line = stream.readline(MAX_JSONL_LINE + 1)
    require(line != '', 'truncated batch certificate')
    require(len(line) <= MAX_JSONL_LINE and line.endswith('\n'), 'oversized/truncated JSONL record')
    record = json.loads(line)
    require(type(record) is dict and type(record.get('kind')) is str, 'JSONL record schema')
    return record


def check_field_record(record, batch_index, D, d, us, ws, surviving):
    require(set(record) == {'kind', 'batch_index', 'D', 'd', 'U_count', 'W_count',
                            'pair_tests', 'survivors', 'rejects'}, 'field record schema')
    require(record['kind'] == 'field' and record['batch_index'] == batch_index
            and record['D'] == D and record['d'] == d, 'field coverage/order')
    for key in ('U_count', 'W_count', 'pair_tests', 'survivors'):
        require(integer(record[key]) and record[key] >= 0, 'field counter type')
    require(record['U_count'] == len(us) and record['W_count'] == len(ws)
            and record['pair_tests'] == len(us) * len(ws) and record['survivors'] == surviving,
            'field domain/survivor counts')
    reject = record['rejects']
    require(type(reject) is dict and set(reject) == set(FAILURES)
            and all(integer(v) and v >= 0 for v in reject.values()), 'rejection counter schema')
    require(sum(reject.values()) + surviving == record['pair_tests'], 'field pair partition')


def combine(args, fixture_names):
    require(len(args.part) == len(args.part_sha256) == 3, 'combine needs exactly three bound parts')
    paths = [(ROOT / p).resolve() for p in args.part]
    require(len(set(paths)) == 3, 'duplicate part path')
    all_fields = fields()
    total, seen = blank_totals(), set()
    bound = []
    for path, expected in zip(paths, args.part_sha256):
        require(len(expected) == 64 and digest(path) == expected, 'part SHA256 mismatch')
        bound.append({'path': str(path.relative_to(ROOT)), 'sha256': expected,
                      'bytes': path.stat().st_size})
    with output_file(args.output) as output:
        emit(output, {'kind': 'header', 'schema': 1, 'unit': 'combine', 'scope': SCOPE,
                      'source_sha256': args.expected_source_sha256, 'fixed_pins': PINS,
                      'discriminants': list(EXPECTED_DISCRIMINANTS), 'parts': bound})
        for path, part_binding in zip(paths, bound):
            subtotal = blank_totals()
            with path.open(encoding='utf-8') as stream:
                header = read_record(stream)
                index = header.get('batch_index')
                require(integer(index) and index in (0, 1, 2) and index not in seen,
                        'duplicate/invalid batch owner')
                require(header == batch_header(index, args.expected_source_sha256), 'batch header/pins')
                seen.add(index)
                emit(output, {'kind': 'part_binding', 'batch_index': index, **part_binding})
                for D, d in all_fields[12 * index:12 * (index + 1)]:
                    us, ws = candidates(D, d, 'U'), candidates(D, d, 'W')
                    surviving, last = 0, (-1, -1)
                    record = read_record(stream)
                    while record['kind'] == 'survivor':
                        require(record.get('batch_index') == index and record.get('D') == D,
                                'survivor field/owner')
                        iu, iw = record.get('index_u'), record.get('index_w')
                        require(integer(iu) and integer(iw) and 0 <= iu < len(us)
                                and 0 <= iw < len(ws) and (iu, iw) > last,
                                'survivor domain/order/uniqueness')
                        root = decode_q(record.get('positive_sqrt_R'), d)
                        expected = survivor_record(index, D, iu, iw, us[iu], ws[iw], root)
                        require(record == expected, 'exact survivor algebra/witness mismatch')
                        emit(output, record)
                        surviving += 1
                        last = (iu, iw)
                        record = read_record(stream)
                    check_field_record(record, index, D, d, us, ws, surviving)
                    emit(output, record)
                    add_field(subtotal, record)
                    add_field(total, record)
                summary = read_record(stream)
                expected_summary = {'kind': 'summary', 'unit': 'batch', 'batch_index': index,
                                    'status': 'PASS', 'scope': SCOPE,
                                    'arithmetic_fixtures': fixture_names, **subtotal}
                require(summary == expected_summary and subtotal['fields'] == 12, 'terminal batch coverage')
                require(stream.read(1) == '', 'extra records after terminal summary')
            require(digest(path) == part_binding['sha256'], 'part changed during combine')
        require(seen == {0, 1, 2} and total['fields'] == 36 and total['pair_tests'] <= 4896,
                'complete disjoint three-batch domain')
        validate_pins(args.expected_source_sha256)
        summary = {'kind': 'summary', 'unit': 'combine', 'status': 'PASS', 'scope': SCOPE,
                   'arithmetic_fixtures': fixture_names, 'batches': 3, **total}
        emit(output, summary)
    print(json.dumps(summary, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    sub = parser.add_subparsers(dest='unit', required=True)
    one = sub.add_parser('batch')
    one.add_argument('--index', type=int, choices=(0, 1, 2), required=True)
    one.add_argument('--output', required=True)
    union = sub.add_parser('combine')
    union.add_argument('--part', action='append', required=True)
    union.add_argument('--part-sha256', action='append', required=True)
    union.add_argument('--output', required=True)
    args = parser.parse_args()
    validate_pins(args.expected_source_sha256)
    fixture_names = arithmetic_fixtures()
    if args.unit == 'batch':
        batch(args, fixture_names)
    else:
        combine(args, fixture_names)


if __name__ == '__main__':
    main()
