#!/usr/bin/env python3
"""Prepared CONDITIONAL necessary census: TWO isolated quadratic A17 ends.

UNEXECUTED. Root-owned guarded jobs only; no repository imports.
The pinned endpoint geometric reduction is NOT accepted by this preparation.
Full25-field integral-basis domain; three9/8/8-field units and a combine.
At most12 models/field and3600 ordered pairs. Orthogonal pairs are allowed.
No output is an example, whole-A17 exclusion or dimension theorem.
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
    'results/arithmetic17_whole_a17_short_endpoints_isolation_and_golden_exclusion_candidate_20261022.md':
        'b80c27b00dc7ecdbdf7201327d620e10b372e35f1047964c9ae204c2d6996558',
    'results/arithmetic_uniform_classical_weak_helper_obstruction_20261021.md':
        '535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905',
    'results/arithmetic17_a17_two_isolated_quadratic_endpoints_exact_census_contract_20261022.md':
        'ba2c8f72de389bea30b8a766e62b7ad11957f72a10c0c178183dd5820e69f220',
    'scripts/check_d17_one_isolated_quadratic_fork_20261022.py':
        '73c44f9435269b21c859acf8b35bfb85927c43f5d4daf4e37cb18ab86d81d5ed',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
EXPECTED_DISCRIMINANTS = (
    5, 8, 12, 13, 17, 21, 24, 28, 29, 33, 37, 40,
    41, 44, 53, 56, 57, 60, 61, 65, 69, 73, 76, 77, 85,
)
BATCH_RANGES = ((0, 9), (9, 17), (17, 25))
LOW = F(36, 17)
HIGH = F(3969, 425)
FAILURES = (
    'endpoint_own_not_native_or_dotted', 'time_not_field_square',
    'physical_g_positive', 'g_nonintegral', 'L_nonintegral',
    'nonzero_cycle_outside_other_open_interval', 'mutual_not_native_or_dotted',
)
SCOPE = 'conditional_two_isolated_quadratic_A17_endpoints_only'
MAX_LINE = 65536
MAX_FILE = 16 * 1024 * 1024


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def valid_sha(value):
    return (type(value) is str and len(value) == 64
            and all(c in '0123456789abcdef' for c in value))


def validate_pins(source_sha):
    require(valid_sha(source_sha), 'invalid expected source SHA256')
    require(digest(SOURCE) == source_sha, 'source SHA256 mismatch')
    for rel, expected in PINS.items():
        require(valid_sha(expected) and digest(ROOT / rel) == expected,
                'input SHA256 mismatch: ' + rel)


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
        require(n != 0, 'nonzero irrational pair cannot have zero norm')
        return ((self.a > 0) - (self.a < 0)) * ((n > 0) - (n < 0))

    def integral(self):
        if self.d % 4 == 1:
            return ((2 * self.b).denominator == 1
                    and (self.a - self.b).denominator == 1)
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
    """Complete K-square test, returning the PHYSICAL nonnegative root."""
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


def product_class(value):
    """Native finite value (including orthogonal0), dotted, or inadmissible."""
    if value.sign() < 0:
        return None
    if value.sign() == 0:
        return 'native2'
    if (value - 4).sign() > 0:
        return 'dotted'
    if (value - 4).sign() == 0:
        return None
    alphabet = [(Q(value.d, 1), 'native3'), (Q(value.d, 2), 'native4'),
                (Q(value.d, 3), 'native6')]
    if value.d == 5:
        alphabet.extend([(Q(5, F(3, 2), F(1, 2)), 'native5'),
                         (Q(5, F(5, 2), F(1, 2)), 'native10')])
    if value.d == 2:
        alphabet.append((Q(2, 2, 1), 'native8'))
    if value.d == 3:
        alphabet.append((Q(3, 2, 1), 'native12'))
    return next((label for item, label in alphabet if value == item), None)


def other_cycle_admissible(g, L):
    if g.sign() == 0:
        require(L.sign() == 0, 'zero g has zero cycle')
        return True
    return L.conjugate().sign() > 0 and (4 - L.conjugate()).sign() > 0


def arithmetic_fixtures():
    """Prepared fixtures only; no execution credit until a root-owned unit."""
    require(Q(2, 2, -1).sign() > 0 and Q(2, 1, -1).sign() < 0,
            'mixed signs')
    require(Q(2, -2, 1).sign() < 0 and Q(2, -1, 1).sign() > 0,
            'opposite mixed signs')
    require(field_sqrt(Q(2, 0)) == Q(2, 0), 'zero square')
    require(field_sqrt(Q(2, 4)) == Q(2, 2), 'rational square')
    require(field_sqrt(Q(2, 2)) == Q(2, 0, 1), 'pure-radical square')
    require(field_sqrt(Q(2, F(4, 3))) is None, 'rational nonsquare')
    require(field_sqrt(Q(2, -1)) is None, 'negative square')
    for root in (Q(2, 2, 1), Q(2, 1, 1), Q(2, -1, 1)):
        found = field_sqrt(root * root)
        require(found is not None and found * found == root * root
                and found.sign() > 0, 'both mixed norm-sign branches')
    require(Q(5, F(1, 2), F(1, 2)).integral(), 'full half basis')
    require(not Q(5, F(1, 2), 0).integral(), 'nonintegral half scalar')
    require(not Q(2, F(3, 2), 1).integral(), 'nonintegral sqrt coefficient')
    zero = Q(2, 0)
    require(product_class(zero) == 'native2'
            and other_cycle_admissible(zero, zero), 'orthogonal pair retained')
    require(product_class(Q(2, 4)) is None, 'compact parallel excluded')
    omega77 = Q(77, F(1, 2), F(1, 2))
    # At b1 the conjugate strip contains THREE a values, not merely two.
    require(all((a + omega77.conjugate()).sign() > 0
                and (LOW - a - omega77.conjugate()).sign() > 0
                for a in (4, 5, 6)), 'three-a conjugate strip')
    return ['exact_signs', 'zero', 'rational_square', 'pure_radical',
            'rational_nonsquare', 'negative', 'both_mixed_norm_signs',
            'full_integer_bases', 'orthogonal_pair', 'compact_parallel',
            'three_a_conjugate_strip']


def squarefree(d):
    return all(d % (p * p) != 0 for p in range(2, isqrt(d) + 1))


def fields():
    result = []
    for d in range(2, 88):
        if squarefree(d):
            D = d if d % 4 == 1 else 4 * d
            if D < 88:
                result.append((D, d))
    result.sort()
    require(tuple(D for D, _ in result) == EXPECTED_DISCRIMINANTS,
            'complete field discriminant domain')
    return result


def candidates(D, d):
    require(D == (d if d % 4 == 1 else 4 * d), 'fundamental discriminant')
    half = d % 4 == 1
    omega = Q(d, F(1, 2), F(1, 2)) if half else Q(d, 0, 1)
    answer = []
    for b in range(1, 5):
        if b * b * D >= HIGH * HIGH:
            continue
        floor = (isqrt(b * b * d) - b) // 2 if half else isqrt(b * b * d)
        for a in (floor + 1, floor + 2, floor + 3):
            U = a + b * omega
            other = U.conjugate()
            if ((U - LOW).sign() > 0 and (HIGH - U).sign() > 0
                    and other.sign() > 0 and (LOW - other).sign() > 0):
                require(U.integral() and U.b > 0, 'full-ring primitive generator')
                answer.append((a, b, U))
    require(len(answer) <= 12 and len(set(answer)) == len(answer),
            'finite model bound/uniqueness')
    return answer


def model_record(item):
    a, b, U = item
    return {'coefficients': [a, b], 'U': U.encode(),
            'own_class': product_class(U)}


def pair_record(index, D, i, j, first, second):
    U1, U2 = first[2], second[2]
    record = {'kind': 'pair', 'batch_index': index, 'D': D,
              'index_first': i, 'index_second': j}
    own1, own2 = product_class(U1), product_class(U2)
    if own1 is None or own2 is None:
        return {**record, 'failure': 'endpoint_own_not_native_or_dotted'}
    F1, F2 = U1 * (17 * U1 - 36), U2 * (17 * U2 - 36)
    R = F1 * F2
    require(R.sign() > 0 and R.conjugate().sign() > 0,
            'positive radicand at both places')
    root = field_sqrt(R)
    if root is None:
        return {**record, 'failure': 'time_not_field_square'}
    g = (U1 * U2 - root) / 18
    if g.sign() > 0:
        return {**record, 'failure': 'physical_g_positive'}
    if not g.integral():
        return {**record, 'failure': 'g_nonintegral'}
    L = g * g / (U1 * U2)
    if not L.integral():
        return {**record, 'failure': 'L_nonintegral'}
    if not other_cycle_admissible(g, L):
        return {**record, 'failure': 'nonzero_cycle_outside_other_open_interval'}
    label = product_class(L)
    if label is None:
        return {**record, 'failure': 'mutual_not_native_or_dotted'}
    return {**record, 'failure': None, 'U_first': U1.encode(),
            'U_second': U2.encode(), 'F_first': F1.encode(),
            'F_second': F2.encode(), 'R': R.encode(),
            'positive_sqrt_R': root.encode(), 'g': g.encode(), 'L': L.encode(),
            'first_own_class': own1, 'second_own_class': own2,
            'mutual_class': label, 'orthogonal': g.sign() == 0}


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


class Writer:
    def __init__(self, stream):
        self.stream = stream
        self.bytes = 0

    def emit(self, record):
        line = (canonical(record) + '\n').encode('utf-8')
        require(len(line) <= MAX_LINE, 'JSONL record byte limit')
        require(self.bytes + len(line) <= MAX_FILE, 'certificate byte limit')
        self.stream.write(line)
        self.stream.flush()
        self.bytes += len(line)


def result_path(argument):
    path = (ROOT / argument).resolve()
    require((ROOT / 'results').resolve() in path.parents,
            'certificate path must be under results/')
    require(path.parent.is_dir(), 'certificate parent must already exist')
    return path


def output_file(argument):
    return result_path(argument).open('xb')


def blank_totals():
    return {'fields': 0, 'models': 0, 'pair_tests': 0, 'survivors': 0,
            'orthogonal_survivors': 0, 'rejects': {key: 0 for key in FAILURES}}


def add_field(total, record):
    total['fields'] += 1
    total['models'] += record['models']
    for key in ('pair_tests', 'survivors', 'orthogonal_survivors'):
        total[key] += record[key]
    for key in FAILURES:
        total['rejects'][key] += record['rejects'][key]


def selected_fields(index):
    low, high = BATCH_RANGES[index]
    return fields()[low:high]


def batch_header(index, source_sha):
    return {'kind': 'header', 'schema': 1, 'unit': 'batch', 'batch_index': index,
            'discriminants': [D for D, _ in selected_fields(index)],
            'source_sha256': source_sha, 'fixed_pins': PINS, 'scope': SCOPE}


def field_header(index, D, d, models):
    return {'kind': 'models', 'batch_index': index, 'D': D, 'd': d,
            'models': [model_record(item) for item in models]}


def field_records(index, D, d, models):
    rejects = {key: 0 for key in FAILURES}
    survivors = orthogonal = 0
    yield field_header(index, D, d, models)
    for i, first in enumerate(models):
        for j, second in enumerate(models):
            record = pair_record(index, D, i, j, first, second)
            failure = record['failure']
            if failure is None:
                survivors += 1
                orthogonal += int(record['orthogonal'])
            else:
                require(failure in rejects, 'unknown failure stage')
                rejects[failure] += 1
            yield record
    require(sum(rejects.values()) + survivors == len(models) ** 2,
            'complete first-failure partition')
    yield {'kind': 'field_summary', 'batch_index': index, 'D': D, 'd': d,
           'models': len(models), 'pair_tests': len(models) ** 2,
           'survivors': survivors, 'orthogonal_survivors': orthogonal,
           'rejects': rejects}


def batch_summary(index, fixtures, totals):
    return {'kind': 'summary', 'unit': 'batch', 'batch_index': index,
            'status': 'PASS', 'scope': SCOPE, 'arithmetic_fixtures': fixtures,
            **totals}


def batch(args, fixtures):
    totals = blank_totals()
    with output_file(args.output) as stream:
        writer = Writer(stream)
        writer.emit(batch_header(args.index, args.expected_source_sha256))
        for D, d in selected_fields(args.index):
            models = candidates(D, d)
            for record in field_records(args.index, D, d, models):
                writer.emit(record)
                if record['kind'] == 'field_summary':
                    add_field(totals, record)
        validate_pins(args.expected_source_sha256)
        summary = batch_summary(args.index, fixtures, totals)
        writer.emit(summary)
    print(canonical(summary))


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON object key')
        result[key] = value
    return result


def invalid_json_number(value):
    raise ValueError('noninteger/nonfinite JSON number is forbidden: ' + value)


def read_record(stream):
    line = stream.readline(MAX_LINE + 1)
    require(line and len(line) <= MAX_LINE and line.endswith(b'\n'),
            'missing/oversized/truncated JSONL record')
    record = json.loads(line, object_pairs_hook=unique_object,
                        parse_float=invalid_json_number, parse_constant=invalid_json_number)
    require(type(record) is dict and type(record.get('kind')) is str,
            'JSONL record schema')
    return record


def match_record(actual, expected):
    # Canonical equality distinguishes JSON true/false from integer1/0.
    require(canonical(actual) == canonical(expected), 'exact record/coverage mismatch')


def combine(args, fixtures):
    require(len(args.part) == len(args.part_sha256) == 3,
            'combine requires exactly three bound parts')
    paths = [result_path(p) for p in args.part]
    require(len(set(paths)) == 3, 'duplicate part path')
    bound = []
    for path, expected in zip(paths, args.part_sha256):
        require(valid_sha(expected) and 0 < path.stat().st_size <= MAX_FILE,
                'invalid part hash/byte size')
        require(digest(path) == expected, 'part SHA256 mismatch')
        bound.append({'path': str(path.relative_to(ROOT)), 'sha256': expected,
                      'bytes': path.stat().st_size})
    total, seen = blank_totals(), set()
    with output_file(args.output) as stream:
        writer = Writer(stream)
        writer.emit({'kind': 'header', 'schema': 1, 'unit': 'combine', 'scope': SCOPE,
                     'source_sha256': args.expected_source_sha256, 'fixed_pins': PINS,
                     'discriminants': list(EXPECTED_DISCRIMINANTS), 'parts': bound})
        for path, binding in zip(paths, bound):
            subtotal = blank_totals()
            with path.open('rb') as part:
                header = read_record(part)
                index = header.get('batch_index')
                require(type(index) is int and index in (0, 1, 2) and index not in seen,
                        'duplicate/invalid batch owner')
                match_record(header, batch_header(index, args.expected_source_sha256))
                seen.add(index)
                writer.emit({'kind': 'part_binding', 'batch_index': index, **binding})
                for D, d in selected_fields(index):
                    models = candidates(D, d)
                    for expected in field_records(index, D, d, models):
                        match_record(read_record(part), expected)
                        writer.emit(expected)
                        if expected['kind'] == 'field_summary':
                            add_field(subtotal, expected)
                            add_field(total, expected)
                summary = batch_summary(index, fixtures, subtotal)
                match_record(read_record(part), summary)
                require(part.read(1) == b'', 'extra bytes after terminal batch summary')
            require(digest(path) == binding['sha256']
                    and path.stat().st_size == binding['bytes'], 'part changed during combine')
        require(seen == {0, 1, 2} and total['fields'] == 25
                and total['models'] <= 300 and total['pair_tests'] <= 3600,
                'complete25-field/three-batch conjunction')
        validate_pins(args.expected_source_sha256)
        summary = {'kind': 'summary', 'unit': 'combine', 'status': 'PASS', 'scope': SCOPE,
                   'arithmetic_fixtures': fixtures, 'batches': 3, **total}
        writer.emit(summary)
    print(canonical(summary))


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
    fixtures = arithmetic_fixtures()
    if args.unit == 'batch':
        batch(args, fixtures)
    else:
        combine(args, fixtures)


if __name__ == '__main__':
    main()
