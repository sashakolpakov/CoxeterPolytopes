#!/usr/bin/env python3
"""Prepared CONDITIONAL A17 finite-cubic-anchor / other-endpoint census.

UNEXECUTED. Root-owned guarded jobs only, no repository imports.
Full original O_K=Z[anchor] basis;2432 coefficient triples in TWO fields.
Complete bounded cubic K-square recovery, not a norm-only test.
Two separate units and a separate every-record replay/combine.
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
    'results/arithmetic17_a17_finite_cubic_endpoint_exact_census_contract_20261022.md':
        '74b213461d250df03e601921744fa2831ea404d816e612562402dae3b8a0ec21',
    'scripts/check_a17_two_isolated_quadratic_endpoints_20261022.py':
        '89139eb11d31259c49d079a64bcef39d80e0cf855fe1e22221fd983f0d335a4d',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
FIELDS = {
    'native7': {
        'polynomial': (-1, 6, -5, 1),
        'brackets_hundredths': ((19, 20), (155, 156), (324, 325)),
        'endpoint_values_millionths': ((-33641, 8000), (11375, -11584), (-35776, 15625)),
        'basis_traces': (3, 5, 13), 'discriminant': 49,
        'coefficient_limits': (4, 9, 3), 'coefficient_trials': 1197,
    },
    'native18': {
        'polynomial': (-3, 9, -6, 1),
        'brackets_hundredths': ((46, 47), (165, 166), (387, 388)),
        'endpoint_values_millionths': ((-32264, 8423), (7125, -19304), (-70797, 4672)),
        'basis_traces': (3, 6, 18), 'discriminant': 81,
        'coefficient_limits': (6, 9, 2), 'coefficient_trials': 1235,
    },
}
FIELD_ORDER = ('native7', 'native18')
LOW, HIGH = F(36, 17), F(3969, 425)
FAILURES = (
    'physical_V_below_strip', 'physical_V_above_strip',
    'other_V_outside_strip', 'own_not_native_or_dotted',
    'time_not_field_square', 'physical_g_positive', 'g_nonintegral',
    'L_nonintegral', 'nonzero_cycle_outside_other_open_intervals',
    'mutual_not_native_or_dotted',
)
SCOPE = 'conditional_A17_native7_or18_isolated_anchor_and_other_isolated_end_only'
MAX_LINE, MAX_FILE = 65536, 16 * 1024 * 1024


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
    require(valid_sha(source_sha) and digest(SOURCE) == source_sha, 'source SHA256 mismatch')
    for rel, expected in PINS.items():
        require(valid_sha(expected) and digest(ROOT / rel) == expected,
                'input SHA256 mismatch: ' + rel)


def polynomial_value(key, x):
    value = F(0)
    for coefficient in reversed(FIELDS[key]['polynomial']):
        value = value * x + coefficient
    return value


def determinant(matrix):
    a, b, c = matrix[0]
    d, e, f = matrix[1]
    g, h, i = matrix[2]
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g)


@dataclass(frozen=True)
class C:
    key: str
    coefficients: tuple

    def __post_init__(self):
        require(self.key in FIELDS and len(self.coefficients) == 3, 'cubic field schema')
        object.__setattr__(self, 'coefficients', tuple(F(x) for x in self.coefficients))

    def lift(self, other):
        if isinstance(other, C):
            require(self.key == other.key, 'different cubic fields')
            return other
        require(isinstance(other, (int, F)), 'non-exact cubic scalar')
        return C(self.key, (other, 0, 0))

    def __add__(self, other):
        other = self.lift(other)
        return C(self.key, tuple(a + b for a, b in zip(self.coefficients, other.coefficients)))

    __radd__ = __add__

    def __neg__(self):
        return C(self.key, tuple(-x for x in self.coefficients))

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        product = [F(0)] * 5
        for i, a in enumerate(self.coefficients):
            for j, b in enumerate(other.coefficients):
                product[i + j] += a * b
        polynomial = FIELDS[self.key]['polynomial']
        for power in (4, 3):
            coefficient = product[power]
            for j in range(3):
                product[power - 3 + j] -= coefficient * polynomial[j]
            product[power] = F(0)
        return C(self.key, tuple(product[:3]))

    __rmul__ = __mul__

    def scaled(self, scalar):
        scalar = F(scalar)
        return C(self.key, tuple(x * scalar for x in self.coefficients))

    def trace(self):
        return sum((a * b for a, b in zip(self.coefficients,
                                         FIELDS[self.key]['basis_traces'])), F(0))

    def e2(self):
        return (self.trace() ** 2 - (self * self).trace()) / 2

    def norm(self):
        columns = []
        for i in range(3):
            basis = tuple(int(i == j) for j in range(3))
            columns.append((self * C(self.key, basis)).coefficients)
        return determinant([[columns[j][i] for j in range(3)] for i in range(3)])

    def inverse(self):
        n = self.norm()
        require(n != 0, 'division by zero cubic element')
        return (self * self - self.trace() * self + self.e2()).scaled(1 / n)

    def __truediv__(self, other):
        if isinstance(other, (int, F)):
            require(other != 0, 'division by zero rational')
            return self.scaled(F(1) / other)
        return self * self.lift(other).inverse()

    def zero(self):
        return all(x == 0 for x in self.coefficients)

    def rational(self):
        return self.coefficients[1:] == (F(0), F(0))

    def integral(self):
        return all(x.denominator == 1 for x in self.coefficients)

    def sign(self, place=2):
        require(place in (0, 1, 2), 'invalid cubic embedding')
        a, b, c = self.coefficients
        if b == c == 0:
            return (a > 0) - (a < 0)
        low, high = (F(x, 100) for x in FIELDS[self.key]['brackets_hundredths'][place])
        left_positive = polynomial_value(self.key, low) > 0
        while True:
            lower = a + min(b * low, b * high) + min(c * low * low, c * high * high)
            upper = a + max(b * low, b * high) + max(c * low * low, c * high * high)
            if lower > 0:
                return 1
            if upper < 0:
                return -1
            middle = (low + high) / 2
            value = polynomial_value(self.key, middle)
            require(value != 0, 'rational root in irreducible cubic')
            if (value > 0) == left_positive:
                low = middle
            else:
                high = middle

    def encode(self):
        return {'field': self.key,
                'coefficients': [[x.numerator, x.denominator] for x in self.coefficients]}


def anchor(key):
    return C(key, (0, 1, 0))


def rational_sqrt(value):
    value = F(value)
    if value < 0:
        return None
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    return F(n, d) if n * n == value.numerator and d * d == value.denominator else None


def bounded_field_sqrt(R):
    """Complete K-square test for the contract's INTEGRAL bounded radicands."""
    require(R.integral(), 'square-test radicand must be integral in the full basis')
    if R.rational():
        root = rational_sqrt(R.coefficients[0])
        return None if root is None else C(R.key, (root, 0, 0))
    if any(R.sign(place) <= 0 for place in (0, 1, 2)):
        return None
    A, B, N = R.trace(), R.e2(), R.norm()
    require(all(x.denominator == 1 for x in (A, B, N)) and N > 0,
            'integral radicand invariants')
    n = isqrt(N.numerator)
    if n * n != N:
        return None
    for signed_norm in (n, -n):
        for t in range(-454, 455):
            z = (F(t * t) - A) / 2
            if z.denominator != 1 or z * z - 2 * signed_norm * t != B:
                continue
            denominator = R + z
            require(not denominator.zero(), 'nonrational R plus rational z is nonzero')
            root = (t * R + signed_norm) / denominator
            if root * root == R:
                if root.sign() < 0:
                    root = -root
                require(root.sign() > 0 and root.integral(), 'positive integral cubic square root')
                return root
    return None


def product_class(value):
    if value.sign() < 0:
        return None
    if value.zero():
        return 'native2'
    if (value - 4).sign() > 0:
        return 'dotted'
    if (value - 4).sign() == 0:
        return None
    U = anchor(value.key)
    alphabet = [(C(value.key, (1, 0, 0)), 'native3'),
                (C(value.key, (2, 0, 0)), 'native4'),
                (C(value.key, (3, 0, 0)), 'native6')]
    if value.key == 'native7':
        alphabet.extend([(U, 'native7'), (U * U - 3 * U + 3, 'native14')])
    else:
        alphabet.extend([(U, 'native18'), (U * U - 4 * U + 4, 'native9')])
    return next((label for item, label in alphabet if value == item), None)


def other_cycle_admissible(g, L):
    if g.zero():
        require(L.zero(), 'zero g has zero cycle')
        return True
    return all(L.sign(place) > 0 and (4 - L).sign(place) > 0 for place in (0, 1))


def arithmetic_fixtures():
    """Prepared only. Actual execution/credit belongs to guarded root units."""
    for key in FIELD_ORDER:
        config = FIELDS[key]
        polynomial = config['polynomial']
        constant, linear, quadratic, _ = polynomial
        discriminant = (quadratic * quadratic * linear * linear - 4 * linear ** 3
                        - 4 * quadratic ** 3 * constant - 27 * constant ** 2
                        + 18 * quadratic * linear * constant)
        require(discriminant == config['discriminant'], 'exact polynomial discriminant')
        for bounds, values in zip(config['brackets_hundredths'],
                                  config['endpoint_values_millionths']):
            require(values[0] * values[1] < 0, 'root-isolating endpoint signs')
            for bound, expected in zip(bounds, values):
                require(polynomial_value(key, F(bound, 100)) == F(expected, 1000000),
                        'exact root bracket endpoint')
        U = anchor(key)
        require((U * U * U + quadratic * U * U + linear * U + constant).zero(),
                'cubic product reduction')
        require(U * U.inverse() == C(key, (1, 0, 0)), 'exact inverse')
        require(U.norm() == -constant and U.trace() == -quadratic, 'trace and norm')
        require(tuple((U - 1).sign(i) for i in (0, 1, 2)) == (-1, 1, 1)
                and tuple((U - 2).sign(i) for i in (0, 1, 2)) == (-1, -1, 1),
                'three separate embedding signs')
        require(U.integral() and not C(key, (F(1, 2), 0, 0)).integral(), 'full original basis')
        zero = C(key, (0, 0, 0))
        require(bounded_field_sqrt(zero) == zero, 'zero square')
        require(bounded_field_sqrt(C(key, (4, 0, 0))) == C(key, (2, 0, 0)),
                'rational square')
        require(bounded_field_sqrt(C(key, (2, 0, 0))) is None, 'rational nonsquare in odd field')
        for root in (U, U - 1, 1 - U):
            found = bounded_field_sqrt(root * root)
            require(found is not None and found * found == root * root and found.sign() > 0,
                    'mixed cubic roots and BOTH norm signs')
        require((U - 1).norm() == -1, 'negative-norm square-root fixture')
        require(product_class(zero) == 'native2' and other_cycle_admissible(zero, zero),
                'orthogonal pair retained')
        require(product_class(C(key, (4, 0, 0))) is None, 'compact parallel excluded')
        extra = U * U - 3 * U + 3 if key == 'native7' else U * U - 4 * U + 4
        require(product_class(extra) == ('native14' if key == 'native7' else 'native9'),
                'complete same-field mutual finite alphabet')
    return ['cubic_reduction', 'trace_norm_inverse', 'exact_three_root_brackets',
            'all_three_signs', 'full_integer_basis', 'zero_rational_squares',
            'odd_field_rational_nonsquare', 'both_signed_norm_cubic_roots',
            'orthogonal_pair', 'parallel_rejection', 'native9_and14_mutual_values']


def trial_record(key, a, b, c):
    V = C(key, (a, b, c))
    record = {'kind': 'trial', 'field': key, 'coefficients': [a, b, c], 'in_strip': False}
    if (V - LOW).sign() <= 0:
        return {**record, 'failure': 'physical_V_below_strip'}
    if (HIGH - V).sign() <= 0:
        return {**record, 'failure': 'physical_V_above_strip'}
    if any(V.sign(place) <= 0 or (LOW - V).sign(place) <= 0 for place in (0, 1)):
        return {**record, 'failure': 'other_V_outside_strip'}
    record['in_strip'] = True
    own_class = product_class(V)
    record['own_class'] = own_class
    if own_class is None:
        return {**record, 'failure': 'own_not_native_or_dotted'}
    U = anchor(key)
    Fanchor, Fother = U * (17 * U - 36), V * (17 * V - 36)
    R = Fanchor * Fother
    require(all(R.sign(place) > 0 for place in (0, 1, 2)), 'positive time radicand at all places')
    root = bounded_field_sqrt(R)
    if root is None:
        return {**record, 'failure': 'time_not_field_square'}
    g = (U * V - root) / 18
    if g.sign() > 0:
        return {**record, 'failure': 'physical_g_positive'}
    if not g.integral():
        return {**record, 'failure': 'g_nonintegral'}
    L = g * g / (U * V)
    if not L.integral():
        return {**record, 'failure': 'L_nonintegral'}
    if not other_cycle_admissible(g, L):
        return {**record, 'failure': 'nonzero_cycle_outside_other_open_intervals'}
    label = product_class(L)
    if label is None:
        return {**record, 'failure': 'mutual_not_native_or_dotted'}
    return {**record, 'failure': None, 'U': U.encode(), 'V': V.encode(),
            'F_anchor': Fanchor.encode(), 'F_other': Fother.encode(), 'R': R.encode(),
            'positive_sqrt_R': root.encode(), 'sqrt_trace': int(root.trace()),
            'sqrt_norm': int(root.norm()), 'g': g.encode(), 'L': L.encode(),
            'mutual_class': label, 'orthogonal': g.zero()}


def trials(key):
    la, lb, lc = FIELDS[key]['coefficient_limits']
    for a in range(-la, la + 1):
        for b in range(-lb, lb + 1):
            for c in range(-lc, lc + 1):
                yield trial_record(key, a, b, c)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


class Writer:
    def __init__(self, stream):
        self.stream, self.bytes = stream, 0

    def emit(self, record):
        line = (canonical(record) + '\n').encode('utf-8')
        require(len(line) <= MAX_LINE and self.bytes + len(line) <= MAX_FILE,
                'certificate record/output byte limit')
        self.stream.write(line)
        self.stream.flush()
        self.bytes += len(line)


def result_path(argument):
    path = (ROOT / argument).resolve()
    require((ROOT / 'results').resolve() in path.parents and path.parent.is_dir(),
            'certificate must be under an existing results/ directory')
    return path


def blank_totals():
    return {'fields': 0, 'coefficient_trials': 0, 'strip_models': 0, 'pair_tests': 0,
            'survivors': 0, 'orthogonal_survivors': 0,
            'rejects': {key: 0 for key in FAILURES}}


def count_trial(total, record):
    total['coefficient_trials'] += 1
    total['strip_models'] += int(record['in_strip'])
    total['pair_tests'] += int(record['in_strip'] and record.get('own_class') is not None)
    failure = record['failure']
    if failure is None:
        total['survivors'] += 1
        total['orthogonal_survivors'] += int(record['orthogonal'])
    else:
        require(failure in total['rejects'], 'unknown failure stage')
        total['rejects'][failure] += 1


def header(key, source_sha):
    return {'kind': 'header', 'schema': 1, 'unit': 'batch', 'field': key,
            'configuration': FIELDS[key], 'source_sha256': source_sha,
            'fixed_pins': PINS, 'scope': SCOPE}


def summary(key, fixtures, total):
    return {'kind': 'summary', 'unit': 'batch', 'field': key, 'status': 'PASS',
            'scope': SCOPE, 'arithmetic_fixtures': fixtures, **total}


def check_totals(total, expected):
    require(total['coefficient_trials'] == expected
            and sum(total['rejects'].values()) + total['survivors'] == expected,
            'complete scalar-triple partition')


def batch(args, fixtures):
    total = blank_totals()
    total['fields'] = 1
    with result_path(args.output).open('xb') as stream:
        writer = Writer(stream)
        writer.emit(header(args.field, args.expected_source_sha256))
        for record in trials(args.field):
            writer.emit(record)
            count_trial(total, record)
        check_totals(total, FIELDS[args.field]['coefficient_trials'])
        validate_pins(args.expected_source_sha256)
        final = summary(args.field, fixtures, total)
        writer.emit(final)
    print(canonical(final))


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON object key')
        result[key] = value
    return result


def invalid_number(value):
    raise ValueError('noninteger/nonfinite JSON number forbidden: ' + value)


def read_record(stream):
    line = stream.readline(MAX_LINE + 1)
    require(line and len(line) <= MAX_LINE and line.endswith(b'\n'),
            'missing/oversized/truncated JSONL record')
    record = json.loads(line, object_pairs_hook=unique_object,
                        parse_float=invalid_number, parse_constant=invalid_number)
    require(type(record) is dict and type(record.get('kind')) is str, 'JSONL schema')
    return record


def match_record(actual, expected):
    require(canonical(actual) == canonical(expected), 'exact record/coverage mismatch')


def combine(args, fixtures):
    require(len(args.part) == len(args.part_sha256) == 2, 'exactly TWO bound parts required')
    paths = [result_path(p) for p in args.part]
    require(len(set(paths)) == 2, 'duplicate part path')
    bound = []
    for path, expected in zip(paths, args.part_sha256):
        require(valid_sha(expected) and 0 < path.stat().st_size <= MAX_FILE,
                'part SHA/size schema')
        require(digest(path) == expected, 'part SHA256 mismatch')
        bound.append({'path': str(path.relative_to(ROOT)), 'sha256': expected,
                      'bytes': path.stat().st_size})
    total, seen = blank_totals(), set()
    with result_path(args.output).open('xb') as stream:
        writer = Writer(stream)
        writer.emit({'kind': 'header', 'schema': 1, 'unit': 'combine', 'scope': SCOPE,
                     'source_sha256': args.expected_source_sha256,
                     'fixed_pins': PINS, 'fields': list(FIELD_ORDER), 'parts': bound})
        for path, binding in zip(paths, bound):
            subtotal = blank_totals()
            subtotal['fields'] = 1
            with path.open('rb') as part:
                actual_header = read_record(part)
                key = actual_header.get('field')
                require(type(key) is str and key in FIELD_ORDER and key not in seen,
                        'duplicate/invalid field owner')
                match_record(actual_header, header(key, args.expected_source_sha256))
                seen.add(key)
                total['fields'] += 1
                writer.emit({'kind': 'part_binding', 'field': key, **binding})
                for expected in trials(key):
                    match_record(read_record(part), expected)
                    writer.emit(expected)
                    count_trial(subtotal, expected)
                    count_trial(total, expected)
                check_totals(subtotal, FIELDS[key]['coefficient_trials'])
                match_record(read_record(part), summary(key, fixtures, subtotal))
                require(part.read(1) == b'', 'extra bytes after terminal summary')
            require(digest(path) == binding['sha256']
                    and path.stat().st_size == binding['bytes'], 'part changed during combine')
        require(seen == set(FIELD_ORDER) and total['fields'] == 2, 'complete two-field conjunction')
        check_totals(total, 2432)
        validate_pins(args.expected_source_sha256)
        final = {'kind': 'summary', 'unit': 'combine', 'status': 'PASS', 'scope': SCOPE,
                 'arithmetic_fixtures': fixtures, 'batches': 2, **total}
        writer.emit(final)
    print(canonical(final))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    sub = parser.add_subparsers(dest='unit', required=True)
    one = sub.add_parser('batch')
    one.add_argument('--field', choices=FIELD_ORDER, required=True)
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
