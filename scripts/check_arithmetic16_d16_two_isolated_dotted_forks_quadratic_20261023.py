#!/usr/bin/env python3
"""UNEXECUTED preparation: conditional D16 quadratic head/two-DOT-fork census.

Standard-library exact arithmetic only; no repository imports or children.
Three parts and a SEPARATE complete same-implementation byte replay.
No geometric reduction, realized example or whole-D16 result is asserted.
Root alone owns any future unchanged resource-guarded launch.
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
RESULTS = (ROOT / 'results').resolve()
PROTOCOL = 'd16_quadratic_two_dotted_forks_v1'
SCOPE = 'conditional_actual_quadratic_head_and_both_isolated_DOT_forks_only'
PARTS = ((0, 40), (40, 80), (80, 119))
MAX_LINE = 2 * 1024 * 1024
MAX_FILE = 64 * 1024 * 1024
PINS = {
    'results/arithmetic16_d16_two_isolated_dotted_forks_exact_quadratic_census_contract_20261023.md':
        '596d1bfce60ae0ad2f82fe0c154f373dd1898eca359b2d03861dcf62ac5b4b70',
    'scripts/check_arithmetic16_b16_isolated_quadratic_head_short_pairs_20261023.py':
        '452587a42fcc03219cc77dabfa1ff0f691c81617c776943276b39e822d8ff42f',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
TOTAL_KEYS = ('fields', 'heads', 'forks', 'head_fork_pairs',
              'ordered_fork_pairs', 'ordered_triples', 'survivors')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def valid_sha(value):
    return (type(value) is str and len(value) == 64
            and all(c in '0123456789abcdef' for c in value))


def digest(path):
    answer = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            answer.update(block)
    return answer.hexdigest()


def validate_pins(expected):
    require(valid_sha(expected), 'invalid expected worker SHA256')
    require(digest(SOURCE) == expected, 'worker source SHA256 mismatch')
    require(len(PINS) == 6, 'complete contract/helper/guard pin count')
    for rel, sha in PINS.items():
        require(valid_sha(sha) and digest(ROOT / rel) == sha,
                'immutable input SHA256 mismatch: ' + rel)


def canonical(record):
    return (json.dumps(record, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True, allow_nan=False) + '\n').encode('ascii')


@dataclass(frozen=True)
class Q:
    d: int
    a: F
    b: F = F(0)

    def __post_init__(self):
        require(type(self.d) is int and self.d > 1, 'invalid field radicand')
        require(isinstance(self.a, (int, F)) and isinstance(self.b, (int, F)),
                'non-exact quadratic coefficient')
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
        product = self * other.conjugate()
        return Q(self.d, product.a / n, product.b / n)

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
        return {'D': self.d,
                'rational': [self.a.numerator, self.a.denominator],
                'radical': [self.b.numerator, self.b.denominator]}


def rational_sqrt(value):
    value = F(value)
    if value < 0:
        return None
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    if n * n == value.numerator and d * d == value.denominator:
        return F(n, d)
    return None


def field_sqrt(value):
    """Complete K-square test; physical nonnegative root, BOTH norm signs."""
    if value.b == 0:
        root = rational_sqrt(value.a)
        if root is not None:
            return Q(value.d, root)
        root = rational_sqrt(value.a / value.d)
        if root is not None:
            answer = Q(value.d, 0, root)
            require(answer * answer == value, 'pure-radical square identity')
            return answer
        return None
    n = rational_sqrt(value.norm())
    if n is None:
        return None
    for signed_norm in (n, -n):
        a = rational_sqrt((value.a + signed_norm) / 2)
        if a is not None and a != 0:
            answer = Q(value.d, a, value.b / (2 * a))
            if answer * answer == value:
                if answer.sign() < 0:
                    answer = -answer
                require(answer.sign() > 0, 'mixed square root sign')
                return answer
    return None


def product_class(value):
    """Entire quadratic native alphabet; zero allowed, compact parallel4 not."""
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


def squarefree(d):
    return all(d % (p * p) != 0 for p in range(2, isqrt(d) + 1))


def field_list():
    answer = []
    for d in range(2, 119):
        if squarefree(d) and (d % 4 == 1 or d <= 29):
            delta = d if d % 4 == 1 else 4 * d
            require(delta <= 118, 'fundamental discriminant bound')
            answer.append((delta, d))
    return sorted(answer)


def models(d, head):
    answer = []
    maximum_a, maximum_gap = (12, 118) if head else (60, 3660)
    physical_max = F(98, 9) if head else F(121, 2)
    other_max = F(2) if head else F(1, 2)
    for a in range(5, maximum_a + 1):
        for b in range(1, isqrt(maximum_gap // d) + 1):
            if d % 4 == 1:
                if (a - b) % 2:
                    continue
            elif a % 2 or b % 2:
                continue
            x = Q(d, F(a, 2), F(b, 2))
            if ((x - 4).sign() <= 0 or (physical_max - x).sign() <= 0
                    or x.conjugate().sign() <= 0
                    or (other_max - x.conjugate()).sign() <= 0):
                continue
            require(x.integral() and x.b > 0, 'full-ring generator invariant')
            time = x * (x - 2) if head else 2 * x * (2 * x - 1)
            require(time.sign() > 0 and time.conjugate().sign() < 0,
                    'raw time-sign invariant')
            answer.append({'a': a, 'b': b, 'x': x, 'time': time})
    require(len(answer) <= (28 if head else 42), 'proved model-count bound')
    return answer


def public_model(model):
    return {'trace_a': model['a'], 'gap_b': model['b'],
            'value': model['x'].encode(), 'time': model['time'].encode()}


def pair_result(left, right, head_fork):
    product = left['x'] * right['x']
    quarter = product / 4
    if not quarter.integral():
        return {'status': 'quarter_product_nonintegral'}
    radicand = left['time'] * right['time']
    require(radicand.sign() > 0 and radicand.conjugate().sign() > 0,
            'time-product total positivity')
    root = field_sqrt(radicand)
    if root is None:
        return {'status': 'time_product_not_K_square'}
    require(root.sign() > 0 and root * root == radicand,
            'physical time-product square root')
    spatial = product / 2 if head_fork else 7 * product / 2
    g = root - spatial
    if g.sign() < 0:
        return {'status': 'physical_g_negative'}
    if not g.integral():
        return {'status': 'path_cycle_g_nonintegral'}
    parameter = g * g / product
    identity_rhs = ((15 * product - 32 * right['x'] - 8 * left['x'] + 16) / 4
                    if head_fork else 15 * product / 4 - 8 * (left['x'] + right['x']) + 4)
    require(parameter + (g if head_fork else 7 * g) == identity_rhs,
            'exact raw pair normalization identity')
    if not parameter.integral():
        return {'status': 'two_cycle_L_nonintegral'}
    if (parameter.conjugate().sign() < 0
            or (4 - parameter.conjugate()).sign() <= 0):
        return {'status': 'other_L_outside_nonnegative_open_upper'}
    label = product_class(parameter)
    if label is None:
        return {'status': 'physical_L_not_native_or_dotted'}
    return {'status': 'retained_pair', 'quarter': quarter, 'sqrt_time': root,
            'g': g, 'L': parameter, 'label': label}


def public_pair(result):
    if result['status'] != 'retained_pair':
        return {'status': result['status']}
    return {'status': 'retained_pair', 'quarter': result['quarter'].encode(),
            'sqrt_time': result['sqrt_time'].encode(),
            'g': result['g'].encode(), 'L': result['L'].encode(),
            'label': result['label']}


def bump(counts, label):
    counts[label] = counts.get(label, 0) + 1


def zero_totals():
    return {key: 0 for key in TOTAL_KEYS}


def add_totals(total, incoming):
    require(set(incoming) == set(TOTAL_KEYS), 'complete total schema')
    for key in TOTAL_KEYS:
        require(type(incoming[key]) is int and incoming[key] >= 0,
                'nonnegative integer total')
        total[key] += incoming[key]


def field_records(delta, d):
    heads, forks = models(d, True), models(d, False)
    hf, ff = {}, {}
    hf_counts, ff_counts = {}, {}
    pairs_digest = hashlib.sha256()
    for h, head in enumerate(heads):
        for i, fork in enumerate(forks):
            result = pair_result(head, fork, True)
            hf[h, i] = result
            bump(hf_counts, result['status'])
            pairs_digest.update(canonical({'kind': 'head_fork', 'indices': [h, i],
                                           'certificate': public_pair(result)}))
    for i, left in enumerate(forks):
        for j, right in enumerate(forks):
            result = pair_result(left, right, False)
            ff[i, j] = result
            bump(ff_counts, result['status'])
            pairs_digest.update(canonical({'kind': 'ordered_fork_pair',
                                           'indices': [i, j],
                                           'certificate': public_pair(result)}))
    require(len(hf) + len(ff) <= 2940, 'proved pair-cache bound')
    yield {'kind': 'field_begin', 'D': d, 'Delta': delta,
           'heads': [public_model(x) for x in heads],
           'forks': [public_model(x) for x in forks]}
    counts, triples_digest = {}, hashlib.sha256()
    for h, head in enumerate(heads):
        for i, left in enumerate(forks):
            for j, right in enumerate(forks):
                p, q, r = hf[h, i], hf[h, j], ff[i, j]
                if p['status'] != 'retained_pair':
                    status = 'head_fork_1/' + p['status']
                elif q['status'] != 'retained_pair':
                    status = 'head_fork_2/' + q['status']
                elif r['status'] != 'retained_pair':
                    status = 'fork_pair/' + r['status']
                else:
                    triangle = (p['g'] * r['g'] * q['g']
                                / (head['x'] * left['x'] * right['x']))
                    status = 'survivor' if triangle.integral() else 'three_cycle_nonintegral'
                bump(counts, status)
                if status == 'survivor':
                    record = {'kind': 'survivor', 'D': d, 'indices': [h, i, j],
                              'head_fork_1': public_pair(p),
                              'head_fork_2': public_pair(q),
                              'ordered_fork_pair': public_pair(r),
                              'three_cycle': triangle.encode()}
                    triples_digest.update(canonical(record))
                    yield record
                else:
                    triples_digest.update(canonical({'kind': 'triple_outcome',
                        'D': d, 'indices': [h, i, j], 'status': status}))
    number = len(heads) * len(forks) * len(forks)
    require(number <= 49392 and sum(counts.values()) == number,
            'complete ordered-triple accounting')
    require(sum(hf_counts.values()) == len(heads) * len(forks)
            and sum(ff_counts.values()) == len(forks) * len(forks),
            'complete pair accounting')
    totals = {'fields': 1, 'heads': len(heads), 'forks': len(forks),
              'head_fork_pairs': len(hf), 'ordered_fork_pairs': len(ff),
              'ordered_triples': number, 'survivors': counts.get('survivor', 0)}
    yield {'kind': 'field_end', 'D': d, 'Delta': delta, 'totals': totals,
           'head_fork_counts': hf_counts, 'ordered_fork_counts': ff_counts,
           'triple_first_outcome_counts': counts,
           'every_pair_sha256': pairs_digest.hexdigest(),
           'every_ordered_triple_sha256': triples_digest.hexdigest()}


def part_records(index, expected):
    validate_pins(expected)
    lower, upper = PARTS[index]
    selected = [(delta, d) for delta, d in field_list() if lower <= delta < upper]
    yield {'kind': 'batch_header', 'protocol': PROTOCOL, 'scope': SCOPE,
           'worker_sha256': expected, 'pins': PINS, 'part_index': index,
           'Delta_interval': [lower, upper],
           'fields': [{'D': d, 'Delta': delta} for delta, d in selected],
           'geometric_reduction': 'EXTERNAL_UNACCEPTED_AT_PREPARATION',
           'ordered_forks': True, 'equal_forks_included': True,
           'independent_implementation': False, 'realization_claim': False}
    totals = zero_totals()
    for delta, d in selected:
        for record in field_records(delta, d):
            if record['kind'] == 'field_end':
                add_totals(totals, record['totals'])
            yield record
    require(totals['fields'] == len(selected), 'complete part fields')
    validate_pins(expected)
    yield {'kind': 'batch_summary', 'protocol': PROTOCOL, 'part_index': index,
           'status': 'ENUMERATION_COMPLETE', 'totals': totals,
           'geometry_conclusion': False, 'independent_implementation': False}


def results_path(argument, new):
    original = ROOT / argument
    require(not original.is_symlink(), 'direct results/ paths required')
    path = original.resolve()
    require(path.parent == RESULTS, 'path must be a flat results/ file')
    if new:
        require(not path.exists() and not path.is_symlink(), 'output must be NEW')
    else:
        require(path.is_file() and 0 < path.stat().st_size <= MAX_FILE,
                'missing/empty/oversized part file')
    return path


class Writer:
    def __init__(self, path):
        self.stream = path.open('xb')
        self.bytes = 0
        self.records = 0

    def write(self, record):
        raw = canonical(record)
        require(len(raw) <= MAX_LINE, 'record exceeds cap; no truncation')
        require(self.bytes + len(raw) <= MAX_FILE, 'file exceeds cap; no truncation')
        self.stream.write(raw)
        self.bytes += len(raw)
        self.records += 1

    def close(self):
        self.stream.close()


def run_batch(args):
    path = results_path(args.output, True)
    validate_pins(args.expected_source_sha256)
    writer = Writer(path)
    summary = None
    try:
        for record in part_records(args.index, args.expected_source_sha256):
            writer.write(record)
            if record['kind'] == 'batch_summary':
                summary = record
    finally:
        writer.close()
    require(summary is not None, 'missing complete part summary')
    validate_pins(args.expected_source_sha256)
    print(canonical({'status': 'ENUMERATION_COMPLETE', 'output': str(path.relative_to(ROOT)),
                     'sha256': digest(path), 'bytes': writer.bytes,
                     'records': writer.records, 'totals': summary['totals']}).decode('ascii'), end='')


def run_combine(args):
    require(len(args.part) == 3 and len(args.part_sha256) == 3,
            'exactly three parts and hashes in index order required')
    paths = [results_path(name, False) for name in args.part]
    require(len(set(paths)) == 3, 'duplicate part paths')
    for path, sha in zip(paths, args.part_sha256):
        require(valid_sha(sha) and digest(path) == sha, 'part SHA256 mismatch')
    path_out = results_path(args.output, True)
    validate_pins(args.expected_source_sha256)
    writer = Writer(path_out)
    totals, summaries = zero_totals(), []
    try:
        writer.write({'kind': 'combine_header', 'protocol': PROTOCOL, 'scope': SCOPE,
                      'worker_sha256': args.expected_source_sha256, 'pins': PINS,
                      'parts': [{'path': str(p.relative_to(ROOT)), 'sha256': sha}
                                for p, sha in zip(paths, args.part_sha256)],
                      'operation': 'SAME_IMPLEMENTATION_EVERY_PAIR_AND_TRIPLE_REPLAY',
                      'independent_implementation': False, 'geometry_conclusion': False})
        for index, (path, sha) in enumerate(zip(paths, args.part_sha256)):
            summary = None
            with path.open('rb') as stream:
                for expected in part_records(index, args.expected_source_sha256):
                    raw = stream.readline(MAX_LINE + 1)
                    require(0 < len(raw) <= MAX_LINE and raw.endswith(b'\n'),
                            'missing/oversized/partial record')
                    require(raw == canonical(expected), 'noncanonical or unequal replay record')
                    writer.write(expected)
                    if expected['kind'] == 'batch_summary':
                        summary = expected
                require(stream.read(1) == b'', 'unexpected extra records or trailing bytes')
            require(digest(path) == sha, 'part changed during replay')
            require(summary is not None, 'missing authentic part summary')
            add_totals(totals, summary['totals'])
            summaries.append(summary)
        require(totals['fields'] == len(field_list()), 'complete disjoint field union')
        validate_pins(args.expected_source_sha256)
        writer.write({'kind': 'combine_summary', 'protocol': PROTOCOL,
                      'status': 'SAME_IMPLEMENTATION_REPLAY_COMPLETE',
                      'enumeration_status': 'ENUMERATION_COMPLETE', 'totals': totals,
                      'part_summaries': summaries, 'independent_implementation': False,
                      'geometry_conclusion': False, 'realization_claim': False})
    finally:
        writer.close()
    validate_pins(args.expected_source_sha256)
    print(canonical({'status': 'SAME_IMPLEMENTATION_REPLAY_COMPLETE',
                     'output': str(path_out.relative_to(ROOT)), 'sha256': digest(path_out),
                     'bytes': writer.bytes, 'records': writer.records,
                     'totals': totals, 'independent_implementation': False}).decode('ascii'), end='')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    commands = parser.add_subparsers(dest='command', required=True)
    batch = commands.add_parser('batch')
    batch.add_argument('--index', type=int, choices=range(3), required=True)
    batch.add_argument('--output', required=True)
    combine = commands.add_parser('combine')
    combine.add_argument('--part', action='append', required=True)
    combine.add_argument('--part-sha256', action='append', required=True)
    combine.add_argument('--output', required=True)
    args = parser.parse_args()
    if args.command == 'batch':
        run_batch(args)
    else:
        run_combine(args)


if __name__ == '__main__':
    main()
