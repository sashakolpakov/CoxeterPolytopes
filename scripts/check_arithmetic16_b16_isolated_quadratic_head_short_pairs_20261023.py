#!/usr/bin/env python3
"""Prepared necessary census: isolated quadratic B16 head/short pairs.

UNEXECUTED at preparation. Root-owned unchanged guarded jobs only.
All quadratic maximal orders, fundamental discriminant<144; three parts
and a SEPARATE authenticating every-pair combine. Local exact helpers
adapt the pinned older A17 source, without repository imports or its
geometry. No survivor is an example; no cubic or whole-B16 claim.
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
PINS = {
    'results/arithmetic16_whole_b16_complete_short_white_and_isolated_dotted_field_reduction_candidate_20261023.md':
        '4fb202752c01f8ace7f38763403ef93dc13f937ddf1ee31cc656cdf80c27a514',
    'results/arithmetic_dimension16_fresh_complete_necessary_forest_frontier_candidate_20261022.md':
        'a8cc202cc364f438e733c5c2536185aa923bc56ee2d3edc9e245a3c8d7e20e3c',
    'results/arithmetic16_whole_b16_long_head_and_near_short_three_scale_obstruction_candidate_20261023.md':
        'b11bfefeb6babec483f4cf5cb5da9dc6a71caf52c967c45af8b7f6bd100e72a4',
    'results/review_arithmetic16_whole_b16_head_near_short_three_scale_root_acceptance_20261023.md':
        '7ceec2f7c3bf9a0a0e0943f48b66afe2317f774b72e97f6e6779e2c1ac77beea',
    'results/arithmetic_uniform_classical_weak_helper_obstruction_20261021.md':
        '535ac0d0fa92d28ccab2976fa147a4ccff0bd3bfcd1e04ca4e134c6ae6857905',
    'results/compact_classical_generic_prefix_energy_20261020.md':
        'f637d7294d17e37581330ca508c914c23a08ce259845614b055e45815d2790ef',
    'results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md':
        '3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84',
    'results/compact_dimension19_and_arithmetic_dimension17_global_composition_20261021.md':
        '0ac802108e0ef8a9ca821c3091fca82f7bf9331a0037639d93888f2fabb37a75',
    'results/arithmetic16_b16_isolated_quadratic_head_short_pair_exact_census_contract_20261023.md':
        '9fb3fa877e82e82bc26a32c10fb37261e6f32cd43c73ef0c5f6fc446c846dc3c',
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
SCOPE = 'conditional_quadratic_B16_isolated_head_and_dotted_short_only'
BATCH_RANGES = ((0, 48), (48, 96), (96, 144))
FAILURES = (
    'head_own_not_native_or_dotted', 'half_product_nonintegral',
    'time_not_field_square', 'physical_g_negative', 'g_nonintegral',
    'L_nonintegral', 'nonphysical_cycle_outside_open_interval',
    'mutual_not_native_or_dotted',
)
MAX_LINE = 2 * 1024 * 1024
MAX_FILE = 64 * 1024 * 1024


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
    require(valid_sha(source_sha), 'invalid externally expected source SHA256')
    require(digest(SOURCE) == source_sha, 'worker source SHA256 mismatch')
    require(len(PINS) == 14, 'complete input/guard pin count')
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
        return {'D': self.d,
                'a': [self.a.numerator, self.a.denominator],
                'b': [self.b.numerator, self.b.denominator]}


def rational_sqrt(value):
    value = F(value)
    if value < 0:
        return None
    n, d = isqrt(value.numerator), isqrt(value.denominator)
    return F(n, d) if n * n == value.numerator and d * d == value.denominator else None


def field_sqrt(value):
    """Complete K-square test; PHYSICAL nonnegative root, both norm signs."""
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
            root = Q(value.d, a, value.b / (2 * a))
            if root * root == value:
                if root.sign() < 0:
                    root = -root
                require(root.sign() > 0, 'mixed square root sign')
                return root
    return None


def product_class(value):
    """4cos^2(pi/m), zero, or dotted; compact parallel4 is excluded."""
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


def fields():
    answer = []
    for d in range(2, 144):
        if squarefree(d):
            delta = d if d % 4 == 1 else 4 * d
            if delta < 144:
                answer.append((delta, d))
    answer.sort()
    require(len(set(answer)) == len(answer)
            and len({delta for delta, _ in answer}) == len(answer),
            'unique complete fundamental discriminants')
    return answer


def selected_fields(index):
    require(type(index) is int and 0 <= index < 3, 'batch selector')
    low, high = BATCH_RANGES[index]
    return [(delta, d) for delta, d in fields() if low <= delta < high]


def models(d, head):
    omega = Q(d, F(1, 2), F(1, 2)) if d % 4 == 1 else Q(d, 0, 1)
    low, high, other_high = (2, 12, F(2)) if head else (4, 64, F(1, 4))
    max_a, max_b, max_models = (7, 5, 10) if head else (32, 28, 28)
    answer = []
    for b in range(1, max_b + 1):
        for a in range(1, max_a + 1):
            value = a + b * omega
            other = value.conjugate()
            if ((value - low).sign() > 0 and (high - value).sign() > 0
                    and other.sign() > 0 and (other_high - other).sign() > 0):
                require(value.integral() and value.b > 0, 'full-ring generator')
                answer.append((a, b, value))
    require(len(answer) <= max_models and len(set(answer)) == len(answer),
            'proved per-field model bound/uniqueness')
    return answer


def model_record(model):
    a, b, value = model
    return {'coefficients': [a, b], 'value': value.encode(),
            'own_class': product_class(value)}


def pair_record(i, j, head, short):
    U, V = head[2], short[2]
    base = {'index_U': i, 'index_V': j}
    if product_class(U) is None:
        return {**base, 'outcome': FAILURES[0]}
    if not (U * V / 2).integral():
        return {**base, 'outcome': FAILURES[1]}
    R = (U - 2) * V / (4 * U * (4 * V - 1))
    require(R.sign() > 0 and R.conjugate().sign() > 0,
            'necessary radicand positive at both places')
    k = field_sqrt(R)
    if k is None:
        return {**base, 'outcome': FAILURES[2]}
    require(k.sign() > 0 and k * k == R, 'actual positive-time K-square root')
    G = 4 * U * k * (4 * V - 1) - U * V
    L = G * G / (2 * U * V)
    require(L + G == (15 * U - 32) * V / 2 - 2 * U + 4,
            'fresh even-rank head/short cyclic identity')
    if G.sign() < 0:
        return {**base, 'outcome': FAILURES[3]}
    if not G.integral():
        return {**base, 'outcome': FAILURES[4]}
    if not L.integral():
        return {**base, 'outcome': FAILURES[5]}
    if G.sign() == 0:
        require(L.sign() == 0, 'orthogonal pair has zero two-cycle')
    elif not (L.conjugate().sign() > 0 and (4 - L.conjugate()).sign() > 0):
        return {**base, 'outcome': FAILURES[6]}
    label = product_class(L)
    if label is None:
        return {**base, 'outcome': FAILURES[7]}
    return {**base, 'outcome': 'survivor', 'R': R.encode(), 'k': k.encode(),
            'G': G.encode(), 'L': L.encode(), 'mutual_class': label}


def canonical(record):
    return (json.dumps(record, sort_keys=True, separators=(',', ':'),
                       ensure_ascii=True) + '\n').encode('ascii')


def empty_counts():
    return {key: 0 for key in (*FAILURES, 'survivor')}


def census_field(index, delta, d):
    require((delta, d) in selected_fields(index), 'field belongs to selector')
    heads, shorts = models(d, True), models(d, False)
    counts, survivors = empty_counts(), []
    pair_hash = hashlib.sha256()
    for i, head in enumerate(heads):
        for j, short in enumerate(shorts):
            record = pair_record(i, j, head, short)
            require(record['outcome'] in counts, 'exhaustive pair outcome')
            counts[record['outcome']] += 1
            pair_hash.update(canonical(record))
            if record['outcome'] == 'survivor':
                survivors.append(record)
    pairs = len(heads) * len(shorts)
    require(pairs <= 280 and sum(counts.values()) == pairs
            and counts['survivor'] == len(survivors), 'all ordered pairs accounted')
    return {'kind': 'field', 'batch_index': index,
            'fundamental_discriminant': delta, 'squarefree_D': d,
            'head_models': [model_record(item) for item in heads],
            'short_models': [model_record(item) for item in shorts],
            'pair_count': pairs, 'counts': counts,
            'pair_outcome_sha256': pair_hash.hexdigest(), 'survivors': survivors}


def empty_totals():
    return {'fields': 0, 'head_models': 0, 'short_models': 0,
            'ordered_pairs': 0, 'counts': empty_counts()}


def accumulate(totals, record):
    totals['fields'] += 1
    totals['head_models'] += len(record['head_models'])
    totals['short_models'] += len(record['short_models'])
    totals['ordered_pairs'] += record['pair_count']
    for key in totals['counts']:
        totals['counts'][key] += record['counts'][key]
    require(sum(totals['counts'].values()) == totals['ordered_pairs'],
            'terminal pair accounting')


def result_path(value, existing):
    path = Path(value)
    path = (ROOT / path).resolve() if not path.is_absolute() else path.resolve()
    require(path.is_relative_to(RESULTS) and path.suffix == '.jsonl',
            'exclusive result path must be JSONL inside results')
    require(path.parent.is_dir(), 'output parent must already exist')
    require(path.is_file() if existing else not path.exists(),
            'missing input or nonexclusive output')
    if existing:
        require(path.stat().st_size <= MAX_FILE, 'oversized part file')
    return path


def relpath(path):
    return str(path.relative_to(ROOT))


def header(index, path, source_sha):
    return {'kind': 'b16_isolated_quadratic_pair_part', 'version': 1,
            'scope': SCOPE, 'batch_index': index,
            'worker': relpath(SOURCE), 'worker_sha256': source_sha,
            'pins': PINS, 'output': relpath(path),
            'field_selector': [list(item) for item in selected_fields(index)]}


def footer(index, totals, fields_sha):
    return {'kind': 'complete', 'status': 'PASS', 'scope': SCOPE,
            'batch_index': index, 'totals': totals, 'fields_sha256': fields_sha}


class Writer:
    def __init__(self, path):
        self.stream = path.open('xb')
        self.size = 0

    def write(self, record):
        raw = canonical(record)
        require(len(raw) <= MAX_LINE, 'oversized record; no truncation allowed')
        require(self.size + len(raw) <= MAX_FILE,
                'storage limit; no survivor truncation allowed')
        self.stream.write(raw)
        self.size += len(raw)
        return raw

    def close(self):
        self.stream.close()


def read_record(stream):
    raw = stream.readline(MAX_LINE + 1)
    require(raw and len(raw) <= MAX_LINE and raw.endswith(b'\n'),
            'missing/truncated/oversized JSONL record')
    item = json.loads(raw)
    require(type(item) is dict and canonical(item) == raw, 'noncanonical record')
    return item


def same_record(left, right):
    return canonical(left) == canonical(right)


def batch(args):
    validate_pins(args.expected_source_sha256)
    path = result_path(args.output, False)
    require(0 <= args.index < 3, 'valid disjoint batch index')
    writer = Writer(path)
    totals, field_hash = empty_totals(), hashlib.sha256()
    try:
        writer.write(header(args.index, path, args.expected_source_sha256))
        for delta, d in selected_fields(args.index):
            record = census_field(args.index, delta, d)
            field_hash.update(writer.write(record))
            accumulate(totals, record)
        validate_pins(args.expected_source_sha256)
        writer.write(footer(args.index, totals, field_hash.hexdigest()))
    finally:
        writer.close()
    return path, totals


def combine(args):
    validate_pins(args.expected_source_sha256)
    require(len(args.part) == 3 and len(args.part_sha256) == 3,
            'exactly three ordered authenticated parts required')
    paths = [result_path(value, True) for value in args.part]
    require(len(set(paths)) == 3, 'duplicate part path')
    for path, expected in zip(paths, args.part_sha256):
        require(valid_sha(expected) and digest(path) == expected,
                'part SHA256 mismatch')
    path = result_path(args.output, False)
    writer = Writer(path)
    totals, field_hash = empty_totals(), hashlib.sha256()
    try:
        writer.write({'kind': 'b16_isolated_quadratic_pair_combined', 'version': 1,
                      'scope': SCOPE, 'worker': relpath(SOURCE),
                      'worker_sha256': args.expected_source_sha256, 'pins': PINS,
                      'output': relpath(path),
                      'full_field_domain': [list(item) for item in fields()],
                      'parts': [{'batch_index': i, 'path': relpath(p), 'sha256': s}
                                for i, (p, s) in enumerate(zip(paths, args.part_sha256))]})
        seen = []
        for index, part in enumerate(paths):
            part_totals, part_hash = empty_totals(), hashlib.sha256()
            with part.open('rb') as stream:
                require(same_record(read_record(stream),
                                    header(index, part, args.expected_source_sha256)),
                        'part header/source/pins/selector/order mismatch')
                for delta, d in selected_fields(index):
                    expected = census_field(index, delta, d)
                    require(same_record(read_record(stream), expected),
                            'full field/every-pair certificate mismatch')
                    raw = writer.write(expected)
                    part_hash.update(raw)
                    field_hash.update(raw)
                    accumulate(part_totals, expected)
                    accumulate(totals, expected)
                    seen.append((delta, d))
                require(same_record(read_record(stream),
                                    footer(index, part_totals, part_hash.hexdigest())),
                        'part terminal summary/counters mismatch')
                require(stream.read(1) == b'', 'extra or duplicated trailing records')
            require(digest(part) == args.part_sha256[index], 'part changed during replay')
        require(seen == fields() and len(set(seen)) == len(seen),
                'complete disjoint all-field union')
        validate_pins(args.expected_source_sha256)
        writer.write(footer(None, totals, field_hash.hexdigest()))
    finally:
        writer.close()
    return path, totals


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    sub = parser.add_subparsers(dest='mode', required=True)
    part = sub.add_parser('batch')
    part.add_argument('--index', type=int, choices=range(3), required=True)
    part.add_argument('--output', required=True)
    whole = sub.add_parser('combine')
    whole.add_argument('--part', action='append', required=True)
    whole.add_argument('--part-sha256', action='append', required=True)
    whole.add_argument('--output', required=True)
    args = parser.parse_args()
    path, totals = batch(args) if args.mode == 'batch' else combine(args)
    print(json.dumps({'status': 'PASS', 'scope': SCOPE, 'mode': args.mode,
                      'output': relpath(path), 'sha256': digest(path),
                      'totals': totals}, sort_keys=True))


if __name__ == '__main__':
    main()
