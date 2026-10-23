#!/usr/bin/env python3
"""Prepared exact A18 BOTH-isolated-finite endpoint-pair census.

UNEXECUTED; root-owned guarded jobs only. No repository imports.
Sixteen (smaller own m, mutual k) batches, followed by an authenticated
coverage-only combine. Survivors are retained; no empty result is assumed.
"""

import argparse
import hashlib
import json
from math import isfinite, isqrt, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = Path(__file__).resolve()
SOURCE_REL = 'scripts/check_compact18_a18_isolated_finite_endpoint_pairs_20261023.py'
PINS = {
    'results/compact18_whole_a18_two_isolated_finite_endpoints_cyclotomic_census_contract_20261023.md':
        '1eb16218407a12b91b05f8a984995c5ad7b43630388b2d51fb75d49c64c00a57',
    'results/compact_farthest_vertex_mixed_replacement_inequalities_20261020.md':
        '3732b82614e13dfe2f7a3c66ab6c7fe67e517cbaeb77be0ca533f559ddf6ed84',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
SCOPE = 'general_compact18_whole_A18_both_actual_endpoint_caps_isolated_and_finite_ONLY'
OWN_M = tuple(range(7, 15))
MUTUAL_K = (3, 4)
UNITS = tuple((m, k) for m in OWN_M for k in MUTUAL_K)
DEGREES = dict(zip(OWN_M, (3, 4, 3, 4, 5, 4, 6, 6)))
MAX_N, MAX_L, MAX_PHI = 105, 11760, 2304
MAX_BITS = 4 * MAX_L + 128
MAX_LINE, MAX_FILE, MAX_LEDGER = 65536, 8 * 1024**2, 65536
FIXTURES = ('totient_and_even_cutoff', 'exact_cyclotomic_divisor_recursion',
            'monic_exact_and_nonexact_division', 'cyclic_integer_convolution',
            'zero_and_nonzero_cyclotomic_remainders')


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
    require(valid_sha(source_sha) and digest(SOURCE) == source_sha,
            'current checker SHA256 mismatch')
    for rel, expected in PINS.items():
        require(digest(ROOT / rel) == expected, 'current input SHA256 mismatch: ' + rel)


def result_path(argument):
    path = (ROOT / argument).resolve()
    require((ROOT / 'results').resolve() in path.parents and path.parent.is_dir(),
            'file must be inside an existing results/ directory')
    return path


def relative(path):
    return str(path.relative_to(ROOT))


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False)


def integer(value):
    require(type(value) is int and value.bit_length() <= MAX_BITS,
            'exact coefficient type/bit bound')
    return value


def phi(n):
    require(type(n) is int and n >= 1, 'positive totient argument')
    remaining, answer, p = n, n, 2
    while p * p <= remaining:
        if remaining % p == 0:
            answer -= answer // p
            while remaining % p == 0:
                remaining //= p
        p += 1
    if remaining > 1:
        answer -= answer // remaining
    return answer


def divisors(n):
    low, high = [], []
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            low.append(d)
            if d * d != n:
                high.append(n // d)
    return low + high[::-1]


def trim(poly):
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def divide_monic(dividend, divisor, exact):
    require(divisor and divisor[-1] == 1 and len(divisor) >= 2,
            'positive-degree monic divisor')
    require(1 <= len(dividend) <= MAX_L + 1, 'bounded dividend degree')
    remainder = trim(list(dividend))
    quotient = [0] * max(1, len(remainder) - len(divisor) + 1)
    terms = tuple((i, c) for i, c in enumerate(divisor[:-1]) if c)
    while len(remainder) >= len(divisor) and remainder != [0]:
        offset = len(remainder) - len(divisor)
        lead = integer(remainder[-1])
        quotient[offset] = lead
        for i, coefficient in terms:
            remainder[offset + i] = integer(remainder[offset + i] - lead * coefficient)
        remainder[-1] = 0
        trim(remainder)
    if exact:
        require(remainder == [0], 'nonexact claimed cyclotomic division')
    return tuple(trim(quotient)), tuple(remainder)


def cyclotomic(n):
    require(1 <= n <= MAX_L, 'bounded cyclotomic order')
    cache = {}
    for d in divisors(n):
        poly = (-1,) + (0,) * (d - 1) + (1,)
        for e in divisors(d)[:-1]:
            poly, _ = divide_monic(poly, cache[e], True)
        require(len(poly) == phi(d) + 1 and poly[-1] == 1 and poly[0] in (-1, 1),
                'cyclotomic degree/monicity/constant')
        cache[d] = poly
    require(sum(len(poly) - 1 for poly in cache.values()) == n,
            'one-order divisor cache degree identity')
    return cache[n]


def cyclic_multiply(left, right, order):
    result = {}
    for i, a in left.items():
        for j, b in right.items():
            exponent = (i + j) % order
            result[exponent] = integer(result.get(exponent, 0) + a * b)
    return {i: a for i, a in result.items() if a}


def sparse_expression(m, n, k, order):
    def cosine(label):
        require(order % (2 * label) == 0, 'exact cosine exponent')
        exponent = order // (2 * label)
        return {exponent: 1, (-exponent) % order: 1}
    a, b, c = (cosine(label) for label in (m, n, k))
    a2 = cyclic_multiply(a, a, order)
    b2 = cyclic_multiply(b, b, order)
    c2 = cyclic_multiply(c, c, order)
    a2b2 = cyclic_multiply(a2, b2, order)
    abc = cyclic_multiply(cyclic_multiply(a, b, order), c, order)
    result = {0: 76}
    for scale, poly in ((17, a2b2), (-36, a2), (-36, b2), (-2, abc), (-19, c2)):
        for exponent, coefficient in poly.items():
            result[exponent] = integer(result.get(exponent, 0) + scale * coefficient)
    return {i: a for i, a in result.items() if a}


def signed_hex(value):
    return ('-' if value < 0 else '+') + format(abs(value), 'x')


def polynomial_hash(poly):
    h = hashlib.sha256((str(len(poly)) + '\n').encode('ascii'))
    for coefficient in poly:
        h.update((signed_hex(coefficient) + '\n').encode('ascii'))
    return h.hexdigest()


def domain(m, k):
    require((m, k) in UNITS, 'one of the sixteen fixed units required')
    bound = 4 * DEGREES[m] * (1 if k == 3 else 2)
    return tuple(n for n in range(m, MAX_N + 1) if phi(2 * n) <= bound)


def metadata(m, n, k):
    require(n in domain(m, k), 'pair outside complete domain')
    order = lcm(2 * m, 2 * n, 2 * k)
    degree = phi(order)
    require(order <= MAX_L and degree <= MAX_PHI, 'proved order/degree bounds')
    return {'kind': 'pair', 'm': m, 'n': n, 'k': k,
            'anchor_degree': DEGREES[m], 'mutual_degree': 1 if k == 3 else 2,
            'other_totient': phi(2 * n), 'order': order, 'cyclotomic_degree': degree}


def pair_record(m, n, k):
    record = metadata(m, n, k)
    order = record['order']
    modulus = cyclotomic(order)
    sparse = sparse_expression(m, n, k, order)
    require(sum(sparse.values()) == -32, 'independent E(2,2,2) augmentation')
    dense = [0] * (max(sparse, default=0) + 1)
    for exponent, coefficient in sparse.items():
        require(abs(coefficient) < 1024, 'sparse input coefficient bound')
        dense[exponent] = coefficient
    _, remainder = divide_monic(dense, modulus, False)
    nonzero = [(i, a) for i, a in enumerate(remainder) if a]
    record.update({'result': 'EXCLUDED' if nonzero else 'SURVIVOR',
                   'cyclotomic_sha256': polynomial_hash(modulus),
                   'remainder_sha256': polynomial_hash(remainder),
                   'remainder_degree': len(remainder) - 1 if nonzero else -1,
                   'remainder_nonzero_terms': len(nonzero),
                   'max_coefficient_bits': max(a.bit_length() for a in remainder),
                   'witness_index': nonzero[0][0] if nonzero else None,
                   'witness_coefficient_hex': signed_hex(nonzero[0][1]) if nonzero else None})
    return record


def arithmetic_fixtures():
    require(tuple(phi(2 * m) // 2 for m in OWN_M) == tuple(DEGREES[m] for m in OWN_M),
            'anchor degree sequence')
    require(phi(210) == 48 and phi(420) == 96
            and all(phi(2 * n) > 48 for n in range(106, 2305)), 'complete even cutoff fixture')
    known = {1: (-1, 1), 2: (1, 1), 3: (1, 1, 1), 4: (1, 0, 1),
             6: (1, -1, 1), 8: (1, 0, 0, 0, 1), 12: (1, 0, -1, 0, 1),
             15: (1, -1, 0, 1, -1, 1, 0, -1, 1)}
    require(all(cyclotomic(n) == poly for n, poly in known.items()), 'exact cyclotomic fixtures')
    require(divide_monic((-1, 0, 0, 1), (-1, 1), True)[0] == (1, 1, 1),
            'exact monic quotient fixture')
    try:
        divide_monic((1, 0, 1), (-1, 1), True)
    except ValueError:
        pass
    else:
        raise ValueError('nonexact division fixture was accepted')
    require(cyclic_multiply({0: 1, 3: -1}, {1: 2}, 4) == {1: 2, 0: -2},
            'cyclic negative-coefficient wrap fixture')
    zero = sparse_expression(3, 3, 3, 6)
    dense = [zero.get(i, 0) for i in range(6)]
    require(divide_monic(dense, known[6], False)[1] == (0,), 'zero remainder fixture')
    require(divide_monic((0, 1), known[4], False)[1] == (0, 1), 'nonzero remainder fixture')
    return list(FIXTURES)


class Writer:
    def __init__(self, stream):
        self.stream, self.bytes = stream, 0

    def emit(self, record):
        line = (canonical(record) + '\n').encode('utf-8')
        require(len(line) <= MAX_LINE and self.bytes + len(line) <= MAX_FILE,
                'bounded certificate line/file size')
        self.stream.write(line)
        self.stream.flush()
        self.bytes += len(line)


def batch_header(m, k, source_sha):
    return {'kind': 'header', 'schema': 1, 'unit': 'batch', 'scope': SCOPE,
            'm': m, 'k': k, 'n_domain': list(domain(m, k)),
            'source_sha256': source_sha, 'fixed_pins': PINS,
            'limits': {'max_n': MAX_N, 'max_order': MAX_L, 'max_phi': MAX_PHI,
                       'max_coefficient_bits': MAX_BITS}}


def totals():
    return {'pairs': 0, 'excluded_nonzero': 0, 'survivors': 0}


def count(total, record):
    total['pairs'] += 1
    total['excluded_nonzero' if record['result'] == 'EXCLUDED' else 'survivors'] += 1


def batch_summary(m, k, total):
    require(total['pairs'] == len(domain(m, k))
            and total['excluded_nonzero'] + total['survivors'] == total['pairs'],
            'complete exact pair partition')
    return {'kind': 'summary', 'unit': 'batch', 'status': 'PASS', 'scope': SCOPE,
            'm': m, 'k': k, 'arithmetic_fixtures': list(FIXTURES), **total}


def seal(path, source_sha, summary):
    return {'kind': 'output_seal', 'status': 'PASS', 'scope': SCOPE,
            'source_sha256': source_sha, 'output': {'path': relative(path),
            'bytes': path.stat().st_size, 'sha256': digest(path)}, 'summary': summary}


def batch(args):
    path = result_path(args.output)
    total = totals()
    with path.open('xb') as stream:
        writer = Writer(stream)
        writer.emit(batch_header(args.m, args.k, args.expected_source_sha256))
        for n in domain(args.m, args.k):
            record = pair_record(args.m, n, args.k)
            writer.emit(record)
            count(total, record)
        final = batch_summary(args.m, args.k, total)
        validate_pins(args.expected_source_sha256)
        writer.emit(final)
    validate_pins(args.expected_source_sha256)
    print(canonical(seal(path, args.expected_source_sha256, final)))


def unique_object(pairs):
    result = {}
    for key, value in pairs:
        require(key not in result, 'duplicate JSON key')
        result[key] = value
    return result


def invalid_number(value):
    raise ValueError('noninteger/nonfinite certificate number: ' + value)


def read_record(stream):
    line = stream.readline(MAX_LINE + 1)
    require(line and len(line) <= MAX_LINE and line.endswith(b'\n'),
            'missing/oversized/truncated JSONL record')
    record = json.loads(line, object_pairs_hook=unique_object,
                        parse_float=invalid_number, parse_constant=invalid_number)
    require(type(record) is dict, 'JSONL object required')
    return record


def match(actual, expected):
    require(canonical(actual) == canonical(expected), 'exact record/schema/domain mismatch')


def validate_pair(record, m, n, k):
    fixed = metadata(m, n, k)
    extra = {'result', 'cyclotomic_sha256', 'remainder_sha256', 'remainder_degree',
             'remainder_nonzero_terms', 'max_coefficient_bits', 'witness_index',
             'witness_coefficient_hex'}
    require(set(record) == set(fixed) | extra, 'exact pair record fields')
    match({key: record[key] for key in fixed}, fixed)
    require(valid_sha(record['cyclotomic_sha256']) and valid_sha(record['remainder_sha256']),
            'pair polynomial hashes')
    degree, terms, bits = (record[key] for key in
                           ('remainder_degree', 'remainder_nonzero_terms', 'max_coefficient_bits'))
    require(all(type(x) is int for x in (degree, terms, bits)), 'integer remainder schema')
    if record['result'] == 'SURVIVOR':
        require((degree, terms, bits, record['witness_index'], record['witness_coefficient_hex'])
                == (-1, 0, 0, None, None)
                and record['remainder_sha256'] == polynomial_hash((0,)), 'zero remainder schema')
    else:
        require(record['result'] == 'EXCLUDED' and 0 <= degree < fixed['cyclotomic_degree']
                and 1 <= terms <= degree + 1 and 1 <= bits <= MAX_BITS, 'nonzero remainder schema')
        index, encoded = record['witness_index'], record['witness_coefficient_hex']
        require(type(index) is int and 0 <= index <= degree and type(encoded) is str
                and 2 <= len(encoded) <= MAX_BITS // 4 + 3
                and encoded[0] in '+-' and all(c in '0123456789abcdef' for c in encoded[1:]),
                'exact coefficient witness schema')
        coefficient = int(encoded, 16)
        require(coefficient != 0 and signed_hex(coefficient) == encoded
                and coefficient.bit_length() <= bits
                and terms <= degree - index + 1, 'canonical nonzero coefficient witness')


def load_ledger(path):
    require(0 < path.stat().st_size <= MAX_LEDGER, 'bounded ledger file')
    with path.open('rb') as stream:
        return json.loads(stream.read(), object_pairs_hook=unique_object,
                          parse_constant=invalid_number)


def validate_ledger(ledger, path, m, k, source_sha, summary):
    command = [SOURCE_REL, '--expected-source-sha256', source_sha, 'batch',
               '--m', str(m), '--k', str(k), '--output', relative(path)]
    require(type(ledger) is dict and ledger.get('status') == 'PASS_BOUNDED_RESEARCH_CHECK'
            and ledger.get('command') == command and ledger.get('checker_sha256') == source_sha
            and ledger.get('runner_sha256') == PINS['scripts/run_bounded_research_check_20261018.py']
            and ledger.get('lock_owner') == 'wrapper'
            and ledger.get('wall_limit_seconds') == 59
            and ledger.get('rss_limit_bytes') == 1024**3 and ledger.get('threads') == 1,
            'actual bounded runner/source/command/lock/caps')
    elapsed, peak = ledger.get('elapsed_seconds'), ledger.get('peak_process_tree_rss_bytes')
    require(type(elapsed) in (int, float) and isfinite(elapsed) and 0 <= elapsed <= 59
            and type(peak) is int and 0 < peak < 1024**3, 'actual bounded elapsed/RSS')
    resources = ledger.get('initial_resources')
    require(type(resources) is dict, 'actual initial resource object')
    for key, floor in (('free_disk_gib', 40), ('available_memory_gib', 6)):
        value = resources.get(key)
        require(type(value) in (int, float) and isfinite(value) and value >= floor,
                'actual resource admission: ' + key)
    stdout = canonical(seal(path, source_sha, summary)) + '\n'
    require(ledger.get('stdout') == stdout
            and ledger.get('stdout_bytes') == len(stdout.encode('utf-8')),
            'actual stdout seal/output hash/terminal summary')


def combine(args):
    require(len(args.part) == len(args.part_sha256) == len(args.ledger)
            == len(args.ledger_sha256) == len(UNITS), 'exactly sixteen complete bindings required')
    paths = [result_path(p) for p in args.part]
    ledgers = [result_path(p) for p in args.ledger]
    output = result_path(args.output)
    require(len(set(paths + ledgers + [output])) == 33, 'duplicate input/output path')
    bindings = []
    for path, sha, ledger, ledger_sha in zip(paths, args.part_sha256, ledgers, args.ledger_sha256):
        require(valid_sha(sha) and valid_sha(ledger_sha)
                and 0 < path.stat().st_size <= MAX_FILE
                and 0 < ledger.stat().st_size <= MAX_LEDGER
                and digest(path) == sha and digest(ledger) == ledger_sha, 'bound part/ledger hashes')
        bindings.append({'path': relative(path), 'bytes': path.stat().st_size, 'sha256': sha,
                         'ledger': {'path': relative(ledger), 'bytes': ledger.stat().st_size,
                                    'sha256': ledger_sha}})
    total, seen = totals(), set()
    with output.open('xb') as stream:
        writer = Writer(stream)
        writer.emit({'kind': 'header', 'schema': 1, 'unit': 'combine', 'scope': SCOPE,
                     'source_sha256': args.expected_source_sha256, 'fixed_pins': PINS,
                     'meaning': 'authenticated_coverage_of_sixteen_actual_batches_no_math_rerun'})
        for path, ledger_path, binding in zip(paths, ledgers, bindings):
            subtotal, survivors = totals(), []
            with path.open('rb') as part:
                header = read_record(part)
                m, k = header.get('m'), header.get('k')
                require(type(m) is int and type(k) is int and (m, k) in UNITS
                        and (m, k) not in seen, 'invalid/duplicate batch owner')
                match(header, batch_header(m, k, args.expected_source_sha256))
                for n in domain(m, k):
                    record = read_record(part)
                    validate_pair(record, m, n, k)
                    count(subtotal, record)
                    if record['result'] == 'SURVIVOR':
                        survivors.append(record)
                final = batch_summary(m, k, subtotal)
                match(read_record(part), final)
                require(part.read(1) == b'', 'extra bytes after batch summary')
            validate_pins(args.expected_source_sha256)
            validate_ledger(load_ledger(ledger_path), path, m, k,
                            args.expected_source_sha256, final)
            require(digest(path) == binding['sha256'] and path.stat().st_size == binding['bytes']
                    and digest(ledger_path) == binding['ledger']['sha256']
                    and ledger_path.stat().st_size == binding['ledger']['bytes'],
                    'part/ledger changed during serial read')
            seen.add((m, k))
            writer.emit({'kind': 'part_binding', 'm': m, 'k': k, **binding, 'summary': final})
            for record in survivors:
                writer.emit(record)
            for key in total:
                total[key] += subtotal[key]
        require(seen == set(UNITS)
                and total['pairs'] == sum(len(domain(m, k)) for m, k in UNITS),
                'complete sixteen-unit domain conjunction')
        validate_pins(args.expected_source_sha256)
        final = {'kind': 'summary', 'unit': 'combine', 'status': 'PASS', 'scope': SCOPE,
                 'coverage_only': True, 'batches': 16, **total}
        writer.emit(final)
    validate_pins(args.expected_source_sha256)
    print(canonical(seal(output, args.expected_source_sha256, final)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256', required=True)
    sub = parser.add_subparsers(dest='unit', required=True)
    one = sub.add_parser('batch')
    one.add_argument('--m', type=int, choices=OWN_M, required=True)
    one.add_argument('--k', type=int, choices=MUTUAL_K, required=True)
    one.add_argument('--output', required=True)
    union = sub.add_parser('combine')
    union.add_argument('--part', action='append', required=True)
    union.add_argument('--part-sha256', action='append', required=True)
    union.add_argument('--ledger', action='append', required=True)
    union.add_argument('--ledger-sha256', action='append', required=True)
    union.add_argument('--output', required=True)
    args = parser.parse_args()
    validate_pins(args.expected_source_sha256)
    arithmetic_fixtures()
    if args.unit == 'batch':
        batch(args)
    else:
        combine(args)


if __name__ == '__main__':
    main()
