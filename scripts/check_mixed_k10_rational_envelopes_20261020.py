#!/usr/bin/env python3
"""Exact Bernstein check of three written k10 envelopes; no geometry census.

Owns no flock: run only inside the parent's bounded serial wrapper.
--write exclusively creates a certificate; --verify is read-only.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path
import signal

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / 'results/right_angled_mixed_k10_rational_envelopes_20261020.json'
PINS = {
    'results/right_angled_mixed_k10_subset_moment_reduction_20261020.md':
        '8c2051b3be0f692e5097eedd81cf494f15313037d36b0f6625ca78cd01aabd34',
    'results/right_angled_mixed_k10_rational_envelopes_20261020.md':
        '7fdcad68b9b4ae113c482e9e17266a6f5e9aea0d21409eb7db332ec915627e45',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(*polys):
    out = [F(0)] * max(map(len, polys))
    for p in polys:
        for i, value in enumerate(p):
            out[i] += value
    return trim(out)


def scale(p, value):
    return trim([F(value) * a for a in p])


def mul(p, q):
    out = [F(0)] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return trim(out)


def power(p, n):
    out = [F(1)]
    for _ in range(n):
        out = mul(out, p)
    return out


def rational(value):
    return [value.numerator, value.denominator]


def bernstein(p, left, right):
    """Power expansion after v=left+(right-left)u, then exact Bernstein."""
    n = len(p) - 1
    width = right - left
    affine = [sum((p[l] * comb(l, j) * left ** (l - j)
                   for l in range(j, n + 1)), F(0)) * width ** j
              for j in range(n + 1)]
    return [sum((affine[j] * F(comb(i, j), comb(n, j))
                 for j in range(i + 1)), F(0)) for i in range(n + 1)]


def constants():
    # Exact rational implications of the square-root comparisons in the proof.
    assert F(227, 100) ** 2 > 5 and F(47, 20) ** 2 > 5
    assert F(29, 12) ** 2 > 5  # v_min > -1/6.
    assert 139 ** 2 > 3 * 80 ** 2  # global a > 11/20.
    assert F(23, 5) ** 2 > 21  # a(0) > 16/25.
    kneg = 4 / (F(21, 20) ** 5 * F(16, 25))
    kneg += 1 / (F(21, 20) ** 4 * F(16, 25) ** 2)
    kpos = 4 / (F(11, 10) ** 5 * F(11, 20))
    kpos += 1 / (F(11, 10) ** 4 * F(11, 20) ** 2)
    assert kneg == F(7220000000, 1045529856) < 7
    assert kpos == F(12000000, 1771561) < 7
    assert F(43, 300) - F(32, 45) / 15 == F(259, 2700) > 0
    assert F(32, 45) < F(17, 20) ** 2
    assert F(8, 15) * F(11, 10) ** 5 < 1
    assert F(21, 10) * F(509, 4000) < F(13, 25) ** 2
    assert F(223, 200) + F(13, 25) < F(41, 25)
    assert F(100, 37) > F(41, 25) ** 2
    assert F(37, 100) ** 3 < F(113, 500) ** 2
    assert F(41, 25) - F(131, 100) > 0
    return {'Knegative': rational(kneg), 'Kpositive': rational(kpos),
            'quartic_joint_lower': rational(F(259, 2700))}


def polynomials():
    mu = [F(11, 10), F(3, 10)]
    var = [F(1, 10), F(3, 5), F(-11, 10)]
    v2 = power(var, 2)
    neg = add(scale(power(mu, 5), F(223, 250)),
              scale(power(mu, 3), -1),
              scale(mul(var, mu), F(-7, 10)),
              scale(var, F(-119, 3000)), scale(v2, F(-119, 120)),
              scale(mul(v2, power(mu, 5)), F(-1813, 2700)))
    zero = add(scale(power(mu, 10), F(223, 250)),
               scale(power(mu, 8), -1),
               scale(mul(var, power(mu, 6)), F(-7, 10)),
               scale(mul(v2, power(mu, 10)), F(-301, 300)),
               scale(var, F(-7, 60)))
    gap = add([F(41, 25)], scale(mu, -1))
    corr = add(mul(gap, zero),
               scale(mul(var, power(mu, 10)), F(-791, 30000)))
    assert [len(p) - 1 for p in [neg, zero, corr]] == [9, 14, 15]
    return [('negative', neg, -140, 0), ('positive_zero', zero, 0, 50),
            ('positive_corrected', corr, 50, 700)]


def certify(name, p, start, stop):
    h = hashlib.sha256()
    floors = []
    minimum = None
    position = None
    total = 0
    for i in range(start, stop):
        coefficients = bernstein(p, F(i, 1000), F(i + 1, 1000))
        for j, value in enumerate(coefficients):
            assert value > 0, (name, i, j, rational(value))
            h.update(f'{value.numerator}/{value.denominator}\n'.encode('ascii'))
            if minimum is None or value < minimum:
                minimum, position = value, [i, j]
            total += 1
        floor = min(coefficients) * 10 ** 18
        floor = floor.numerator // floor.denominator
        assert floor > 0
        floors.append(floor)
    return {'id': name, 'degree': len(p) - 1,
            'power_coefficients': [rational(a) for a in p],
            'grid_denominator': 1000, 'start': start, 'stop': stop,
            'intervals': stop - start, 'bernstein_coefficients': total,
            'minimum': rational(minimum), 'minimum_position': position,
            'all_coefficients_sha256': h.hexdigest(),
            'interval_lower_bound_denominator': 10 ** 18,
            'interval_lower_bound_numerators': floors}


def run():
    pins = []
    for relative, expected in PINS.items():
        path = ROOT / relative
        actual = digest(path)
        assert actual == expected, relative
        pins.append({'path': relative, 'sha256': actual,
                     'bytes': path.stat().st_size})
    source = Path(__file__).resolve()
    result = {'status': 'PASS', 'scope': 'three rational k10 envelope polynomials',
              'geometry_checked': False, 'dimension_claim': False,
              'owns_lock': False, 'inputs': pins,
              'checker_sha256': digest(source), 'constants': constants(),
              'certificates': [certify(*row) for row in polynomials()]}
    return json.loads(json.dumps(result, sort_keys=True))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write', action='store_true')
    mode.add_argument('--verify', action='store_true')
    parser.add_argument('--certificate', type=Path, default=DEFAULT)
    args = parser.parse_args()
    signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(TimeoutError()))
    signal.alarm(45)
    data = run()
    if args.write:
        with args.certificate.open('x', encoding='utf-8') as handle:
            json.dump(data, handle, sort_keys=True, separators=(',', ':'))
            handle.write('\n')
    else:
        with args.certificate.open(encoding='utf-8') as handle:
            assert data == json.load(handle), 'certificate mismatch'
    rows = data['certificates']
    print(json.dumps({'status': 'PASS', 'mode': 'write' if args.write else 'verify',
                      'intervals': sum(r['intervals'] for r in rows),
                      'coefficients': sum(r['bernstein_coefficients'] for r in rows),
                      'certificate': str(args.certificate)},
                     sort_keys=True))


if __name__ == '__main__':
    main()
