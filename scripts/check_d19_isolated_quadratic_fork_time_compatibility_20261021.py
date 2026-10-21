#!/usr/bin/env python3
"""Exact finite necessary compatibility census; root-owned guarded execution only.

No scientific packages, floating arithmetic, external data tables or imports of
repository workers. This file is source preparation until a guarded run occurs.
An independent mathematical inventory note is a REQUIRED hashed input, supplied
by --inventory/--inventory-sha. Hashing a note does not prove its mathematics.

Domain and termination:
* Fork B has minimal polynomial X^2-T X+N, T=5..53, N=1..22. Keep precisely
  nonsquare positive discriminants and 4<B<53, 0<Bbar<8/19. These bounds follow
  from 59e6: trace<53+8/19 and norm<53*8/19. All 1078 integer rows are visited.
* Isolated-prefix B_i has trace=5..16, norm=1..6, with 4<B_i<16 and
  0<B_i_bar<2/i, i=5..16. All 864 integer/index rows are visited. Other prefix
  candidates are the explicit finite silver/golden lists below, including the
  fresh D3 silver dotted-fork row16. Golden right-fork rows12/13 are already
  geometrically strong; the norm-pass/square-fail regression12 below is separate
  from this census. A surviving pair does not prove that a wall exists.
* A discriminant D is written uniquely as q^2*d with d>1 squarefree by bounded
  trial division. This gives the physically embedded root (T+q*sqrt(d))/2 and
  exact same-field tests. No numerical root identification or field census.
* Arithmetic is in Q(sqrt(d)), with rational pairs, positive sqrt(d), and exact
  signs by comparing a^2 with d*b^2. A ratio a+b*sqrt(d), b!=0, is a field square
  iff its norm is a nonnegative rational square n^2 and one of (a+n)/2,(a-n)/2
  is a nonzero rational square x^2 with y=b/(2x) and (x+y*sqrt(d))^2 equal to the
  ratio. For b=0, test a and a/d for rational squares. Every positive result
  includes its exact verified square root. Every negative branch is exhaustive.

Actual geometry, completeness of the supplied weak inventory and the same-vector
E16(.48)>1.009327 certificate remain explicit mathematical inputs. The output is
a finite necessary compatibility result, not an ambient polytope construction,
not a whole-D19 exclusion and not by itself a new dimension bound. All candidate
and accepted sources stay read-only; the required JSON destination must be new.
"""

import argparse
from dataclasses import dataclass
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PINS = {
    'results/arithmetic_d19_isolated_dotted_fork_continuation_20261021.md':
        '59e6c68bd60bce749ac7da60e0ba80ad82f91e0aa431dbfb297369125871fa5e',
    'results/arithmetic_d19_golden_near_fork_time_anchor_exclusion_20261021.md':
        '635ea0aa09efcfbfdd097c625bf79c341ffda724594c326d5db654e36db50f7c',
    'results/arithmetic_d19_silver_centre_actual_prefix_time_exclusion_20261021.md':
        '7d9f4facaa3c622ba41a03d570a54b39b9ac0c3fadc3d3c13104798fa4914fcc',
}


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def rational_square_root(value):
    value = F(value)
    if value < 0:
        return None
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    if a*a == value.numerator and b*b == value.denominator:
        return F(a, b)
    return None


def squarefree_decomposition(value):
    original, factor, scale, squarefree = value, 2, 1, 1
    while factor*factor <= value:
        exponent = 0
        while value % factor == 0:
            value //= factor
            exponent += 1
        scale *= factor**(exponent//2)
        if exponent % 2:
            squarefree *= factor
        factor += 1
    if value > 1:
        squarefree *= value
    assert squarefree*scale*scale == original and squarefree > 1
    return squarefree, scale


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
            assert self.d == other.d
            return other
        return Q(self.d, F(other))

    def __add__(self, other):
        other = self.lift(other)
        return Q(self.d, self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Q(self.d, -self.a, -self.b)

    def __sub__(self, other):
        return self + (-self.lift(other))

    def __rsub__(self, other):
        return self.lift(other) - self

    def __mul__(self, other):
        other = self.lift(other)
        return Q(self.d, self.a*other.a+self.d*self.b*other.b,
                 self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def conjugate(self):
        return Q(self.d, self.a, -self.b)

    def norm(self):
        return self.a*self.a-self.d*self.b*self.b

    def __truediv__(self, other):
        other = self.lift(other)
        n = other.norm()
        assert n != 0
        product = self*other.conjugate()
        return Q(self.d, product.a/n, product.b/n)

    def sign(self):
        if self.b == 0:
            return (self.a > 0)-(self.a < 0)
        if self.a == 0:
            return (self.b > 0)-(self.b < 0)
        if self.a*self.b > 0:
            return (self.a > 0)-(self.a < 0)
        delta = self.norm()
        assert delta != 0
        return ((self.a > 0)-(self.a < 0))*((delta > 0)-(delta < 0))

    def encode(self):
        return {'squarefree_d': self.d,
                'rational': [self.a.numerator, self.a.denominator],
                'radical': [self.b.numerator, self.b.denominator]}


def positive_root(trace, norm):
    discriminant = trace*trace-4*norm
    if discriminant <= 0 or isqrt(discriminant)**2 == discriminant:
        return None
    d, scale = squarefree_decomposition(discriminant)
    root = Q(d, F(trace, 2), F(scale, 2))
    assert root*root-trace*root+norm == Q(d, 0)
    return root


def square_in_field(value):
    if value.b == 0:
        rational = rational_square_root(value.a)
        if rational is not None:
            witness = Q(value.d, rational)
            assert witness*witness == value
            return witness, 'rational_root'
        radical = rational_square_root(value.a/value.d)
        if radical is not None:
            witness = Q(value.d, 0, radical)
            assert witness*witness == value
            return witness, 'pure_radical_root'
        return None, 'rational_and_radical_nonsquares'
    n = rational_square_root(value.norm())
    if n is None:
        return None, 'norm_not_a_rational_square'
    for signed_norm in (n, -n):
        x = rational_square_root((value.a+signed_norm)/2)
        if x is not None and x != 0:
            witness = Q(value.d, x, value.b/(2*x))
            if witness*witness == value:
                return witness, 'mixed_root'
    return None, 'both_rational_half_sums_nonsquares'


def arithmetic_regressions():
    # A genuine same-norm-class pair: rational norms alone do not exclude it.
    phi = Q(5, F(1, 2), F(1, 2))
    fork_b = 2+3*phi
    fork_time = fork_b*(19*fork_b-8)/4
    prefix_time = (59+96*phi)/4
    ratio = prefix_time/fork_time
    assert ratio == Q(5, F(49, 142), F(-17, 426))
    assert ratio.norm() == F(1, 9)
    assert square_in_field(ratio) == (None, 'both_rational_half_sums_nonsquares')
    assert square_in_field(Q(2, F(4, 3)))[0] is None
    assert square_in_field(Q(2, 2))[0] == Q(2, 0, 1)
    mixed = Q(2, 1, 1)
    witness, _ = square_in_field(mixed*mixed)
    assert witness is not None and witness*witness == mixed*mixed
    assert square_in_field(Q(2, -1))[0] is None
    return ['phi4_fork_vs_native5_fork12_norm_pass_field_square_fail',
            'rational_nonsquare_in_quadratic_field',
            'pure_radical_square', 'mixed_square', 'negative_nonsquare']


def anchors():
    output = []
    for trace in range(5, 54):
        for norm in range(1, 23):
            b = positive_root(trace, norm)
            if b is None:
                continue
            if not ((b-4).sign() > 0 and (53-b).sign() > 0
                    and b.conjugate().sign() > 0
                    and (F(8, 19)-b.conjugate()).sign() > 0):
                continue
            time = b*(19*b-8)/4
            assert time.sign() > 0 and time.conjugate().sign() < 0
            output.append({'trace': trace, 'norm': norm, 'B': b, 'time': time})
    return output


def prefix_candidates():
    output = []

    def add(family, index, time, parameter=None):
        assert 5 <= index <= 16
        assert time.sign() > 0 and time.conjugate().sign() < 0
        output.append({'family': family, 'index': index, 'time': time,
                       'parameter': parameter})

    add('silver_near_finite4', 5, Q(2, 14, 10))
    for index in range(5, 11):
        add('silver_far_dotted3', index, Q(2, 3*index+1, 2*index+2))
    add('silver_fork_dotted3_D3', 16, Q(2, F(267, 4), 48))
    phi = Q(5, F(1, 2), F(1, 2))
    for index in (5, 6):
        add('golden_far_D_position', index, index+(index+4)*phi)
    add('golden_near_D_position', 5, 6+11*phi)
    add('golden_far_helper3', 5, 4+7*phi)
    for index in range(14, 17):
        add('golden_right_fork3', index, (4*index+11+8*index*phi)/4)
    for index in range(5, 17):
        for trace in range(5, 17):
            for norm in range(1, 7):
                b = positive_root(trace, norm)
                if b is None:
                    continue
                if ((b-4).sign() > 0 and (16-b).sign() > 0
                        and b.conjugate().sign() > 0
                        and (F(2, index)-b.conjugate()).sign() > 0):
                    add('isolated_dotted', index, b*(index*b-2),
                        {'trace': trace, 'norm': norm, 'B': b.encode()})
    return output


def census():
    regressions = arithmetic_regressions()
    anchor_rows, prefix_rows = anchors(), prefix_candidates()
    survivors, zero_match, square_tests, same_field_pairs = [], [], 0, 0
    by_family, square_results = {}, {}
    for row in prefix_rows:
        by_family[row['family']] = by_family.get(row['family'], 0)+1
    for anchor in anchor_rows:
        matches = []
        for prefix in prefix_rows:
            if prefix['time'].d != anchor['time'].d:
                continue
            same_field_pairs += 1
            ratio = prefix['time']/anchor['time']
            assert ratio.sign() > 0 and ratio.conjugate().sign() > 0
            square_tests += 1
            witness, reason = square_in_field(ratio)
            square_results[reason] = square_results.get(reason, 0)+1
            if witness is not None:
                matches.append({'family': prefix['family'], 'index': prefix['index'],
                    'prefix_parameter': prefix['parameter'],
                    'prefix_time': prefix['time'].encode(), 'ratio': ratio.encode(),
                    'verified_square_root': witness.encode(), 'branch': reason})
        base = {'trace': anchor['trace'], 'norm': anchor['norm'],
                'B': anchor['B'].encode(), 'fork_time': anchor['time'].encode()}
        if matches:
            survivors.append(dict(base, compatible_weak_candidates=matches))
        else:
            zero_match.append(base)
    return {
        'status': 'COMPLETE_EXACT_NECESSARY_COMPATIBILITY_CENSUS',
        'scope': 'Algebraic local-data compatibility only; no wall or polytope realization.',
        'executed_exact_arithmetic_regressions': regressions,
        'anchor_integer_domain_rows': 49*22,
        'isolated_prefix_integer_index_domain_rows': 12*12*6,
        'anchor_models': len(anchor_rows), 'prefix_candidate_rows': len(prefix_rows),
        'prefix_rows_by_family': by_family,
        'same_field_pairs': same_field_pairs, 'exact_field_square_tests': square_tests,
        'field_square_result_counts': square_results,
        'no_compatible_pair': not survivors,
        'remaining_anchor_models': len(survivors),
        'remaining_compatible_pairs': sum(len(r['compatible_weak_candidates']) for r in survivors),
        'remaining_necessary_models_not_geometric_realizations': survivors,
        'zero_match_anchor_models': zero_match,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inventory', required=True,
                        help='Independent complete weak-inventory note under results/')
    parser.add_argument('--inventory-sha', required=True)
    parser.add_argument('--out', required=True, help='New JSON file directly under results/')
    args = parser.parse_args()
    inventory = (ROOT/args.inventory).resolve()
    assert inventory.parent == ROOT/'results' and inventory.suffix == '.md'
    assert args.inventory == str(inventory.relative_to(ROOT))
    assert len(args.inventory_sha) == 64 and all(c in '0123456789abcdef' for c in args.inventory_sha)
    destination = (ROOT/args.out).resolve()
    assert destination.parent == ROOT/'results' and destination.suffix == '.json'
    assert not destination.exists(), 'Output must be new'
    pins = dict(PINS)
    assert args.inventory not in pins
    pins[args.inventory] = args.inventory_sha
    for name, pin in pins.items():
        assert digest(ROOT/name) == pin, 'Input changed: '+name
    original = digest(Path(__file__))
    output = census()
    output['checker_sha256'] = original
    output['mathematical_input_pins'] = pins
    assert digest(Path(__file__)) == original
    for name, pin in pins.items():
        assert digest(ROOT/name) == pin, 'Input changed: '+name
    with destination.open('x') as stream:
        json.dump(output, stream, indent=2, sort_keys=True)
        stream.write('\n')
    summary = {key: value for key, value in output.items()
               if key not in ('remaining_necessary_models_not_geometric_realizations',
                              'zero_match_anchor_models', 'mathematical_input_pins')}
    print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    main()
