#!/usr/bin/env python3
"""Exact, streaming prime2099/k6 transverse-edge certificate.

Prepared for source review; not executed during preparation. Source preparation
is not verification. Run only through the pinned parent-lock bounded wrapper.
No scientific imports, corpus loading, Pell cutoff or row-order winner rule.
The actual starting germ, complete finite simple-profile theorem and integral
arrangement/germ principle are explicit PINNED mathematical inputs, not proved
by this worker. A period conclusion additionally requires a certified unique
simple germ with exactly the displayed upper mirror.

The certificate is JSONL: header, each admitted mirror, final result. A missing
final result or a failure record is not a completed certificate. Storage is
constant in the scalar-domain size: only counters and at most five minimum
ties are retained. Every input/source pin is checked before and after.
"""

import argparse
import fcntl
from fractions import Fraction as Q
import hashlib
import json
from math import isqrt
import os
from pathlib import Path
import resource
import signal
import sys


ROOT = Path(__file__).resolve().parents[1]
SPEC = 'results/compact_h9_prime2099_first_edge_checker_spec_20261021.md'
PINS = {
    'results/compact_h9_next_host_after_prime811_continuation_20261021.md':
        '1fe722198c81088723975a88e2ab5a1e721f94f130ddc0376554070c13b34589',
    'scripts/check_compact_h9_prime811_first_edge_20261021.py':
        'b72f3274cd5cb953f424d24bf358e22b355cc7ff005c02a0edcff4c6597174d0',
    'results/compact_h9_prime811_literal_host_obstruction_20261021.md':
        '6796ff613364d043538f3a98b5770c5053557334ab347cf787c78589e27597d9',
    'results/golden_e8_binary_family_full_base_facet_20261020.md':
        'ea5d96e331ac91cb9efc457db2ddaa24bac3341e9e6781038150765895675afd',
    'results/golden_e8_binary_family_uniform_first_edge_20261020.md':
        '71cfa11bcadf5fffae017d712d168eeec44eecf7d5b74f595bcebca7ad92af2c',
    'results/golden_e8_binary_family_discriminant_collar_and_finite_cutoff_20261020.md':
        '841a4b7f9e75cdd49b6384920439d998cc401f00bbf75aed92444892c851d425',
    'results/compact_h9_golden_k4_transverse_edge_and_period_candidate_20261020.md':
        '7d12748c606ecc9aad469afd9fbf3b4e6a1709c70c2d1b90c7658a14743fb092',
    'results/bugaenko_facet_binary_extension_and_collar_audit_20261020.md':
        '774e87f7abb33680a06879a17a5be9f06e2218bf12a95023bcff5c7046ba1227',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}


class F:
    """Q(phi), phi^2=phi+1; exact rational coefficients and physical signs."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a, self.b = Q(a), Q(b)

    @staticmethod
    def cast(value):
        return value if isinstance(value, F) else F(value)

    def __add__(self, other):
        other = F.cast(other)
        return F(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b)

    def __sub__(self, other):
        return self + -F.cast(other)

    def __rsub__(self, other):
        return F.cast(other) + -self

    def __mul__(self, other):
        other = F.cast(other)
        return F(self.a*other.a+self.b*other.b,
                 self.a*other.b+self.b*other.a+self.b*other.b)

    __rmul__ = __mul__

    def conjugate(self):
        return F(self.a+self.b, -self.b)

    def norm(self):
        return self.a*self.a+self.a*self.b-self.b*self.b

    def __truediv__(self, other):
        other = F.cast(other)
        n = other.norm()
        assert n, 'division by zero'
        v = self*other.conjugate()
        return F(v.a/n, v.b/n)

    def __rtruediv__(self, other):
        return F.cast(other)/self

    def __eq__(self, other):
        other = F.cast(other)
        return self.a == other.a and self.b == other.b

    def sign(self):
        # Twice the physical element is u+v*sqrt(5).
        u, v = 2*self.a+self.b, self.b
        if not v:
            return (u > 0)-(u < 0)
        if not u or (u > 0) == (v > 0):
            return (v > 0)-(v < 0)
        difference = u*u-5*v*v
        assert difference, 'sqrt(5) is irrational'
        return ((u > 0)-(u < 0))*((difference > 0)-(difference < 0))

    def integral(self):
        return self.a.denominator == self.b.denominator == 1

    def encode(self):
        return [str(self.a), str(self.b)]


ZERO, ONE, PHI = F(), F(1), F(0, 1)
U, A, B, D = F(5, 8), F(8, 13), F(148, 233), F(233, 381)
XB, VB = F(679607, 1099627), F(37729, 61046)
TAU_B = PHI*XB/VB
EDGES = ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7))


def vector(values):
    assert len(values) == 10
    return tuple(F.cast(x) for x in values)


def add(v, w):
    return tuple(x+y for x, y in zip(v, w))


def scale(k, v):
    return tuple(k*x for x in v)


def subtract(v, w):
    return add(v, scale(-1, w))


def dot(v, w):
    value = -2*PHI*v[0]*w[0]+A*(v[0]*w[1]+v[1]*w[0])+2*v[1]*w[1]
    value += sum((2*v[2+j]*w[2+j] for j in range(8)), ZERO)
    for i, j in EDGES:
        value -= v[2+i]*w[2+j]+v[2+j]*w[2+i]
    return value


NATIVE = tuple(vector([int(i == j) for i in range(10)]) for j in range(10))
C = scale(-1, NATIVE[1])
P1 = vector([2, -A]+[0]*8)
OM7 = vector([0, 0, 2, 3, 4, 6, 5, 4, 3, 2])
OM6 = vector([0, 0, 4, 6, 8, 12, 10, 8, 6, 3])
RETAINED = tuple(scale(-1, NATIVE[2+j]) for j in range(7))+(C,)
RELEASED = scale(-1, NATIVE[9])
MB = add(scale(XB, P1), scale(B*VB, OM7))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(65536), b''):
            h.update(chunk)
    return h.hexdigest()


def check_pins(pins):
    for name, pin in pins.items():
        assert digest(ROOT/name) == pin, 'changed pinned input: '+name


def parent_lock():
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 1 and os.getppid() == owner, 'bounded parent wrapper required'
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise AssertionError('parent must retain the shared exclusive lock')
    return owner


def integer_sqrt(q):
    if q < 0 or q.denominator != 1:
        return None
    root = isqrt(q.numerator)
    return root if root*root == q else None


def positive_integral_square_root(w):
    """Complete O-square test, including both norm signs and trace zero."""
    if not w.integral() or w.sign() < 0:
        return None
    if w == ZERO:
        return ZERO
    root_norm = integer_sqrt(w.norm())
    if root_norm is None:
        return None
    trace = 2*w.a+w.b
    for n in (root_norm, -root_norm):
        t = integer_sqrt(trace+2*n)
        b = integer_sqrt((trace-2*n)/5)
        if t is None or b is None:
            continue
        for tt in (t, -t):
            for bb in (b, -b):
                if (tt-bb) % 2:
                    continue
                v = F((tt-bb)//2, bb)
                if v*v == w and v.sign() > 0:
                    return v
    return None


def scalar_domain(max_b, upper, conjugate_bound):
    """Complete thin strip; at most four trial A per B, constant storage."""
    for b in range(max_b+1):
        a0 = (isqrt(5*b*b)-b)//2
        for a in range(a0-1, a0+3):
            z = F(a, b)
            if z.sign() <= 0 or (F(upper)-z).sign() <= 0:
                continue
            conjugate = z.conjugate()
            if (F(conjugate_bound)-conjugate).sign() <= 0:
                continue
            if (F(conjugate_bound)+conjugate).sign() <= 0:
                continue
            yield z


def even(z):
    return z.integral() and int(z.a) % 2 == int(z.b) % 2 == 0


def matrix_solve(matrix, rhs):
    n = len(rhs)
    rows = [list(row)+[value] for row, value in zip(matrix, rhs)]
    for col in range(n):
        pivot = next((i for i in range(col, n) if rows[i][col] != ZERO), None)
        assert pivot is not None, 'singular exact basis'
        rows[col], rows[pivot] = rows[pivot], rows[col]
        divisor = rows[col][col]
        rows[col] = [x/divisor for x in rows[col]]
        for i in range(n):
            if i != col:
                multiplier = rows[i][col]
                rows[i] = [x-multiplier*y for x, y in zip(rows[i], rows[col])]
    return tuple(row[-1] for row in rows)


def positive_pivots(matrix):
    rows = [list(row) for row in matrix]
    pivots = []
    for j in range(len(rows)):
        p = rows[j][j]
        if p.sign() <= 0:
            return None
        pivots.append(p)
        for i in range(j+1, len(rows)):
            for k in range(j+1, len(rows)):
                rows[i][k] -= rows[i][j]*rows[j][k]/p
    return pivots


def linear_combination(coefficients, vectors):
    value = vector([0]*10)
    for coefficient, v in zip(coefficients, vectors):
        value = add(value, scale(coefficient, v))
    return value


def assert_profile(kind, normal, d):
    assert all(x.integral() for x in normal)
    assert dot(normal, normal) == 2*d
    # r_m(z)=z-(z.m/d)m preserves the ORIGINAL integral lattice.
    assert all((dot(e, normal)/d).integral() for e in NATIVE)
    products = tuple(dot(normal, wall) for wall in RETAINED)
    target = [ZERO]*8
    if kind == 'U1':
        target[7] = -ONE
    elif kind == 'Uphi':
        target[7] = -PHI
    elif kind == 'E8':
        target[6] = -ONE
    assert products == tuple(target)
    assert dot(normal, P1).sign() < 0
    assert dot(normal, RELEASED).sign() <= 0


def encode(value):
    if isinstance(value, F):
        return value.encode()
    if isinstance(value, dict):
        return {k: encode(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [encode(v) for v in value]
    return value


class Writer:
    MAX_BYTES = 256*1024**2

    def __init__(self, stream):
        self.stream, self.bytes, self.records = stream, 0, 0

    def write(self, record):
        payload = (json.dumps(encode(record), sort_keys=True, separators=(',', ':'))+'\n').encode()
        self.bytes += len(payload)
        assert self.bytes <= self.MAX_BYTES, 'certificate byte limit'
        self.stream.write(payload)
        self.records += 1


class Census:
    def __init__(self, writer, owner):
        self.writer, self.owner = writer, owner
        self.counts = {kind: {'strip_pairs': 0, 'parity_pass': 0,
                             'field_square': 0, 'admitted_mirrors': 0}
                       for kind in ('U0', 'U1', 'Uphi', 'E8', 'B')}
        self.minimum, self.ties, self.serial = None, [], 0
        self.upper_admitted = False

    def touch(self, kind):
        self.counts[kind]['strip_pairs'] += 1
        if sum(row['strip_pairs'] for row in self.counts.values()) % 256 == 0:
            assert os.getppid() == self.owner, 'bounded owner exited'

    def admit(self, kind, scalar, square, normal, d, tau):
        self.counts[kind]['field_square'] += 1
        assert_profile(kind, normal, d)
        assert tau.sign() > 0
        assert dot(normal, add(P1, scale(tau, OM7))) == ZERO
        if (TAU_B-tau).sign() < 0:
            return
        if kind == 'B' and normal == MB and tau == TAU_B:
            self.upper_admitted = True
        self.counts[kind]['admitted_mirrors'] += 1
        record = {'record': 'mirror', 'id': self.serial, 'profile': kind,
                  'scalar': scalar, 'positive_square_root': square,
                  'normal': normal, 'half_norm': d, 'tau': tau}
        self.serial += 1
        self.writer.write(record)
        if self.minimum is None or (tau-self.minimum).sign() < 0:
            self.minimum, self.ties = tau, [record]
        elif tau == self.minimum:
            self.ties.append(record)
        assert len(self.ties) <= 5, 'at most one mirror per profile at a fixed tau'

    def run(self):
        # E8 is exactly the even-X0/odd-Y0 subset of the U0 norm-1 row.
        for x0 in scalar_domain(50000, 110000, Q(2, 3)):
            self.touch('U0')
            self.counts['U0']['parity_pass'] += 1
            e8_scalar = even(x0)
            if e8_scalar:
                self.touch('E8')
            y0 = positive_integral_square_root(ONE+D*x0*x0)
            if y0 is None:
                continue
            self.admit('U0', x0, y0, add(scale(x0, P1), scale(y0, OM7)),
                       ONE, D*x0/y0)
            if e8_scalar and even(y0-3) and (y0-3).sign() >= 0:
                self.counts['E8']['parity_pass'] += 1
                x, y = x0/2, (y0-3)/2
                normal = add(add(scale(x, P1), OM6), scale(y, OM7))
                self.admit('E8', x, y0, normal, ONE, D*x0/y0)
                coefficients = (2, 3, 4, 6, 5, 4, 3)
                difference = subtract(add(scale(x0, P1), scale(y0, OM7)), scale(2, normal))
                assert difference == linear_combination(coefficients, RETAINED[:7])
        for kind, eta, maximum, upper, bound in (
                ('U1', ONE, 87000, 190000, Q(9, 8)),
                ('Uphi', PHI, 60000, 130000, Q(5, 4))):
            for u in scalar_domain(maximum, upper, bound):
                self.touch(kind)
                if not even(eta-A*u):
                    continue
                self.counts[kind]['parity_pass'] += 1
                y = positive_integral_square_root((D*u*u+4-eta*eta)/4)
                if y is None:
                    continue
                normal = add(add(scale(u/2, P1), scale(-eta/2, C)), scale(y, OM7))
                self.admit(kind, u, y, normal, ONE, D*u/(2*y))
        for v in scalar_domain(63000, 137000, Q(251, 500)):
            self.touch('B')
            self.counts['B']['parity_pass'] += 1
            x = positive_integral_square_root((B*v*v-1)/PHI)
            if x is None:
                continue
            normal = add(scale(x, P1), scale(B*v, OM7))
            self.admit('B', v, x, normal, B, PHI*x/v)
        maxima = {'U0': 100002, 'E8': 25001, 'U1': 261003, 'Uphi': 180003, 'B': 126002}
        for kind, maximum in maxima.items():
            assert self.counts[kind]['strip_pairs'] <= maximum
        assert sum(row['strip_pairs'] for row in self.counts.values()) <= 692011
        assert self.minimum is not None and self.upper_admitted, 'explicit B upper missing'
        return self.minimum, self.ties


def endpoint_germ(tau, ties):
    point = add(P1, scale(tau, OM7))
    assert dot(point, point).sign() < 0 and point[0].sign() > 0
    assert all(dot(point, wall) == ZERO for wall in RETAINED)
    proposals, accepted = [], []
    for q in ties:
        basis = RETAINED+(q['normal'],)
        gram = tuple(tuple(dot(v, w) for w in basis) for v in basis)
        pivots = positive_pivots(gram)
        assert pivots is not None, 'finite endpoint tangent Gram must be positive'
        witnesses, contains_all = [], True
        for other in ties:
            coefficients = matrix_solve(gram, tuple(dot(v, other['normal']) for v in basis))
            assert linear_combination(coefficients, basis) == other['normal']
            nonnegative = all(x.sign() >= 0 for x in coefficients)
            contains_all = contains_all and nonnegative
            witnesses.append({'mirror': other['id'], 'coefficients': coefficients,
                              'all_nonnegative': nonnegative})
        proposals.append({'simple_mirror': q['id'], 'positive_gram_pivots': pivots,
                          'all_ties_in_simple_cone': contains_all, 'witnesses': witnesses})
        if contains_all:
            accepted.append(q)
    return point, proposals, accepted


def apply_matrix(columns, v):
    return linear_combination(v, columns)


def determinant(columns):
    rows = [[columns[j][i] for j in range(10)] for i in range(10)]
    value = ONE
    for j in range(10):
        pivot = next((i for i in range(j, 10) if rows[i][j] != ZERO), None)
        assert pivot is not None
        if pivot != j:
            rows[pivot], rows[j] = rows[j], rows[pivot]
            value = -value
        p = rows[j][j]
        value *= p
        for i in range(j+1, 10):
            multiplier = rows[i][j]/p
            for k in range(j+1, 10):
                rows[i][k] -= multiplier*rows[j][k]
    return value


def period_certificate(winner):
    if winner['profile'] != 'B' or winner['normal'] != MB or winner['tau'] != TAU_B:
        return {'certified': False, 'reason': 'actual first simple is not the displayed B upper'}
    r = vector([U, 2]+[0]*8)
    g = add(NATIVE[0], scale(PHI, OM7))
    old_point = add(NATIVE[0], OM7)
    new_point = add(scale(VB/PHI, P1), scale(XB, OM7))
    tx = add(add(scale(PHI*VB, P1), scale(-PHI, C)), scale(PHI*PHI*XB, OM7))
    ty = scale(Q(1, 2), subtract(MB, scale(U, tx)))
    tom = add(add(scale(PHI, C), scale(-VB, P1)), scale(-PHI*XB, OM7))
    te7 = add(NATIVE[9], scale(Q(1, 2), subtract(tom, OM7)))
    columns = (tx, ty)+NATIVE[2:9]+(te7,)
    assert all(x.integral() for col in columns for x in col)
    assert determinant(columns) == ONE
    for i in range(10):
        for j in range(10):
            assert dot(columns[i], columns[j]) == dot(NATIVE[i], NATIVE[j])
    assert apply_matrix(columns, r) == MB
    assert apply_matrix(columns, g) == C
    assert apply_matrix(columns, old_point) == new_point
    assert new_point[0].sign() > 0 and dot(new_point, new_point) == -2/PHI
    old_germ = RETAINED[:7]+(r, g)
    new_germ = RETAINED[:7]+(MB, C)
    assert tuple(apply_matrix(columns, v) for v in old_germ) == new_germ
    trace = sum((columns[j][j] for j in range(10)), ZERO)
    assert trace == F(5059064, 8185727) and (trace-10).sign() > 0
    return {'certified': True, 'matrix_columns': columns, 'determinant': ONE,
            'trace': trace, 'old_point': old_point, 'new_point': new_point,
            'all_original_basis_entries_integral': True,
            'integral_inverse_reason': 'determinant1 over O implies integral adjugate inverse',
            'actual_oriented_germ_matches_conditional_on_pinned_geometry': True,
            'conditional_conclusion': 'infinite-order full-chamber stabilizer; literal k6 nonreflectivity',
            'external_inputs': ['actual complete old and starting germs',
                'complete first-simple classification and local finiteness',
                'integral full-arrangement/oriented-germ principle',
                'arithmetic lattice and nonphysical-definiteness inputs for infinite volume']}


def fixed_checks():
    assert B.norm() == 2099 and D == PHI*B and A == PHI*U
    assert B*VB*VB-PHI*XB*XB == ONE
    assert dot(P1, P1) == -2*D
    assert dot(OM7, OM7) == 2 and dot(OM6, OM6) == 6 and dot(OM6, OM7) == 3
    for j in range(8):
        assert dot(OM7, NATIVE[2+j]) == int(j == 7)
        assert dot(OM6, NATIVE[2+j]) == int(j == 6)
    assert_profile('B', MB, B)
    assert dot(MB, add(P1, scale(TAU_B, OM7))) == ZERO
    # Strict physical/conjugate domain constants from the pinned derivation.
    assert (F(110000)*F(110000)*B-XB*XB).sign() > 0
    assert (F(190000)*F(190000)*B-3*XB*XB).sign() > 0
    assert (F(130000)*F(130000)*B-(3-PHI)*XB*XB).sign() > 0
    assert (F(137000)-VB).sign() > 0
    assert (Q(4, 9)*(-D.conjugate())-1).sign() > 0
    assert (Q(81, 64)*(-D.conjugate())-3).sign() > 0
    assert (Q(25, 16)*(-D.conjugate())-(3-PHI).conjugate()).sign() > 0
    assert (Q(63001, 250000)*B.conjugate()-1).sign() > 0
    # Exact sufficient coefficient maxima: z<upper, z'>-bound implies
    # Bcoef<(upper+bound)/sqrt(5). The quotient must be below max_b+1.
    for maximum, upper, bound in ((50000, 110000, Q(2, 3)),
            (25000, 55000, Q(1, 3)), (87000, 190000, Q(9, 8)),
            (60000, 130000, Q(5, 4)), (63000, 137000, Q(251, 500))):
        assert (upper+bound)**2 < 5*(maximum+1)**2


def main():
    if not __debug__:
        raise RuntimeError('assertions must be enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--expected-source-sha', required=True)
    parser.add_argument('--spec-sha', required=True)
    parser.add_argument('--out', required=True, help='NEW .jsonl file directly under results/')
    args = parser.parse_args()
    assert sys.dont_write_bytecode, 'invoke through the -B bounded wrapper'
    owner = parent_lock()
    signal.alarm(45)
    resource.setrlimit(resource.RLIMIT_CPU, (45, 45))
    own_path = Path(__file__).resolve()
    assert own_path == ROOT/'scripts/check_compact_h9_prime2099_first_edge_20261021.py'
    pins = dict(PINS)
    pins[str(own_path.relative_to(ROOT))] = args.expected_source_sha
    pins[SPEC] = args.spec_sha
    assert all(len(v) == 64 and all(c in '0123456789abcdef' for c in v) for v in pins.values())
    check_pins(pins)
    destination = (ROOT/args.out).resolve()
    assert destination.parent == ROOT/'results' and destination.suffix == '.jsonl'
    assert not destination.exists(), 'certificate destination must be new'
    fixed_checks()
    with destination.open('xb') as stream:
        writer = Writer(stream)
        writer.write({'record': 'header', 'status': 'INCOMPLETE_UNTIL_FINAL_RESULT',
            'encoding': '[a,b] means a+b*phi, rational strings; original basis X0,Y0,e0,...,e7',
            'pins': pins, 'scalar_pair_upper_bound': 692011,
            'scope': 'Literal k6 next-edge census; geometric inputs remain explicit pinned theorems'})
        try:
            census = Census(writer, owner)
            tau, ties = census.run()
            point, proposals, accepted = endpoint_germ(tau, ties)
            unique = len(accepted) == 1
            period = period_certificate(accepted[0]) if unique else {
                'certified': False, 'reason': 'simple cone not uniquely certified'}
            if not unique:
                status = 'COMPLETE_SCALAR_DOMAIN_UNRESOLVED_SIMPLE_GERM'
            elif period['certified']:
                status = 'COMPLETE_FIRST_EDGE_AND_CONDITIONAL_LITERAL_K6_OBSTRUCTION'
            else:
                status = 'COMPLETE_FIRST_EDGE_WITH_NEW_ACTUAL_GERM_NO_HOST_EXCLUSION'
            check_pins(pins)
            assert os.getppid() == owner, 'bounded owner exited'
            writer.write({'record': 'result', 'complete': True, 'status': status,
                'counts': census.counts, 'admitted_mirrors': census.serial,
                'minimum_tau': tau, 'endpoint': point, 'minimum_ties': ties,
                'germ_cone_proposals': proposals, 'unique_simple_germ': unique,
                'simple_mirror_id': accepted[0]['id'] if unique else None,
                'actual_firstness_conditional_on_pinned_geometry': unique,
                'period': period, 'pins_rechecked_after': True,
                'no_new_example_or_form_class_exclusion': True})
            stream.flush()
            os.fsync(stream.fileno())
        except BaseException as error:
            writer.write({'record': 'failure', 'complete': False,
                          'exception': type(error).__name__, 'message': str(error)})
            stream.flush()
            raise
    signal.alarm(0)
    print(json.dumps({'status': status, 'scalar_pairs': sum(
        row['strip_pairs'] for row in census.counts.values()),
        'admitted_mirrors': census.serial, 'minimum_ties': len(ties),
        'unique_simple_germ': unique, 'period_certified': period['certified'],
        'output': str(destination.relative_to(ROOT)), 'output_bytes': destination.stat().st_size,
        'output_sha256': digest(destination)}, sort_keys=True))


if __name__ == '__main__':
    main()
