#!/usr/bin/env python3
"""Exact nine-model norm79 A7+c interval certificate; no host conclusion.

Prepared without execution. Use only the pinned parent-lock bounded runner.
The complete model/box theorem and actual starting germ are mathematical
inputs. This stdlib-only worker has no local imports or subprocesses.
JSONL keeps every reconstructed integral mirror, its known-wall decision,
the exact minimum and ALL ties. Only a successful final record is complete.
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
SPEC = 'results/compact_h9_prime79_a7_c_interval_checker_spec_20261021.md'
PINS = {
    'results/compact_h9_prime79_a7_c_finite_upper_and_complete_interval_contract_20261021.md':
        'b4de91853f6fd625831d591554cb47f7a61400fde49c233017ef905e8b24c3fa',
    'results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md':
        '33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3',
    'data/known_examples/bugaenko_compact_h8.json':
        'c17ac49e443918e61e541c9de35c7f29da9120e22bd9883678cc8091a2f4fc5c',
    'data/known_examples/bugaenko_vinberg_transcript.json':
        'c4178fc4678bacb48c430b2cec18e16577078824fdd6dea39375e63ed6ae015a',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/build_illustrated_ra_manuscript_20261018.py':
        '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py':
        'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


class F:
    """Exact Q(phi), phi^2=phi+1, in the ORIGINAL integral basis."""
    __slots__ = ('a', 'b')

    def __init__(self, a=0, b=0):
        self.a, self.b = Q(a), Q(b)

    @staticmethod
    def cast(x):
        return x if isinstance(x, F) else F(x)

    def __add__(self, x):
        x = F.cast(x)
        return F(self.a+x.a, self.b+x.b)

    __radd__ = __add__

    def __neg__(self):
        return F(-self.a, -self.b)

    def __sub__(self, x):
        return self+-F.cast(x)

    def __rsub__(self, x):
        return F.cast(x)+-self

    def __mul__(self, x):
        x = F.cast(x)
        return F(self.a*x.a+self.b*x.b,
                 self.a*x.b+self.b*x.a+self.b*x.b)

    __rmul__ = __mul__

    def conjugate(self):
        return F(self.a+self.b, -self.b)

    def norm(self):
        return self.a*self.a+self.a*self.b-self.b*self.b

    def __truediv__(self, x):
        x = F.cast(x)
        n = x.norm()
        require(n != 0, 'field division by zero')
        y = self*x.conjugate()
        return F(y.a/n, y.b/n)

    def __rtruediv__(self, x):
        return F.cast(x)/self

    def __eq__(self, x):
        x = F.cast(x)
        return self.a == x.a and self.b == x.b

    def sign(self):
        u, v = 2*self.a+self.b, self.b  # 2*self=u+v*sqrt5
        if not v:
            return (u > 0)-(u < 0)
        if not u or (u > 0) == (v > 0):
            return (v > 0)-(v < 0)
        d = u*u-5*v*v
        require(d != 0, 'irrational sqrt5 became rational')
        return ((u > 0)-(u < 0))*((d > 0)-(d < 0))

    def integral(self):
        return self.a.denominator == self.b.denominator == 1


ZERO, ONE, PHI = F(), F(1), F(0, 1)
COUPLING, BINARY_B, DISC = F(3, 3), F(13, 18), F(18, 31)
EDGES = ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7))
PATH = (0, 2, 3, 4, 5, 6, 7)
USTAR, ZSTAR, YSTAR = F(1206, 1953), F(19925, 32240), F(2490, 4030)
TAUSTAR = DISC*USTAR/ZSTAR
NATIVE = ((F(1), 3), (F(2), 4), (PHI*PHI, 5), (F(3), 6), (F(2)+PHI, 10))


def vec(xs):
    require(len(xs) == 10, 'vector length')
    return tuple(F.cast(x) for x in xs)


def add(v, w):
    return tuple(x+y for x, y in zip(v, w))


def scale(a, v):
    return tuple(a*x for x in v)


def dot(v, w):
    value = -2*PHI*v[0]*w[0]+COUPLING*(v[0]*w[1]+v[1]*w[0])+2*v[1]*w[1]
    value += sum((2*v[2+i]*w[2+i] for i in range(8)), ZERO)
    for i, j in EDGES:
        value -= v[2+i]*w[2+j]+v[2+j]*w[2+i]
    return value


UNIT = tuple(vec([int(i == j) for i in range(10)]) for j in range(10))
C = scale(-1, UNIT[1])
P1 = vec([2, -COUPLING]+[0]*8)


def integer_sqrt(q):
    if q < 0 or q.denominator != 1:
        return None
    r = isqrt(q.numerator)
    return r if r*r == q else None


def positive_o_sqrt(w):
    """COMPLETE O-square test: both root-norm signs, trace0, exact parity."""
    if not w.integral() or w.sign() < 0:
        return None
    if w == ZERO:
        return ZERO
    n0 = integer_sqrt(w.norm())
    if n0 is None:
        return None
    tr = 2*w.a+w.b
    for n in (n0, -n0):
        t, b = integer_sqrt(tr+2*n), integer_sqrt((tr-2*n)/5)
        if t is None or b is None:
            continue
        for tt in (t, -t):
            for bb in (b, -b):
                if (tt-bb) % 2:
                    continue
                r = F((tt-bb)//2, bb)
                if r*r == w and r.sign() > 0:
                    return r
    return None


def solve(matrix, rhs):
    rows = [list(row)+[x] for row, x in zip(matrix, rhs)]
    n = len(rhs)
    for j in range(n):
        pivot = next((i for i in range(j, n) if rows[i][j] != ZERO), None)
        require(pivot is not None, 'singular exact system')
        rows[j], rows[pivot] = rows[pivot], rows[j]
        d = rows[j][j]
        rows[j] = [x/d for x in rows[j]]
        for i in range(n):
            if i != j:
                m = rows[i][j]
                rows[i] = [x-m*y for x, y in zip(rows[i], rows[j])]
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


def combination(cs, vs):
    out = vec([0]*10)
    for a, v in zip(cs, vs):
        out = add(out, scale(a, v))
    return out


def covector(h):
    n = dot(h, h)
    return tuple(2*dot(h, u)/n for u in UNIT)


def encode(x):
    if isinstance(x, F):
        return [str(x.a), str(x.b)]
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, dict):
        return {k: encode(v) for k, v in x.items()}
    if isinstance(x, (tuple, list)):
        return [encode(v) for v in x]
    return x


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda: stream.read(65536), b''):
            h.update(b)
    return h.hexdigest()


def check_pins(pins):
    for p, sha in pins.items():
        require(digest(ROOT/p) == sha, 'changed pinned input: '+p)


def parent_lock():
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    require(owner > 1 and os.getppid() == owner, 'bounded parent wrapper required')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise RuntimeError('parent must retain shared exclusive lock')
    return owner


def setup():
    cartan = tuple(tuple(dot(UNIT[2+i], UNIT[2+j]) for j in range(8)) for i in range(8))
    omega = {j: vec([0, 0]+list(solve(cartan, [F(int(i == j)) for i in range(8)])))
             for j in (0, 1, 2, 3, 5, 6, 7)}
    q = omega[1]
    require(q == vec([0, 0, 5, 8, 10, 15, 12, 9, 6, 3]), 'omega1 calibration')
    expected = {0: (5, 4, 9), 2: (10, 14, 4), 3: (15, 30, 1),
                5: (9, 12, 1), 6: (6, 6, 4), 7: (3, 2, 9)}
    for j, (i, n, r) in expected.items():
        require(dot(q, omega[j]) == i and dot(omega[j], omega[j]) == n, 'E8 inverse row')
        require(i*i-8*n+16 == r, 'E8 norm residual')
    data = json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text())
    base = []
    for row in data['simple_roots']:
        v = [F(x['a'], x['b']) for x in row]
        base.append(vec([v[0], ZERO]+v[1:]))
    require(len(base) == 11, 'eleven base facets required')
    for i, v in enumerate(base):
        for j, w in enumerate(base):
            s = data['wall_gram'][i][j]
            require(dot(v, w) == F(s['a'], s['b']), 'base Gram calibration')
    r = vec([3*PHI, 2]+[0]*8)
    m = add(vec([PHI-1, -PHI]+[0]*8), scale(PHI*PHI, omega[7]))
    n = add(add(scale(PHI/2, P1), scale(-Q(1, 2), C)),
            add(omega[7], scale(PHI*PHI*PHI, omega[0])))
    known = [('base'+str(i), v) for i, v in enumerate(base)]+[('r', r), ('c', C), ('m', m), ('n', n)]
    for name, v in known:
        require(all(x.integral() for x in v), 'known wall original integrality: '+name)
        require(dot(v, v) == (2*BINARY_B if name == 'r' else 2), 'known wall norm: '+name)
        require(all(x.integral() for x in covector(v)), 'known wall reflection: '+name)
        require(dot(P1, v).sign() <= 0, 'known orientation: '+name)
    require(any(v == scale(-1, UNIT[3]) for _, v in known), 'released wall must be known')
    retained = tuple(scale(-1, UNIT[2+j]) for j in PATH)+(C,)
    require(all(dot(P1, v) == ZERO and dot(q, v) == ZERO for v in retained), 'retained edge')
    upper = add(add(scale(USTAR/2, P1), scale(-Q(1, 2), C)),
                add(omega[0], scale(YSTAR, q)))
    require(dot(P1, P1) == -2*DISC and dot(q, q) == 8, 'edge plane')
    require(ZSTAR*ZSTAR-4*DISC*USTAR*USTAR == 5 and 8*YSTAR+5 == ZSTAR, 'upper identities')
    require(all(x.integral() for x in upper) and dot(upper, upper) == 2, 'upper original lattice/norm')
    require(upper[:2] == (F(1206, 1953), F(-4738, -7668)), 'upper original binary coordinates')
    require(all(x.integral() for x in covector(upper)), 'upper reflection')
    require(dot(P1, upper) == -DISC*USTAR and dot(q, upper) == ZSTAR, 'upper collision')
    require(TAUSTAR.sign() > 0 and (DISC-4*TAUSTAR*TAUSTAR).sign() > 0, 'upper is finite')
    require(positive_o_sqrt(F(0)) == ZERO and positive_o_sqrt(F(4)) == 2, 'square zero/rational')
    require(positive_o_sqrt(F(5)) == F(-1, 2), 'square trace-zero branch')
    require(positive_o_sqrt(PHI*PHI) == PHI and positive_o_sqrt(F(3)) is None, 'square not norm-only')
    return omega, q, retained, known, upper


# name, kind, attachment, residual, coefficient maximum, radius, physical cap
MODELS = (
    ('c-only3', 'c', None, 3, 3050, Q(161, 100), 6800),
    ('A8-0', 'extension', 0, 9, 1350, Q(7, 10), 3000),
    ('A8-7', 'extension', 7, 9, 1350, Q(7, 10), 3000),
    ('D8-2', 'extension', 2, 4, 900, Q(93, 200), 2000),
    ('D8-6', 'extension', 6, 4, 900, Q(93, 200), 2000),
    ('E8-3', 'extension', 3, 1, 450, Q(233, 1000), 1000),
    ('E8-5', 'extension', 5, 1, 450, Q(233, 1000), 1000),
    ('A9-0', 'bridge', 0, 5, 1970, Q(26, 25), 4400),
    ('A9-7', 'bridge', 7, 5, 1970, Q(26, 25), 4400),
)


def known_decision(h, known):
    labels = []
    for name, w in known:
        if h == w:  # A known first wall is not rejected against itself.
            labels.append([name, 'same mirror'])
            continue
        ip = dot(h, w)
        if ip.sign() > 0:
            return False, {'wall': name, 'reason': 'positive_product', 'product': ip}
        if ip == ZERO:
            labels.append([name, 2])
            continue
        angle = 4*ip*ip/(dot(h, h)*dot(w, w))
        require(angle.integral(), 'integral reflections gave nonintegral two-cycle')
        gap = (angle-4).sign()
        if gap == 0:
            return False, {'wall': name, 'reason': 'parallel_anisotropic', 'product': ip}
        if gap > 0:
            labels.append([name, 'dotted'])
            continue
        label = next((m for a, m in NATIVE if angle == a), None)
        if label is None:
            return False, {'wall': name, 'reason': 'nonnative_finite_angle', 'angular_invariant': angle}
        labels.append([name, label])
    return True, {'known_wall_labels': labels}


def census(write, omega, q, retained, known, owner):
    counts, minimum, ties, serial = {}, None, [], 0
    for name, kind, j, residual, max_b, radius, cap in MODELS:
        c = dict.fromkeys(('raw_trials', 'strip_controls', 'exact_box_controls',
                          'field_squares', 'congruences', 'integral_mirrors', 'admitted'), 0)
        counts[name] = c
        for b in range(max_b+1):
            require(os.getppid() == owner, 'bounded owner exited')
            a0 = (isqrt(5*b*b)-b)//2
            offsets = (-1, 0, 1, 2) if radius >= 1 else (0, 1)
            for offset in offsets:
                c['raw_trials'] += 1
                a, t = a0+offset, F(a0+offset, b)
                if t.sign() <= 0 or (t-cap).sign() >= 0:
                    continue
                tc = t.conjugate()
                if (tc-radius).sign() >= 0 or (tc+radius).sign() <= 0:
                    continue
                c['strip_controls'] += 1
                physical = (Q(12, 5) if kind == 'c' else Q(residual, 20) if kind == 'extension' else Q(1))
                conjugate_factor = 1 if kind == 'c' else 16 if kind == 'extension' else 4
                conjugate_rhs = 3 if kind == 'c' else residual if kind == 'extension' else 5
                if (t*t-physical*USTAR*USTAR).sign() > 0:
                    continue
                if (conjugate_factor*(-DISC.conjugate())*tc*tc-conjugate_rhs).sign() > 0:
                    continue
                c['exact_box_controls'] += 1
                radicand = DISC*t*t+3 if kind == 'c' else F(residual)+16*DISC*t*t if kind == 'extension' else F(5)+4*DISC*t*t
                z = positive_o_sqrt(radicand)
                if z is None:
                    continue
                c['field_squares'] += 1
                y = z/4 if kind == 'c' else (z-dot(q, omega[j]))/8
                if not y.integral() or y.sign() < 0 or (kind == 'c' and y == ZERO):
                    continue
                if kind != 'extension' and not ((1-COUPLING*t)/2).integral():
                    continue
                c['congruences'] += 1
                h = add(scale(t if kind == 'extension' else t/2, P1), scale(y, q))
                if kind != 'extension':
                    h = add(h, scale(-Q(1, 2), C))
                if kind != 'c':
                    h = add(h, omega[j])
                require(all(x.integral() for x in h) and dot(h, h) == 2, 'reconstructed original norm-two mirror')
                cv = covector(h)
                require(all(x.integral() for x in cv), 'original reflection covector')
                expected = tuple(F(-1 if k == j else 0) for k in PATH)+(F(0 if kind == 'extension' else -1),)
                require(tuple(dot(h, w) for w in retained) == expected, 'retained profile mismatch')
                tau = -dot(P1, h)/dot(q, h)
                require(tau.sign() > 0 and (tau-TAUSTAR).sign() <= 0, 'collision outside exact closed interval')
                c['integral_mirrors'] += 1
                serial += 1
                okay, decision = known_decision(h, known)
                record = {'record': 'mirror', 'id': serial, 'model': name,
                          'control_coefficients': [a, b], 'control': t, 'positive_o_square_root': z,
                          'Y': y, 'normal': h, 'norm': F(2), 'reflection_covector': cv,
                          'primitive_by_even_norm_two': True, 'tau': tau,
                          'known_wall_admitted': okay, 'known_wall_decision': decision}
                write(record)
                if not okay:
                    continue
                c['admitted'] += 1
                tie = {'id': serial, 'model': name, 'normal': h, 'tau': tau}
                comparison = -1 if minimum is None else (tau-minimum).sign()
                if comparison < 0:
                    minimum, ties = tau, [tie]
                elif comparison == 0:
                    ties.append(tie)
                require(len(ties) <= len(MODELS), 'more than one minimum per monotone model')
    require(sum(c['raw_trials'] for c in counts.values()) == 38784, 'incomplete raw domain')
    require(minimum is not None, 'no candidate despite actual finite-upper theorem')
    return counts, minimum, ties, serial


def germ(tau, ties, q, retained):
    point = add(P1, scale(tau, q))
    require(dot(point, point).sign() < 0 and point[0].sign() > 0, 'endpoint not future finite')
    require(all(dot(point, v) == ZERO for v in retained), 'endpoint retained equations')
    proposals, accepted = [], []
    for candidate in ties:
        basis = retained+(candidate['normal'],)
        gram = tuple(tuple(dot(v, w) for w in basis) for v in basis)
        pivots = positive_pivots(gram)
        require(pivots is not None, 'finite candidate tangent Gram not positive')
        witnesses, okay = [], True
        for other in ties:
            require(dot(point, other['normal']) == ZERO, 'not an exact tie')
            cs = solve(gram, tuple(dot(v, other['normal']) for v in basis))
            require(combination(cs, basis) == other['normal'], 'cone reconstruction failed')
            nonnegative = all(x.sign() >= 0 for x in cs)
            witnesses.append({'id': other['id'], 'coefficients': cs, 'nonnegative': nonnegative})
            okay = okay and nonnegative
        proposals.append({'candidate_id': candidate['id'], 'gram': gram,
                          'positive_pivots': pivots, 'all_ties_in_simple_cone': okay,
                          'witnesses': witnesses})
        if okay:
            accepted.append(candidate)
    unique = bool(accepted) and all(v['normal'] == accepted[0]['normal'] for v in accepted)
    return point, proposals, accepted, unique


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--expected-source-sha', required=True)
    parser.add_argument('--spec-sha', required=True)
    parser.add_argument('--out', required=True, help='NEW .jsonl directly under results/')
    args = parser.parse_args()
    require(sys.dont_write_bytecode, 'use -B bounded wrapper')
    owner = parent_lock()
    signal.alarm(45)
    resource.setrlimit(resource.RLIMIT_CPU, (45, 45))
    own = Path(__file__).resolve()
    require(own == ROOT/'scripts/check_compact_h9_prime79_a7_c_interval_20261021.py', 'worker path')
    pins = dict(PINS)
    pins[str(own.relative_to(ROOT))], pins[SPEC] = args.expected_source_sha, args.spec_sha
    require(all(len(s) == 64 and all(x in '0123456789abcdef' for x in s) for s in pins.values()), 'invalid SHA')
    check_pins(pins)
    target = (ROOT/args.out).resolve()
    require(target.parent == ROOT/'results' and target.suffix == '.jsonl', 'new results JSONL path required')
    with target.open('xb') as stream:
        def write(record):
            stream.write((json.dumps(encode(record), sort_keys=True)+'\n').encode())
            stream.flush()
        write({'record': 'header', 'status': 'INCOMPLETE_UNTIL_FINAL_RESULT', 'pins': pins,
               'encoding': '[a,b] means a+b*phi with rational strings; original X,Y,e0,...,e7 basis',
               'raw_scalar_trials': 38784, 'models': [m[0] for m in MODELS],
               'excluded_profiles_by_pinned_proof': ['independent norm2', 'independent norm2b', 'c-only5'],
               'scope': 'Literal norm79 A7+c interval only; no period, nonreflectivity or example'})
        try:
            omega, q, retained, known, upper = setup()
            write({'record': 'setup', 'upper_normal': upper, 'upper_tau': TAUSTAR,
                   'known_walls': [{'name': name, 'normal': v} for name, v in known],
                   'retained': retained, 'square_and_inverse_fixtures_passed': True})
            counts, tau, ties, serial = census(write, omega, q, retained, known, owner)
            point, proposals, accepted, unique = germ(tau, ties, q, retained)
            check_pins(pins)
            require(os.getppid() == owner, 'bounded owner exited before completion')
            status = ('COMPLETE_INTERVAL_WITH_UNIQUE_ACTUAL_GERM_AT_PINNED_GEOMETRIC_BOUNDARY'
                      if unique else 'COMPLETE_NECESSARY_INTERVAL_UNRESOLVED_SIMPLE_GERM')
            write({'record': 'result', 'complete': True, 'status': status, 'counts': counts,
                   'integral_mirrors': serial, 'minimum_tau': tau, 'endpoint': point,
                   'minimum_ties': ties, 'germ_cone_proposals': proposals,
                   'unique_simple_germ': unique, 'simple_mirror_ids': [v['id'] for v in accepted] if unique else [],
                   'actual_firstness_conditional_on_pinned_geometry': unique,
                   'pins_rechecked_after': True, 'period_attempted': False,
                   'no_host_or_form_class_exclusion': True, 'no_new_example_or_bound': True})
            os.fsync(stream.fileno())
        except BaseException as error:
            write({'record': 'failure', 'complete': False, 'exception': type(error).__name__, 'message': str(error)})
            os.fsync(stream.fileno())
            raise
    signal.alarm(0)
    print(json.dumps({'status': status, 'raw_scalar_trials': sum(v['raw_trials'] for v in counts.values()),
                      'integral_mirrors': serial, 'minimum_ties': len(ties), 'unique_simple_germ': unique,
                      'output': str(target.relative_to(ROOT)), 'output_bytes': target.stat().st_size,
                      'output_sha256': digest(target)}, sort_keys=True))


if __name__ == '__main__':
    main()
