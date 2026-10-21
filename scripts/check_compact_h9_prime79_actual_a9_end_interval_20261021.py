#!/usr/bin/env python3
"""Prepared exact norm79 asymmetric-A9 END interval checker, not executed.

Standalone stdlib only. Bulk loops use integer pairs, with Fraction geometry
only after a complete O-square test. A selected profile is only a profile
certificate; all six profiles are required before a next-germ conclusion.
The actual starting corner remains conditional on the pinned geometric and
previous-certificate boundary. No chamber-period/host/example claim is made.
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
SPEC = 'results/compact_h9_prime79_actual_a9_end_interval_checker_spec_20261021.md'
START = 'results/prime79_a7_c_interval_20261021.jsonl'
START_LEDGER = 'results/prime79_a7_c_interval_execution_20261021.json'
PINS = {
    'results/compact_h9_prime79_actual_a9_end_edge_upper_and_complete_contract_candidate_20261021.md':
        '615d2915b1f95ffe5dc4d0ec634840b67991a7afb7b620aaca63d498bb10c2d3',
    'results/compact_h9_nonunit_coupling_prime79_full_base_and_first_transverse_20261021.md':
        '33f6cb9dd20ba91e7e7974586075309f2f437564f2e8dd978a16594349e48da3',
    'results/compact_h9_prime79_a7_c_finite_upper_and_complete_interval_contract_20261021.md':
        'b4de91853f6fd625831d591554cb47f7a61400fde49c233017ef905e8b24c3fa',
    'results/compact_h9_prime79_a9_bridge_integral_reversal_and_next_edge_candidate_20261021.md':
        '68a8a49491e4aa856ae67eacdbf367f411b1fc9c97911546b5961c980dc141a0',
    START: 'b5adb33eaa2c21dd2d9fc4f564c1a319914fd197911243b1bde83a4263d06050',
    START_LEDGER: 'ad2b9772be35bc956ecb20a761f6397fcd02d6e739a3d6fac2db791cf7769af4',
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


# Integer O=Z[phi] arithmetic for the COMPLETE bulk domain; no floats.
def ia(x, y):
    return x[0]+y[0], x[1]+y[1]


def ine(x):
    return -x[0], -x[1]


def isu(x, y):
    return ia(x, ine(y))


def im(x, y):
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def isc(n, x):
    return n*x[0], n*x[1]


def ic(x):
    return x[0]+x[1], -x[1]


def ins(x):
    a, b = x
    return a*a+a*b-b*b


def isign(x):
    u, v = 2*x[0]+x[1], x[1]
    if not v:
        return (u > 0)-(u < 0)
    if not u or (u > 0) == (v > 0):
        return (v > 0)-(v < 0)
    d = u*u-5*v*v
    require(d != 0, 'irrational sqrt5 became rational')
    return ((u > 0)-(u < 0))*((d > 0)-(d < 0))


def exact_isqrt(n):
    if n < 0:
        return None
    r = isqrt(n)
    return r if r*r == n else None


def positive_integer_o_sqrt(w):
    """Complete O-square test: BOTH norm signs, trace0, parity, resquaring."""
    if isign(w) < 0:
        return None
    if w == (0, 0):
        return (0, 0)
    root_norm = exact_isqrt(ins(w))
    if root_norm is None:
        return None
    trace = 2*w[0]+w[1]
    for n in (root_norm, -root_norm):
        t = exact_isqrt(trace+2*n)
        b_square = trace-2*n
        if t is None or b_square % 5:
            continue
        b = exact_isqrt(b_square//5)
        if b is None:
            continue
        for tt in (t, -t):
            for bb in (b, -b):
                if (tt-bb) % 2:
                    continue
                r = ((tt-bb)//2, bb)
                if im(r, r) == w and isign(r) > 0:
                    return r
    return None


def strip_floor(b):
    """EXACT floor(b*(phi-1)), including the required negative b values."""
    if b >= 0:
        return (isqrt(5*b*b)-b)//2
    k = -b
    return -((isqrt(5*k*k)-k)//2)-1


class F:
    """Exact Q(phi), confined to setup and integer-square survivors."""
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
        return F(self.a*x.a+self.b*x.b, self.a*x.b+self.b*x.a+self.b*x.b)

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
        # Clear only rational denominators; this is the SAME exact integer sign.
        den = self.a.denominator*self.b.denominator
        return isign((int(self.a*den), int(self.b*den)))

    def integral(self):
        return self.a.denominator == self.b.denominator == 1


ZERO, ONE, PHI = F(), F(1), F(0, 1)
COUPLING, BINARY_B, DISC = F(3, 3), F(13, 18), F(18, 31)
UOBS, ZOBS, YOBS = F(16, 25), F(259, 416), F(32, 52)
PELL_X, PELL_Z = F(8911, 14418), F(1080, 1746)
IDISC, IB, IPZ = (18, 31), (13, 18), (1080, 1746)
EDGES = ((0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (5, 6), (6, 7))
PATH = (0, 2, 3, 4, 5, 6, 7)
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


def combination(cs, vs):
    out = vec([0]*10)
    for a, v in zip(cs, vs):
        out = add(out, scale(a, v))
    return out


def covector(h):
    n = dot(h, h)
    return tuple(2*dot(h, u)/n for u in UNIT)


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


def decode_vector(xs):
    require(isinstance(xs, list) and len(xs) == 10, 'encoded vector length')
    return vec([F(*x) for x in xs])


def digest(p):
    h = hashlib.sha256()
    with p.open('rb') as stream:
        for b in iter(lambda: stream.read(65536), b''):
            h.update(b)
    return h.hexdigest()


def check_pins(pins):
    for p, sha in pins.items():
        target = (ROOT/p).resolve()
        require(target.is_relative_to(ROOT) and len(sha) == 64, 'invalid input binding')
        require(digest(target) == sha, 'changed pinned input: '+p)


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


def starting_records(pins):
    require((ROOT/START).stat().st_size == 11033, 'starting certificate size')
    with (ROOT/START).open() as stream:
        rows = [json.loads(line) for line in stream]
    require(len(rows) == 5 and rows[0]['record'] == 'header', 'starting JSONL schema')
    nested = rows[0]['pins']
    require(isinstance(nested, dict), 'starting nested pins')
    for path, sha in nested.items():
        require(path not in pins or pins[path] == sha, 'conflicting inherited pin: '+path)
        pins[path] = sha
    check_pins(pins)  # Includes EVERY inherited source/resource pin, not only outer inputs.
    setup = rows[1]
    result = rows[-1]
    require(setup['record'] == 'setup' and result['record'] == 'result', 'starting records')
    require(result['complete'] is True and result['unique_simple_germ'] is True,
            'starting germ not complete and unique')
    require(result['status'] == 'COMPLETE_INTERVAL_WITH_UNIQUE_ACTUAL_GERM_AT_PINNED_GEOMETRIC_BOUNDARY',
            'starting result status')
    require(result['actual_firstness_conditional_on_pinned_geometry'] is True,
            'starting geometric boundary missing')
    require(result['integral_mirrors'] == 2 and result['simple_mirror_ids'] == [2], 'starting totals')
    require(len(result['minimum_ties']) == 1 and result['minimum_ties'][0]['id'] == 2,
            'starting minimum/tie coverage')
    mirror = next(row for row in rows if row.get('record') == 'mirror' and row.get('id') == 2)
    require(mirror['model'] == 'A9-7' and mirror['known_wall_admitted'] is True, 'starting model')
    require(F(*mirror['control']) == UOBS and F(*mirror['Y']) == YOBS
            and F(*mirror['positive_o_square_root']) == ZOBS, 'starting exact controls')
    require(F(*result['minimum_tau']) == F(Q(1589, 1769), Q(3531, 1769)), 'starting collision')
    require(decode_vector(result['minimum_ties'][0]['normal']) == decode_vector(mirror['normal']),
            'starting minimum normal')
    ledger = json.loads((ROOT/START_LEDGER).read_text())
    require(ledger['status'] == 'PASS_BOUNDED_RESEARCH_CHECK'
            and ledger['wall_limit_seconds'] == 59 and ledger['rss_limit_bytes'] == 1073741824
            and ledger['threads'] == 1 and ledger['lock_owner'] == 'wrapper', 'starting bounded ledger')
    require(ledger['initial_resources']['free_disk_gib'] >= 40
            and ledger['initial_resources']['available_memory_gib'] >= 6, 'starting admission')
    command = ledger['command']
    require(command[0] == 'scripts/check_compact_h9_prime79_a7_c_interval_20261021.py'
            and command[-2:] == ['--out', START], 'starting command output')
    require(command[command.index('--expected-source-sha')+1] == nested[command[0]]
            and ledger['checker_sha256'] == nested[command[0]], 'starting worker binding')
    require(command[command.index('--spec-sha')+1]
            == nested['results/compact_h9_prime79_a7_c_interval_checker_spec_20261021.md']
            and ledger['runner_sha256'] == pins['scripts/run_bounded_research_check_20261018.py'],
            'starting contract/resource binding')
    summary = json.loads(ledger['stdout'])
    require(summary['output_sha256'] == pins[START] and summary['output_bytes'] == 11033
            and summary['raw_scalar_trials'] == 38784 and summary['integral_mirrors'] == 2
            and summary['minimum_ties'] == 1 and summary['unique_simple_germ'] is True,
            'starting stdout/certificate binding')
    return setup, mirror, result


def setup(end, pins):
    old_setup, observed, old_result = starting_records(pins)
    s = decode_vector(observed['normal'])
    cartan = tuple(tuple(dot(UNIT[2+i], UNIT[2+j]) for j in range(8)) for i in range(8))
    omega1 = vec([0, 0]+list(solve(cartan, [F(int(i == 1)) for i in range(8)])))
    omega7 = vec([0, 0]+list(solve(cartan, [F(int(i == 7)) for i in range(8)])))
    require(omega1 == vec([0, 0, 5, 8, 10, 15, 12, 9, 6, 3]), 'E8 omega1')
    expected_s = add(add(scale(UOBS/2, P1), scale(-Q(1, 2), C)), add(omega7, scale(YOBS, omega1)))
    require(s == expected_s and 8*YOBS+3 == ZOBS and ZOBS*ZOBS-4*DISC*UOBS*UOBS == 5,
            'observed original normal and norm equation')
    beta = tuple(scale(-1, UNIT[2+j]) for j in PATH)+(s, C)
    require(all(dot(v, w) == F(2 if i == j else -1 if abs(i-j) == 1 else 0)
                for i, v in enumerate(beta) for j, w in enumerate(beta)), 'complete ordered A9 Gram')
    p = add(scale(ZOBS, P1), scale(DISC*UOBS, omega1))
    require(dot(p, p) == -10*DISC and p[0].sign() > 0
            and all(dot(p, v) == ZERO for v in beta), 'observed finite future corner')
    require(decode_vector(old_result['endpoint']) == scale(ONE/ZOBS, p), 'starting endpoint ray')
    co = tuple(F(Q(min(end+1, k+1)*(10-max(end+1, k+1)), 10)) for k in range(9))
    om = combination(co, beta)
    ell = scale(10, om)
    require(all(x.integral() for x in ell) and dot(om, om) == F(Q(9, 10))
            and dot(p, om) == ZERO, 'actual end dual/glue')
    require(all(dot(om, v) == F(int(k == end)) for k, v in enumerate(beta)), 'dual pairings')
    if end == 0:
        mu = combination(tuple(F(Q(7-k, 8)) for k in range(7)), beta[:7])
        literal_ell = add(add(scale(UOBS, P1), scale(ZOBS/4, omega1)), scale(10, mu))
    else:
        literal_ell = add(add(scale(4*UOBS, P1), scale(5, C)), scale(ZOBS, omega1))
    require(ell == literal_ell, 'literal original-module end-dual formula')
    require(PELL_X*PELL_X-DISC*PELL_Z*PELL_Z == 1
            and (PELL_X-1)/9 == F(990, 1602) and PELL_Z/3 == F(360, 582), 'Pell/glue identities')

    def boost(w, direction):
        v, k = dot(w, p)/DISC, dot(w, ell)
        aa = -(PELL_X-1)*v+direction*(PELL_Z/3)*k
        bb = -direction*DISC*(PELL_Z/3)*v+((PELL_X-1)/9)*k
        return add(w, scale(Q(1, 10), add(scale(aa, p), scale(bb, ell))))

    images = tuple(boost(v, 1) for v in UNIT)
    inverse = tuple(boost(v, -1) for v in UNIT)
    require(all(x.integral() for v in images+inverse for x in v), 'boost original-module integrality')
    require(all(dot(v, w) == dot(UNIT[i], UNIT[j])
                for i, v in enumerate(images) for j, w in enumerate(images)), 'boost isometry')
    require(all(boost(v, -1) == UNIT[i] for i, v in enumerate(images)), 'boost inverse')
    retained = tuple(v for k, v in enumerate(beta) if k != end)
    require(all(boost(v, 1) == v for v in retained), 'boost retained A8')
    upper = scale(-1, boost(beta[end], -1))
    upper_t = (10*DISC/3)*PELL_Z/PELL_X
    require(dot(p, upper) == -(10*DISC/3)*PELL_Z and dot(om, upper) == -PELL_X,
            'finite-upper pairings')
    require(dot(upper, upper) == 2 and all(x.integral() for x in upper)
            and all(x.integral() for x in covector(upper)), 'integral upper reflector')
    require(upper_t.sign() > 0 and dot(add(p, scale(-upper_t, om)),
            add(p, scale(-upper_t, om))).sign() < 0, 'upper before light cone')
    known = [(row['name'], decode_vector(row['normal'])) for row in old_setup['known_walls']]+[('observed_s', s)]
    require(len(known) == 16, 'complete inherited known-wall list')
    for name, wall in known:
        require(all(x.integral() for x in wall) and all(x.integral() for x in covector(wall)),
                'known original wall: '+name)
        require(dot(wall, wall) == (2*BINARY_B if name == 'r' else 2), 'known norm: '+name)
        require(dot(p, wall).sign() <= 0, 'known orientation at actual corner: '+name)
    gram = tuple(tuple(dot(v, w) for w in retained) for v in retained)
    require(all(gram[i][j] == F(2 if i == j else -1 if abs(i-j) == 1 else 0)
                for i in range(8) for j in range(8)), 'retained A8 order')
    projections = {}
    for name, attach, *_ in MODELS:
        products = tuple(F(-int(k == attach)) if attach is not None else ZERO for k in range(8))
        projections[name] = combination(solve(gram, products), retained)
    # Algorithm fixtures, not a substitute for source review or complete enumeration.
    for w, want in (((0, 0), (0, 0)), ((4, 0), (2, 0)), ((5, 0), (-1, 2)),
                    ((1, 1), (0, 1)), ((3, 0), None)):
        require(positive_integer_o_sqrt(w) == want, 'integer-pair square fixture')
    require(strip_floor(-2) == -2 and strip_floor(-1) == -1
            and strip_floor(0) == 0 and strip_floor(1) == 0 and strip_floor(2) == 1,
            'negative/positive strip floors')
    return p, om, ell, retained, known, upper, upper_t, projections


# name, A8 attachment, norm multiplier, integer-box multiplier, base RHS,
# min/max phi coefficient, radius numerator/denominator, physical cap, offsets.
# Exact boxes are d*A0^2<=base*PELL_Z^2 and (-D')*d*A0'^2<=base'.
MODELS = (
    ('independent2', None, (1, 0), 1, (20, 0), -1, 7830, 104, 25, 17500, -4, 5),
    ('independent2b', None, IB, 1, (260, 360), -2, 50770, 57, 10, 113500, -5, 6),
    ('A9-left', 0, (1, 0), 9, (100, 0), -1, 5830, 31, 10, 13020, -3, 4),
    ('A9-right', 7, (1, 0), 9, (100, 0), -1, 5830, 31, 10, 13020, -3, 4),
    ('D9-left', 1, (1, 0), 9, (40, 0), 0, 3690, 49, 25, 8240, -1, 2),
    ('D9-right', 6, (1, 0), 9, (40, 0), 0, 3690, 49, 25, 8240, -1, 2),
)


def known_decision(h, known):
    labels = []
    for name, w in known:
        if h == w:
            labels.append([name, 'same mirror'])
            continue
        ip = dot(h, w)
        if ip.sign() > 0:
            return False, {'wall': name, 'reason': 'positive_product', 'product': ip}
        if ip == ZERO:
            labels.append([name, 2])
            continue
        angle = 4*ip*ip/(dot(h, h)*dot(w, w))
        require(angle.integral(), 'integral reflectors gave nonintegral two-cycle')
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


def census(write, models, geometry, owner):
    p, om, ell, retained, known, _, upper_t, projections = geometry
    counts, minimum, ties, serial = {}, None, [], 0
    iz2, minus_dc = im(IPZ, IPZ), ine(ic(IDISC))
    for model in models:
        name, attach, norm_d, factor, base, lo_b, hi_b, rn, rd, cap, lo_a, hi_a = model
        ct = dict.fromkeys(('raw_trials', 'strip_controls', 'exact_box_controls', 'field_squares',
                            'original_integral', 'integral_mirrors', 'admitted'), 0)
        counts[name] = ct
        rhs = isc(9, base) if factor == 1 else base  # 90Delta=9base/factor.
        physical_rhs = im(base, iz2)
        conjugate_rhs = ic(base)
        for b in range(lo_b, hi_b+1):
            require(os.getppid() == owner, 'bounded owner exited')
            floor = strip_floor(b)
            for offset in range(lo_a, hi_a+1):
                ct['raw_trials'] += 1
                a = floor+offset
                control = (a, b)
                if isign(control) <= 0 or isign(isu(control, (cap, 0))) >= 0:
                    continue
                cc = ic(control)
                if isign(isu(isc(rd, cc), (rn, 0))) >= 0 or isign(ia(isc(rd, cc), (rn, 0))) <= 0:
                    continue
                ct['strip_controls'] += 1
                control2 = im(control, control)
                if isign(isu(isc(factor, control2), physical_rhs)) > 0:
                    continue
                if isign(isu(im(minus_dc, isc(factor, ic(control2))), conjugate_rhs)) > 0:
                    continue
                ct['exact_box_controls'] += 1
                square = ia(isc(9, im(IDISC, control2)), rhs)
                e0 = positive_integer_o_sqrt(square)
                if e0 is None:
                    continue
                ct['field_squares'] += 1
                aa, ee = F(*control), F(*e0)
                h = add(projections[name], add(scale(aa/10, p), scale(-ee/9, om)))
                if not all(x.integral() for x in h):
                    continue
                ct['original_integral'] += 1
                expected_norm = 2*F(*norm_d)
                require(dot(h, h) == expected_norm and -dot(h, p)/DISC == aa and -dot(h, ell) == ee,
                        'reconstructed exact norm/controls')
                cv = covector(h)
                if not all(x.integral() for x in cv):
                    continue
                expected = tuple(F(-int(k == attach)) if attach is not None else ZERO for k in range(8))
                require(tuple(dot(h, w) for w in retained) == expected, 'retained profile mismatch')
                tau = 10*DISC*aa/ee
                require(tau.sign() > 0 and (tau-upper_t).sign() <= 0
                        and dot(add(p, scale(-tau, om)), add(p, scale(-tau, om))).sign() < 0,
                        'collision outside finite CLOSED interval')
                ct['integral_mirrors'] += 1
                serial += 1
                okay, decision = known_decision(h, known)
                write({'record': 'mirror', 'id': serial, 'model': name,
                       'A0_coefficients': list(control), 'E0_coefficients': list(e0),
                       'A0': aa, 'E0': ee, 'normal': h, 'norm': expected_norm,
                       'reflection_covector': cv, 'primitive_by_even_squarefree_norm': True,
                       't': tau, 'known_wall_admitted': okay, 'known_wall_decision': decision})
                if not okay:
                    continue
                ct['admitted'] += 1
                tie = {'id': serial, 'model': name, 'normal': h, 't': tau}
                comparison = -1 if minimum is None else (tau-minimum).sign()
                if comparison < 0:
                    minimum, ties = tau, [tie]
                elif comparison == 0:
                    ties.append(tie)
                require(len(ties) <= len(models), 'more than one minimum per monotone profile')
    expected_trials = sum((m[6]-m[5]+1)*(m[11]-m[10]+1) for m in models)
    require(sum(c['raw_trials'] for c in counts.values()) == expected_trials, 'incomplete scalar domain')
    if len(models) == 6:
        require(expected_trials == 810436 and minimum is not None, 'complete interval lacks a candidate')
    return counts, minimum, ties, serial


def germ(tau, ties, geometry):
    p, om, _, retained, _, _, _, _ = geometry
    point = add(p, scale(-tau, om))
    require(dot(point, point).sign() < 0 and dot(point, p).sign() < 0, 'endpoint not finite/future')
    require(all(dot(point, v) == ZERO for v in retained), 'endpoint retained equations')
    proposals, accepted = [], []
    for candidate in ties:
        basis = retained+(candidate['normal'],)
        gram = tuple(tuple(dot(v, w) for w in basis) for v in basis)
        pivots = positive_pivots(gram)
        require(pivots is not None, 'candidate tangent Gram not positive definite')
        witnesses, okay = [], True
        for other in ties:
            require(dot(point, other['normal']) == ZERO, 'inexact collision tie')
            cs = solve(gram, tuple(dot(v, other['normal']) for v in basis))
            require(combination(cs, basis) == other['normal'], 'tie cone reconstruction')
            nonnegative = all(x.sign() >= 0 for x in cs)
            witnesses.append({'id': other['id'], 'coefficients': cs, 'nonnegative': nonnegative})
            okay = okay and nonnegative
        proposals.append({'candidate_id': candidate['id'], 'gram': gram, 'positive_pivots': pivots,
                          'all_ties_in_simple_cone': okay, 'witnesses': witnesses})
        if okay:
            accepted.append(candidate)
    unique = bool(accepted) and all(v['normal'] == accepted[0]['normal'] for v in accepted)
    return point, proposals, accepted, unique


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--end', type=int, choices=(0, 8), required=True)
    parser.add_argument('--profile', choices=('all',)+tuple(m[0] for m in MODELS), default='all')
    parser.add_argument('--expected-source-sha', required=True)
    parser.add_argument('--spec-sha', required=True)
    parser.add_argument('--out', required=True, help='NEW JSONL directly under results/')
    args = parser.parse_args()
    require(sys.dont_write_bytecode, 'use -B bounded wrapper')
    owner = parent_lock()
    signal.alarm(45)
    resource.setrlimit(resource.RLIMIT_CPU, (45, 45))
    own = Path(__file__).resolve()
    require(own == ROOT/'scripts/check_compact_h9_prime79_actual_a9_end_interval_20261021.py', 'worker path')
    pins = dict(PINS)
    pins[str(own.relative_to(ROOT))], pins[SPEC] = args.expected_source_sha, args.spec_sha
    require(all(len(s) == 64 and all(x in '0123456789abcdef' for x in s) for s in pins.values()), 'invalid SHA')
    check_pins(pins)
    models = tuple(m for m in MODELS if args.profile == 'all' or m[0] == args.profile)
    raw_trials = sum((m[6]-m[5]+1)*(m[11]-m[10]+1) for m in models)
    target = (ROOT/args.out).resolve()
    require(target.parent == ROOT/'results' and target.suffix == '.jsonl', 'new results JSONL path required')
    with target.open('xb') as stream:
        def write(record):
            stream.write((json.dumps(encode(record), sort_keys=True)+'\n').encode())
            stream.flush()
        try:
            geometry = setup(args.end, pins)
            write({'record': 'header', 'status': 'INCOMPLETE_UNTIL_FINAL_RESULT', 'pins': pins,
                   'end': args.end, 'profile_selection': args.profile, 'models': [m[0] for m in models],
                   'raw_scalar_trials': raw_trials,
                   'encoding': '[a,b]=a+b*phi with rational strings; ORIGINAL X,Y,e0,...,e7 coordinates',
                   'scope': 'Literal norm79 observed A9-7 END interval; actual start conditional on pinned geometry',
                   'no_chamber_period_or_host_exclusion': True, 'no_new_example_or_bound': True})
            p, om, ell, retained, known, upper, upper_t, _ = geometry
            write({'record': 'setup', 'corner': p, 'dual': om, 'integral_dual10': ell,
                   'upper_normal': upper, 'upper_t': upper_t,
                   'upper_is_arrangement_mirror_only': True,
                   'retained': retained, 'known_walls': [{'name': name, 'normal': v} for name, v in known],
                   'boost_original_integrality_and_inverse_fixtures_passed': True,
                   'square_and_negative_strip_floor_fixtures_passed': True})
            counts, minimum, ties, serial = census(write, models, geometry, owner)
            unique, proposals, accepted, point = False, [], [], None
            if args.profile == 'all':
                point, proposals, accepted, unique = germ(minimum, ties, geometry)
                status = ('COMPLETE_END_INTERVAL_WITH_UNIQUE_ACTUAL_GERM_AT_PINNED_GEOMETRIC_BOUNDARY'
                          if unique else 'COMPLETE_NECESSARY_END_INTERVAL_UNRESOLVED_SIMPLE_GERM')
            else:
                status = 'COMPLETE_NECESSARY_PROFILE_INTERVAL_ONLY'
            check_pins(pins)
            require(os.getppid() == owner, 'bounded owner exited before completion')
            write({'record': 'result', 'complete': True, 'status': status, 'end': args.end,
                   'profile_selection': args.profile, 'all_six_profiles_complete': args.profile == 'all',
                   'counts': counts, 'integral_mirrors': serial, 'minimum_t': minimum,
                   'endpoint': point, 'minimum_ties': ties, 'germ_cone_proposals': proposals,
                   'unique_simple_germ': unique,
                   'simple_mirror_ids': [v['id'] for v in accepted] if unique else [],
                   'actual_firstness_conditional_on_pinned_geometry': unique,
                   'pins_rechecked_after_including_all_starting_nested_pins': True,
                   'period_attempted': False, 'no_host_or_form_class_exclusion': True,
                   'no_new_example_or_bound': True})
            os.fsync(stream.fileno())
        except BaseException as error:
            write({'record': 'failure', 'complete': False, 'exception': type(error).__name__, 'message': str(error)})
            os.fsync(stream.fileno())
            raise
    signal.alarm(0)
    print(json.dumps({'status': status, 'end': args.end, 'profile_selection': args.profile,
                      'raw_scalar_trials': sum(v['raw_trials'] for v in counts.values()),
                      'integral_mirrors': serial, 'minimum_ties': len(ties), 'unique_simple_germ': unique,
                      'output': str(target.relative_to(ROOT)), 'output_bytes': target.stat().st_size,
                      'output_sha256': digest(target)}, sort_keys=True))


if __name__ == '__main__':
    main()
