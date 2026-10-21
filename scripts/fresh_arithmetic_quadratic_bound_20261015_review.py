#!/usr/bin/env python3
"""Independent Fraction Qsqrt5 profile and nested-radical cycle review.

No SymPy and no producer imports. All real square-root membership tests are
algebraic equations in Qsqrt5 or Qsqrt2sqrt5, with exact rational signs.
"""
import argparse
from collections import Counter, deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd, isqrt, lcm
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_quadratic_bound_20261015_review'
PRIMARY = 'fresh_arithmetic_quadratic_bound_20261015'


def rational_sqrt(x):
    x = F(x)
    if x < 0:
        return None
    a, b = isqrt(x.numerator), isqrt(x.denominator)
    return F(a, b) if a*a == x.numerator and b*b == x.denominator else None


class Q5:
    def __init__(self, a=0, b=0):
        if isinstance(a, Q5):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = F(a), F(b)
    def __add__(self, y):
        y = Q5(y)
        return Q5(self.a+y.a, self.b+y.b)
    __radd__ = __add__
    def __neg__(self):
        return Q5(-self.a, -self.b)
    def __sub__(self, y):
        return self+-Q5(y)
    def __rsub__(self, y):
        return Q5(y)+-self
    def __mul__(self, y):
        y = Q5(y)
        return Q5(self.a*y.a+5*self.b*y.b, self.a*y.b+self.b*y.a)
    __rmul__ = __mul__
    def __truediv__(self, y):
        y = Q5(y)
        n = y.a*y.a-5*y.b*y.b
        assert n
        z = self*Q5(y.a, -y.b)
        return Q5(z.a/n, z.b/n)
    def __eq__(self, y):
        y = Q5(y)
        return (self.a, self.b) == (y.a, y.b)
    def __bool__(self):
        return bool(self.a or self.b)
    def sign(self):
        sa, sb = (self.a > 0)-(self.a < 0), (self.b > 0)-(self.b < 0)
        if not sa:
            return sb
        if not sb or sa == sb:
            return sa
        q = self.a*self.a-5*self.b*self.b
        return sa*((q > 0)-(q < 0))
    def conjugate(self):
        return Q5(self.a, -self.b)
    def encoded(self):
        return [str(self.a), str(self.b)]
    def integral(self):
        # O_Qsqrt5=Z[(1+sqrt5)/2].
        return (self.a-self.b).denominator == 1 and (2*self.b).denominator == 1


def sqrt_q5(z):
    z = Q5(z)
    if z.sign() < 0:
        return None
    candidates = []
    if not z.b:
        u = rational_sqrt(z.a)
        if u is not None:
            candidates.append(Q5(u))
        v = rational_sqrt(z.a/5)
        if v is not None:
            candidates.append(Q5(0, v))
    else:
        q = rational_sqrt(z.a*z.a-5*z.b*z.b)
        if q is not None:
            for signed in (q, -q):
                u = rational_sqrt((z.a+signed)/2)
                if u:
                    candidates.append(Q5(u, z.b/(2*u)))
    for x in candidates:
        if x*x == z:
            return x if x.sign() >= 0 else -x
    return None


class Biquad:
    """u+sqrt2*v, u,v in independently implemented Qsqrt5."""
    def __init__(self, u=0, v=0):
        if isinstance(u, Biquad):
            self.u, self.v = u.u, u.v
        else:
            self.u, self.v = Q5(u), Q5(v)
    def __add__(self, y):
        y = Biquad(y)
        return Biquad(self.u+y.u, self.v+y.v)
    __radd__ = __add__
    def __neg__(self):
        return Biquad(-self.u, -self.v)
    def __sub__(self, y):
        return self+-Biquad(y)
    def __mul__(self, y):
        y = Biquad(y)
        return Biquad(self.u*y.u+2*self.v*y.v, self.u*y.v+self.v*y.u)
    __rmul__ = __mul__
    def __bool__(self):
        return bool(self.u or self.v)
    def __eq__(self, y):
        y = Biquad(y)
        return self.u == y.u and self.v == y.v
    def sign(self):
        a, b = self.u.sign(), self.v.sign()
        if not a:
            return b
        if not b or a == b:
            return a
        return a*(self.u*self.u-2*self.v*self.v).sign()
    def encoded(self):
        return [self.u.encoded(), self.v.encoded()]


def sqrt_biquad(z):
    z = Biquad(z)
    assert z.sign() > 0
    candidates = []
    if not z.v:
        u, v = sqrt_q5(z.u), sqrt_q5(z.u/2)
        if u is not None:
            candidates.append(Biquad(u))
        if v is not None:
            candidates.append(Biquad(0, v))
    else:
        discriminant = sqrt_q5(z.u*z.u-2*z.v*z.v)
        if discriminant is not None:
            for d in (discriminant, -discriminant):
                u = sqrt_q5((z.u+d)/2)
                if u:
                    candidates.append(Biquad(u, z.v/(2*u)))
    for x in candidates:
        if x*x == z:
            return x if x.sign() > 0 else -x
    return None


def decode_phi(v):
    a, b = map(F, v)
    return Q5(a+b/2, b/2)


def decode_biquad(v):
    a, b, c, d = map(F, v)
    return Biquad(Q5(a, c), Q5(b, d))


def gram(n):
    a = [[F(2*(i == j)) for j in range(n)] for i in range(n)]
    for i, j in [(i, i+1) for i in range(n-2)]+[(2, n-1)]:
        a[i][j] = a[j][i] = F(-1)
    return a


def inverse(a):
    n = len(a)
    x = [list(row)+[F(i == j) for j in range(n)] for i, row in enumerate(a)]
    for j in range(n):
        pivot = x[j][j]
        assert pivot > 0
        x[j] = [v/pivot for v in x[j]]
        for i in range(n):
            if i != j:
                p = x[i][j]
                x[i] = [u-p*v for u, v in zip(x[i], x[j])]
    result = [row[n:] for row in x]
    assert all(sum(a[i][k]*result[k][j] for k in range(n)) == (i == j)
               for i in range(n) for j in range(n))
    return result


def quadratic_form(j, u, v):
    return sum(j[a][b]*u[a]*v[b] for a in range(len(j)) for b in range(len(j)))


def profile_review(record):
    n = int(record['type'][1:])
    j = inverse(gram(n))
    assert [[str(v) for v in row] for row in j] == record['inverse']
    minimum = min(j[i][i] for i in range(n))
    assert minimum == (F(3, 2) if n == 7 else F(2))
    phi = Q5(F(1, 2), F(1, 2))
    pure_margins = [minimum*q-2 for q in [Q5(2), Q5(3), 3-phi]]
    assert all(q.sign() > 0 for q in pure_margins)
    assert all(x == decode_phi(y) for x, y in zip(pure_margins, record['pure_class_exclusion_margins']))
    counts = Counter(total=0, conjugate_not_positive=0, physical_not_indefinite=0, accepted=0)
    survivors, subcritical = [], []
    pair_indices = {p: [] for p in combinations(range(n), 2)}
    for labels in product((0, 1, 2), repeat=n):
        counts['total'] += 1
        u = [F(0) if x == 0 else F(1) if x == 1 else F(1, 2) for x in labels]
        v = [F(x == 2, 2) for x in labels]
        cost = Q5(quadratic_form(j, u, u)+5*quadratic_form(j, v, v),
                  2*quadratic_form(j, u, v))
        if (2-cost.conjugate()).sign() <= 0:
            counts['conjugate_not_positive'] += 1
            continue
        if (cost-2).sign() <= 0:
            counts['physical_not_indefinite'] += 1
            if any(labels):
                subcritical.append(dict(labels=labels, cost=cost.encoded()))
            continue
        counts['accepted'] += 1
        dual = [Q5(sum(j[i][k]*u[k] for k in range(n)),
                   sum(j[i][k]*v[k] for k in range(n))) for i in range(n)]
        deleted = [cost-dual[i]*dual[i]/j[i][i] for i in range(n)]
        finite = [i for i, q in enumerate(deleted) if (2-q).sign() > 0]
        contacts = [[i, 3 if x == 1 else 5] for i, x in enumerate(labels) if x]
        index = len(survivors)
        saved = record['profiles'][index]
        assert saved['contacts'] == contacts
        assert cost == decode_phi(saved['cost']) and cost.conjugate() == decode_phi(saved['conjugate_cost'])
        assert finite == saved['finite_deletions']
        assert all(x == decode_phi(y) for x, y in zip(deleted, saved['deletion_costs']))
        for pair in combinations(finite, 2):
            pair_indices[pair].append(index)
        survivors.append(dict(contacts=contacts, cost=cost.encoded(), finite_deletions=finite))
    assert dict(counts) == record['counts']
    expected_subcritical = [dict(contacts=[], cost=Q5().encoded())]
    expected_subcritical += [dict(contacts=[[i, 3 if x == 1 else 5]
                                           for i, x in enumerate(r['labels']) if x],
                                 cost=r['cost']) for r in subcritical]
    assert expected_subcritical == [dict(contacts=r['contacts'], cost=decode_phi(r['cost']).encoded())
                                    for r in record['subcritical_profiles']]
    for row in record['marked_pairs']:
        assert pair_indices[tuple(row['pair'])] == row['profiles']
    assert sum(bool(x) for x in pair_indices.values()) == record['triangle_covered_pairs'] == 10
    return dict(type=record['type'], counts=dict(counts), accepted_profiles=survivors,
                covered_pairs=[p for p, ix in pair_indices.items() if ix],
                nonzero_subcritical=subcritical, pure_class_margins=[x.encoded() for x in pure_margins])


def primitive_minpoly(z):
    if not z.b:
        coefficients = [-z.a, F(1)]
    else:
        coefficients = [z.a*z.a-5*z.b*z.b, -2*z.a, F(1)]
    denominator = lcm(*(x.denominator for x in coefficients))
    values = [int(x*denominator) for x in coefficients]
    common = gcd(*values)
    values = [x//common for x in values]
    if values[-1] < 0:
        values = [-x for x in values]
    return values


def old_tree_path(n, start, target):
    a = gram(n)
    queue = deque([(start, [start])])
    seen = {start}
    while queue:
        v, path = queue.popleft()
        if v == target:
            return path
        for w in range(n):
            if a[v][w] == -1 and w not in seen:
                seen.add(w)
                queue.append((w, path+[w]))
    raise AssertionError('disconnected old tree')


def quad_review(primary):
    old = json.loads((ROOT/'results/fresh_compact29_small_faces_20261015.json').read_text())
    phi = Q5(F(1, 2), F(1, 2))
    magnitude = {2: Biquad(), 3: Biquad(1), 4: Biquad(0, 1), 5: Biquad(phi)}
    rows = []
    for data in old['cases']:
        rec, proof = data['recipe'], data['proof']
        if not rec['type'].startswith('E') or rec['kind'] != 'quadrilateral':
            continue
        n = int(rec['type'][1:])
        old_a, old_b = rec['marked']
        a, b = map(decode_biquad, proof['dotted_magnitude'])
        D = decode_biquad(proof['positive_radicand'])
        assert D.sign() > 0
        if rec['cd_label'] == 2:
            A, B = a*a+b*b*D, 2*a*b
            cycle = [old_b, n+1]
            kind = 'dotted two-cycle'
        else:
            if rec['ac_label'] != 2:
                k, label = old_a, rec['ac_label']
            else:
                k, label = rec['c_contacts'][0]
            factor = magnitude[rec['cd_label']]*magnitude[label]
            assert factor
            A, B = factor*a, factor*b
            path = old_tree_path(n, k, old_b)
            cycle = [old_b, n+1, n]+path[:-1]
            assert len(cycle) == len(set(cycle))
            kind = 'fundamental cycle through dotted edge'
        square_root = sqrt_biquad(D) if B else None
        obstruction, value = None, None
        if B and square_root is None:
            obstruction = 'cyclic product outside Q(sqrt5)'
            method = 'nontrivial quadratic extension over Q(sqrt2,sqrt5)'
        else:
            z = A+B*square_root if B else A
            if z.v:
                obstruction = 'cyclic product outside Q(sqrt5)'
                method = 'nonzero sqrt2 coefficient over Q(sqrt5)'
            else:
                value = z.u
                assert not value.integral()
                obstruction = 'cyclic product in Q(sqrt5) but not an algebraic integer'
                method = 'exact O_Qsqrt5 coefficient test'
        saved = primary['literal_prior_quad_obstructions'][len(rows)]
        assert (rec['type'], rec['marked'], kind, obstruction) == (
            saved['type'], saved['marked'], saved['cycle_kind'], saved['obstruction'])
        polynomial = primitive_minpoly(value) if value is not None else None
        if polynomial:
            assert polynomial[-1] != 1
        rows.append(dict(type=rec['type'], marked=rec['marked'], literal_cycle=cycle,
                         expression_constant=A.encoded(), expression_radical_coefficient=B.encoded(),
                         radicand=D.encoded(), square_root_in_biquad=square_root.encoded() if square_root else None,
                         obstruction=obstruction, independent_method=method,
                         value_in_Qsqrt5=value.encoded() if value else None,
                         primitive_minpoly_low_to_high=polynomial))
    assert len(rows) == 23
    assert sum('outside' in r['obstruction'] for r in rows) == 20
    return rows


def spectator_checks(records):
    j = inverse(gram(7))
    assert [i for i in range(7) if j[i][i] < 2] == [5]
    assert j[5][5] == F(3, 2)
    assert all(j[i][i]+j[k][k]+2*j[i][k] >= 3 for i, k in combinations(range(7), 2))
    sizes = []
    for omitted in range(7):
        retained = set(range(7))-{omitted}
        if 5 not in retained:
            sizes.append(0)
        else:
            sizes.append(sum(k in retained and omitted not in old_tree_path(7, 5, k) for k in retained))
    assert sizes == [6, 5, 3, 2, 1, 0, 6]
    near = [i for i, n in enumerate(sizes) if n <= 2]
    assert near == [3, 4, 5]
    assert all(tuple(p) in set(map(tuple, records[0]['covered_pairs'])) for p in combinations(near, 2))
    phi = Q5(F(1, 2), F(1, 2))
    full_cost = Q5(F(3, 2))+phi*phi/2
    assert full_cost == Q5(F(9, 4), F(1, 4))
    assert (full_cost-2).sign() > 0 and (2-full_cost.conjugate()).sign() > 0
    finite = []
    for i in range(7):
        cost_after_deletion = full_cost-j[i][5]*j[i][5]/j[i][i]
        if (2-cost_after_deletion).sign() > 0:
            finite.append(i)
    assert finite == [3, 4, 5]
    return dict(retained_component_sizes=sizes, forced_rescue_nodes=near,
                explicit_A1_rescue=dict(E7_contact=[5, 3], A1_contact_label=5,
                                       total_cost=full_cost.encoded(), conjugate_cost=full_cost.conjugate().encoded(),
                                       finite_E7_replacements=finite),
                conclusion='Every spectator rescue uses only the already covered marked pairs among3,4,5.')


def run():
    primary = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    records = [profile_review(r) for r in primary['records']]
    quads = quad_review(primary)
    spectator = spectator_checks(records)
    saved_spectator = primary['spectator_lemma']
    assert spectator['retained_component_sizes'] == saved_spectator['E7_leaf5_component_sizes']
    assert spectator['forced_rescue_nodes'] == saved_spectator['possible_marked_nodes']
    rescue = spectator['explicit_A1_rescue']
    assert Q5(*rescue['total_cost']) == decode_phi(saved_spectator['explicit_A1_rescue_total_cost'])
    assert Q5(*rescue['conjugate_cost']) == decode_phi(saved_spectator['rescue_conjugate_cost'])
    assert rescue['finite_E7_replacements'] == saved_spectator['rescued_finite_deletions']
    # phi(m)<=4 forces prime factors2,3,5 and exponents at most3,1,1.
    # This finite divisor enumeration is the consequence of that analytic
    # reduction, not a finite-label cutoff.
    cyclotomic_candidates = []
    for a, b, c in product(range(4), range(2), range(2)):
        m = 2**a*3**b*5**c
        totient = (2**(a-1) if a else 1)*(2 if b else 1)*(4 if c else 1)
        if m >= 2 and totient <= 4:
            cyclotomic_candidates.append(m)
    assert sorted(cyclotomic_candidates) == [2,3,4,5,6,8,10,12]
    # Independent squareclass checks supporting the four complete classes.
    phi = Q5(F(1, 2), F(1, 2))
    for q in [Q5(2), Q5(3), Q5(6)]:
        assert sqrt_q5(q) is None
    nonsquares = []
    for m in [1, 2, 3, 6]:
        q = m*(2+phi)
        norm = q.a*q.a-5*q.b*q.b
        assert rational_sqrt(norm) is None and sqrt_q5(q) is None
        nonsquares.append(str(norm))
    paths = [f'scripts/{TAG}.py', f'scripts/{PRIMARY}.py', f'results/{PRIMARY}.json',
             'scripts/fresh_compact29_small_faces_20261015.py',
             'results/fresh_compact29_small_faces_20261015.json']
    for p, h in primary['sources'].items():
        assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    return dict(status='PASS_INDEPENDENT_FRACTION_QSQRT5_TRIANGLES_SPECTATORS_AND_23_LITERAL_QUAD_OBSTRUCTIONS',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                scope=dict(complete_finite_triangle_profiles=True, orthogonal_spectators_allowed=True,
                           literal_quad_recipes_only=23, exhaustive_quad_classification=False,
                           dimension_bound_claim=False, producer_imports=False, symbolic_library=False),
                records=records, literal_quadrilaterals=quads, spectator_lemma=spectator,
                real_cyclotomic_degree_candidates=sorted(cyclotomic_candidates),
                complete_finite_label_alphabet=[2,3,4,5,6,10],
                squareclass_norms=nonsquares,
                proof_audit=['All finite labels reduce analytically to2,3,4,5,6,10.',
                             'Connected E tree cycle products require one of four squareclasses.',
                             'Conjugate-positive Schur cost excludes pure4/6/10 classes.',
                             'Every3/5 profile and all finite replacement inequalities are independently reconstructed.',
                             'The unique subcritical E7 profile can be rescued, but finite H3/H4 geometry confines its marks to3,4,5.',
                             'The old quadrilateral obstruction checks only the23 literal recipes; no all-quad conclusion follows.'])


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
        print('PASS saved replay')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
        print(result['status'])
    print(json.dumps(dict(profile_counts=[r['counts']['accepted'] for r in result['records']],
                          covered_pairs=[len(r['covered_pairs']) for r in result['records']],
                          literal_quads=len(result['literal_quadrilaterals']))))
