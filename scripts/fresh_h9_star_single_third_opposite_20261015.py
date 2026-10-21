#!/usr/bin/env python3
"""A single-own third leaf and opposite single anchors: complete rank bounds."""
import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_single_third_opposite_20261015'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH = 'results/fresh_h9_star_single_third_20261014_search_batch_013.json'
S = [3, 4, 6, 7, 9, 10, 12, 13]
CORE = sorted(S + [11, 14, 17])
R2, R3, R5, R10 = s.sqrt(2), s.sqrt(3), s.sqrt(5), s.sqrt(10)


def sha(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def enc(expr):
    expr = s.expand(expr)
    b, c, d = [expr.coeff(r) for r in (R2, R5, R10)]
    a = s.expand(expr - b*R2 - c*R5 - d*R10)
    assert all(x.is_Rational for x in (a, b, c, d)), expr
    return [str(x) for x in (a, b, c, d)]


def sign5(a, b):
    if a == 0:
        return (b > 0) - (b < 0)
    if b == 0:
        return (a > 0) - (a < 0)
    if (a > 0) == (b > 0):
        return (a > 0) - (a < 0)
    v = a*a - 5*b*b
    return ((a > 0) - (a < 0)) * ((v > 0) - (v < 0))


def sign(expr):
    a, b, c, d = map(Fraction, enc(expr))
    sa, sb = sign5(a, c), sign5(b, d)
    if not sa:
        return sb
    if not sb or sa == sb:
        return sa
    # (a+c*sqrt5)^2 - 2*(b+d*sqrt5)^2.
    difference = sign5(a*a + 5*c*c - 2*b*b - 10*d*d, 2*a*c - 4*b*d)
    return sa*difference


def run():
    seen = {}

    def read(path):
        if path in seen:
            return seen[path]
        data = json.loads((ROOT / path).read_text())
        seen[path] = data
        for p, h in data.get('sources', {}).items():
            assert sha(p) == h, p
            if p.endswith('.json'):
                read(p)
        return data

    checkpoint_path = 'results/polar_extensions_and_discriminant55_obstruction_2026_10_15_audit.json'
    checkpoint = read(checkpoint_path)
    assert checkpoint['status'] == 'PASS_RECURSIVE_SOURCES_FROZEN_PREDECESSOR_AND_LOCAL_LINKS'
    geom, batch = read(GEOM)['representatives'][0]['geometry'], read(BATCH)
    assert batch['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert batch['cumulative_additional_nodes'] == 26951
    vertices = [set(v['incident']) for v in geom['vertices']]
    minimal = {frozenset(m) for m in geom['minimal_nonfaces']}
    assert any(set(S) <= v for v in vertices)
    assert all(any(set(p) <= v for v in vertices) for p in combinations(CORE, 2))

    def face(nodes):
        return dict(support=sorted(nodes), actual_vertex=next(i for i, v in enumerate(vertices)
                                                            if set(nodes) <= v))

    guards = []
    cos = {3: s.Rational(1, 2), 4: R2/2, 6: R3/2}

    def guard(nodes, edges, kernel, kind, conditions):
        G = s.eye(len(nodes))
        for i, j, m in edges:
            G[i, j] = G[j, i] = -cos[m]
        w = s.Matrix(kernel)
        assert all(t.is_positive for t in w)
        assert all(s.simplify(t) == 0 for t in G*w)
        guards.append(dict(kind=kind, ordered_support=nodes, edges=edges,
                           positive_kernel=list(map(str, w)),
                           gram=[[str(t) for t in row] for row in G.tolist()],
                           conditional=conditions, **face(nodes)))

    def path(nodes, labels, kernel, kind, conditions):
        guard(nodes, [[i, i+1, m] for i, m in enumerate(labels)], kernel, kind, conditions)

    def star(nodes, kind, conditions):
        guard(nodes, [[0, 1, 4], [0, 2, 3], [0, 3, 3]], [2, R2, 1, 1], kind, conditions)

    def fork(nodes, kind, conditions):
        guard(nodes, [[0, 1, 3], [0, 2, 3], [0, 3, 3], [3, 4, 4]],
              [2, 1, 1, 2, R2], kind, conditions)

    for i, j in combinations(range(4), 2):
        for a, b in product(S[2*i:2*i+2], S[2*j:2*j+2]):
            aa = next(x for x in S[2*i:2*i+2] if x != a)
            bb = next(x for x in S[2*j:2*j+2] if x != b)
            path([aa, a, b, bb], [4, 3, 4], [1, R2, R2, 1], 'cross_block_right',
                 dict(trial_nonright=[a, b]))
    for anchor in (11, 17):
        for pair in ((6, 7), (9, 10)):
            for contact in pair:
                other = next(x for x in pair if x != contact)
                path([other, contact, anchor], [6, 3], [R3, 2, 1], 'high_middle_anchor_right',
                     dict(trial_nonright=[contact, anchor]))
        for contact in (3, 4):
            path([7-contact, contact, anchor], [6, 3], [R3, 2, 1], 'owner_foreign_requires_internal4or5',
                 dict(trial_nonright=[contact, anchor], trial_owner_internal_at_least=6))
        path([4, 14, anchor], [4, 4], [1, R2, 1], 'third_anchor_label_at_most_three',
             dict(trial_label_at_least_four=[14, anchor]))
        guard([4, 14, anchor], [[0, 1, 3], [1, 2, 3], [0, 2, 3]], [1, 1, 1],
              'nonright_owner_foreign_forces_own_third_right',
              dict(owner_foreign=[4, anchor], trial_nonright=[14, anchor]))
    path([13, 11, 17], [3, 6], [1, 2, R3], 'mutual_anchor_label_at_most_five',
         dict(trial_mutual_label_at_least=6))
    for root in (3, 4):
        star([root, 7-root, 11, 17], 'same_owner_two_foreign_forbidden',
             dict(foreign_contacts=[[root, 11], [root, 17]]))
    for r, t in ((3, 4), (4, 3)):
        path([13, 11, r, t, 17], [3, 3, 4, 3], [1, 2, 3, 2*R2, R2],
             'different_owner_two_foreign_forbidden', dict(foreign_contacts=[[r, 11], [t, 17]]))
    star([14, 4, 11, 17], 'both_third_contacts_forbidden',
         dict(third_contacts=[[14, 11], [14, 17]]))
    ownroot = {11: 13, 17: 12}
    for anchor, other in ((11, 17), (17, 11)):
        fork([anchor, ownroot[anchor], other, 14, 4], 'mutual_and_third_contact_forbidden',
             dict(third_contact=[14, anchor], trial_mutual_nonright=[11, 17]))
        for root in (3, 4):
            fork([anchor, ownroot[anchor], other, root, 7-root], 'mutual_and_owner_foreign_forbidden',
                 dict(foreign_contact=[root, anchor], trial_mutual_nonright=[11, 17]))
        path([anchor, 4, 14, other, ownroot[other]], [3, 4, 3, 3], [R2, 2*R2, 3, 2, 1],
             'nonright_owner_foreign_forces_other_third_right',
             dict(foreign_contact=[4, anchor], trial_third_nonright=[14, other]))
    assert len(guards) == 54

    budgets = [dict(leaf=a, removed_own_root=25-ownroot[a], retained_own_root=ownroot[a],
                    own_cost_lower='1/4', **face((set(S)-{25-ownroot[a]})|{a})) for a in (11, 17)]
    budgets.append(dict(leaf=14, removed_own_root=3, retained_own_root=4,
                        own_cost_lower='(3+sqrt(5))/8 > 1/2', **face((set(S)-{3})|{14})))
    assert len(budgets) == 3 and all(len(r['support']) == 8 for r in budgets)
    assert sign((3+R5)/8-s.Rational(1, 2)) > 0

    # Complete discrete contact domain, before any rank calculation.
    literal_patterns = []
    for fx, fy, dx, dy, mu in product((None, 3, 4), (None, 3, 4), (2, 3), (2, 3), (2, 3, 4, 5)):
        if fx is not None and fy is not None:
            continue
        if dx == dy == 3:
            continue
        if mu > 2 and (fx is not None or fy is not None or dx == 3 or dy == 3):
            continue
        if (fx == 4 or fy == 4) and (dx == 3 or dy == 3):
            continue
        for owner in ([None] if fx is None and fy is None else [4, 5]):
            literal_patterns.append([fx, fy, dx, dy, mu, owner])
    assert len(literal_patterns) == 22

    # Exact rational bounds on cos(pi/7), and a safe label8 lower bound.
    u = s.symbols('u')
    cubic = 8*u**3 + 4*u**2 - 4*u - 1
    assert s.expand(s.chebyshevt(7, u)-1-(u-1)*cubic**2) == 0
    lower_u = s.Rational(31, 50)
    upper_u = 2*s.Rational(901, 1000)**2-1
    assert cubic.subs(u, lower_u) < 0 < cubic.subs(u, upper_u)
    assert s.diff(cubic, u).subs(u, s.Rational(1, 2)) > 0
    assert s.diff(cubic, u, 2).subs(u, s.Rational(1, 2)) > 0
    assert (1+lower_u)/2 == s.Rational(81, 100)
    assert sign((2+R2)/4-s.Rational(23, 25)**2) > 0
    assert s.Rational(3, 4) > s.Rational(6, 7)**2

    beta, c, t, x, y, v = s.symbols('beta c t x y v', real=True)
    q = 1-c*c
    D = 1/(1-t*t)
    mu_cos = {2: 0, 3: s.Rational(1, 2), 4: R2/2, 5: (1+R5)/4}
    canonical = [(None, 0, 0, m, None) for m in (2, 3, 4, 5)]
    canonical += [(None, 0, 1, 2, None), (None, 1, 0, 2, None)]
    canonical += [(0, a, b, 2, owner) for a, b in ((0, 0), (0, 1), (1, 0)) for owner in (4, 5)]
    canonical += [(1, 0, 0, 2, owner) for owner in (4, 5)]
    assert len(canonical) == 14
    # Naming the unique foreign anchor first preserves its opposite own profile
    # after swapping the two last roots; all literal role supports were checked.
    reduced = set()
    for fx, fy, dx, dy, mu, owner in literal_patterns:
        if fy is not None:
            fx, fy, dx, dy = fy, fx, dy, dx
        reduced.add((None if fx is None else fx-3, int(dx == 3), int(dy == 3), mu, owner))
    assert reduced == set(canonical)

    cases = []
    positive = negative = 0
    for foreign, dx, dy, mutual, owner in canonical:
        F = 0 if foreign is None else 1/(4*q)
        T = 0 if foreign is None else beta*(c if foreign == 0 else 1)/(2*q)
        R = s.Matrix([[1-D/4-F, -mu_cos[mutual]-t*D/4, -s.Rational(dx, 2)-T],
                      [-mu_cos[mutual]-t*D/4, 1-D/4, -s.Rational(dy, 2)],
                      [-s.Rational(dx, 2)-T, -s.Rational(dy, 2), 1-beta*beta/q]])
        P = s.cancel(q*(1-t*t)*R.det())
        assert not s.denom(P).free_symbols
        assert s.Poly(P, t).degree() <= 2
        exceptional = (foreign, dx, dy, mutual, owner) == (0, 0, 1, 2, 4)
        branches = ['third_ge6', 'third5_last_ge8', 'third5_last7'] if exceptional else ['uniform']
        for branch in branches:
            lo, hi, target = s.Rational(9, 10), s.Integer(1), 1
            beta_sub = (1+R5)/4+x
            if branch == 'third_ge6':
                beta_sub = s.Rational(6, 7)+x
            elif branch == 'third5_last_ge8':
                beta_sub = (1+R5)/4
                lo = s.Rational(23, 25)
            elif branch == 'third5_last7':
                beta_sub = (1+R5)/4
                hi, target = s.Rational(901, 1000), -1
            c_sub = R2/2+y if owner is None else R2/2 if owner == 4 else (1+R5)/4
            Q = s.Poly(s.expand(P.subs({beta: beta_sub, c: c_sub, t: lo+(hi-lo)*v})), x, y, v)
            grouped = {}
            for (a, b, k), value in Q.terms():
                grouped.setdefault((a, b), {})[k] = value
            coefficients = []
            for powers, row in sorted(grouped.items()):
                for k in range(3):
                    value = s.expand(sum(row.get(j, 0)*s.binomial(k, j)/s.binomial(2, j)
                                         for j in range(k+1)))
                    assert sign(value) == target
                    coefficients.append(dict(nonnegative_power=list(powers), bernstein_index=k,
                                             coefficient=enc(value)))
            assert (0, 0) in grouped
            if target == 1:
                positive += len(coefficients)
            else:
                negative += len(coefficients)
            cases.append(dict(foreign_owner_root=foreign, third_anchor_cosines=[str(s.Rational(dx, 2)), str(s.Rational(dy, 2))],
                              mutual_anchor_label=mutual, owner_internal=owner,
                              branch=branch, target_sign=target, scaled_determinant=str(s.expand(P)),
                              substitutions=dict(beta=str(beta_sub), c=str(c_sub), t=str(lo+(hi-lo)*v)),
                              tail_interval=[str(lo), str(hi)], nonnegative_variables=['x', 'y'],
                              bernstein_degree=2, coefficients=coefficients))
    assert len(cases) == 16 and positive == 165 and negative == 3

    families, orientations = [], []
    for owner_reverse, last_reverse in product((False, True), repeat=2):
        perm = list(range(18))
        if owner_reverse:
            perm[3], perm[4] = 4, 3
        if last_reverse:
            perm[12], perm[13] = 13, 12
        assert {frozenset(perm[i] for i in v) for v in vertices} == set(map(frozenset, vertices))
        assert {frozenset(perm[i] for i in m) for m in minimal} == minimal
        fixed = [[3, 14, 2], [11, 12, 2], [11, 13, 3], [12, 17, 3], [13, 17, 2]]
        bounds = [[3, 4, 4], [6, 7, 7], [9, 10, 7], [12, 13, 7], [4, 14, 5]]
        image = {k: sorted([sorted((perm[a], perm[b]))+[m] for a, b, m in rows])
                 for k, rows in [('fixed', fixed), ('bounds', bounds)]}
        families.append(dict(id='single_third_opposite_tail'+('_owner_swap' if owner_reverse else '')
                             +('_last_swap' if last_reverse else ''),
                             origin='single_third_opposite_anchor_rank_three', **image))
        orientations.append(perm)
    labels = {tuple(r[:2]): r[2] for r in batch['solution']['pair_labels']}
    assert any(all(labels[a, b] == m for a, b, m in f['fixed']) and
               all(labels[a, b] >= m for a, b, m in f['bounds']) for f in families)
    faces = {tuple(sorted(v & set(CORE))) for v in vertices}
    maximal = sorted(f for f in faces if not any(set(f) < set(g) for g in faces))
    assert len(maximal) == 9
    paths = [GEOM, BATCH, checkpoint_path, f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_LABEL_SINGLE_THIRD_OPPOSITE_ANCHOR_CORE', case=0,
                sources={p: sha(p) for p in paths}, rank_core=CORE,
                complete_maximal_core_faces=[list(f) for f in maximal],
                positive_kernel_guards=guards, projection_budgets=budgets,
                literal_patterns=literal_patterns, canonical_patterns=[list(p) for p in canonical],
                rational_tail=dict(cubic=str(cubic), chebyshev_identity='T7(u)-1=(u-1)*cubic(u)^2',
                                   lower_u=str(lower_u), lower_value=str(cubic.subs(u, lower_u)),
                                   upper_u=str(upper_u), upper_value=str(cubic.subs(u, upper_u)),
                                   conclusion='9/10<cos(pi/7)<901/1000; cos(pi/m)>=cos(pi/8)>23/25 for m>=8',
                                   third_ge6_lower='cos(pi/k)>=sqrt3/2>6/7 for k>=6'),
                residual_order=['anchor_with_possible_foreign', 'other_anchor', 'third'],
                residual_matrix=[['1-D/4-F', '-mu-t*D/4', '-delta_x-T'],
                                 ['-mu-t*D/4', '1-D/4', '-delta_y'],
                                 ['-delta_x-T', '-delta_y', '1-beta^2/(1-c^2)']],
                continuous_cases=cases, positive_coefficients=positive, negative_coefficients=negative,
                exceptional_partition='Only owner4, foreign at right ownerroot, delta_x=0,delta_y=1/2,mu=0 '
                                      'needs a split: thirdlabel>=6; third5,last>=8; third5,last7. '
                                      'The last piece is strictly negative; no false continuous positivity claim.',
                rank_obstruction='Multiply the3x3 residual determinant by(1-c^2)*(1-t^2)>0. '
                                 'Every discrete label case belongs to a certified strict-sign polynomial '
                                 'region. Three residuals cannot have nonzero Gram determinant in the '
                                 'two-dimensional complement of the actual PD8 support in Lorentz dimension10.',
                literal_orientations=orientations, families=families, covered_saved_cursor=26951,
                scope='Star0 owner internal>=4, other three internal labels>=7, third singleown(2,>=5), '
                      'anchors opposite singleownpi3 contacts. All finite integer labels in these ranges '
                      'and every noncore contact/distance are unrestricted. No whole-incidence exclusion.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    out = ROOT / f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], len(data['continuous_cases']), data['positive_coefficients'], data['negative_coefficients'])


