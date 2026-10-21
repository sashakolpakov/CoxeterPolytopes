#!/usr/bin/env python3
"""A single-own third leaf and shared single anchors: eleven rank polynomials."""
import argparse
import json
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_star_single_third_20261014'
GEOM = 'results/fresh_h9_tree_cube_wedges_20261014.json'
BATCH = 'results/fresh_h9_star_third34_20261014_search_batch_012.json'
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

    geom, batch = read(GEOM)['representatives'][0]['geometry'], read(BATCH)
    assert batch['status'] == 'FULL_MANDATORY_BRIDGE_SURVIVOR'
    assert batch['cumulative_additional_nodes'] == 26241
    vertices = [set(v['incident']) for v in geom['vertices']]
    minimal = {frozenset(m) for m in geom['minimal_nonfaces']}
    assert any(set(S) <= v for v in vertices)
    assert all(any(set(p) <= v for v in vertices) for p in combinations(CORE, 2))

    def face(nodes):
        nodes = sorted(nodes)
        return dict(support=nodes, actual_vertex=next(i for i, v in enumerate(vertices)
                                                    if set(nodes) <= v))

    guards = []
    cosine = {3: s.Rational(1, 2), 4: R2/2, 6: R3/2}

    def path(nodes, labels, kernel, kind, conditional=None):
        G = s.eye(len(nodes))
        for i, m in enumerate(labels):
            G[i, i+1] = G[i+1, i] = -cosine[m]
        w = s.Matrix(kernel)
        assert all(t.is_positive for t in w)
        assert all(s.simplify(t) == 0 for t in G*w)
        guards.append(dict(kind=kind, ordered_support=nodes, labels=labels,
                           cycle=False, positive_kernel=list(map(str, w)),
                           gram=[[str(t) for t in row] for row in G.tolist()],
                           conditional=conditional, **face(nodes)))

    def cycle(nodes, kind, conditional=None):
        G = s.eye(len(nodes))
        for i in range(len(nodes)):
            j = (i+1) % len(nodes)
            G[i, j] = G[j, i] = -s.Rational(1, 2)
        assert G*s.ones(len(nodes), 1) == s.zeros(len(nodes), 1)
        guards.append(dict(kind=kind, ordered_support=nodes, labels=[3]*len(nodes),
                           cycle=True, positive_kernel=['1']*len(nodes),
                           gram=[[str(t) for t in row] for row in G.tolist()],
                           conditional=conditional, **face(nodes)))

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
                path([other, contact, anchor], [6, 3], [R3, 2, 1],
                     'high_middle_forces_anchor_right', dict(trial_nonright=[contact, anchor]))
        for contact in (3, 4):
            other = 7-contact
            path([other, contact, anchor], [6, 3], [R3, 2, 1],
                 'owner_foreign_requires_internal_at_most_five',
                 dict(trial_nonright=[contact, anchor], trial_owner_internal_at_least=6))
        path([4, 14, anchor], [4, 4], [1, R2, 1], 'third_anchor_label_at_most_three',
             dict(trial_label_at_least_four=[14, anchor]))
        cycle([4, 14, anchor], 'nonright_own_root_foreign_forces_third_right',
              dict(owner_foreign=[4, anchor], trial_nonright=[14, anchor]))
    cycle([11, 13, 17], 'shared_anchors_are_right', dict(trial_nonright=[11, 17]))
    for r, t in product((3, 4), repeat=2):
        nodes = [r, 11, 13, 17] if r == t else [r, 11, 13, 17, t]
        cycle(nodes, 'at_most_one_owner_foreign_across_both_anchors',
              dict(two_foreign_contacts=[[r, 11], [t, 17]]))
    cycle([14, 11, 13, 17], 'third_cannot_contact_both_anchors',
          dict(two_third_contacts=[[14, 11], [14, 17]]))
    for anchor, other in ((11, 17), (17, 11)):
        cycle([4, 14, other, 13, anchor], 'nonright_own_root_foreign_forces_other_third_right',
              dict(owner_foreign=[4, anchor], trial_nonright=[14, other]))
    assert len(guards) == 48

    budgets = [dict(leaf=a, removed_own_root=12, retained_own_root=13,
                    own_cost_lower='1/4', **face((set(S)-{12})|{a})) for a in (11, 17)]
    budgets.append(dict(leaf=14, removed_own_root=3, retained_own_root=4,
                        own_cost_lower='(3+sqrt(5))/8 > 1/2',
                        **face((set(S)-{3})|{14})))
    assert all(len(r['support']) == 8 for r in budgets)
    assert sign((3+R5)/8-s.Rational(1, 2)) > 0

    # A rigorous label-independent tail bound; no decimal cosine estimate.
    u = s.symbols('u')
    cubic = 8*u**3 + 4*u**2 - 4*u - 1
    assert s.expand(s.chebyshevt(7, u)-1-(u-1)*cubic**2) == 0
    at = s.Rational(31, 50)
    assert cubic.subs(u, at) < 0
    assert s.diff(cubic, u).subs(u, s.Rational(1, 2)) > 0
    assert s.diff(cubic, u, 2).subs(u, s.Rational(1, 2)) > 0
    assert (1+at)/2 == s.Rational(81, 100)
    assert 1/(1-s.Rational(81, 100)) == s.Rational(100, 19)

    b, c, D, x, y, z = s.symbols('beta c D x y z', real=True)
    den = 1-c*c
    cases = []
    total_coefficients = 0
    patterns = []
    for foreign in (None, 0, 1):
        for dx, dy in product((s.Rational(0), s.Rational(1, 2)), repeat=2):
            if dx and dy:
                continue
            if foreign == 1 and (dx or dy):
                continue
            patterns.append((foreign, dx, dy))
    assert len(patterns) == 7
    for foreign, dx, dy in patterns:
        F = 0 if foreign is None else 1/(4*den)
        T = 0 if foreign is None else b*(c if foreign == 0 else 1)/(2*den)
        G = s.Matrix([[1-D/4-F, -D/4, -dx-T],
                      [-D/4, 1-D/4, -dy],
                      [-dx-T, -dy, 1-b*b/den]])
        P = s.cancel(den*G.det())
        assert s.Poly(P, D).degree() == 1
        owner_values = [None] if foreign is None else [4, 5]
        for internal in owner_values:
            replacement = {b: (1+R5)/4+x, D: s.Rational(100, 19)+z}
            if internal is None:
                replacement[c] = R2/2+y
                variables = [x, y, z]
            else:
                replacement[c] = R2/2 if internal == 4 else (1+R5)/4
                variables = [x, z]
            Q = s.Poly(s.expand(P.subs(replacement)), *variables)
            terms = [dict(exponents=list(mon), coefficient=enc(coef))
                     for mon, coef in Q.terms()]
            assert all(sign(coef) > 0 for mon, coef in Q.terms())
            assert sign(Q.coeff_monomial(1)) > 0
            assert len(terms) == (10 if internal is None else 6)
            total_coefficients += len(terms)
            cases.append(dict(foreign_owner_root=foreign, third_anchor_cosines=[str(dx), str(dy)],
                              owner_internal=internal if internal is not None else 'all finite >=4',
                              variable_order=list(map(str, variables)),
                              scaled_determinant=str(s.expand(P)),
                              shifted_terms=terms, shifted_constant=enc(Q.coeff_monomial(1))))
    assert len(cases) == 11 and total_coefficients == 78

    families, orientations = [], []
    for owner_reverse, last_reverse in product((False, True), repeat=2):
        perm = list(range(18))
        if owner_reverse:
            perm[3], perm[4] = 4, 3
        if last_reverse:
            perm[12], perm[13] = 13, 12
        assert {frozenset(perm[i] for i in v) for v in vertices} == set(map(frozenset, vertices))
        assert {frozenset(perm[i] for i in m) for m in minimal} == minimal
        fixed = [[3, 14, 2], [11, 12, 2], [11, 13, 3], [12, 17, 2], [13, 17, 3]]
        bounds = [[3, 4, 4], [6, 7, 7], [9, 10, 7], [12, 13, 7], [4, 14, 5]]
        image = {k: sorted([sorted((perm[a], perm[b]))+[m] for a, b, m in rows])
                 for k, rows in [('fixed', fixed), ('bounds', bounds)]}
        families.append(dict(id='single_third_shared_tail'+('_owner_swap' if owner_reverse else '')
                             +('_last_swap' if last_reverse else ''),
                             origin='single_third_shared_anchor_rank_three', **image))
        orientations.append(perm)
    labels = {tuple(r[:2]): r[2] for r in batch['solution']['pair_labels']}
    assert any(all(labels[a, b] == m for a, b, m in f['fixed'])
               and all(labels[a, b] >= m for a, b, m in f['bounds']) for f in families)
    faces = {tuple(sorted(v & set(CORE))) for v in vertices}
    maximal = sorted(f for f in faces if not any(set(f) < set(g) for g in faces))
    assert len(maximal) == 9
    paths = [GEOM, BATCH, f'scripts/{TAG}.py']
    return dict(status='PASS_ALL_LABEL_SINGLE_THIRD_SHARED_ANCHOR_CORE', case=0,
                sources={p: sha(p) for p in paths}, rank_core=CORE,
                complete_maximal_core_faces=[list(f) for f in maximal],
                positive_kernel_guards=guards, projection_budgets=budgets,
                tail=dict(cubic=str(cubic), argument='u=cos(2*pi/7)>1/2; the cubic is increasing on[1/2,1]',
                          rational_test='31/50', value_at_test=str(cubic.subs(u, at)),
                          consequence='cos(pi/m)^2>81/100 and D=1/sin(pi/m)^2>100/19 for every finite m>=7'),
                domain='Blocks are mutually orthogonal. Third is right to all foreign blocks. '
                       'Anchors are mutually right and right to both high middle blocks. Each anchor '
                       'has at most one foreign pi/3 contact; together they have at most one owner '
                       'contact. Such a contact forces owner internal4or5. Third-anchor labels are2or3; '
                       'both cannot be3, and if the foreign contact hits the nonright third ownroot, '
                       'both third-anchor contacts are right.',
                residual_order=['anchor_with_possible_foreign', 'other_anchor', 'third'],
                residual_matrix=[['1-D/4-F', '-D/4', '-delta_x-T'],
                                 ['-D/4', '1-D/4', '-delta_y'],
                                 ['-delta_x-T', '-delta_y', '1-beta^2/(1-c^2)']],
                owner_projection='F=0,T=0 if no foreign; otherwise F=1/[4(1-c^2)], '
                                 'T=beta*c/[2(1-c^2)] at right third ownroot, '
                                 'T=beta/[2(1-c^2)] at nonright ownroot.',
                residual_patterns=7, continuous_cases=cases, positive_coefficient_count=78,
                rank_obstruction='Every scaled3x3 determinant is strictly positive after nonnegative '
                                 'shifts. Its denominator1-c^2 is positive. Three residuals cannot '
                                 'have rank3 in the two-dimensional orthogonal complement of PD8 '
                                 'inside the ten-dimensional Lorentz space.',
                literal_orientations=orientations, families=families,
                covered_saved_cursor=26241,
                scope='Star0 owner internal>=4, other three internal labels>=7, third has exactly '
                      'one right owncontact and the other label>=5; both anchors share a single '
                      'ownpi/3 contact. All finite integer labels in these ranges and all other '
                      'contacts/distances are unrestricted. No noncore label hypothesis; '
                      'not a whole-incidence exclusion.')


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
    print(data['status'], len(data['continuous_cases']), data['positive_coefficient_count'])
