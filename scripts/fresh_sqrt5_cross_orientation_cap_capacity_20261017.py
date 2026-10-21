#!/usr/bin/env python3
"""Tighten two-orientation cap capacities using exact cross two-cycles."""
import argparse
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path

import fresh_sqrt5_cross_coface_cap_capacity_20261017_review as V

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_cross_orientation_cap_capacity_20261017'
TAGS = dict(old='fresh_sqrt5_orbit_integral_cap_capacity_20261017',
    incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
    model='fresh_sqrt5_cross_coface_cap_full_model_20261017')
K = V.K


def integral(z):
    basis = (z.a-z.b).denominator == 1 and (2*z.b).denominator == 1
    polynomial = (2*z.a).denominator == 1 and (z.a*z.a-5*z.b*z.b).denominator == 1
    assert basis == polynomial
    return polynomial


def read(tag):
    d = json.loads((ROOT/f'results/{tag}.json').read_text())
    for obj in (d, d.get('prepared', {})):
        for p, h in obj.get('sources', {}).items():
            assert sha256((ROOT/p).read_bytes()).hexdigest() == h, p
    return d


def run():
    D = {k: read(t) for k,t in TAGS.items()}
    old, I, M = (D[k] for k in ('old', 'incidence', 'model'))
    assert old['status'] == 'PASS_EXACT_ORBIT_AND_INTEGRALITY_CAP_CAPACITIES'
    dual = next(c for c in M['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    alpha = {i:F(a) for i,a in dual['face_multipliers']}
    occurrences = {o['base_index']:o for o in old['occurrences']}
    columns = deepcopy(old['columns'])
    required = defaultdict(F)
    cache = {}
    checks = []
    for c in columns:
        c['inherited_capacity'] = c['capacity']
        if c['capacity'] == 2:
            rows = [occurrences[i] for i in c['source_bases']]
            assert rows and all(o['orientation_orbit'] == c['orientation_orbit'] for o in rows)
            assert all(o['core'] == c['core'] and o['dotted_squared_product'] == c['dotted_squared_product'] for o in rows)
            ori = c['orientation_orbit']
            assert len(ori) == 2 and ori[0] != ori[1]
            g = I['cores'][c['core']]
            ports = [[i for i,e in enumerate(x[0]) if e[0] == 3] for x in ori]
            assert all(len(x[0]) == g['nodes'] for x in ori)
            assert all(e[0] in (2,3) for x in ori for e in x[0])
            overlap = sorted(set(ori[0][1]) & set(ori[1][1]))
            tau = K(*c['dotted_squared_product'])
            h = K(*c['repeated_wall_product'])
            assert h == 2-tau/2 and h.sign() < 0
            if all(e[2][0] in (3,5) for e in g['edges']):
                if c['core'] not in cache:
                    G, scale = V.physical_core(V.core_graph(g))
                    assert scale is None
                    J, det = V.inverse(G)
                    assert det.sign() < 0
                    n = len(G)
                    assert all(sum((G[i][k]*J[k][j] for k in range(n)), K()) == int(i == j)
                               for i in range(n) for j in range(n))
                    cache[c['core']] = J
                J = cache[c['core']]
                norms = [sum((J[i][j] for i in p for j in p), K()) for p in ports]
                cross = sum((J[i][j] for i in ports[0] for j in ports[1]), K())
                cross_square = cross*cross
                method = 'physical_Qsqrt5_core_inverse'
                product = cross.out()
            else:
                assert g['nodes'] == 2 and len(g['edges']) == 1
                e = g['edges'][0][2]
                assert e[0] == 0 and ports in ([[0],[1]], [[1],[0]])
                d = K(*e[1:])
                assert (d-4).sign() > 0
                # G=[[2,-sqrt(d)],[-sqrt(d),2]]. Its inverse is
                # [[2,sqrt(d)],[sqrt(d),2]]/(4-d). Avoid adjoining sqrt(d).
                norms = [K(2)/(4-d)]*2
                cross_square = d/((d-4)*(d-4))
                method = 'two_node_dotted_symbolic_inverse'
                product = None
            assert norms == [h,h]
            nonintegral = not integral(cross_square)
            tightened = bool(overlap) and nonintegral
            checks.append(dict(column=c['index'], core=c['core'], ports=ports,
                common_base_vertices=overlap, method=method,
                projection_norms=[z.out() for z in norms],
                cross_product=product, cross_two_cycle=cross_square.out(),
                cross_two_cycle_trace=str(2*cross_square.a),
                cross_two_cycle_norm=str(cross_square.a*cross_square.a-5*cross_square.b*cross_square.b),
                nonintegral_two_cycle=nonintegral, tightened=tightened))
            if tightened:
                c['capacity'] = 1
        co = c['positive_face_coefficients']
        load = sum((m*alpha.get(fi,F()) for fi,m in co), F())
        mass = alpha.get(c['simplex_face'],F())
        assert load == F(c['load']) and mass == F(c['initial_simplex_mass'])
        c['slack'] = str(c['capacity']*mass-load)
        required[c['simplex_face']] = max(required[c['simplex_face']], load/c['capacity'])
    strengthened = [r['column'] for r in checks if r['tightened']]
    violations = [c['index'] for c in columns if F(c['slack']) < 0]
    new_violations = sorted(set(violations)-set(old['violating_columns']))
    assert len(checks) == 243 and len(strengthened) == 215
    assert len(violations) == 42 and len(new_violations) == 23
    assert set(old['violating_columns']) <= set(violations)
    assert len(columns) == 815
    paths = [f'scripts/{TAG}.py', str(Path(V.__file__).relative_to(ROOT))]
    paths += [f'results/{t}.json' for t in TAGS.values()]
    return dict(status='PASS_EXACT_CROSS_ORIENTATION_CAP_CAPACITIES',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
        counts=dict(columns=len(columns), inherited_capacity_two=len(checks),
            tightened=len(strengthened), positive_columns=sum(F(c['load'])>0 for c in columns),
            tightened_positive_columns=sum(F(columns[i]['load'])>0 for i in strengthened),
            violations=len(violations), new_violations=len(new_violations)),
        capacity_histogram=sorted(Counter(c['capacity'] for c in columns).items()),
        pair_checks=checks, strengthened_columns=strengthened,
        columns=columns, violating_columns=violations, new_violating_columns=new_violations,
        required_simplex_masses=[[i,str(alpha.get(i,F())),str(a)] for i,a in sorted(required.items())],
        input_bound=dual['bound'], new_repaired_bound=None, optimizer_invocations=0,
        scope='Retain all815 inherited selected cap columns; tighten215 capacities from2 to1. Within an orientation use the inherited obstruction. Distinct orientations at one coface have intersecting base subsets; at different orthogonal A1 cofaces they have a nonintegral arithmetic two-cycle. No optimum or dimension bound is claimed.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify',action='store_true'); a = p.parse_args()
    r = json.loads(json.dumps(run())); out = ROOT/f'results/{TAG}.json'
    if a.verify: assert r == json.loads(out.read_text())
    else: out.write_text(json.dumps(r,indent=2)+'\n')
    print(r['status'],r['counts'],flush=True)
