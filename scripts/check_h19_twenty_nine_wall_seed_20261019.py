#!/usr/bin/env python3
"""Fixed29-wall algebra and one64-vertex diagnostic walk; no root search.

Coordinator only: bounded parent-lock wrapper, --write or --verify.
Only pure helpers are imported. No old inputs/check/run body is called.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import combinations
import json
from pathlib import Path
import signal
from search_h19_omega2_finite_profiles_20261019 import (
    Q, native as N, decode, dot, mv, digest, inverse, pair_kind)
from check_arithmetic_dn_opposed_fork_models_20261019 import inertia
from check_h19_twenty_eight_wall_seed_20261019 import add, encoded, bounded_walk

ROOT = Path(__file__).resolve().parents[1]
CERT = ROOT/'results/h19_twenty_nine_wall_seed_20261019.json'
PINS = {
    'data/candidates/h19_vinberg_prefix.json':
    '91c16b581e5763d6a283b6818f0e087e941fbb987b6b05a7e283febf1e32d0f1',
    'results/h19_twenty_six_wall_seed_20261019.json':
    'e310723817ce78bce5c8ed1c11f2a210925d6edc3c6fcb8faf443d5458fc011b',
    'results/h19_twenty_seven_wall_seed_20261019.json':
    'b5b529b485a01b9315be6f31bbb122c170652273d3932b1eecece85680b5ce7d',
    'results/h19_twenty_eight_wall_seed_20261019.json':
    '3996be94b07f251685a15650e50321456a4f3ed5548e0ad927d94abba1525b83',
    'results/h19_omega3_finite_profiles_20261019.json':
    'aa85c4c7c5373f66f124702563a55f9c88a7b7a28fe962f78db5bbe01af54c1d',
    'results/h19_next_edge_short_finite_profiles_20261019.json':
    '72a0873fbd6df501162444ef8608c572e51ef218fbac3cfe1254679a30313e0c',
    'scripts/check_h19_twenty_eight_wall_seed_20261019.py':
    'f9ccdb1dbb5fe6a810bf456be56b3f15d1b83be69fd0d080931713eae3738600',
    'scripts/search_h19_next_edge_short_finite_profiles_20261019.py':
    '23f496f097e35b2da90999ebebee78907122ec86fc1f62cc67cc1e0c83f5db57',
    'scripts/search_h19_omega2_finite_profiles_20261019.py':
    '02f82e5ddfb72486593a761d8a1ca86fef8e5ae3ffaaa72e901dccf129967b3e',
    'scripts/check_arithmetic_dn_opposed_fork_models_20261019.py':
    '185d20648368696429b675c1b6c023d30bb1be262ff83d3622e4f48801e1bae9',
    'results/h19_next_edge_short_finite_profiles_completion_20261019.md':
    'ad7795420d28c1efb5a9425ee4844429b85a6176955597550448c1c258529f2b',
}


def inputs():
    for name, expected in PINS.items():
        assert digest(ROOT/name) == expected, name
    read = lambda stem: json.loads((ROOT/'results'/(stem+'_20261019.json')).read_text())
    old, middle, prior, hfound, found = map(read, ('h19_twenty_six_wall_seed',
        'h19_twenty_seven_wall_seed', 'h19_twenty_eight_wall_seed',
        'h19_omega3_finite_profiles', 'h19_next_edge_short_finite_profiles'))
    assert [old['walls'], middle['walls'], prior['walls']] == [26, 27, 28]
    assert prior['status'] == 'PASS_EXACT_28_WALL_AND_BOUNDED_EDGE_DIAGNOSTIC'
    assert found['status'] == 'PASS_BOUNDED_NEXT_EDGE_11_SHORT_PROFILE_SEARCH'
    assert prior['checker_sha256'] == PINS['scripts/check_h19_twenty_eight_wall_seed_20261019.py']
    assert found['source_sha256'] == PINS['scripts/search_h19_next_edge_short_finite_profiles_20261019.py']
    for record in (old, middle, prior, hfound, found):
        assert record['status'].startswith('PASS')
        for name, expected in record['sources'].items():
            assert digest(ROOT/name) == expected, name
    phi = N(0, 1)
    G = [[Q() for _ in range(10)] for _ in range(10)]
    G[0][0], G[0][1], G[1][0], G[1][1] = -2*phi, phi, phi, Q(2)
    for i in range(8):
        G[i+2][i+2] = Q(2)
    for i, j in [(0,1),(1,2),(2,3),(3,4),(4,5),(5,6),(2,7)]:
        G[i+2][j+2] = G[j+2][i+2] = Q(-1)
    assert inertia(G) == [9,1,0] and inertia(G, True) == [10,0,0]
    gi = inverse(G, Q); ci = [row[2:] for row in gi[2:]]
    prefix = json.loads((ROOT/'data/candidates/h19_vinberg_prefix.json').read_text())
    roots = {f'old{i}': [N(v['a'],v['b']) for v in row]
        for i, row in enumerate(prefix['output']['roots_in_vinberg_order']) if i not in (10,20,21)}
    for name, cv in old['new_covectors'].items():
        roots[name] = mv(gi, [decode(v) for v in cv])
    roots['q'] = [decode(v) for v in middle['new_root']]
    assert len(hfound['survivors']) == 1
    roots['h'] = [decode(v) for v in hfound['survivors'][0]['coordinates']]
    assert len(roots) == 28
    assert len(found['survivors']) == 2
    small = next(s for s in found['survivors'] if s['x_phi_coefficients'] == [6,10])
    large = next(s for s in found['survivors'] if s['x_phi_coefficients'] == [12,20])
    assert small['profile'] == 'next_edge_profile_05' and large['profile'] == 'next_edge_profile_02'
    roots['j'] = [decode(v) for v in small['coordinates']]
    cov = {n: mv(G,v) for n,v in roots.items()}
    norms = {n: dot(v,cov[n]) for n,v in roots.items()}
    assert roots['j'][0] == N(6,10)
    assert cov['j'] == [N(-28,-45),Q(),N(5,7),2*phi,Q(),Q(),Q(1),Q(),Q(),Q()]
    assert cov['j'] == [decode(v) for v in small['covector']]
    for n, root in roots.items():
        expected = Q(2) if n in ('q','h','j') else decode(old['norms'][n])
        assert norms[n] == expected and expected.sign() > 0
        assert all(v.integral() for v in root)
        assert all((2*v/expected).integral() for v in cov[n])
    p = [Q(2),-phi]+[Q()]*8
    c = add(p,[Q(),Q(F(1,10))]+[sum(row,Q())/10 for row in ci])
    return phi,G,ci,roots,cov,norms,p,c,prior,small,large


def check():
    phi,G,ci,roots,cov,norms,p,c,prior,small,large = inputs()
    evaluate = lambda v: {n: dot(v,cv) for n,cv in cov.items()}
    square = lambda v: dot(v,mv(G,v))
    cp = evaluate(c)
    assert square(c).sign() < 0 and dot(p,mv(G,c)).sign() < 0
    assert all(v.sign() < 0 for v in cp.values())
    assert cp['j'] == -(123+198*phi)/5
    categories,pairs,newproducts,oldpairs = Counter(),[],{},[]
    for a,b in combinations(roots,2):
        v = dot(roots[a],cov[b])
        if norms[a] == 2: kind = pair_kind(v,norms[b],phi)
        elif norms[b] == 2: kind = pair_kind(v,norms[a],phi)
        else:
            assert v == 0 or (v.sign() < 0 and (v*v-norms[a]*norms[b]).sign() > 0)
            kind = 'right' if v == 0 else 'ultraparallel'
        assert not kind.startswith('reject'), (a,b,kind)
        row = [a,b,v.encoded(),kind]; pairs.append(row); categories[kind] += 1
        if b == 'j': newproducts[a] = [v.encoded(),kind]
        else: oldpairs.append(row)
    pair_hash = lambda rows: hashlib.sha256(json.dumps(rows,sort_keys=True,separators=(',',':')).encode()).hexdigest()
    assert len(pairs) == 406 and len(oldpairs) == 378
    assert pair_hash(oldpairs) == prior['wall_pairs_sha256']
    assert newproducts == small['wall_products']
    for n,root in roots.items():
        point = add(c,root,-cp[n]/norms[n]); values = evaluate(point)
        assert square(point).sign() < 0 and dot(c,mv(G,point)).sign() < 0
        assert values[n] == 0 and all(v.sign() < 0 for m,v in values.items() if m != n)
    u = [Q(),Q(),Q(),Q(),Q(1),Q(2),Q(3),Q(2),Q(1),Q()]
    J = [decode(v) for v in large['coordinates']]
    assert J == add([2*v for v in roots['j']],u,-1)
    assert square(u) == 6 and dot(roots['j'],mv(G,u)) == 3
    assert dot(J,cov['j']) == 1
    vertices,retention = {},{}

    def vertex(name,point):
        values,sq = evaluate(point),square(point)
        active = [n for n,v in values.items() if v == 0]
        assert all(v.sign() <= 0 for v in values.values()) and len(active) == 9
        assert sq.sign() < 0 and dot(p,mv(G,point)).sign() < 0
        assert inertia([[dot(roots[a],cov[b]) for b in active] for a in active]) == [9,0,0]
        assert values['r'].sign() < 0 and values['w'].sign() < 0
        vertices[name] = dict(point=encoded(point),active=active,square=sq.encoded(),j_product=values['j'].encoded())
        return active

    for name,record in prior['finite_vertex_witnesses'].items():
        point = [decode(v) for v in record['point']]
        v = dot(point,cov['j']); retained = v.sign() <= 0
        retention[name] = dict(retained=retained,j_product=v.encoded())
        if retained: vertex(name,point)
    V = [Q(1),-phi/2]+[phi*phi*phi*v/4 for v in ci[0]]
    direction = [Q(),Q()]+[b-F(7,4)*a for a,b in zip(ci[0],ci[1])]
    tau = decode(small['endpoint_parameter']); endpoint = add(V,direction,tau)
    assert tau == (-3+4*phi)/19 and (Q(F(1,5))-tau).sign() > 0
    assert endpoint == [decode(v) for v in small['endpoint']]
    assert square(endpoint) == decode(small['endpoint_square'])
    active = vertex('p_j',endpoint)
    # Explicit E7+2A1 simple graph: each arm is listed from the branch outward.
    assert set(active) == {'old0','z','j'} | {f'old{i}' for i in range(3,9)}
    edges = {frozenset((a,b)) for a,b in combinations(active,2) if dot(roots[a],cov[b]) != 0}
    expected_edges = [('old5','j'),('old5','old6'),('old6','old7'),
                      ('old5','old4'),('old4','old3'),('old3','old8')]
    assert edges == {frozenset(e) for e in expected_edges}
    assert all(dot(roots[a],cov[b]) == -1 for a,b in expected_edges)
    caps = {n:dict(product=newproducts[n][0],kind=newproducts[n][1],
        inherited_full_facet_preserved=newproducts[n][1]=='ultraparallel') for n in ('r','w')}
    assert all(record['inherited_full_facet_preserved'] for record in caps.values())
    assert len({tuple(tuple(v) for v in r['point']) for r in vertices.values()}) == len(vertices)
    # Recheck the nine original p-rays; do not infer cone completeness from them.
    initial,pv = [],evaluate(p)
    for i,direction in enumerate([[Q(),Q(1)]+[Q()]*8]+[[Q(),Q()]+row for row in ci]):
        slopes = evaluate(direction)
        assert slopes[f'old{i}'].sign() < 0
        assert all(slopes[f'old{k}'] == 0 for k in range(9) if k != i)
        hits = [(-pv[n]/v,n) for n,v in slopes.items() if pv[n].sign() < 0 and v.sign() > 0]
        assert hits
        first = hits[0][0]
        for value,_ in hits[1:]:
            if (value-first).sign() < 0: first = value
        point = add(p,direction,first); values,sq = evaluate(point),square(point)
        assert first.sign() > 0 and all(v.sign() <= 0 for v in values.values())
        assert sq.sign() < 0
        names = [n for n,v in values.items() if v == 0]
        assert len(names) == 9
        assert inertia([[dot(roots[a],cov[b]) for b in names] for a in names]) == [9,0,0]
        initial.append(dict(freed_wall=f'old{i}',parameter=first.encoded(),
            first_hit=[n for v,n in hits if v == first],square=sq.encoded()))
    walk = bounded_walk(G,roots,cov,p,64)
    return dict(status='PASS_EXACT_29_WALL_AND_BOUNDED_EDGE_DIAGNOSTIC',
        checker_sha256=digest(Path(__file__)),sources=PINS,walls=29,wall_pairs=406,
        wall_pairs_sha256=pair_hash(pairs),wall_pair_categories=dict(categories),
        new_root=encoded(roots['j']),new_covector=encoded(cov['j']),new_wall_products=newproducts,
        reflection_integrality_checked=True,actual_facet_projection_witnesses=29,
        strict_center=encoded(c),strict_center_products={n:v.encoded() for n,v in cp.items()},
        alternate_survivor_product=Q(1).encoded(),prior_named_vertex_retention=retention,
        finite_vertex_witnesses=vertices,new_endpoint_type='E7+2A1',old_cap_facets=caps,
        initial_p_edge_first_hits=initial,initial_p_finite_edges=9,bounded_edge_walk=walk,
        inherited_finite_vertex_lower_bound=82+len(vertices),global_finite_volume_proved=False,
        inherited_only='Two full disjoint41-vertex compact H8 facets. Preservation follows from '
        'new-wall ultraparallelism and the checked strict facet projection witness; no cap census rerun.',
        scope='Fixed29-wall algebra,29 actual facets,named finite vertices,nine initial edges and '
        'one64-vertex walk. No all-height root census,complete volume conclusion or dimension claim. '
        'Queue closure is diagnostic only. A spacelike escape needs its separate written volume argument.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true'); mode.add_argument('--verify',action='store_true')
    args = parser.parse_args(); signal.alarm(45)
    result = check(); payload = json.dumps(result,sort_keys=True,indent=2)+'\n'
    assert len(payload.encode()) < 1024**2
    if args.verify: assert json.loads(CERT.read_text()) == result, 'Saved certificate differs'
    else:
        with CERT.open('x') as stream: stream.write(payload)
    signal.alarm(0)
    print(json.dumps(dict(status=result['status'],walls=29,finite_named=len(result['finite_vertex_witnesses']),
        walk_stop=result['bounded_edge_walk']['stop'])))


if __name__ == '__main__':
    main()
