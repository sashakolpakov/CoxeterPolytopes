#!/usr/bin/env python3
"""Independent complete marked transitions, free-pair tests and facet incidence."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from math import comb, isqrt

import fresh_sqrt5_marked_flow_20261016_review as M

I, C, K, ROOT = M.I, M.C, M.K, M.ROOT
TAG = 'fresh_sqrt5_marked_short_faces_20261016_review'
PRIMARY = 'results/fresh_sqrt5_marked_short_faces_20261016.json'
TRANSITIONS = M.FLOW


def decode(v):
    return K(*map(F,v))


def code(g, facet=None, port=None):
    remaining, out = set(g), []
    while remaining:
        part = set(M.distances(g,min(remaining)))
        remaining -= part
        if sum(len(g[i]) for i in part) != 2*len(part)-2:
            return None
        out.append(min(M.rooted_tree(g,i,None,facet,port) for i in part))
    return tuple(sorted(out))


def rational_square(q):
    q = F(q)
    return q >= 0 and isqrt(q.numerator)**2 == q.numerator and isqrt(q.denominator)**2 == q.denominator


def run():
    saved, trans = map(I.read,(PRIMARY,TRANSITIONS))
    assert saved['status'] == 'PASS_EXACT_ALL_PARAMETER_MARKED_SHORT_FACES_AND_FREQUENCY_BOUND'
    models = {name:g for name,g in C.models() if len(g)<=21}
    first = {0:{}} | {i+1:{j+1:m for j,m in row.items()} for i,row in models['D20'].items()}
    cases = [('A1+D20',first,0),('D21',models['D21'],19)]
    target_codes = {code(g,e):name for name,g,e in cases}
    finite_codes = {code(g)[0]:name for name,g in models.items()}
    assert len(finite_codes) == 68
    # If the normalized finite cycle lies in the quadratic field, phi(m)<=4.
    # Every prime dividing m therefore has p-1 dividing4; only2,3,5 occur.
    def totient(m):
        out = m
        for p in (2,3,5):
            if m%p == 0:
                out = out//p*(p-1)
        return out
    labels = sorted({2**a*3**b*5**c for a in range(4) for b in range(2) for c in range(2)
                     if 2**a*3**b*5**c >= 2 and totient(2**a*3**b*5**c)<=4})
    assert labels == [2,3,4,5,6,8,10,12]
    assert all(not rational_square(q) and not rational_square(F(q,5)) for q in (2,3))
    assert (2+K(0,1)).norm() == 5 and not rational_square(5)
    # Labels8/12 require sqrt2/sqrt3 inK; labels4/6/10 have x^2 inK
    # but x notinK, contradicting the corank equation when c!=0.
    finite_products = [(2,K()),(3,K(-1)),(5,-K(0,1))]
    aggregate = Counter()
    outputs = []
    for name,g,facet in cases:
        current = next(row for row in saved['type_results'] if row['source_marked_type']==name)
        A = [[K(2 if i==j else -int(j in g[i])) for j in range(21)] for i in range(21)]
        assert len(I.pivots(A)) == 21
        AJ,_ = I.inverse(A)
        choices = {}
        inputs = {row['port']:row for row in trans['port_records'] if row['source_marked_type']==name}
        saved_choices = {row['port']:row['choices'] for row in current['choices_by_port']}
        attachment_count = 0
        for port in sorted(set(g)-{facet}):
            common = {i:{j:m for j,m in row.items() if j!=port} for i,row in g.items() if i!=port}
            nodes = sorted(common)
            S = [[A[i][j] for j in nodes] for i in nodes]
            SJ,_ = I.inverse(S)
            u = [A[port][i] for i in nodes]
            residual = 2-I.dot(u,I.mv(SJ,u))
            exhaustive = {}
            for degree in range(4):
                for attachments in combinations(nodes,degree):
                    aggregate['attachment_subsets'] += 1
                    newg = {i:dict(row) for i,row in common.items()}
                    newg[21] = {}
                    for j in attachments:
                        newg[j][21] = newg[21][j] = 3
                    target = target_codes.get(code(newg,facet))
                    if target is None:
                        continue
                    attachment_count += 1
                    v = [K(-int(i in attachments)) for i in nodes]
                    nr = 2-I.dot(v,I.mv(SJ,v))
                    cross = I.dot(u,I.mv(SJ,v))
                    assert residual>0 and nr>0 and cross>0
                    assert cross*cross >= residual*nr and 2-cross>=0 and (2-cross)*(2-cross)>=residual*nr
                    # A unit-bond path from the old framing root to a new
                    # attachment gives a cyclic product +-g and hence integrality.
                    assert set(M.distances(g,port)) & set(attachments)
                    scalars = []
                    for aa in range(-5,2):
                        for bb in range(-7,0):
                            value = K(aa,bb)
                            h = value-cross
                            if -14<value<=0 and residual*nr-h.conjugate()*h.conjugate()>0:
                                if value < -2 or value == -K(0,1):
                                    assert residual*nr-h*h<0
                                    scalars.append(value)
                                else:
                                    assert -1<value<0
                    exhaustive[attachments] = (target,sorted(scalars))
            original = inputs[port]['complete_marked_reattachments']
            assert set(exhaustive) == {tuple(row['attachments']) for row in original}
            for row in original:
                assert exhaustive[tuple(row['attachments'])] == (row['target_marked_type'],sorted(decode(x['g']) for x in row['all_narrow_admissible_products']))
            rows = []
            for row in saved_choices[port]:
                attachments, scalar = row['attachments'],decode(row['g'])
                assert scalar in exhaustive[tuple(attachments)][1]
                b = [scalar if i==port else K(-int(i in attachments)) for i in range(21)]
                jb = I.mv(AJ,b)
                nr = 2-I.dot(b,jb)
                assert nr<0 and nr.conjugate()>0 and nr == decode(row['residual'])
                rows.append(dict(attachments=attachments,g=scalar,b=b,jb=jb,residual=nr))
            assert len(rows) == sum(len(scalars) for target,scalars in exhaustive.values())
            assert len({(tuple(r['attachments']),r['g']) for r in rows}) == len(rows)
            choices[port] = rows
        assert attachment_count == (43 if name=='A1+D20' else 23)
        pair_map = {tuple(row['ports']):row['edge_choice_pairs'] for row in current['port_pair_records']}
        assert set(pair_map) == set(combinations(sorted(choices),2))
        quads = []
        for i,j in combinations(sorted(choices),2):
            indexed = {tuple(row['choice_indices']):row for row in pair_map[i,j]}
            assert set(indexed) == {(ri,si) for ri in range(len(choices[i])) for si in range(len(choices[j]))}
            for ri,r in enumerate(choices[i]):
                for si,s in enumerate(choices[j]):
                    aggregate['edge_choice_pairs'] += 1
                    row = indexed[ri,si]
                    a,d,c = r['residual'],s['residual'],I.dot(r['b'],s['jb'])
                    assert c>0 and c==decode(row['cross_projection'])
                    determinant = lambda x:a*d-(x-c)*(x-c)
                    assert determinant(-2)==decode(row['determinant_at_minus2']) and determinant(0)==decode(row['determinant_at0'])
                    if determinant(-2)>0 and determinant(0)>0:
                        aggregate['continuous_exclusions'] += 1
                        assert row['status']=='NO_REAL_FINITE_PAIRING_BY_STRICT_CONCAVITY'
                        continue
                    assert row['status']=='EXACT_COMPLETE_FINITE_LABEL_TEST'
                    for (m,x),test in zip(finite_products,row['finite_label_tests']):
                        aggregate['angular_tests'] += 1
                        assert m==test['label'] and determinant(x)==decode(test['determinant'])
                        if determinant(x):
                            continue
                        aggregate['angular_rank_roots'] += 1
                        witness = test['actual_quadrilateral']
                        B = [ar+[r['b'][k]] for k,ar in enumerate(A)]+[r['b']+[K(2)]]
                        ww = [decode(z) for z in witness['extra_root_coefficients']]
                        extra = s['b']+[x]
                        assert I.mv(B,ww)==extra and I.dot(ww,extra)==2
                        rays = [[decode(z) for z in v] for v in witness['cyclic_vertex_rays']]
                        expected = [[21,22],[i,22],[i,j],[j,21]]
                        assert len(rays)==4
                        fourth = None
                        for v,negative in zip(rays,expected):
                            values = I.mv(B,v)+[I.dot(v,extra)]
                            assert all(z<=0 for z in values) and [k for k,z in enumerate(values) if z<0]==negative
                            assert I.dot(v,I.mv(B,v))<0
                            if negative==[i,j]:
                                ids = [k for k in range(23) if k not in negative]
                                full = [ar+[extra[k]] for k,ar in enumerate(B)]+[extra+[K(2)]]
                                fourth = {k:{} for k in ids}
                                for k,l in combinations(ids,2):
                                    assert full[k][l] in (K(),K(-1))
                                    if full[k][l]:fourth[k][l]=fourth[l][k]=3
                        assert all(I.dot(v,I.mv(B,w))<0 for v,w in combinations(rays,2))
                        assert code(fourth,facet) in target_codes
                        quads.append(dict(ports=[i,j],label=m,fourth_type=target_codes[code(fourth,facet)]))
        triangles = []
        for i, options in choices.items():
            for ri,r in enumerate(options):
                if r['g'] != -K(0,1):
                    assert r['g'] < -2
                    continue
                for j in sorted(set(g)-{facet,i}):
                    newg = {k:{l:m for l,m in row.items() if l!=j} for k,row in g.items() if k!=j}
                    newg[21] = {}
                    for k in r['attachments']:
                        if k!=j:newg[k][21]=newg[21][k]=3
                    newg[i][21]=newg[21][i]=5
                    cc = code(newg)
                    if cc is None or any(c not in finite_codes for c in cc):
                        continue
                    types = sorted(finite_codes[c] for c in cc)
                    found = next(t for t in current['triangles_with_two_maximal_vertices'] if t['max_to_max_port']==i and t['choice_index']==ri and t['ports']==sorted([i,j]))
                    assert types==found['third_vertex_type']
                    cost = F()
                    for k,l in combinations(newg,2):
                        if facet in (k,l):continue
                        distance = M.distances(newg,k).get(l)
                        if distance is not None:cost += int(distance<=2)
                    assert cost==35 and all(t!='H4' for t in types)
                    B = [ar+[r['b'][k]] for k,ar in enumerate(A)]+[r['b']+[K(2)]]
                    rays = [[decode(z) for z in v] for v in found['actual_triangle']['future_vertex_rays']]
                    assert len(rays)==3
                    for v,omitted in zip(rays,(21,i,j)):
                        assert I.mv(B,v)==[-K(int(k==omitted)) for k in range(22)]
                        assert I.dot(v,I.mv(B,v))<0
                    assert all(I.dot(v,I.mv(B,w))<0 for v,w in combinations(rays,2))
                    triangles.append(dict(ports=sorted([i,j]),third_types=types,third_baseline='35'))
        assert len(triangles)==len(current['triangles_with_two_maximal_vertices'])
        designated = [18,19] if name=='A1+D20' else [17,18]
        assert all(row['ports']==designated for row in triangles+quads)
        assert len(quads)==(1 if name=='A1+D20' else 0) and len(triangles)==(2 if name=='A1+D20' else 1)
        outputs.append(dict(type=name,quadrilaterals=quads,triangles=triangles,designated_corner=designated))
        print('independent full transition and short-face census',name,flush=True)
    assert aggregate==dict(attachment_subsets=54040,edge_choice_pairs=3028,continuous_exclusions=3010,angular_tests=54,angular_rank_roots=1)
    # Exhaust the finite list of possible short-face H/L counts for the
    # charge inequality, after separately budgeting T2 andQ4 at H corners.
    assert all(2 <= 3-h for h in (0,1))
    assert all(1 <= 4-h for h in (0,1,2))
    assert comb(20,2)==190 and 5*45-190==35
    assert 35-1==34 and 190-35==155
    assert F(37)+F(155,189)==F(7148,189)
    paths = [f'scripts/{TAG}.py',PRIMARY,TRANSITIONS,
             'scripts/fresh_sqrt5_marked_flow_20261016_review.py',
             'results/fresh_sqrt5_marked_flow_20261016_review.json',
             'scripts/fresh_h9_beta131_h300_symmetry_20261016_review.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
             'src/coxeter_polytopes/exact/golden.py']
    return dict(status='PASS_INDEPENDENT_COMPLETE_MARKED_SHORT_FACES_AND_INCIDENCE_BOUND',
                sources={p:I.sha(p) for p in paths},counts=dict(aggregate),type_results=outputs,
                selected_facet_inequality='34H<155L',baseline_average_strict_upper_bound='7148/189',
                scope='Independent complete54040attachment and enlarged scalar-box reconstruction, all3028free-pair tests, finite-angle field gate, unique quadrilateral and three directed triangles with actual future rays; designated-corner counting on the selected narrow facet. Uses frozen independently reviewed finite diagram catalogue and quadratic arithmetic. No new dimension ceiling or finite-volume realization.')


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'])
