#!/usr/bin/env python3
"""Local arithmetic geometry and independent counts for the H3 model barrier.

No payment-model or barrier-producer import. Shares only a standalone exact
Fraction Q(sqrt5) field and Gaussian algebra from the earlier geometry witness.
"""
import argparse
import json
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import fresh_sqrt5_adjacent_discount_witness_20261016 as X

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_multiplicity_geometry_review_20261016'
PRIMARY = 'fresh_sqrt5_h3_multiplicity_barrier_20261016'
K, Z, O = X.K, X.Z, X.O


def pieces(ids, bonds):
    graph = {i: {} for i in ids}
    for a, b, m in bonds:
        if a in graph and b in graph:
            graph[a][b] = graph[b][a] = m
    unseen = set(ids); answer = []
    while unseen:
        todo = [min(unseen)]; nodes = set()
        while todo:
            u = todo.pop()
            if u in nodes: continue
            nodes.add(u); todo.extend(graph[u])
        unseen -= nodes
        answer.append({u: graph[u] for u in sorted(nodes)})
    return answer


def distances(graph, u):
    d = {u: 0}; queue = deque([u])
    while queue:
        a = queue.popleft()
        for b in graph[a]:
            if b not in d: d[b] = d[a]+1; queue.append(b)
    return d


def kind(g):
    n = len(g)
    labels = [m for u in g for v, m in g[u].items() if u < v]
    assert len(labels) == n-1
    if not labels: return 'A1'
    if max(map(len, g.values())) <= 2:
        if all(m == 3 for m in labels): return 'A'+str(n)
        assert sorted(labels) == [3]*(n-2)+[5] and n in (3,4)
        return 'H'+str(n)
    branch = [u for u in g if len(g[u]) == 3]
    assert len(branch) == 1 and all(m == 3 for m in labels)
    arms = sorted(distances(g, branch[0])[u] for u in g if len(g[u]) == 1)
    if arms[:2] == [1,1]: return 'D'+str(n)
    assert arms == [1,2,n-4] and n in (6,7,8)
    return 'E'+str(n)


def run():
    primary = json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert primary['status'] == 'PASS_EXACT_H3_MULTIPLICITY_MODEL_BARRIER'
    row = primary['face']['row']; n = row['nodes']; assert n == 10
    bonds = [(i,j,m) for (i,j),m in zip(combinations(range(n),2),row['labels']) if m != 2]
    assert bonds == [(0,6,5),(1,9,5),(2,3,3),(2,5,3),(3,7,3),(3,8,3),(4,5,3),(4,6,3),(8,9,3)]
    phi = K(F(1,2),F(1,2))
    G = [[K(2*(i == j)) for j in range(n)] for i in range(n)]
    for i,j,m in bonds: G[i][j] = G[j][i] = -{3:O,5:phi}[m]
    assert X.inertia(G) == [8,1,1]
    assert X.inertia([[x.conj() for x in r] for r in G]) == [9,0,1]
    ids = [i for i in range(n) if i != 1]; Q = X.sub(G,ids); Qi = X.inverse(Q)
    normals = [X.mv(Qi,[G[j][i] for j in ids]) for i in range(n)]
    assert [[X.pair(Q,a,b) for b in normals] for a in normals] == G
    H = X.mv(Qi,[-O]*9); hp = [X.pair(Q,H,v) for v in normals]; hq = X.pair(Q,H,H)
    assert hq == K(381,F(-421,2)) and hq.sign() < 0 and max(x.sign() for x in hp) < 0
    white = list(range(4,n)); assert row['white'] == white
    assert X.inertia(X.sub(G,white)) == [6,0,0]
    vertices = []; patterns = []
    for a,b in product((0,2),(1,3)):
        S = white+[a,b]; SG = X.sub(G,S); assert X.inertia(SG) == [8,0,0]
        c = X.mv(X.inverse(SG),[hp[i] for i in S])
        v = [H[j]-sum((t*normals[i][j] for t,i in zip(c,S)),Z) for j in range(9)]
        p = [X.pair(Q,v,r) for r in normals]
        assert all(not p[i] for i in S) and all(p[i].sign() < 0 for i in range(n) if i not in S)
        assert X.pair(Q,v,v).sign() < 0 and X.pair(Q,v,H).sign() < 0
        vertices.append(dict(black=[a,b],active=sorted(S),norm=X.pair(Q,v,v).enc(),point=[x.enc() for x in v],products=[x.enc() for x in p]))
        pcs = pieces(S,bonds); left = next(g for g in pcs if a in g); right = next(g for g in pcs if b in g)
        patterns.append(dict(black=[a,b],forest=sorted(kind(g) for g in pcs),
                             black_components=[kind(left),kind(right)],
                             baseline=int(b in left and distances(left,a)[b] <= 2),
                             adjacent_discount=int(b in left and left[a].get(b) == 3 and len(left) >= 6),
                             H3_count=sum(kind(g) == 'H3' for g in pcs)+4,
                             A1_count=sum(kind(g) == 'A1' for g in pcs)+1))
    assert [p['baseline'] for p in patterns] == [0,0,0,1]
    assert [p['adjacent_discount'] for p in patterns] == [0,0,0,1]
    paying = patterns[2]; assert paying['black'] == [2,1] and paying['H3_count'] == 5 and paying['A1_count'] == 2
    pcs = pieces(white+[2,1],bonds); h3 = next(g for g in pcs if 1 in g); a4 = next(g for g in pcs if 2 in g)
    assert kind(h3) == 'H3' and list(h3[1].values()) == [5]
    assert kind(a4) == 'A4' and len(a4[2]) == 1 and len(a4[4]) == 2
    cone = []; accepted = []
    for black in combinations(range(4),2):
        ker = X.kernel([X.mv(Q,normals[i]) for i in white+list(black)]); assert len(ker) == 1
        signs = []
        for s in (-1,1):
            v = [s*x for x in ker[0]]; p = [X.pair(Q,v,r) for r in normals]
            if max(x.sign() for x in p) <= 0:
                assert X.pair(Q,v,v).sign() < 0 and X.pair(Q,v,H).sign() < 0
                signs.append(s); accepted.append(sorted(black))
        cone.append(dict(black=list(black),feasible_signs=signs))
    assert sorted(accepted) == sorted(sorted(v['black']) for v in vertices)
    # The literal common spectator context A1+4H3 has rank13, not just a
    # formal count vector. It gives an admissible local chamber inH21.
    spectator = [[K(2*(i == j)) for j in range(13)] for i in range(13)]
    for i in (1,4,7,10):
        spectator[i][i+1] = spectator[i+1][i] = -phi
        spectator[i+1][i+2] = spectator[i+2][i+1] = -O
    assert X.inertia(spectator) == [13,0,0]
    assert X.inertia([[x.conj() for x in r] for r in spectator]) == [13,0,0]
    sh = [-x/K(1024) for x in X.mv(X.inverse(spectator),[O]*13)]
    assert (hq+X.pair(spectator,sh,sh)).sign() < 0
    big = [[Z]*23 for _ in range(23)]
    for i in range(10):
        for j in range(10): big[i][j] = G[i][j]
    for i in range(13):
        for j in range(13): big[10+i][10+j] = spectator[i][j]
    assert X.inertia(big) == [21,1,1] and X.inertia([[x.conj() for x in r] for r in big]) == [22,0,1]
    assert all((2*big[i][j]/big[j][j]).integral() for i in range(23) for j in range(23))
    # Independent actual root-pair recounts, without the payment compiler.
    d20 = {i:{} for i in range(20)}
    for a,b in [(0,2),(1,2)]+[(i,i+1) for i in range(2,19)]: d20[a][b] = d20[b][a] = 3
    assert kind(d20) == 'D20'
    count1 = sum(distances(d20,a)[b] == 1 for a,b in combinations(d20,2))
    count2 = sum(distances(d20,a)[b] == 2 for a,b in combinations(d20,2))
    assert count1 == count2 == 19
    a4path = {i:{} for i in range(4)}
    for i in range(3): a4path[i][i+1] = a4path[i+1][i] = 3
    a4base = sum(distances(a4path,a)[b] <= 2 for a,b in combinations(a4path,2))
    assert a4base == 5 and sum(len(g) == 1 for g in a4path.values()) == 2
    base_second = 3+4*a4base; beta_second = F(4*2*3,3*1)
    assert base_second == 23 and beta_second == 8
    beta_face = F(1,3*paying['H3_count']); assert beta_face == F(1,15)
    beta = primary['bonus_index']; face = {j:F(z) for j,z in primary['face']['profile']}
    assert face == {-1:F(1),1:F(-1),beta:beta_face}
    vd, vw = [{j:F(z) for j,z in rec['profile']} for rec in primary['vertices']]
    assert vd[-1] == 38 and vd[1] == -count1 and vd[0] == -count2 and vd.get(beta,0) == 0
    assert vw[-1] == base_second and vw[beta] == beta_second and vw.get(0,0) == vw.get(1,0) == 0
    assert all(z >= 0 for p in (vd,vw) for j,z in p.items() if j >= 2)
    lower = F(120,139)*38+F(19,139)*23
    assert F(120,139)*(-19)+F(19,139)*120 == 0
    assert lower == F(4997,139) > 35 and F(primary['lower_bound']) == lower
    paths = [f'scripts/{TAG}.py',f'scripts/{X.TAG}.py',f'results/{PRIMARY}.json',f'scripts/{PRIMARY}.py',
             'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py',
             'scripts/fresh_sqrt5_h3_center_e6_orbit_exempt_20261016.py']
    return dict(status='PASS_INDEPENDENT_H3_MULTIPLICITY_COUNTS_AND_LOCAL_GEOMETRY',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                field='Fraction Q(sqrt5), shared standalone geometry helper; no payment compiler import',
                geometry=dict(class_index=7661,gram=X.enc(G),kernel=[[x.enc() for x in v] for v in X.kernel(G)],
                              basis_indices=ids,center_norm=hq.enc(),vertices=vertices,complete_face_cone=cone,
                              spectator_gram=X.enc(spectator),full_inertias=[[21,1,1],[22,0,1]],integral_cartan_entries=529),
                corners=patterns,role=dict(H3_black=1,H3_role='5-edge leaf',source_black=2,source_type='A4 endpoint',facet_mark=4),
                independent_counts=dict(face_actual_H3=5,face_beta=str(beta_face),vertex_actual_H3=1,
                                        A4_endpoints=8,vertex_beta=str(beta_second),D20_adjacent=count1,D20_distance2=count2),
                lower_bound=str(lower),
                scope='Exact local arithmetic compact quadrilateral plus independent literal normalization/scalar proof. Sparse absence of other face bonuses is checked against the frozen full profile and audited directly in the source definitions, not through a second complete model compiler. No compact ambient realization or new dimension bound. Actual-count splitting removes this coupling only; successor feasibility and full coverage remain open.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(); ap.add_argument('--verify',action='store_true'); args = ap.parse_args()
    data = json.loads(json.dumps(run())); path = ROOT/f'results/{TAG}.json'
    if args.verify: assert data == json.loads(path.read_text()); print('PASS saved independent multiplicity review')
    else: path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['lower_bound'])
