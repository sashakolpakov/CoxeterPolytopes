#!/usr/bin/env python3
"""Two exact local arithmetic quadrilaterals force an E6 discount repair.
Physical calculations reuse only the new standalone Fraction field helper.
"""
import argparse,json
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as X
K,F,Z,O=X.K,X.F,X.Z,X.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_e6_discount_witness_20261016'


def certify(ci,row):
    n=row['nodes'];phi=K(F(1,2),F(1,2));G=[[K(2*(i==j))for j in range(n)]for i in range(n)]
    if ci==4266:G[0][0]=2*(phi+2)
    bonds=[]
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        if m!=2:
            G[i][j]=G[j][i]=-{3:O,5:phi,10:phi+2}[m];bonds.append((i,j,m))
    assert X.inertia(G)==[n-2,1,1]and X.inertia([[x.conj()for x in r]for r in G])==[n-1,0,1]
    omitted=0 if ci==4266 else 1;ids=[i for i in range(n)if i!=omitted];Q=X.sub(G,ids);Qi=X.inverse(Q)
    N=[]
    for i in range(n):
        N.append(X.mv(Qi,[G[j][i]for j in ids]))
    assert[[X.pair(Q,a,b)for b in N]for a in N]==G
    H=X.mv(Qi,[-O]*(n-1));hp=[X.pair(Q,H,v)for v in N];hq=X.pair(Q,H,H)
    assert hq.sign()<0 and max(x.sign()for x in hp)<0
    white=row['white'];assert X.inertia(X.sub(G,white))==[n-4,0,0]
    vertices=[]
    for a,b in product((0,2),(1,3)):
        S=white+[a,b];SG=X.sub(G,S);assert X.inertia(SG)==[n-2,0,0]
        c=X.mv(X.inverse(SG),[hp[i]for i in S])
        v=[H[j]-sum((d*N[i][j]for d,i in zip(c,S)),Z)for j in range(n-1)]
        p=[X.pair(Q,v,r)for r in N]
        assert all(not p[i]for i in S)and all(p[i].sign()<0 for i in range(n)if i not in S)
        assert X.pair(Q,v,v).sign()<0 and X.pair(Q,v,H).sign()<0
        vertices.append(dict(black=[a,b],active=sorted(S),coordinates=[x.enc()for x in v],norm=X.pair(Q,v,v).enc(),products=[x.enc()for x in p]))
    tests=[];accepted=[]
    for black in combinations(range(4),2):
        S=white+list(black);ker=X.kernel([X.mv(Q,N[i])for i in S]);assert len(ker)==1
        orientations=[]
        for sign in(-1,1):
            v=[sign*x for x in ker[0]];p=[X.pair(Q,v,r)for r in N]
            if max(x.sign()for x in p)<=0:
                assert X.pair(Q,v,v).sign()<0 and X.pair(Q,v,H).sign()<0
                orientations.append(sign);accepted.append(list(black))
        tests.append(dict(black=list(black),feasible_orientations=orientations))
    assert sorted(accepted)==sorted([sorted(r['black'])for r in vertices])and len(accepted)==4
    # A single large A component realizes zero reserved spectator counts.
    r=23-n;S=[[K(2*(i==j))for j in range(r)]for i in range(r)]
    for i in range(r-1):S[i][i+1]=S[i+1][i]=-O
    assert X.inertia(S)==[r,0,0]
    hs=[-x/K(1024)for x in X.mv(X.inverse(S),[O]*r)];fullnorm=hq+X.pair(S,hs,hs)
    assert fullnorm.sign()<0
    big=[[Z]*23 for _ in range(23)]
    for i in range(n):
        for j in range(n):big[i][j]=G[i][j]
    for i in range(r):
        for j in range(r):big[n+i][n+j]=S[i][j]
    assert X.inertia(big)==[21,1,1]and X.inertia([[x.conj()for x in t]for t in big])==[22,0,1]
    assert all((2*big[i][j]/big[j][j]).integral()for i in range(23)for j in range(23))
    import fresh_sqrt5_marked_h3_center_cross_20261016 as M
    edges=bonds+[(i,i+1,3)for i in range(n,22)];fullwhite=list(range(4,23));profiles=[]
    for mark in fullwhite:
        compiled=[M.compile_corner(fullwhite+list(p),edges,p,[mark])for p in product((0,2),(1,3))]
        p,tables=M.add_bonus(compiled,(0,)*7)
        for table in tables:p=M.plus(p,table[mark])
        assert p=={-1:F(1),1:F(-1)}
        profiles.append(dict(facet_root=mark,profile=[[-1,'1'],[1,'-1']]))
    return dict(class_index=ci,row=row,scaled_gram=X.enc(G),basis_indices=ids,basis_gram=X.enc(Q),
                root_coordinates=X.enc(N),kernel=[[x.enc()for x in z]for z in X.kernel(G)],
                core_inertias=dict(physical=[n-2,1,1],conjugate=[n-1,0,1]),
                strict_center=[x.enc()for x in H],center_norm=hq.enc(),center_products=[x.enc()for x in hp],
                white=white,vertices=vertices,complete_face_cone=tests,
                spectator_forest=['A'+str(r)],spectator_gram=X.enc(S),full_center_norm=fullnorm.enc(),
                full_inertias=dict(physical=[21,1,1],conjugate=[22,0,1]),
                integral_cartan_entries=529,all_facet_mark_profiles=profiles,
                necessary_payment_inequality='1-delta_adjacent3 >= 1, hence delta_adjacent3=0 under nonnegative discount')


def run():
    data=json.loads((ROOT/'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json').read_text())
    rows={r['class_index']:r['row']for r in data['records']}
    records=[certify(ci,rows[ci])for ci in(4266,7746)]
    paths=[f'scripts/{TAG}.py',f'scripts/{X.TAG}.py',f'results/{X.TAG}.json',
           'results/fresh_sqrt5_adjacent_discount_frontier_20261016.json',
           'scripts/fresh_sqrt5_marked_h3_center_cross_20261016.py']
    return dict(status='PASS_EXACT_E6_DISCOUNT_LOCAL_OBSTRUCTION',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
                field='a+b*sqrt(5)',records=records,
                scope='Two actual compact quadrilateral faces of local admissible integral-reflection chambers. The current compensation parametrization identically pays1-delta for every white facet mark in the displayed rank21 spectator context. This forces delta0 if that local face domain is to be covered. It is not a global nonexistence theorem, nor a compact ambient realization claim; additional arithmetic exclusions, new bonus channels or orbit-specific discounts remain possible.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text());print('PASS saved E6 witness replay')
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],[r['class_index']for r in d['records']])
