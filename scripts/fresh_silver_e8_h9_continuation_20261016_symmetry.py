#!/usr/bin/env python3
"""Exact infinite-order stabilizer of the full silver E8 reflection chamber."""
import argparse,json
from hashlib import sha256
from pathlib import Path
from fractions import Fraction as F
import fresh_silver_e8_h9_20261015 as A

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_silver_e8_h9_continuation_20261016_symmetry'
PREFIX='results/fresh_silver_e8_h9_continuation_20261016.json'
CORNERS='results/fresh_silver_e8_h9_continuation_20261016_h100_corners.json'
S,H=A.S,A.H
Z,O=S(0),S(1)
SOURCE=list(range(9))
TARGETS=([4,3,2,1,0,17,11,7,8],[0,1,2,41,6,5,4,7,8])
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,S):return x.to_json()
    if isinstance(x,(tuple,list)):return [enc(v)for v in x]
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    return x
def trans(m):return tuple(zip(*m))
def mm(a,b):return tuple(tuple(H.dot(r,c)for c in zip(*b))for r in a)
def pd(m):
    m=[list(r)for r in m];p=[]
    for k in range(len(m)):
        q=m[k][k];assert q.sign()>0;p.append(q)
        for i in range(k+1,len(m)):
            for j in range(k+1,len(m)):m[i][j]-=m[i][k]*m[k][j]/q
    return p

def run():
    d=json.loads((ROOT/PREFIX).read_text());census=json.loads((ROOT/CORNERS).read_text())
    for row in(d,census):
        for p,h in row['sources'].items():assert sha(p)==h,p
    assert d['status']=='PASS_COMPLETE_BOUNDED_FULL_NORM_VINBERG_PREFIX'
    assert d['height']==S(100).to_json() and d['time_controls']==29 and d['wall_count']==48
    G=tuple(tuple(map(S.from_json,r))for r in d['gram'])
    N=tuple(tuple(map(S.from_json,r))for r in d['all_walls'])
    C=tuple(map(S.from_json,d['strict_interior']));base=(O,)+(Z,)*9
    assert G==A.GRAM and A.silver_determinant(G)==-A.EPS
    assert A.silver_inertia(G).positive==9 and A.silver_inertia(G).negative==1
    conjugate_pivots=pd([[x.conjugate()for x in r]for r in G])
    cov=[H.mv(G,r)for r in N];gram=[[H.dot(r,v)for v in cov]for r in N]
    assert H.pair(C,G,C).sign()<0 and H.pair(C,G,base).sign()<0
    assert all(H.dot(C,v).sign()<0 for v in cov)
    assert all(gram[i][j].sign()<=0 for i in range(48)for j in range(i))
    assert all((2*x/gram[i][i]).is_integral for i in range(48)for x in cov[i])
    assert [gram[i][i]for i in range(9)]==[S(2)]*8+[O]
    facets=[]
    for i,n in enumerate(N):
        P=tuple(x-H.dot(C,cov[i])/gram[i][i]*y for x,y in zip(C,n))
        assert H.pair(P,G,P).sign()<0 and H.pair(P,G,base).sign()<0
        assert all(H.dot(P,v).sign()==(0 if i==j else -1)for j,v in enumerate(cov))
        facets.append(P)
    SG=tuple(tuple(gram[i][j]for j in SOURCE)for i in SOURCE)
    pivots=pd(SG);J=H.inverse(SG);assert all(x.sign()>=0 for r in J for x in r)
    F0=trans([N[i]for i in SOURCE]+[base]);F0i=H.inverse(F0)
    corner_rows=[];maps=[]
    for ids in (SOURCE,*TARGETS):
        assert tuple(tuple(gram[i][j]for j in ids)for i in ids)==SG
        spatial=H.mv(H.inverse(tuple(tuple(cov[i][j]for j in range(1,10))for i in ids)),tuple(-cov[i][0]for i in ids))
        V=(O,*spatial);q=H.pair(V,G,V);vals=[H.dot(V,v)for v in cov]
        assert q.sign()<0 and H.pair(V,G,base).sign()<0
        assert [i for i,x in enumerate(vals)if not x]==sorted(ids) and all(x.sign()<=0 for x in vals)
        coefficient=H.pair(base,G,V)/q;b=H.mv(J,tuple(H.dot(base,cov[i])for i in ids))
        assert coefficient.sign()>0 and all(x.sign()<=0 for x in b)
        assert base==tuple(coefficient*V[k]+sum((b[j]*N[ids[j]][k]for j in range(9)),Z)for k in range(10))
        corner_rows.append(dict(ordered_indices=ids,point=V,norm=q,wall_pairings=vals,
            basepoint_decomposition=dict(vertex_coefficient=coefficient,normal_coefficients=b)))
        if ids==SOURCE:continue
        D=trans([N[i]for i in ids]+[V]);scale=A.silver_determinant(F0)/A.silver_determinant(D)
        if scale.sign()<0:scale=-scale
        assert H.pair(base,G,base)==scale*scale*q
        T=mm(trans([N[i]for i in ids]+[tuple(scale*x for x in V)]),F0i);TI=H.inverse(T)
        assert all(x.is_integral for r in(*T,*TI)for x in r)
        assert mm(mm(trans(T),G),T)==G
        assert all(H.mv(T,N[i])==N[j]for i,j in zip(SOURCE,ids))
        assert H.mv(T,base)==tuple(scale*x for x in V)
        maps.append(dict(source_indices=SOURCE,target_indices=ids,positive_corner_scale=scale,
                         matrix=T,inverse=TI,determinant=A.silver_determinant(T)))
    assert [m['positive_corner_scale']for m in maps]==[S(11,8),S(29,20)]
    U=mm(maps[1]['matrix'],maps[0]['matrix']);UI=H.inverse(U);T=mm(U,U);TI=H.inverse(T)
    assert all(x.is_integral for r in(*U,*UI,*T,*TI)for x in r)
    assert mm(mm(trans(T),G),T)==G and A.silver_determinant(T)==1
    assert H.pair(H.mv(T,base),G,base).sign()<0
    traceU=sum((U[i][i]for i in range(10)),Z);traceT=sum((T[i][i]for i in range(10)),Z)
    assert traceU==S(4,1) and traceT==S(10,4) and (traceT-10).sign()>0
    oldE=tuple(map(S.from_json,d['geometry']['hyperideal_witness']['vector']))
    eps=S(F(1,2))
    while True:
        E=tuple(x+eps*y for x,y in zip(oldE,C));q=H.pair(E,G,E)
        if q.sign()>0 and H.pair(E,G,base).sign()<0 and all(H.dot(E,v).sign()<0 for v in cov):break
        eps/=2
    paths=[f'scripts/{TAG}.py',PREFIX,CORNERS,'scripts/fresh_silver_e8_h9_20261015.py',
           'scripts/fresh_silver_f4_h9_20261013.py','src/coxeter_polytopes/exact/silver.py']
    return enc(dict(status='PASS_EXACT_INFINITE_ORDER_SILVER_E8_FULL_CHAMBER_STABILIZER',
        sources={p:sha(p)for p in paths},host='<-1-sqrt2>+E8+<1>',field='Q(sqrt2)',dimension=9,
        determinant=-A.EPS,conjugate_positive_pivots=conjugate_pivots,
        complete_prefix=dict(height=100,controls=29,walls=48,norm_classes=A.NORMS),
        all_pair_gram=gram,strict_center=C,actual_facet_points=facets,
        spherical_gram=SG,spherical_pivots=pivots,inverse_spherical_gram=J,corners=corner_rows,
        generators=maps,word=[1,0,1,0],matrix=T,inverse=TI,matrix_determinant=S(1),
        product_matrix=U,product_trace=traceU,trace=traceT,trace_minus_ten=traceT-10,
        strict_escape=dict(point=E,norm=q,epsilon=eps,wall_pairings=[H.dot(E,v)for v in cov]),
        theorem='The literal silver E8 lattice is nonreflective: its full integral reflection chamber has infinite covolume. No finite-volume Coxeter polytope can have all its wall reflections preserving this lattice.',
        proof_scope='Exact full-norm bounded prefix, primitive norm alphabet and zero sector are explicit inherited proof inputs. The displayed future integral maps carry complete persistent outward corner sectors to each other; their product squared stabilizes the full chamber and has trace10+4sqrt2>10. No exhaustive stabilizer search is claimed. Other lattices or nonintegral wall reflections are not excluded. Independent review is a separate requirement.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved full-chamber symmetry replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
