#!/usr/bin/env python3
"""Literal Bugaenko H8 facet and redundancy of every other current wall."""
import argparse,json,sys
from hashlib import sha256
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from coxeter_polytopes.exact.golden import Golden as K
from coxeter_polytopes.verify.chamber import diagnose_chamber_cone,_pair,_full_rank_kernel,_mat_vec
TAG='fresh_h9_beta131_h100_20261016_facet'
PREFIX='results/fresh_h9_beta131_h100_20261016.json'
GEOMETRY='results/fresh_h9_beta131_h100_20261016_geometry.json'
KNOWN='data/known_examples/bugaenko_compact_h8.json'
IDS=[11,0,2,3,4,5,6,7,9,1,40]
def vec(v):return tuple(map(K.from_json,v))
def mat(m):return tuple(map(vec,m))
def enc(x):
    if isinstance(x,K):return x.to_json()
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d,g,s=[json.loads((ROOT/p).read_text())for p in(PREFIX,GEOMETRY,KNOWN)]
    for p in(d,g):
        for q,h in p['sources'].items():assert sha(q)==h,q
    G,N,H=mat(d['form']),mat(d['roots']),vec(d['basepoint'])
    assert len(N)==41 and N[8]==(K(1),)+(K(),)*8+(K(0,-2),)
    delta=K(1,12);beta=K(12,13);r=N[8];qr=_pair(r,G,r)
    assert qr==2*beta and tuple(_mat_vec(G,r))==(K(),)*9+(-delta,)
    Q=tuple(row[:9]for row in G[:9]);U=tuple(N[i][:9]for i in IDS)
    assert Q==mat(s['lattice_gram']) and U==mat(s['simple_roots'])
    assert all(N[i][9]==0 and _pair(N[i],G,r)==0 for i in IDS)
    assert [i for i,n in enumerate(N)if _pair(n,G,r)==0]==sorted(IDS)
    gram=tuple(tuple(_pair(u,Q,v)for v in U)for u in U)
    assert gram==mat(s['wall_gram'])
    C=vec(s['geometry']['interior_point']);C10=(*C,K())
    assert _pair(C,Q,C)<0 and _pair(C10,G,H)<0
    assert all(_pair(C10,G,n).sign()==(0 if i==8 else -1)for i,n in enumerate(N))
    redundant=[]
    for i,n in enumerate(N):
        if i==8 or i in IDS:continue
        p=_pair(n,G,r);projection=tuple(x-p/qr*y for x,y in zip(n,r))
        q=_pair(projection,G,projection);c=_pair(projection,G,C10)
        assert p<0 and p*p>_pair(n,G,n)*qr and q<0 and c<0
        assert _pair(projection,G,r)==0
        redundant.append(dict(wall=i,pairing_with_facet=p,projection=projection,
                              projected_norm=q,future_center_pairing=c))
    assert len(redundant)==29
    cone=diagnose_chamber_cone(U,Q,C,dimension=8)
    assert len(cone.finite_vertices)==41 and not cone.cusps and not cone.hyperideal_rays and not cone.opposite_sheet_rays
    assert {tuple(x)for x in s['geometry']['vertices']}==set(cone.finite_vertices)
    vertices=[]
    for support in cone.finite_vertices:
        V=_full_rank_kernel([_mat_vec(Q,U[i])for i in support]);assert V is not None
        if _pair(V,Q,C)>0:V=tuple(-x for x in V)
        assert _pair(V,Q,V)<0 and _pair(V,Q,C)<0
        full=(*V,K());active=[i for i,n in enumerate(N)if _pair(full,G,n)==0]
        assert active==sorted([8]+[IDS[i]for i in support])
        assert all(_pair(full,G,n)<=0 for n in N)
        vertices.append(dict(source_support=list(support),full_incident=active,point=full,norm=_pair(full,G,full)))
    assert {tuple(v['full_incident'])for v in vertices}=={tuple(v['incident'])for v in g['finite_vertices']if 8 in v['incident']}
    paths=[f'scripts/{TAG}.py',PREFIX,GEOMETRY,KNOWN,'src/coxeter_polytopes/exact/golden.py',
           'src/coxeter_polytopes/exact/matrix.py','src/coxeter_polytopes/verify/chamber.py']
    return enc(dict(status='PASS_LITERAL_COMPACT_BUGAENKO_H8_FACET_WITH_ALL_OTHER_WALLS_REDUNDANT',
        sources={p:sha(p)for p in paths},facet_wall=8,facet_normal=r,facet_norm=qr,
        source_to_current_wall_map=IDS,restricted_form=Q,source_roots=U,source_gram=gram,
        strict_facet_center=C10,other_wall_redundancy=redundant,
        complete_vertices=vertices,counts=dict(facet_dimension=8,side_facets=11,finite_vertices=41,
            cusps=0,hyperideal=0,other_walls_redundant=29),
        scope='Exact literal identification and complete compact cone of one actual facet of the current41wall prefix. All29 other nonorthogonal walls have strictly future timelike projections and are redundant throughout this facet. This does not prove ambient finite volume.'))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved H8 facet replay')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'])
