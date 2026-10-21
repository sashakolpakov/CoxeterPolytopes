#!/usr/bin/env python3
"""Independent primitive-element arithmetic and complete new-wall finite census.

The complete old 32-wall cone is inherited from its frozen independent review.
No producer arithmetic, proposal engine or chamber engine is imported.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sympy as s
import fresh_h9_h8_all_high_20261014_review as C
import fresh_h9_h8_facet_cover_20261014_line_review as L

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_joint_corner_20261015_bridge_review'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json'
OLD_REVIEW='results/fresh_h9_joint_pair_20261015_two_review.json'
PRIMARY='results/fresh_h9_joint_corner_20261015_bridge.json'
SQ5=s.sqrt(5);RAD=6+s.Rational(16,5)*SQ5
DELTA=s.Rational(1988340828,1990805)+s.Rational(4448872824,9954025)*SQ5
K=s.QQ.algebraic_field(SQ5,s.sqrt(RAD),s.sqrt(DELTA))
Z,O=K.zero,K.one
C.K,C.Z,C.O=K,Z,O
L.BITS=384;L.SCALE=1<<L.BITS
THETA=C.ev(K.ext.as_expr())
BASIS=[K.from_sympy(x)for x in[s.S.One,SQ5,s.sqrt(RAD),SQ5*s.sqrt(RAD)]]
EXTRA=K.from_sympy(s.sqrt(DELTA))

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
@lru_cache(None)
def decode4(t):return sum((K.convert(F(a))*b for a,b in zip(t,BASIS)),Z)
def d4(a):return decode4(tuple(a[0]+a[1]))
def d8(a):return d4(a[0])+EXTRA*d4(a[1])
def v4(a):return list(map(d4,a))
def m4(a):return list(map(v4,a))
def v8(a):return list(map(d8,a))
@lru_cache(None)
def sg(a):
    if not a:return 0
    out=L.I()
    for x in a.to_list():out=out*THETA+L.I(F(str(x)))
    assert out.l>0 or out.h<0,('unresolved sign',a)
    return 1 if out.l>0 else-1
def enc(a):
    if isinstance(a,(list,tuple)):return[enc(x)for x in a]
    return list(map(str,a.to_list()))
def canon(v):return tuple(x/next(a for a in v if a)for x in v)
def pd(G):
    B=[row[:]for row in G];piv=[]
    for k in range(len(B)):
        q=B[k][k];piv.append(q)
        if sg(q)<=0:return False,piv
        for i in range(k+1,len(B)):
            for j in range(k+1,len(B)):B[i][j]-=B[i][k]*B[k][j]/q
    return True,piv
def relation(g,q,r):
    if sg(g)>0:return'obtuse',None,None
    if not g:return'right',2,None
    u=g*g/(q*r)
    if sg(u-O)>=0:return('parallel'if u==O else'dotted'),None,u
    allowed={3:O/4,4:O/2,5:(3*O+BASIS[1])/8,6:3*O/4,10:(5*O+BASIS[1])/8}
    for label,value in allowed.items():
        if u==value:return'finite',label,u
    return'forbidden',None,u
def components(n,edges):
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in edges:
        a,b=find(i),find(j)
        if a!=b:parent[max(a,b)]=min(a,b)
    out={}
    for i in range(n):out.setdefault(find(i),[]).append(i)
    return sorted(out.values())


def run():
    paths=[BASE,OLD_REVIEW,PRIMARY];b,old,d=[json.loads((ROOT/p).read_text())for p in paths]
    for doc in(b,old,d):
        for p,h in doc['sources'].items():assert sha(p)==h,p
    assert old['status']=='PASS_INDEPENDENT_TWO_CAPS_AND_COMPLETE_32_WALL_CONE'
    assert old['sources'][BASE]==sha(BASE)
    assert old['counts']==b['counts']and old['finite_vertices']==90
    assert len(K.mod.to_list())==9 and sg(EXTRA)>0 and EXTRA*EXTRA==K.from_sympy(DELTA)
    assert sg(K.from_sympy(RAD))>0
    assert d4(d['field']['extra_positive_radical_square'])==K.from_sympy(DELTA)
    assert d['field']['absolute_degree']==8 and d['field']['relative_degree']==2
    assert d['field']['base']==b['field']
    Q,N,H=m4(b['form']),m4(b['normals']),v4(b['center'])
    assert len(N)==32 and[sg(x)for x in C.ldl(Q)].count(-1)==1
    pair=lambda a,bb:C.pair(a,bb,Q)
    cov=[C.mv(Q,n)for n in N]
    norm=[pair(n,n)for n in N]
    assert norm==v4(b['norms'])and all(sg(x)>0 for x in norm)
    ids=[0,1,2,5,6,10,29,30,8]
    kernel=C.null_one([cov[i]for i in ids]);assert kernel[9]
    kernel=[x/kernel[9]for x in kernel]
    phi=(O+BASIS[1])/2
    rhs=[-O if i==30 else-phi if i==8 else Z for i in ids]
    part=C.mv(C.inverse([cov[i]for i in ids]+[C.mv(Q,kernel)]),rhs+[Z])
    assert not pair(part,kernel)
    assert (2*O-pair(part,part))/pair(kernel,kernel)==K.from_sympy(DELTA)
    rec=d['cap']['reconstruction']
    assert rec['ids']==ids and part==v4(rec['particular'])and kernel==v4(rec['kernel'])
    assert d4(rec['delta'])==K.from_sympy(DELTA)
    accepted=[];branches=[]
    for sign in(-1,1):
        n=[a+sign*EXTRA*v for a,v in zip(part,kernel)]
        assert pair(n,n)==2*O
        relations=[relation(C.dot(n,row),2*O,q)for row,q in zip(cov,norm)]
        good=sg(pair(n,H))<0 and all(t not in('obtuse','forbidden')for t,_,_ in relations)
        branches.append(dict(sign=sign,center_sign=sg(pair(n,H)),accepted=good,
                             bad_walls=[i for i,(t,_,_)in enumerate(relations)if t in('obtuse','forbidden')]))
        if good:accepted.append((sign,n))
    assert len(accepted)==1
    sign,n=accepted[0];assert sign==d['cap']['branch']and n==v8(d['cap']['normal'])
    assert pair(n,n)==d8(d['cap']['norm'])and pair(n,H)==d8(d['cap']['center_pairing'])
    N.append(n);cov.append(C.mv(Q,n));norm.append(2*O)
    G=[[C.dot(a,row)for row in cov]for a in N]
    assert[row[:32]for row in G[:32]]==m4(b['gram'])
    counts=Counter();new_counts=Counter();neighbors=[]
    for i,j in combinations(range(33),2):
        typ,label,u=relation(G[i][j],norm[i],norm[j]);assert typ not in('obtuse','forbidden')
        counts[typ+(str(label)if label else'')]+=1
        if j==32:
            rr=d['cap']['relations'][i]
            assert rr['wall']==i and(typ,label)==(rr['kind'],rr['label'])
            assert u==(None if rr['invariant']is None else d8(rr['invariant']))
            new_counts[typ+(str(label)if label else'')]+=1
            if typ in('right','finite'):neighbors.append(i)
    assert G[32][:32]==v8(d['cap']['pairings'])and neighbors==sorted(ids)
    assert sg(pair(H,H))<0 and all(sg(C.dot(H,c))<0 for c in cov)
    for i,row in enumerate(d['facet_witnesses']):
        P=[h-Gi/norm[i]*a for h,a in zip(H,N[i])for Gi in[C.dot(H,cov[i])]]
        assert row['wall']==i and P==v8(row['point'])and pair(P,P)==d8(row['norm'])
        assert sg(pair(P,P))<0 and sg(pair(P,H))<0
        assert all(sg(C.dot(P,c))==(0 if j==i else-1)for j,c in enumerate(cov))
    assert len(d['facet_witnesses'])==33
    print('independent octic recipe,528 pairs,33 facets PASS',flush=True)
    # Reconstruct every old finite ray from its literal support; the frozen
    # independent complete cone establishes that this is the full old list.
    vertices=[];preservation=[]
    for index,row in enumerate(b['finite_vertices']):
        active=row['incident'];assert len(active)==9
        V=C.null_one([cov[i]for i in active]);V=[x/-pair(V,H)for x in V]
        saved=v4(row['point']);assert canon(V)==canon(saved)and sg(pair(saved,H))<0
        assert pd([[G[i][j]for j in active]for i in active])[0]
        assert sg(pair(V,V))<0 and pair(saved,saved)==d4(row['norm'])
        vals=[C.dot(V,c)for c in cov]
        assert all(sg(g)==(0 if i in active else-1)for i,g in enumerate(vals))
        value=pair(saved,n);assert value==d8(d['old_finite_vertex_preservation_pairings'][index])
        preservation.append(value);vertices.append((active,V))
    assert len(vertices)==90 and len(d['old_finite_vertex_preservation_pairings'])==90
    census=d['complete_new_cap_finite_census'];assert census['finite_neighbors']==neighbors
    tested=[];new=[]
    for supp in combinations(neighbors,8):
        active=list(supp)+[32];positive,piv=pd([[G[i][j]for j in active]for i in active])
        tested.append(dict(support=active,positive_definite=positive,pivots=enc(piv)))
        if not positive:continue
        V=C.null_one([cov[i]for i in active]);V=[x/-pair(V,H)for x in V]
        assert sg(pair(V,V))<0
        vals=[C.dot(V,c)for c in cov]
        assert all(sg(g)==(0 if i in active else-1)for i,g in enumerate(vals))
        rr=next(x for x in census['vertices']if x['incident']==active)
        assert V==v8(rr['point'])and pair(V,V)==d8(rr['norm'])and vals==v8(rr['pairings'])
        new.append(active);vertices.append((active,V))
    assert[{k:r[k]for k in('support','positive_definite')}for r in tested]==census['candidate_supports']
    assert len(tested)==9 and len(new)==2
    assert new==[[0,1,2,5,6,8,10,29,32],[0,1,2,5,6,10,29,30,32]]
    incidents=[ids for ids,_ in vertices];assert incidents==d['finite_vertex_incidence']
    edges=[]
    for i,j in combinations(range(len(vertices)),2):
        common=sorted(set(incidents[i])&set(incidents[j]))
        if len(common)==8:
            assert pd([[G[a][bb]for bb in common]for a in common])[0]
            edges.append(dict(vertices=[i,j],support=common))
    assert edges==d['finite_edges']
    old_edges=[e for e in edges if max(e['vertices'])<90]
    assert old_edges==b['finite_edges']
    comps=components(92,[r['vertices']for r in edges]);old_comps=components(90,[r['vertices']for r in old_edges])
    assert comps==d['finite_graph_components']and old_comps==d['old_finite_graph_components']
    assert len(comps)==7 and len(old_comps)==8
    joint=next(i for i,ids in enumerate(incidents[:90])if 30 in ids and 31 in ids)
    merged=next(x for x in comps if joint in x)
    assert merged==d['joint_corner_joined_component']==[19,20,21,90,91]
    new_edges=[x for x in edges if max(x['vertices'])>=90]
    assert len(new_edges)==3 and all(sum(i in e['vertices']for e in new_edges)==2 for i in(90,91))
    # Each full spherical corner is simple, so nine directions per vertex.
    germ_count=9*len(vertices)-2*len(edges)
    assert len(vertices)==d['finite_vertex_count']==92
    assert len(edges)==d['finite_edge_count']==253 and germ_count==d['open_finite_germ_count']==322
    E=v4(b['strict_escape']['point']);ee=d['retained_strict_escape']
    assert E==v8(ee['point'])and sg(pair(E,E))>0 and sg(pair(E,H))<0
    assert pair(E,E)==d8(ee['norm'])and[C.dot(E,c)for c in cov]==v8(ee['pairings'])
    assert all(sg(C.dot(E,c))<0 for c in cov)and pair(E,n)==d8(ee['cap_pairing'])
    paths += [f'scripts/{TAG}.py','scripts/fresh_h9_joint_corner_20261015_bridge.py',
              'scripts/fresh_h9_joint_pair_20261015_two_review.py',
              'scripts/fresh_h9_h8_all_high_20261014_review.py',
              'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return dict(status='PASS_INDEPENDENT_OCTIC_BRIDGE_AND_COMPLETE_FINITE_VERTEX_CENSUS',
        sources={p:sha(p)for p in paths},arithmetic='Independent SymPy degree-eight primitive-element field and384-bit outward Fraction sign bounds. No producer arithmetic imported.',
        minimal_polynomial=list(map(str,K.mod.to_list())),radical_square=str(DELTA),
        inherited_complete_old_cone=OLD_REVIEW,branch_checks=branches,new_pair_counts=dict(new_counts),
        all_pair_counts=dict(counts),pair_checks=528,actual_facets=33,facet_sign_checks=1089,
        old_finite_vertices_reconstructed=90,preservation_pairings=enc(preservation),
        complete_new_wall_supports=tested,new_finite_vertices=new,finite_vertices=92,
        finite_edges=253,new_edges=new_edges,open_finite_germs=germ_count,
        old_graph_components=old_comps,graph_components=comps,merged_component=merged,
        strict_escape_norm=enc(pair(E,E)),strict_escape_center_pairing=enc(pair(E,H)),
        strict_escape_pairings=enc([C.dot(E,c)for c in cov]),
        scope='Complete finite-vertex/finite-edge/germ census follows from the inherited independently complete32-wall cone, reconstruction and strict preservation of all90old finite vertices, and all9 possible new-wall finite supports. The two new vertices produce a finite chain joining two old components. No exterior-ray census repeated or claimed. Strict spacelike escape proves infinite volume.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertices'],result['finite_edges'],result['open_finite_germs'],flush=True)
