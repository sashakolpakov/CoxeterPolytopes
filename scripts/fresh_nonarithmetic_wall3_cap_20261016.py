#!/usr/bin/env python3
"""An exact escape-cutting cap and complete finite-corner update on 33 old walls."""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import fresh_h9_joint_corner_20261015_bridge_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_nonarithmetic_wall3_cap_20261016'
BASE='results/fresh_h9_joint_pair_20261015_two_extension.json'
BRIDGE='results/fresh_h9_joint_corner_20261015_bridge.json'
OLD_REVIEW='results/fresh_h9_joint_corner_20261015_bridge_review.json'
TARGET='results/fresh_nonarithmetic_escape_target_20261016.json'
C,K,Z,O=R.C,R.K,R.Z,R.O
D=O+R.BASIS[1]
T=[0,1,2,5,6,10,30]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d['sources'].items():assert sha(f)==h,f
    return d
def coerce(a):
    return a if hasattr(a,'to_list') else K.convert(a)
class X:
    """Formal a+b*sqrt(D), using only base-scalar division."""
    __slots__=('a','b')
    def __init__(self,a=Z,b=Z):
        if isinstance(a,X):self.a,self.b=a.a,a.b
        else:self.a,self.b=coerce(a),coerce(b)
    def __add__(self,o):
        o=X(o);return X(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self):return X(-self.a,-self.b)
    def __sub__(self,o):return self+-X(o)
    def __rsub__(self,o):return X(o)+-self
    def __mul__(self,o):
        o=X(o);return X(self.a*o.a+D*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=X(o);assert not o.b and o.a;return X(self.a/o.a,self.b/o.a)
    def __eq__(self,o):
        o=X(o);return self.a==o.a and self.b==o.b
    def sign(self):
        a,b=R.sg(self.a),R.sg(self.b)
        if not a:return b
        if not b or a==b:return a
        return a*R.sg(self.a*self.a-D*self.b*self.b)
def enc(x):
    if isinstance(x,X):return [R.enc(x.a),R.enc(x.b)]
    if isinstance(x,dict):return {k:enc(v) for k,v in x.items()}
    if isinstance(x,(tuple,list)):return [enc(v) for v in x]
    if hasattr(x,'to_list'):return R.enc(x)
    return x
def dot(x,y):return sum((X(a)*X(b) for a,b in zip(x,y)),X())
def mv(M,x):return [dot(row,x) for row in M]
def pair(x,y,Q):return dot(x,mv(Q,y))
def components(n,edges):
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    for row in edges:
        a,b=map(find,row['vertices'])
        if a!=b:parent[max(a,b)]=min(a,b)
    out={}
    for i in range(n):out.setdefault(find(i),[]).append(i)
    return sorted(out.values())

def run():
    base,bridge,review,target=map(read,(BASE,BRIDGE,OLD_REVIEW,TARGET))
    assert review['status']=='PASS_INDEPENDENT_OCTIC_BRIDGE_AND_COMPLETE_FINITE_VERTEX_CENSUS'
    assert review['sources'][BRIDGE]==sha(BRIDGE)
    assert review['finite_vertices']==92 and review['finite_edges']==253
    Q,N,H=R.m4(base['form']),R.m4(base['normals']),R.v4(base['center'])
    N.append(R.v8(bridge['cap']['normal']));assert len(N)==33
    oldE=R.v4(base['strict_escape']['point'])
    oldpair=lambda a,b:C.pair(a,b,Q)
    oldcov=[C.mv(Q,n) for n in N]
    oldnorms=[oldpair(n,n) for n in N]
    support_gram=[[oldpair(N[i],N[j]) for j in T] for i in T]
    assert R.pd(support_gram)[0]
    support_inverse=C.inverse(support_gram)
    wall3_contacts=[oldpair(N[3],N[i]) for i in T]
    projected3=oldnorms[3]-C.dot(wall3_contacts,C.mv(support_inverse,wall3_contacts))
    assert projected3==(3*O-R.BASIS[1])/4 and R.sg(O-2*projected3)>0
    ids=T+[3,7];rows=[oldcov[i] for i in ids]
    _,pivots=C.eliminate(rows);assert len(pivots)==9
    k=C.null_one(rows)
    p=C.mv(C.inverse(rows+[C.mv(Q,k)]),[Z]*8+[-O,Z])
    factor=(337156*O+149731*R.BASIS[1])/4691
    w=[factor*a for a in k]
    assert oldpair(p,p)==-8*O-4*R.BASIS[1]
    assert oldpair(w,w)==(5*O+3*R.BASIS[1])/2
    assert not oldpair(p,w) and R.sg(D)>0
    assert C.mv(rows,p)==[Z]*8+[-O] and C.mv(rows,w)==[Z]*9
    branches=[]
    for sign in (-1,1):
        n=[X(a,sign*b) for a,b in zip(p,w)]
        assert pair(n,n,Q)==X(2)
        vals=[dot(n,c) for c in oldcov]
        branches.append(dict(sign=sign,center_sign=pair(n,H,Q).sign(),
            escape_sign=pair(n,oldE,Q).sign(),old_pairing_signs=[x.sign() for x in vals],
            positive_walls=[i for i,x in enumerate(vals) if x.sign()>0]))
    assert not branches[0]['positive_walls'] and branches[0]['center_sign']==-1
    assert branches[0]['escape_sign']==1 and branches[1]['positive_walls']
    r=[X(a,-b) for a,b in zip(p,w)]
    normals=[[X(a) for a in n] for n in N]+[r]
    cov=[mv(Q,n) for n in normals]
    norms=[dot(n,c) for n,c in zip(normals,cov)]
    gram=[[dot(n,c) for c in cov] for n in normals]
    relations=[];neighbors=[];counts=Counter()
    for i in range(33):
        g=gram[33][i];assert g.sign()<=0
        cycle=(4*g*g)/(norms[33]*norms[i])
        if g==X():kind,label='right',2
        elif i==7:
            assert g==X(-1) and norms[i]==X(2) and cycle==X(1)
            kind,label='finite',3
        else:
            assert (cycle-4).sign()>0
            kind,label='dotted',None
        if kind!='dotted':neighbors.append(i)
        counts[kind]+=1
        relations.append(dict(wall=i,pairing=g,cartan_cycle=cycle,kind=kind,label=label))
    assert neighbors==sorted(ids) and counts==dict(right=8,finite=1,dotted=24)
    assert all(gram[i][j].sign()<=0 for i in range(34) for j in range(i))
    hp=[dot(H,c) for c in cov]
    assert all(x.sign()<0 for x in hp) and pair(H,H,Q).sign()<0
    assert dot(oldE,cov[33]).sign()>0
    facets=[]
    for i,n in enumerate(normals):
        coefficient=hp[i]/norms[i]
        V=[X(h)-coefficient*a for h,a in zip(H,n)]
        values=[dot(V,c) for c in cov];norm=pair(V,V,Q)
        assert norm.sign()<0 and pair(V,H,Q).sign()<0
        assert [v.sign() for v in values]==[0 if j==i else -1 for j in range(34)]
        facets.append(dict(wall=i,point=V,norm=norm,pairings=values))
    print('34 actual facets;33 new Coxeter pairs;old escape cut PASS',flush=True)
    # The frozen independent33-wall census supplies all old finite corners.
    old_vertices=[dict(incident=v['incident'],point=R.v4(v['point'])) for v in base['finite_vertices']]
    old_vertices += [dict(incident=v['incident'],point=R.v8(v['point']))
                     for v in bridge['complete_new_cap_finite_census']['vertices']]
    assert [v['incident'] for v in old_vertices]==bridge['finite_vertex_incidence']
    assert len(old_vertices)==92
    preservation=[];vertices=[]
    for v in old_vertices:
        value=dot(v['point'],cov[33]);assert value.sign()<0
        preservation.append(value)
        vertices.append(dict(incident=v['incident'],point=list(map(X,v['point'])),
                             norm=pair(v['point'],v['point'],Q)))
    candidates=[];new=[]
    for sub in combinations(neighbors,8):
        support=list(sub)+[33]
        G=[[gram[i][j] for j in support] for i in support]
        assert all(not x.b for row in G for x in row)
        GK=[[x.a for x in row] for row in G]
        positive,piv=R.pd(GK)
        candidates.append(dict(support=support,positive_definite=positive,pivots=piv))
        if not positive:continue
        inverse=C.inverse(GK)
        coefficients=[sum((X(c)*hp[j] for c,j in zip(row,support)),X()) for row in inverse]
        assert all(a.sign()<0 for a in coefficients)
        V=[X(H[k])-sum((a*normals[i][k] for a,i in zip(coefficients,support)),X()) for k in range(10)]
        vals=[dot(V,c) for c in cov];norm=pair(V,V,Q)
        assert norm.sign()<0 and pair(V,H,Q)==norm
        assert [v.sign() for v in vals]==[0 if i in support else -1 for i in range(34)]
        row=dict(incident=support,point=V,norm=norm,pairings=vals,
                 active_gram=GK,inverse_active_gram=inverse,projection_coefficients=coefficients)
        new.append(row);vertices.append(dict(incident=support,point=V,norm=norm))
    assert len(candidates)==9 and len(new)==3
    incidence=[v['incident'] for v in vertices]
    edges=[]
    for i,j in combinations(range(95),2):
        common=sorted(set(incidence[i])&set(incidence[j]))
        if len(common)==8:edges.append(dict(vertices=[i,j],support=common))
    assert [e for e in edges if max(e['vertices'])<92]==bridge['finite_edges']
    new_edges=[e for e in edges if max(e['vertices'])>=92]
    assert len(edges)==257 and len(new_edges)==4
    assert new_edges[0]['vertices']==[2,92]
    germs=9*95-2*len(edges);assert germs==341
    print('Complete new-facet census3;total95 finite vertices257 finite edges PASS',flush=True)
    # One literal exterior ray from the independent32-wall cone survives both
    # the old octic bridge and the new cap. No new exterior-ray census is needed.
    escape_ids=[0,1,2,3,5,6,7,12,30]
    assert base['complete_census']['hyperideal_rays'][1]==escape_ids
    V=C.null_one([oldcov[i] for i in escape_ids])
    V=[a/-oldpair(V,H) for a in V]
    assert R.sg(oldpair(V,V))>0 and oldpair(V,H)==-O
    ray_pairs=[dot(V,c) for c in cov]
    assert all(a.sign()<=0 for a in ray_pairs)
    epsilon=K.convert(F(1,8192))
    E=[v+epsilon*h for v,h in zip(V,H)]
    en=oldpair(E,E);ep=[dot(E,c) for c in cov]
    assert R.sg(en)>0 and R.sg(oldpair(E,H))<0 and all(a.sign()<0 for a in ep)
    paths=[f'scripts/{TAG}.py',BASE,BRIDGE,OLD_REVIEW,TARGET,
        'scripts/fresh_h9_joint_corner_20261015_bridge_review.py',
        'scripts/fresh_h9_h8_all_high_20261014_review.py',
        'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
    return enc(dict(status='PASS_EXACT_ESCAPE_CUTTING_WALL3_CAP_AND_FINITE_CENSUS',
        sources={p:sha(p) for p in paths},
        field=dict(base='Frozen octic bridge review primitive element',
            base_encoding='ANP coefficient lists, descending powers of the frozen primitive element',
            extra_radical_square=D,extra_radical_expression='1+sqrt(5)',
            extension_encoding='Each new scalar is [base a,base b] meaning a+b*positive_sqrt(1+sqrt5).',
            degree_scope='At most a quadratic extension of the old octic field; no irreducibility or minimal-degree claim is needed.'),
        form=Q,old_normals=N,target_support=T,projected_wall3_norm=projected3,
        cap=dict(index=33,linear_ids=ids,linear_pivots=pivots,
            particular=p,kernel=k,kernel_rescaling=factor,radical_vector=w,
            particular_norm=oldpair(p,p),radical_vector_norm=oldpair(w,w),cross_pairing=Z,
            normal=r,norm=norms[33],branch=-1,branch_checks=branches,
            center_pairing=hp[33],old_escape_pairing=dot(oldE,cov[33]),relations=relations),
        normals=normals,norms=norms,gram=gram,center=H,new_pair_counts=dict(counts),
        actual_facets=34,facet_witnesses=facets,
        old_finite_vertex_preservation_pairings=preservation,
        complete_new_cap_finite_census=dict(finite_neighbors=neighbors,candidate_supports=candidates,vertices=new),
        finite_vertices=vertices,finite_vertex_incidence=incidence,
        finite_vertex_count=95,finite_edges=edges,finite_edge_count=257,
        new_finite_edges=new_edges,open_finite_germ_count=germs,
        old_finite_graph_components=bridge['finite_graph_components'],
        finite_graph_components=components(95,edges),
        old_escape=dict(point=oldE,norm=oldpair(oldE,oldE),new_cap_pairing=dot(oldE,cov[33])),
        strict_escape=dict(point=E,norm=en,center_pairing=oldpair(E,H),pairings=ep,
            inherited_hyperideal_ray=V,ray_active=escape_ids,ray_index=1,
            ray_pairings=ray_pairs,center_epsilon='1/8192'),
        scope='All34 facets are actual;the new cap cuts the old strict escape and closes one old open finite edge. Complete finite-vertex census follows from independently certified92 old vertices,strict preservation,and all9 new-facet candidate supports. The new3 vertices and4 finite edges are exact. No complete exterior-ray census is claimed. A renewed strict spacelike escape proves the34-wall chamber still has infinite volume. No novelty or finite-volume conclusion follows from wall count.'))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text());print('PASS saved wall3 cap replay')
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertex_count'],result['finite_edge_count'],result['open_finite_germ_count'])
