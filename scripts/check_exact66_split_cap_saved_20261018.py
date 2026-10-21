#!/usr/bin/env python3
"""Independent saved low-cube/profile/norm/native check; no project imports."""
import json,math
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,permutations,product
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_exact66_split_cap_native_pilot_20261018'
d=json.loads((ROOT/f'results/{TAG}.json').read_text())
assert d['status']=='PASS_BOUNDED_EXACT66_SPLIT_CAP_NUMERICAL_PILOT'
for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p

def relation(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return None
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=round(math.pi/math.acos(-x/2))
    if m<3:return None
    theta=lambda n:2*math.cos(math.pi/n)
    spacing=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return m if abs(x+theta(m))<=min(2e-7,spacing/100)else None

@lru_cache(None)
def canonical(flat):
    n=math.isqrt(len(flat));assert n*n==len(flat)
    return min(tuple(flat[p[i]*n+p[j]]for i in range(n)for j in range(i+1,n))for p in permutations(range(n)))

catalog={}
def register(name,n,edges):
    g=[[2]*n for _ in range(n)]
    for i in range(n):g[i][i]=1
    for i,j,m in edges:g[i][j]=g[j][i]=m
    catalog[n,canonical(tuple(sum(g,[])))]=name
for n in range(1,6):register('A'+str(n),n,[(i,i+1,3)for i in range(n-1)])
for n in range(2,6):register('B'+str(n),n,[(i,i+1,4 if i==0 else 3)for i in range(n-1)])
for n in(4,5):register('D'+str(n),n,[(i,i+1,3)for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)])
register('F4',4,[(0,1,3),(1,2,4),(2,3,3)])
for n in(3,4):register('H'+str(n),n,[(i,i+1,5 if i==0 else 3)for i in range(n-1)])
register('I2(5)',2,[(0,1,5)])
def duty(pi,b,m,n):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    if pi:g[0 if pi<=3 else 1][4]=g[4][0 if pi<=3 else 1]=3+(pi-1)%3
    for i,j,l in[(2,3,b),(2,4,m),(3,4,n)]:g[i][j]=g[j][i]=l
    unseen=set(range(5));names=[]
    while unseen:
        todo=[min(unseen)];comp=set()
        while todo:
            i=todo.pop()
            if i in comp:continue
            comp.add(i);todo.extend(j for j in unseen if j!=i and g[i][j]!=2)
        unseen-=comp;vs=sorted(comp)
        name=catalog.get((len(vs),canonical(tuple(g[i][j]for i in vs for j in vs))))
        if name is None:return None
        names.append(name)
    return sorted(names)

z=d['fixed_seed'];old=json.loads((ROOT/'results'/f"{z['stem']}.json").read_text())
assert(z['exact_facets'],z['exact_finite_vertices'])==(66,205)
assert old['physical']['summary']['total_walls']==66 and old['physical']['summary']['complete_finite_vertices']==205
assert old['physical']['deleted_original_walls']==[]
controls=[r for r in old['numerical']['accepted']if r['profile']==0 and r['labels']==[2,2,2]]
assert len(controls)==1 and z['wall77_control']==controls[0]
expected={**old['numerical']['snapshot']['old_normals'],'77':controls[0]['normal']}
assert z['projected_normals']==expected
nd={int(i):np.array(v)for i,v in z['projected_normals'].items()};assert len(nd)==60
original=sorted(set(nd)-{75,76,77});assert len(original)==57
index={i:k for k,i in enumerate(original)}
metric=np.array(z['metric']);center=np.array(z['center']);profiles=z['profiles']
assert z['metric']==np.diag([-1,1,1,1]).tolist()
assert profiles==old['numerical']['snapshot']['profiles']and len(profiles)==7
assert z['center']==old['numerical']['snapshot']['center']
assert z['profile_types']=={'75':0,'76':0,'77':0}
invA2=np.array([[2,1],[1,2]])/3
for pi,pr in enumerate(profiles):
    cov=np.zeros(2)
    if pi:cov[0 if pi<=3 else 1]=-2*math.cos(math.pi/(3+(pi-1)%3))
    assert np.max(abs(cov-np.array(pr['contacts'])))<1e-12
    assert abs(pr['square']-(2-cov@invA2@cov))<1e-12

assert len(d['banks'])==2
for bank,new,frame,bonds,eligible,branches in zip(d['banks'],[78,79],[[30,52,47],[47,76,30]],[[4,2],[2,2]],[38,96],[76,192],strict=True):
    assert(bank['new_wall'],bank['frame'],bank['frame_bonds'])==(new,frame,bonds)
    expect={(pi,ms)for pi,ms in product(range(7),product(range(2,6),repeat=3))};seen=set();ok=roots=0;hist=Counter()
    assert len(bank['recipes'])==448
    for r in bank['recipes']:
        pi=r['profile'];ms=tuple(r['labels']);assert(pi,ms)in expect and(pi,ms)not in seen;seen.add((pi,ms))
        duties=[duty(pi,b,ms[j],ms[j+1])for j,b in enumerate(bonds)]
        saved=[None if x is None else sorted({'I2(3)':'A2','I2(4)':'B2'}.get(t,t)for t in x)for x in r['duties']]
        assert saved==duties
        if any(x is None for x in duties):assert r['status']=='finite_duty_rejected'and r['branches']==[];continue
        ok+=1;assert r['status']==r['norm']['status']=='two_real_norm_branches'
        assert len(r['branches'])==2 and [b['branch']for b in r['branches']]==[-1,1]
        ff=np.array([nd[i]for i in frame]);gram=ff@metric@ff.T;rhs=-2*np.cos(np.pi/np.array(ms))
        assert abs(np.linalg.det(gram))>1e-7 and np.max(abs(gram-np.array(r['norm']['gram'])))<1e-8
        normals=[]
        for b in r['branches']:
            v=np.array(b['normal']);normals.append(v);roots+=1
            assert abs(v@metric@v-profiles[pi]['square'])<2e-8
            assert np.max(abs(ff@metric@v-rhs))<2e-8
            vals={i:float(v@metric@w)+(profiles[pi]['offsets'][index[i]]if i in index else 0)for i,w in nd.items()}
            assert len(b['products'])==60 and[i for i,x in b['products']]==sorted(nd)
            assert all(abs(vals[i]-x)<2e-8 for i,x in b['products'])
            bad=[i for i in nd if relation(vals[i])is None];hist.update(bad)
            assert b['blockers']==[[i,dict(b['products'])[i]]for i in bad]
            cp=float(v@metric@center)+profiles[pi]['center_offset']
            assert abs(cp-b['center_product'])<2e-8
            assert b['native']==(not bad and cp<0)==False
        assert np.linalg.norm(normals[0]-normals[1])>1e-5
    assert seen==expect and(ok,roots)==(eligible,branches)
    assert bank['native_survivors']==[]
    assert bank['blocking_histogram']=={str(i):n for i,n in hist.items()}
    assert bank['summary']==dict(complete_recipes=448,finite_eligible=eligible,norm_branches=branches,native_survivors=0)

tri=[[30,52,78],[52,47,78],[47,78,79],[47,76,79],[76,30,79],[30,78,79]]
t=d['topology'];assert t['triangles']==tri and t['boundary']==[30,52,47,76]
edges=Counter(tuple(sorted(e))for f in tri for e in combinations(f,2));vertices=set(sum(tri,[]))
assert(len(vertices),len(edges),len(tri))==(t['vertices'],t['edges'],t['faces'])==(6,11,6)
assert set(edges.values())=={1,2}and len(vertices)-len(edges)+len(tri)==1
assert {e for e,n in edges.items()if n==1}=={tuple(sorted(e))for e in[(30,52),(52,47),(47,76),(76,30)]}
for v in vertices:
    links=[set(f)-{v}for f in tri if v in f];degree=Counter(x for e in links for x in e)
    assert sorted(degree.values())==([2]*len(degree)if v in(78,79)else[1,1]+[2]*(len(degree)-2))
assert(30,47)not in edges
assert d['joins']==d['cones']==[]
assert d['summary']==dict(K_native=0,L_native=0,Cartesian_pairs=0,accepted_joins=0,cones_checked=0)
assert d['deleted_original_walls']==[]and d['old_finite_vertex_debit']==0
assert d['scope']==dict(numerical_new_normals=True,profiles=list(range(7)),labels=[2,3,4,5],two_independent_banks=True,all_declared_pairs_checked=True,high_label_families_untested=True,exact_new_normal_certified=False,new_finite_volume_example=False)
print('PASS896 recipes/268 numerical norm roots/all retained gates; empty join; zero deletion')
