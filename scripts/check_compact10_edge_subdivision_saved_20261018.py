#!/usr/bin/env python3
"""Independent saved-domain/interface/debit checker; no project imports."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact10_edge_subdivision_20261018'
M=np.diag([-1.,1.,1.,1.])


def relation(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return 'parallel_or_unresolved'
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=round(np.pi/np.arccos(-x/2))
    if m<3:return None
    theta=lambda n:2*np.cos(np.pi/n)
    spacing=min(abs(theta(m)-theta(m-1)),abs(theta(m+1)-theta(m)))
    return int(m) if abs(x+theta(m))<=min(2e-7,spacing/100) else None


def native_interval(pri):
    dec=lambda x:tuple(map(F,x))
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    h=list(map(dec,pri['common_center']));g=[[dec(x) for x in r] for r in pri['computational_gram']]
    ns={};gates={r['wall']:r for r in pri['native_gates']}
    for row in pri['facet_witnesses']:
        hp=dec(pri['center_products'][str(row['wall'])]);assert hp[1]<0
        ns[row['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,row['point']))]
    def pair(i,j):
        if i==j:return F(2),F(2)
        if 66 in(i,j):
            r=gates[j if i==66 else i]
            return (F(0),F(0)) if r['label']==2 else dec(r['product'])
        q=(F(0),F(0))
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);q=(q[0]+z[0],q[1]+z[1])
        return q
    return pair


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='BOUNDED_FINITE_COMPACT10_EDGE_SUBDIVISION'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    pair=native_interval(pri);profiles=[[0.,0.]]
    for node in(0,1):
        for m in(3,4,5):
            p=[0.,0.];p[node]=-2*np.cos(np.pi/m);profiles.append(p)
    J=np.array([[2.,1.],[1.,2.]])/3
    inf=d['old_frame_interface'];assert len(inf['actual_side_profiles'])==18
    for r in inf['actual_side_profiles']:
        p=np.array(profiles[r['profile']]);assert max(abs(p-np.array(r['stored_A2_contacts'])))<2e-5
        assert abs(r['intrinsic_square']-(2-p@J@p))<2e-5
        assert [z['white'] for z in r['original_native_checks']]==[0,1,2,3,5,6]
        for z in r['original_native_checks']:
            q=pair(r['wall'],z['white']);assert list(map(str,q))==z['exact_native_enclosure']
            assert float(q[0])-2e-5<=z['numerical_product']<=float(q[1])+2e-5
    assert {r['wall']:r['profile'] for r in inf['actual_side_profiles']}[10]==1
    for r in inf['all_old_frame_pairs']+d['frame_native_checks']:
        q=pair(*r['walls']);assert list(map(str,q))==r['exact_native_enclosure']
        assert float(q[0])-2e-5<=r['native_product']<=float(q[1])+2e-5
        assert relation(r['native_product'])==r.get('relation',r.get('native_relation'))
    expected=sorted({tuple(sorted(e)) for f in([7,55,33],[10,44,33]) for e in combinations(f,2)})
    assert [r['walls'] for r in d['frame_native_checks']]==list(map(list,expected))
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if 30 in v};assert len(lost)==13 and all(8 not in vs[i] for i in lost)
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,ii in faces.items():
        if not(lost&set(ii)) or len(ii)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i] if j!=k)]==2 for i in ii for k in vs[i] if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i]) for i in ii]))
    debit=d['exact_deletion_debit']
    assert debit['lost_finite_supports']==[list(vs[i]) for i in sorted(lost)]
    assert debit['compact_faces_touched']==touched and Counter(r['dimension'] for r in touched)=={2:10,3:1}
    assert d['deleted_walls']==[30] and d['retained_old_walls']==62
    top=d['topology'];ts=top['triangles'];boundary=[10,33,67,44]
    assert ts==[[10,68,33],[67,68,33],[10,68,44],[67,68,44]] and top['boundary']==boundary
    ee=Counter(tuple(sorted(e)) for f in ts for e in combinations(f,2))
    assert len(set(sum(ts,[])))-len(ee)+len(ts)==1 and (10,67) not in ee and (33,44) not in ee
    assert {e for e,n in ee.items() if n==1}=={tuple(sorted(e)) for e in zip(boundary,boundary[1:]+boundary[:1])}
    assert len(top['links'])==5
    for link in top['links']:
        v=link['vertex'];els=[list(x for x in f if x!=v) for f in ts if v in f];assert link['edges']==els
        deg=Counter(x for e in els for x in e);seen={next(iter(deg))}
        while True:
            fresh=seen|{b for a,b in els if a in seen}|{a for a,b in els if b in seen}
            if fresh==seen:break
            seen=fresh
        assert seen==set(deg) and sorted(deg.values())==([2]*4 if v==68 else [1,1]+[2]*(len(deg)-2))
    assert len(d['banks'])==2 and d['low_labels']==[2,3,4,5]
    counts=Counter();single=[];oldkeys=None
    for bi,b in enumerate(d['banks']):
        assert b['frame']==[[7,55,33],[10,44,33]][bi] and b['duties']==[[[0,1]],[[0,1],[0,2]]][bi]
        assert b['complete_low_recipe_count']==len(b['recipes'])==[64,42][bi]
        assert len({(r['profile'],*r['labels']) for r in b['recipes']})==len(b['recipes'])
        assert len(b['caps'])==1;refs=[]
        for ri,r in enumerate(b['recipes']):
            counts['norm_equations']+=1
            assert r['index']==ri and 0<=r['profile']<7 and len(r['labels'])==3 and all(m in range(2,6) for m in r['labels'])
            assert r['norm']['status']=='two_real_norm_branches' and sorted(br['branch'] for br in r['branches'])==[-1,1]
            for br in r['branches']:
                counts['norm_branches']+=1
                assert all(relation(v) in(None,'parallel_or_unresolved') for _,v in br['blockers'])
                good=not br['blockers'] and br['center_product']<-1e-8 and br['norm_residual']<2e-5
                assert good==('cap' in br or br.get('duplicate',False))
                if 'cap' in br:
                    refs.append(br['cap']);assert b['caps'][br['cap']]['recipe']==ri and b['caps'][br['cap']]['branch']==br['branch']
                if len(br['blockers'])==1:single.append([bi,r['profile'],r['labels'],br['branch'],br['blockers']])
        assert refs==[0];c=b['caps'][0]
        assert c['index']==0 and c['profile']==0 and c['labels']==[[5,2,2],[2,2,2]][bi]
        assert c['restores_old_numerically'] and c['cone_validated'] and c['score'] is not None
        vals=dict(c['native_products']);assert len(vals)==len(c['native_products'])==56 and 30 not in vals
        if oldkeys is None:oldkeys=set(vals)
        assert set(vals)==oldkeys and all(relation(v) not in(None,'parallel_or_unresolved') for v in vals.values())
        assert all(relation(vals[i])==m for i,m in zip(b['frame'],c['labels']))
        v=np.array(c['normal']);assert abs(v@M@v-2)<2e-5 and c['center_product']<-1e-8
    assert counts==Counter(d['counts'])=={'norm_equations':106,'norm_branches':212}
    a,b=(bank['caps'][0] for bank in d['banks']);assert max(abs(np.array(a['normal'])-np.array(b['normal'])))<2e-6
    assert a['changed_incidence_relations']=={'10':2,'44':2} and not a['eligible_for_subdivision']
    p10=np.array(profiles[1]);pr=np.array(profiles[a['profile']])
    ip=dict(a['native_products'])[10]-float(pr@J@p10)
    det=(2-pr@J@pr)*(2-p10@J@p10)-ip*ip
    assert abs(ip-a['projected10_pair']['product'])<2e-7 and abs(det-a['projected10_pair']['determinant'])<2e-7 and det>0
    assert len(single)==2
    assert [(x[0],x[1],x[2],x[3],x[4][0][0]) for x in single]==[(0,0,[5,2,4],1,10),(0,0,[5,2,5],1,10)]
    assert all(-1<x[4][0][1]<0 for x in single)
    assert d['selected_R']==d['pairs']==d['coherent_pairs']==d['completed_connector_pairs']==[]
    assert d['complete_pair_count']==0 and d['pair_histogram']=={} and d['complete_pair_digest']==sha256().hexdigest()
    assert d['required_corners']==ts+[[7,55,67]]
    assert len(d['boundary_ledgers'])==3 and d['deletion_score'] is None and len(d['unresolved'])==1
    assert d['unresolved'][0]['stage']=='delete30' and d['unresolved'][0]['boundary']['warnings']
    for ledger in d['boundary_ledgers']:
        assert not ledger['warnings'] and ledger['PD_status']=='numerically_separated_acute_incidence'
        assert ledger['counts']=={'finite':26,'exterior':6}
        assert ledger['open_ends']==ledger['PD_open_ends']==10
        assert sum(1 if e['kind']=='one_open' else 2 if e['kind']=='two_open' else 0 for e in ledger['edges'])==10
        assert {tuple(r['active']) for r in ledger['rays'] if r['kind']=='finite'}=={tuple(r) for r in ledger['PD_finite_supports']}
        assert all(r['kind']=='finite' for r in ledger['rays'] if 8 in r['active'] or 10 in r['active'])
    print('PASS_COMPACT10_EDGE_SUBDIVISION_SAVED',dict(old_supports=199,debit=[13,10,1],recipes=106,branches=212,banks=[1,1],selected=0,pairs=0,ledgers=3,one_wall_gap_rejections=2))


if __name__=='__main__':run()
