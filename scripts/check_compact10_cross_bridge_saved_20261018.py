#!/usr/bin/env python3
"""Independent saved branch/native/incidence checks; numpy, no project imports.

Does not solve norm equations or reconstruct the geometric cone.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_compact10_cross_bridge_20261018'
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


def original_intervals(pri,pairs):
    dec=lambda x:tuple(map(F,x))
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    h=list(map(dec,pri['common_center']));g=[[dec(q) for q in r] for r in pri['computational_gram']]
    needed=set(sum((list(p) for p in pairs),[]));ns={}
    for r in pri['facet_witnesses']:
        if r['wall'] not in needed:continue
        hp=dec(pri['center_products'][str(r['wall'])]);assert hp[1]<0
        ns[r['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,r['point']))]
    out={}
    for i,j in pairs:
        q=(F(),F())
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);q=q[0]+z[0],q[1]+z[1]
        out[i,j]=list(map(str,q))
    return out


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='BOUNDED_COMPACT10_CHANGED_CROSS_BRIDGE'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    prior=json.loads((ROOT/'results/fresh_exterior_polar_native10_frame_20261018.json').read_text())
    assert d['snapshot']==prior['numerical_snapshot']
    snap=d['snapshot'];old={int(k):np.array(v) for k,v in snap['old_normals'].items()};center=np.array(snap['center'])
    pp=snap['profiles'];assert len(old)==57 and len(pp)==7
    profiles=[np.zeros(2)]
    for node in (0,1):
        for label in (3,4,5):
            v=np.zeros(2);v[node]=-2*np.cos(np.pi/label);profiles.append(v)
    J=np.array([[2.,1.],[1.,2.]])/3
    off=lambda a,b:float(profiles[a]@J@profiles[b])
    for pi,p in enumerate(pp):
        assert max(abs(profiles[pi]-p['contacts']))<2e-7
        assert abs(p['square']-(2-off(pi,pi)))<2e-7
        assert len(p['native_old_offsets'])==57 and set(dict(p['native_old_offsets']))==set(old)
    types={r['wall']:r['profile'] for r in d['actual_side_profiles']}
    assert types=={r['wall']:r['profile'] for r in prior['old_frame_interface']['actual_side_profiles']}
    for r in d['actual_side_profiles']:assert max(abs(profiles[r['profile']]-r['contacts']))<2e-5
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    pairs=list(combinations([8,10,33],2));iv=original_intervals(pri,pairs)
    assert [r['walls'] for r in d['original_frame_native_intervals']]==list(map(list,pairs))
    for r in d['original_frame_native_intervals']:
        i,j=r['walls'];x=float(old[i]@M@old[j])+off(types[i],types[j]);q=iv[i,j]
        assert r['exact_native_enclosure']==q and abs(x-r['native_product'])<2e-7
        assert float(F(q[0]))-2e-5<=x<=float(F(q[1]))+2e-5 and relation(x)==r['relation']
    assert relation(float(old[8]@M@old[33]))=='dotted'
    # Exact independent old-support/compact-certificate deletion debit.
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&{30,44}}
    assert len(lost)==16 and sum(8 in vs[i] for i in lost)==2
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
    assert debit['lost_finite_supports']==[list(vs[i]) for i in sorted(lost)] and debit['compact_faces_touched']==touched
    assert Counter(r['dimension'] for r in touched)=={2:12,3:1}
    assert d['deleted_walls']==[30,44] and d['retained_old_walls']==61 and d['retained_projected_old_walls']==55
    assert d['fixed_white_walls']==[0,1,2,3,5,6]
    # Independent complete low cube using the actual five-root A2 corner Gram.
    bank=d['bank'];assert bank['frame']==[8,10,33] and bank['duties']==[[0,1],[1,2]] and bank['nonedge']==[8,33]
    domain=[]
    for pi in range(7):
        for ms in product(range(2,6),repeat=3):
            good=True
            for a,b in ((0,1),(1,2)):
                gram=2*np.eye(5);gram[0,1]=gram[1,0]=-1
                for k,p in enumerate((types[bank['frame'][a]],types[bank['frame'][b]],pi)):
                    gram[k+2,:2]=profiles[p];gram[:2,k+2]=profiles[p]
                gram[2,4]=gram[4,2]=-2*np.cos(np.pi/ms[a]);gram[3,4]=gram[4,3]=-2*np.cos(np.pi/ms[b])
                good &= np.linalg.eigvalsh(gram).min()>1e-8
            if good:domain.append((pi,list(ms)))
    assert len(domain)==42 and [(r['profile'],r['labels']) for r in bank['recipes']]==domain
    keys=[i for i in old if i not in (30,44)];digest=sha256();hist=Counter();allbranches=[]
    for ri,r in enumerate(bank['recipes']):
        assert r['index']==ri and r['norm']['status']=='two_real_norm_branches'
        assert sorted(b['branch'] for b in r['branches'])==[-1,1]
        pi=r['profile'];offsets=dict(pp[pi]['native_old_offsets'])
        for b in r['branches']:
            v=np.array(b['normal']);vals={i:float(v@M@old[i])+offsets[i] for i in keys}
            assert len(b['native_products'])==55 and set(dict(b['native_products']))==set(keys)
            assert max(abs(vals[i]-x) for i,x in b['native_products'])<2e-7
            bad=[i for i in keys if relation(vals[i]) in (None,'parallel_or_unresolved')]
            assert [i for i,_ in b['blockers']]==bad
            assert max(abs(vals[i]-x) for i,x in b['blockers'])<2e-7
            assert {7,55}<=set(bad) and 'cap' not in b
            assert abs(v@M@v-pp[pi]['square'])<2e-5 and b['norm_residual']<2e-5
            assert abs(b['center_product']-(float(v@M@center)+pp[pi]['center_offset']))<2e-7
            assert all(relation(vals[i])==m for i,m in zip(bank['frame'],r['labels']))
            digest.update(json.dumps([ri,b['branch'],pi,r['labels'],b['native_products'],b['center_product']],separators=(',',':')).encode())
            hist[str(tuple(i for i,_ in b['blockers']))]+=1;allbranches.append(b)
    assert len(allbranches)==84 and digest.hexdigest()==bank['complete_branch_digest'] and dict(hist)==bank['branch_histogram']
    shared=sorted(set.intersection(*(set(i for i,_ in b['blockers']) for b in allbranches)))
    assert shared==d['bounded_obstruction']['shared_native_blockers']==[7,55]
    for i in shared:
        vals=[dict(b['native_products'])[i] for b in allbranches]
        assert d['bounded_obstruction']['native_product_ranges'][str(i)]==[min(vals),max(vals)]
        assert d['bounded_obstruction']['positive_products'][str(i)]==sum(x>2e-7 for x in vals)
    assert not bank['caps'] and not d['joins'] and not d['accepted'] and d['frozen_R_domain_count']==35
    assert bank['symbolic_high_families_not_tested']==9 and len(d['unresolved'])==1
    assert d['deletion']['score'] is None and d['deletion']['cone_diagnostic']==d['unresolved']
    assert d['unresolved'][0]['stage']=='delete30_44' and d['unresolved'][0]['boundary']['warnings']
    assert len(d['boundary_ledgers'])==1
    for rec,walls in ((d['control'],[8,10,11,30,33,44,55,63,64,65,66]),(d['deletion'],[8,10,11,33,55,63,64,65,66])):
        bd=d['boundary_ledgers'][rec['score']] if rec['score'] is not None else rec['cone_diagnostic'][0]['boundary']
        assert len(rec['whole_faces'])==len(walls)
        for f,w in zip(rec['whole_faces'],walls):
            rays=[r for r in bd['rays'] if w in r['active']]
            assert f==dict(wall=w,rays=rays,all_rays_finite=bool(rays) and all(r['kind']=='finite' for r in rays))
        assert rec['finite_rays']==bd['counts']['finite'] and rec['open_ends']==bd['open_ends']
        assert rec['exterior_supports']==[r['active'] for r in bd['rays'] if r['kind']=='exterior']
        assert rec['open_edges']==[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')]
    bd=d['boundary_ledgers'][0]
    assert d['control']['open_ends']==10 and d['control']['finite_rays']==26
    assert not bd['warnings'] and bd['PD_status']=='numerically_separated_acute_incidence'
    assert bd['open_ends']==bd['PD_open_ends']==sum(1 if e['kind']=='one_open' else 2 if e['kind']=='two_open' else 0 for e in bd['edges'])
    assert {tuple(r['active']) for r in bd['rays'] if r['kind']=='finite'}=={tuple(r) for r in bd['PD_finite_supports']}
    print('PASS_COMPACT10_CROSS_BRIDGE_SAVED',dict(old_supports=199,debit=[16,12,1],on8=2,recipes=42,branches=84,bridges=0,joins=0))


if __name__=='__main__':run()
