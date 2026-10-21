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
TAG='fresh_moving33_compact10_bridge_20261018'
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
    assert d['status']=='BOUNDED_MOVING33_COMPACT10_BRIDGE'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    snap=d['numerical_snapshot'];old={int(k):np.array(v) for k,v in snap['old_normals'].items()};center=np.array(snap['center'])
    pp=snap['profiles'];assert len(old)==57 and len(pp)==7 and snap['metric']==M.tolist()
    labels=snap['labels'];assert set(labels)==set(old) and len(labels)==57
    profiles=[np.zeros(2)]
    for node in (0,1):
        for label in (3,4,5):
            v=np.zeros(2);v[node]=-2*np.cos(np.pi/label);profiles.append(v)
    J=np.array([[2.,1.],[1.,2.]])/3;off=lambda a,b:float(profiles[a]@J@profiles[b])
    for pi,p in enumerate(pp):
        assert max(abs(profiles[pi]-p['contacts']))<2e-7
        assert abs(p['square']-(2-off(pi,pi)))<2e-7 and len(p['offsets'])==57
    types={r['wall']:r['profile'] for r in d['old_frame_interface']['actual_side_profiles']}
    for r in d['old_frame_interface']['actual_side_profiles']:
        assert max(abs(profiles[r['profile']]-r['stored_A2_contacts']))<2e-5
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    pairs=list(combinations([7,8,10,33,66],2));iv=original_intervals(pri,pairs)
    assert [r['walls'] for r in d['original_native_intervals']]==list(map(list,pairs))
    for r in d['original_native_intervals']:
        i,j=r['walls'];x=float(old[i]@M@old[j])+off(types[i],types[j]);q=iv[i,j]
        assert r['exact_native_enclosure']==q and abs(x-r['native_product'])<2e-7
        assert float(F(q[0]))-2e-5<=x<=float(F(q[1]))+2e-5 and relation(x)==r['relation']
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&{30,33,44}}
    assert len(lost)==18 and sum(8 in vs[i] for i in lost)==2
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
    assert d['deleted_walls']==[30,33,44] and d['retained_old_walls']==60
    base={i:v for i,v in old.items() if i not in d['deleted_walls']}
    def ledger(rec):
        if rec['score'] is None:
            assert len(rec['cone_diagnostic'])==1;bd=rec['cone_diagnostic'][0]['boundary']
            assert rec['cone_diagnostic'][0] in d['unresolved'] and bd['warnings']
        else:
            bd=d['boundary_ledgers'][rec['score']]
            assert not bd['warnings'] and bd['PD_status']=='numerically_separated_acute_incidence'
            assert bd['open_ends']==bd['PD_open_ends']
            assert {tuple(r['active']) for r in bd['rays'] if r['kind']=='finite'}=={tuple(r) for r in bd['PD_finite_supports']}
        assert rec['finite_rays']==bd['counts']['finite'] and rec['open_ends']==bd['open_ends']
        assert rec['exterior_supports']==[r['active'] for r in bd['rays'] if r['kind']=='exterior']
        assert rec['open_edges']==[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')]
        assert len({f['wall'] for f in rec['whole_faces']})==len(rec['whole_faces'])
        for f in rec['whole_faces']:
            rr=[r for r in bd['rays'] if f['wall'] in r['active']]
            assert f==dict(wall=f['wall'],rays=rr,all_rays_finite=bool(rr) and all(r['kind']=='finite' for r in rr))
    totals=Counter();allsecond=[]
    def checkbank(bank,nd,tt,second):
        frame=bank['frame'];duties=bank['duties'];pt=lambda i:tt.get(i,types.get(i,0))
        assert bank['frame_profiles']==[pt(i) for i in frame]
        bonds=bank['duty_bonds'];domain=[]
        for a,b in duties:
            assert relation(float(nd[frame[a]]@M@nd[frame[b]])+off(pt(frame[a]),pt(frame[b])))==bonds[3*a+b]
        for pi in range(7):
            for ms in product(range(2,6),repeat=3):
                good=True
                for a,b in duties:
                    gram=2*np.eye(5);gram[0,1]=gram[1,0]=-1
                    for k,p in enumerate((pt(frame[a]),pt(frame[b]),pi)):
                        gram[k+2,:2]=profiles[p];gram[:2,k+2]=profiles[p]
                    gram[2,3]=gram[3,2]=-2*np.cos(np.pi/bonds[3*a+b])
                    gram[2,4]=gram[4,2]=-2*np.cos(np.pi/ms[a]);gram[3,4]=gram[4,3]=-2*np.cos(np.pi/ms[b])
                    good &= np.linalg.eigvalsh(gram).min()>1e-8
                if good:domain.append((pi,list(ms)))
        assert [(r['profile'],r['labels']) for r in bank['recipes']]==domain
        assert len(domain)==bank['complete_low_recipe_count']
        digest=sha256();hist=Counter();seen=set();capindex=0
        for ri,r in enumerate(bank['recipes']):
            assert r['index']==ri and r['norm']['status']=='two_real_norm_branches'
            assert sorted(b['branch'] for b in r['branches'])==[-1,1]
            pi=r['profile'];offsets=dict(zip(labels,pp[pi]['offsets']))
            for b in r['branches']:
                v=np.array(b['normal']);vals={i:float(v@M@w)+(offsets[i] if i in offsets else off(pi,tt[i])) for i,w in nd.items()}
                assert len(b['native_products'])==len(nd) and set(dict(b['native_products']))==set(nd)
                assert max(abs(vals[i]-x) for i,x in b['native_products'])<2e-7
                bad=[i for i in nd if relation(vals[i]) in (None,'parallel_or_unresolved')]
                assert [i for i,_ in b['blockers']]==bad
                assert all(abs(vals[i]-x)<2e-7 for i,x in b['blockers'])
                assert abs(v@M@v-pp[pi]['square'])<2e-5 and b['norm_residual']<2e-5
                cp=float(v@M@center)+pp[pi]['center_offset'];assert abs(b['center_product']-cp)<2e-7
                assert all(relation(vals[i])==m for i,m in zip(frame,r['labels']))
                digest.update(json.dumps([ri,b['branch'],pi,r['labels'],b['native_products'],b['center_product']],separators=(',',':')).encode())
                hist[str(tuple(bad))]+=1;totals['branches']+=1
                if second:assert 7 in bad and 'cap' not in b;allsecond.append(b)
                if bad or cp>=-1e-8:assert 'cap' not in b;continue
                key=(pi,)+tuple(np.round(v,8))
                if key in seen:assert b['duplicate'] and 'cap' not in b;continue
                seen.add(key);assert b['cap']==capindex
                cap=bank['caps'][capindex];capindex+=1
                assert cap['index']==b['cap'] and cap['recipe']==ri and cap['branch']==b['branch']
                assert cap['profile']==pi and cap['labels']==r['labels'] and cap['normal']==b['normal']
                assert cap['native_products']==b['native_products'] and cap['center_product']==b['center_product']
                assert cap['restores_old33']==bool(pi==types[33] and np.max(abs(v-old[33]))<2e-6)
                assert cap['score'] is None and cap['open_ends']==14 and cap['finite_rays']==16
                ledger(cap)
        assert capindex==len(bank['caps']) and digest.hexdigest()==bank['complete_branch_digest'] and dict(hist)==bank['branch_histogram']
        totals['recipes']+=len(domain)
    first=d['first_bank'];assert first['frame']==[7,10,66] and first['duties']==[[0,1],[0,2]]
    checkbank(first,base,{},False);assert len(first['caps'])==12 and sum(r['restores_old33'] for r in first['caps'])==1
    assert len(d['records'])==12
    for r,a in zip(d['records'],first['caps']):
        assert r['S']==a['index'] and r['S_restores_old33']==a['restores_old33']
        b=r['bank'];assert b['frame']==[8,10,67] and b['duties']==[[0,1],[1,2]]
        checkbank(b,{**base,67:np.array(a['normal'])},{67:a['profile']},True);assert not b['caps']
    assert totals==dict(recipes=301,branches=602)
    assert len(allsecond)==560 and sorted(set.intersection(*(set(i for i,_ in b['blockers']) for b in allsecond)))==[7]
    assert not d['accepted'] and len(d['unresolved'])==13 and len(d['boundary_ledgers'])==1
    assert Counter(u['stage'] for u in d['unresolved'])=={'delete30_33_44':1,'moving33':12}
    ledger(d['control']);ledger(d['deletion'])
    assert d['control']['finite_rays']==26 and d['control']['open_ends']==10
    assert d['deletion']['finite_rays']==14 and d['deletion']['open_ends']==14
    assert d['summary']==dict(first_recipes=21,first_branches=42,first_native=12,moving_S=11,second_banks=12,second_recipes=280,second_branches=560,native_pairs=0,clean_pairs=0,whole_compact10=0,whole_new_faces=0,fewer_than_ten_ends=0,unresolved=13,complete_low_domains=True)
    print('PASS_MOVING33_COMPACT10_BRIDGE_SAVED',dict(old_supports=199,debit=[18,12,1],on8=2,recipes=301,branches=602,native_S=12,bridges=0,provisional_cones=13))


if __name__=='__main__':run()
