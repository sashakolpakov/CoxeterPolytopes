#!/usr/bin/env python3
"""Saved-data audit only: no project imports, root solves, or cone reconstruction."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_joint10_30_native_replacement_20261018'
M=np.diag([-1.,1.,1.,1.])
J=np.array([[2.,1.],[1.,2.]])/3


def relation(x):
    if x < -2-2e-7:return 'dotted'
    if abs(x+2)<=2e-7:return 'parallel_or_unresolved'
    if abs(x)<=2e-7:return 2
    if x>0 or x < -2:return None
    m=round(np.pi/np.arccos(-x/2))
    if m<3:return None
    t=lambda n:2*np.cos(np.pi/n)
    spacing=min(abs(t(m)-t(m-1)),abs(t(m+1)-t(m)))
    return int(m) if abs(x+t(m))<=min(2e-7,spacing/100) else None


def interval_pair(pri):
    dec=lambda q:tuple(map(F,q))
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    h=list(map(dec,pri['common_center']));g=[[dec(q) for q in row] for row in pri['computational_gram']]
    ns={};gates={r['wall']:r for r in pri['native_gates']}
    for row in pri['facet_witnesses']:
        hp=dec(pri['center_products'][str(row['wall'])]);assert hp[1]<0
        ns[row['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,row['point']))]
    def pair(i,j):
        if i==j:return F(2),F(2)
        if 66 in (i,j):
            r=gates[j if i==66 else i]
            return (F(),F()) if r['label']==2 else dec(r['product'])
        out=(F(),F())
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);out=(out[0]+z[0],out[1]+z[1])
        return out
    return pair


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='BOUNDED_JOINT10_30_NATIVE_REPLACEMENT'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    intervals=interval_pair(pri)
    cs=[[0.,0.]]
    for k in (0,1):
        for m in (3,4,5):
            c=[0.,0.];c[k]=-2*np.cos(np.pi/m);cs.append(c)
    inf=d['old_frame_interface'];assert len(inf['actual_side_profiles'])==18
    oldtypes={r['wall']:r['profile'] for r in inf['actual_side_profiles']}
    for r in inf['actual_side_profiles']:
        c=np.array(cs[r['profile']]);assert max(abs(c-np.array(r['stored_A2_contacts'])))<2e-5
        assert abs(r['intrinsic_square']-(2-c@J@c))<2e-5
        assert [z['white'] for z in r['original_native_checks']]==[0,1,2,3,5,6]
        for z in r['original_native_checks']:
            q=intervals(r['wall'],z['white']);assert list(map(str,q))==z['exact_native_enclosure']
            assert float(q[0])-2e-5<=z['numerical_product']<=float(q[1])+2e-5
    assert oldtypes[10]==1
    for r in inf['all_old_frame_pairs']+d['old_frame_native_checks']+d['changed_frame_native_checks']:
        q=intervals(*r['walls']);assert list(map(str,q))==r['exact_native_enclosure']
        assert float(q[0])-2e-5<=r['native_product']<=float(q[1])+2e-5
        assert relation(r['native_product'])==r.get('relation',r.get('native_relation'))
    assert [r['walls'] for r in d['changed_frame_native_checks']]==[[8,44],[8,63],[44,63]]
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&{10,30}};assert len(lost)==75
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
    assert Counter(r['dimension'] for r in touched)=={2:313,3:325,4:244,5:129,6:46,7:10,8:1}
    assert d['lost_old8_finite_supports']==[list(vs[i]) for i in sorted(lost) if 8 in vs[i]] and len(d['lost_old8_finite_supports'])==44
    assert d['deleted_walls']==[10,30] and d['retained_old_wall_count']==61
    assert debit['compact_S6_sides_touched']==[8,10,11,44,55,63,64]
    snap=d['numerical_replay_snapshot'];old={i:np.array(v) for i,v in snap['old_normals']};center=np.array(snap['center'])
    assert len(old)==57 and len(snap['profiles'])==7 and abs(center@M@center+1)<2e-8
    offsets=[]
    for pi,p in enumerate(snap['profiles']):
        assert max(abs(np.array(cs[pi])-np.array(p['contacts'])))<2e-8
        assert abs(p['square']-(2-np.array(cs[pi])@J@np.array(cs[pi])))<2e-8
        oo=dict(p['native_old_offsets']);assert len(oo)==len(p['native_old_offsets'])==57 and set(oo)==set(old);offsets.append(oo)
    base={i:v for i,v in old.items() if i not in (10,30)};assert len(base)==55
    seed=json.loads((ROOT/'results/fresh_compact10_edge_subdivision_20261018.json').read_text())
    sole=[(r,b) for r in seed['banks'][0]['recipes'] for b in r['branches'] if len(b['blockers'])==1 and b['blockers'][0][0]==10 and b['center_product']<-1e-8]
    assert len(sole)==2 and len(d['prefixes'])==len(d['new10_banks'])==3
    counts=Counter();digest=sha256();pairs=[];single=[]
    def corner(pa,pb,pc,mab,mac,mbc):
        g=np.eye(5)*2;g[0,1]=g[1,0]=-1
        for k,pi in enumerate((pa,pb,pc),2):g[:2,k]=g[k,:2]=cs[pi]
        for a,b,m in ((2,3,mab),(2,4,mac),(3,4,mbc)):g[a,b]=g[b,a]=-2*np.cos(np.pi/m) if m>2 else 0
        return np.linalg.eigvalsh(g)[0]>1e-8
    for bi,(p,b) in enumerate(zip(d['prefixes'],d['new10_banks'])):
        assert p['index']==bi and p['profile']==0 and p['labels']==[[5,2,4],[5,2,5],[5,2,2]][bi] and p['branch']==1
        v=np.array(p['normal']);assert abs(v@M@v-2)<2e-5 and p['norm_residual']<2e-5
        assert abs(p['center_product']-(v@M@center+snap['profiles'][0]['center_offset']))<2e-8
        before={i:float(v@M@w)+offsets[0][i] for i,w in old.items() if i!=30}
        assert len(p['before_dropping10_native_products'])==56 and set(dict(p['before_dropping10_native_products']))==set(before)
        assert max(abs(q-before[i]) for i,q in p['before_dropping10_native_products'])<2e-7
        assert [i for i,q in before.items() if relation(q) in(None,'parallel_or_unresolved')]==([] if bi==2 else [10])
        if bi<2:
            sr,sb=sole[bi];assert p['source_recipe']==sr['index'] and p['before_dropping10_blockers']==sb['blockers']
        assert p['native44_relation']==relation(before[44])==(2 if bi==2 else 'dotted')
        assert len(p['retained_native_products'])==55 and set(dict(p['retained_native_products']))==set(base)
        assert all(relation(q) not in(None,'parallel_or_unresolved') for i,q in p['retained_native_products'])
        assert max(abs(q-before[i]) for i,q in p['retained_native_products'])<2e-7
        nd={**base,67:v};types={**oldtypes,67:0}
        assert b['frame']==([44,33,67] if bi==2 else [8,44,63])
        assert b['duties']==([[0,2],[1,2]] if bi==2 else [[0,1],[1,2]])
        assert len(b['frame_native_checks'])==3 and b['symbolic_high_families_not_tested']==10
        bonds={}
        for r,(a,c) in zip(b['frame_native_checks'],combinations(range(3),2)):
            i,j=b['frame'][a],b['frame'][c];assert r['walls']==[i,j]
            q=float(nd[i]@M@nd[j])+np.array(cs[types[i]])@J@np.array(cs[types[j]])
            assert abs(q-r['native_product'])<2e-7 and relation(q)==r['relation'];bonds[a,c]=r['relation']
        recipes=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
            if all(corner(types[b['frame'][a]],types[b['frame'][c]],pi,bonds[a,c],ms[a],ms[c]) for a,c in b['duties'])]
        assert sorted(recipes)==sorted((r['profile'],r['labels']) for r in b['recipes']) and len(recipes)==b['complete_low_recipe_count']==96
        refs=[];one65=[]
        for ri,r in enumerate(b['recipes']):
            assert r['index']==ri and r['norm']['status']=='two_real_norm_branches' and sorted(q['branch'] for q in r['branches'])==[-1,1]
            counts['new10_norm_equations']+=1
            for q in r['branches']:
                counts['new10_norm_branches']+=1;pi=r['profile'];w=np.array(q['normal'])
                assert abs(w@M@w-snap['profiles'][pi]['square'])<2e-5 and q['norm_residual']<2e-5
                assert abs(q['center_product']-(w@M@center+snap['profiles'][pi]['center_offset']))<2e-7
                vals={i:float(w@M@x)+(offsets[pi][i] if i in old else np.array(cs[pi])@J@np.array(cs[types[i]])) for i,x in nd.items()}
                assert len(q['native_products'])==len(vals)==56 and set(dict(q['native_products']))==set(vals)
                assert max(abs(x-vals[i]) for i,x in q['native_products'])<2e-7
                bad=[[i,x] for i,x in q['native_products'] if relation(x) in(None,'parallel_or_unresolved')]
                assert bad==q['blockers']
                assert all(relation(vals[i])==m for i,m in zip(b['frame'],r['labels']))
                good=not bad and q['center_product']<-1e-8
                assert good==('cap' in q or q.get('duplicate',False))
                if len(bad)==1 and bad[0][0]==65:assert bad[0][1]>0;one65.append(bad[0][1])
                if 'cap' in q:
                    refs.append(q['cap']);c=b['caps'][q['cap']];assert c['recipe']==ri and c['branch']==q['branch'] and c['normal']==q['normal']
                    digest.update(json.dumps([bi,c['index'],pi,c['labels'],c['branch'],c['native_products']],separators=(',',':')).encode());pairs.append((bi,c))
        assert refs==([0] if bi==2 else []) and len(b['caps'])==len(refs)
        if bi<2:assert len(one65)==94 and abs(min(one65)-0.14837496801192174)<2e-7 and abs(max(one65)-4.058347906335257)<2e-7
    assert counts==d['counts']=={'new10_norm_equations':288,'new10_norm_branches':576}
    assert d['complete_native_pair_count']==len(pairs)==1 and d['complete_native_pair_digest']==digest.hexdigest()
    assert len(d['coherent_pairs'])==1 and len(d['boundary_ledgers'])==2
    c=d['coherent_pairs'][0];assert (c['prefix'],c['cap'])==(2,0) and c['old30_restored'] and c['old10_restored']
    assert c['all_old_compact_S6_faces_rebuilt'] and c['mandatory_7_33_new10'] and not c['failed_old_S6_corner_restorations']
    assert len(c['whole_face_records'])==9 and [f['wall'] for f in c['whole_face_records']]==[8,11,33,44,55,63,64,67,68]
    assert len(d['required_old_S6_corner_restorations'])==11
    assert d['deletion_score'] is None and len(d['unresolved'])==4
    assert [r['stage'] for r in d['unresolved']]==['delete10_30']+['new30_prefix']*3
    assert all(r['boundary']['warnings'] for r in d['unresolved']) and all(p['score'] is None for p in d['prefixes'])
    for ledger in d['boundary_ledgers']:
        assert not ledger['warnings'] and ledger['PD_status']=='numerically_separated_acute_incidence'
        assert ledger['counts']=={'finite':26,'exterior':6} and ledger['open_ends']==ledger['PD_open_ends']==10
        assert sum(1 if e['kind']=='one_open' else 2 if e['kind']=='two_open' else 0 for e in ledger['edges'])==10
        assert {tuple(r['active']) for r in ledger['rays'] if r['kind']=='finite'}=={tuple(r) for r in ledger['PD_finite_supports']}
    ledger=d['boundary_ledgers'][c['score']]
    for f in c['whole_face_records']:
        rr=[r for r in ledger['rays'] if f['wall'] in r['active']]
        assert f['rays']==rr and f['all_rays_finite']==(bool(rr) and all(r['kind']=='finite' for r in rr))
    print('PASS_JOINT10_30_SAVED',dict(old_supports=199,lost=75,lost_on8=44,recipes=288,branches=576,sole65_per_changed_bank=94,pairs=1,open_ends=10))


if __name__=='__main__':run()
