#!/usr/bin/env python3
"""Independent saved algebra/domain/debit/ledger checks, without producer imports."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_joint10_30_65_compact_pair_20261018'
M=np.diag([-1.,1.,1.,1.]);J=np.array([[2.,1.],[1.,2.]])/3


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


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='BOUNDED_JOINT10_30_65_COMPACT_PAIR'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    old=json.loads((ROOT/'results/fresh_joint10_30_native_replacement_20261018.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&{10,30,65}};prior={i for i,v in enumerate(vs) if set(v)&{10,30}}
    assert len(lost)==80 and len(prior)==75 and sum(8 in vs[i] for i in lost)==46
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
    assert Counter(r['dimension'] for r in touched)=={2:314,3:325,4:244,5:129,6:46,7:10,8:1}
    assert d['additional_lost_supports']==list(map(list,sorted(vs[i] for i in lost-prior))) and len(d['additional_lost_supports'])==5
    assert d['deleted_walls']==[10,30,65] and d['retained_old_walls']==60
    top=d['topology'];ts=[[8,44,68],[44,63,68],[8,68,69],[63,68,69],[8,52,69],[52,64,69],[63,64,69]]
    assert top['triangles']==ts and top['boundary']==[8,44,63,64,52]
    ee=Counter(tuple(sorted(e)) for t in ts for e in combinations(t,2));vv=set(sum(ts,[]))
    assert (len(vv),len(ee),len(ts))==(7,13,7) and len(vv)-len(ee)+len(ts)==1 and set(ee.values())=={1,2}
    bc=top['boundary'];assert {e for e,n in ee.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert len(top['links'])==7 and {r['vertex'] for r in top['links']}==vv
    for r in top['links']:
        v=r['vertex'];es=[list(i for i in t if i!=v) for t in ts if v in t];assert r['edges']==es
        deg=Counter(i for e in es for i in e);seen={next(iter(deg))}
        while True:
            z=seen|{b for a,b in es if a in seen}|{a for a,b in es if b in seen}
            if z==seen:break
            seen=z
        assert seen==set(deg) and sorted(deg.values())==([2]*len(deg) if v in(68,69) else [1,1]+[2]*(len(deg)-2))
    snap=d['numerical_replay_snapshot'];assert snap==old['numerical_replay_snapshot']
    assert d['old_frame_interface']==old['old_frame_interface'] and d['old8_44_63_frame_checks']==old['changed_frame_native_checks']
    normals={i:np.array(v) for i,v in snap['old_normals']};center=np.array(snap['center']);assert len(normals)==57
    profiles=snap['profiles'];assert len(profiles)==7
    cs=[[0.,0.]]
    for k in (0,1):
        for m in (3,4,5):
            c=[0.,0.];c[k]=-2*np.cos(np.pi/m);cs.append(c)
    offsets=[]
    for pi,p in enumerate(profiles):
        assert max(abs(np.array(p['contacts'])-np.array(cs[pi])))<2e-8
        assert abs(p['square']-(2-np.array(cs[pi])@J@np.array(cs[pi])))<2e-8
        off=dict(p['native_old_offsets']);assert len(off)==57 and set(off)==set(normals);offsets.append(off)
    oldtypes={r['wall']:r['profile'] for r in d['old_frame_interface']['actual_side_profiles']}
    base={i:v for i,v in normals.items() if i not in(10,30,65)};assert len(base)==54
    inventory=[];selected=[]
    for bi in(0,1):
        pool=[(q['blockers'][0][1],r['profile'],r['labels'],q['branch'],r['index'])
              for r in old['new10_banks'][bi]['recipes'] for q in r['branches']
              if len(q['blockers'])==1 and q['blockers'][0][0]==65 and q['center_product']<-1e-8]
        assert len(pool)==94;rows=[]
        for pi in range(7):
            v,_,ms,sign,ri=min(r for r in pool if r[1]==pi)
            z=dict(R_prefix=bi,N_profile=pi,source_recipe=ri,branch=sign,labels=ms,old65_product=v,selected=pi in(0,3,6))
            rows.append(z)
            if z['selected']:selected.append(z)
        inventory.append(dict(R_prefix=bi,sole65_branches=94,profile_minima=rows))
    assert d['selection_inventory']==inventory and d['selected_prefixes']==selected and len(d['records'])==len(selected)==6
    def native(v,pi,nd,types):
        return {i:float(v@M@w)+(offsets[pi][i] if i in normals else np.array(cs[pi])@J@np.array(cs[types[i]])) for i,w in nd.items()}
    def check_values(saved,vals):
        assert len(saved)==len(vals) and set(dict(saved))==set(vals)
        assert max(abs(q-vals[i]) for i,q in saved)<2e-7
    def corner(pa,pb,pc,mab,mac,mbc):
        g=np.eye(5)*2;g[0,1]=g[1,0]=-1
        for k,pi in enumerate((pa,pb,pc),2):g[:2,k]=g[k,:2]=cs[pi]
        for a,b,m in ((2,3,mab),(2,4,mac),(3,4,mbc)):g[a,b]=g[b,a]=-2*np.cos(np.pi/m) if m>2 else 0
        return np.linalg.eigvalsh(g)[0]>1e-8
    counts=Counter();digest=sha256();ncaps=0
    for ri,(rec,sel) in enumerate(zip(d['records'],selected)):
        assert rec['index']==ri and rec['selection']==sel and rec['score'] is None and len(rec['native'])==2
        nd=dict(base);types=dict(oldtypes)
        assert [r['wall'] for r in rec['native']]==[67,68]
        for r in rec['native']:nd[r['wall']]=np.array(r['normal']);types[r['wall']]=r['profile']
        assert types[67]==0 and types[68]==sel['N_profile']
        assert max(abs(nd[67]-np.array(old['prefixes'][sel['R_prefix']]['normal'])))<2e-6
        sr=old['new10_banks'][sel['R_prefix']]['recipes'][sel['source_recipe']]
        sb=next(r for r in sr['branches'] if r['branch']==sel['branch'])
        assert max(abs(nd[68]-np.array(sb['normal'])))<2e-6
        for r in rec['native']:
            pi=r['profile'];v=nd[r['wall']];assert abs(v@M@v-profiles[pi]['square'])<2e-5
            vv=native(v,pi,{i:w for i,w in nd.items() if i!=r['wall']},types);check_values(r['native_products'],vv)
            assert all(relation(q) not in(None,'parallel_or_unresolved') for q in vv.values())
            assert abs(r['center_product']-(v@M@center+profiles[pi]['center_offset']))<2e-7 and r['center_product']<-1e-8
        raw=rec['cone_diagnostic'][0]['boundary'];assert raw['counts']=={'finite':13,'exterior':16} and raw['open_ends']==23 and raw['warnings']
        b=rec['J_bank'];assert b['frame']==[8,63,68] and b['duties']==[[0,2],[1,2]] and len(b['frame_native_checks'])==3
        bonds={}
        for r,(a,c) in zip(b['frame_native_checks'],combinations(range(3),2)):
            i,j=b['frame'][a],b['frame'][c];assert r['walls']==[i,j]
            q=float(nd[i]@M@nd[j])+np.array(cs[types[i]])@J@np.array(cs[types[j]])
            assert abs(q-r['native_product'])<2e-7 and relation(q)==r['relation'];bonds[a,c]=r['relation']
        low=[(pi,list(ms)) for pi in range(7) for ms in product(range(2,6),repeat=3)
            if all(corner(types[b['frame'][a]],types[b['frame'][c]],pi,bonds[a,c],ms[a],ms[c]) for a,c in b['duties'])]
        assert sorted(low)==sorted((r['profile'],r['labels']) for r in b['recipes'])
        assert len(low)==b['complete_low_recipe_count']==[96,17,17,96,17,17][ri]
        assert b['symbolic_high_families_not_tested']==[10,9,9,10,9,9][ri]
        refs=[]
        for ki,r in enumerate(b['recipes']):
            assert r['index']==ki and r['norm']['status']=='two_real_norm_branches' and sorted(q['branch'] for q in r['branches'])==[-1,1]
            counts['J_equations']+=1
            for q in r['branches']:
                counts['J_branches']+=1;pi=r['profile'];v=np.array(q['normal'])
                assert abs(v@M@v-profiles[pi]['square'])<2e-5 and q['norm_residual']<2e-5
                vv=native(v,pi,nd,types);check_values(q['native_products'],vv)
                bad=[[i,x] for i,x in q['native_products'] if relation(x) in(None,'parallel_or_unresolved')];assert bad==q['blockers']
                assert all(relation(vv[i])==m for i,m in zip(b['frame'],r['labels']))
                assert abs(q['center_product']-(v@M@center+profiles[pi]['center_offset']))<2e-7
                good=not bad and q['center_product']<-1e-8
                assert good==('cap' in q or q.get('duplicate',False))
                if 'cap' not in q:continue
                refs.append(q['cap']);c=b['caps'][q['cap']];assert c['recipe']==ki and c['branch']==q['branch'] and c['normal']==q['normal']
                assert c['native_products']==q['native_products'] and c['score'] is None
                assert c['failed_disk_corners']==[[8,52,69],[52,64,69],[63,64,69]]
                assert relation(vv[52])==relation(vv[64])=='dotted'
                assert len(c['cone_diagnostic'])==1;diag=c['cone_diagnostic'][0];bd=diag['boundary'];pd=diag['PD_incidence']
                assert bd['counts']=={'finite':15,'exterior':16} and bd['open_ends']==23
                assert {tuple(w['sides']) for w in bd['warnings']}=={(8,33),(11,44)}
                assert pd['status']=='numerically_separated_acute_incidence' and pd['finite']==15 and pd['open_ends']==23
                assert {tuple(r['active']) for r in bd['rays'] if r['kind']=='finite'}=={tuple(r) for r in pd['finite_supports']}
                assert len(c['raw_whole_face_records'])==10 and [f['wall'] for f in c['raw_whole_face_records']]==[8,11,33,44,55,63,64,67,68,69]
                for f in c['raw_whole_face_records']:
                    rr=[r for r in bd['rays'] if f['wall'] in r['active']]
                    assert f['rays']==rr and f['all_rays_finite']==(bool(rr) and all(r['kind']=='finite' for r in rr))
                    if f['wall']==68:assert {tuple(sorted(r['active'])) for r in rr}=={tuple(sorted(t)) for t in ts if 68 in t} and f['all_rays_finite']
                    else:assert not f['all_rays_finite']
                assert c['raw_whole_N'] and not c['raw_whole_J']
                digest.update(json.dumps([ri,c['index'],pi,c['labels'],c['branch'],c['native_products']],separators=(',',':')).encode());ncaps+=1
        assert refs==list(range(len(b['caps']))) and len(refs)==[78,15,15,78,15,15][ri]
    assert counts==d['counts']=={'J_equations':260,'J_branches':520}
    assert ncaps==d['complete_native_triple_count']==216 and digest.hexdigest()==d['complete_native_triple_digest'] and d['coherent_triples']==[]
    ctrl=d['old65_control'];assert ctrl['frame']==[8,63,64] and ctrl['labels']==[2,2,2] and len(ctrl['branches'])==2
    accepted=0
    for c in ctrl['branches']:
        pi=c['profile'];v=np.array(c['normal']);assert pi==oldtypes[65] and abs(v@M@v-profiles[pi]['square'])<2e-5
        vv=native(v,pi,{i:w for i,w in normals.items() if i!=65},oldtypes);check_values(c['native_products'],vv)
        assert [[i,x] for i,x in c['native_products'] if relation(x) in(None,'parallel_or_unresolved')]==c['blockers']
        good=not c['blockers'] and c['center_product']<-1e-8;assert good==(c['score'] is not None)
        if good:assert max(abs(v-normals[65]))<2e-6 and relation(vv[44])==15;accepted+=1
    assert accepted==1 and len(d['boundary_ledgers'])==2
    for bd in d['boundary_ledgers']:
        assert not bd['warnings'] and bd['counts']=={'finite':26,'exterior':6} and bd['open_ends']==bd['PD_open_ends']==10
        assert {tuple(r['active']) for r in bd['rays'] if r['kind']=='finite'}=={tuple(r) for r in bd['PD_finite_supports']}
    assert len(d['unresolved'])==223 and Counter(r['stage'] for r in d['unresolved'])=={'delete10_30_65':1,'selected_prefix':6,'complete_J':216}
    assert all(r['boundary']['warnings'] for r in d['unresolved'])
    print('PASS_JOINT10_30_65_SAVED',dict(old_supports=199,debit=80,increment=5,prefixes=6,recipes=260,branches=520,native_triples=216,whole_N=216,whole_J=0,clean_triples=0))


if __name__=='__main__':run()
