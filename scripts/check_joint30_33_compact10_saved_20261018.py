#!/usr/bin/env python3
"""Saved native-interface, domain and boundary accounting; no project imports."""
from collections import Counter,defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_joint30_33_compact10_first_20261018'
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
    return int(m)if abs(x+theta(m))<=min(2e-7,spacing/100)else None


def exact_interface(pri,saved):
    dec=lambda a:tuple(map(Q,a))
    def plus(a,b):return a[0]+b[0],a[1]+b[1]
    def times(a,b):
        x=[a[0]*b[0],a[0]*b[1],a[1]*b[0],a[1]*b[1]]
        return min(x),max(x)
    h=list(map(dec,pri['common_center']));G=[[dec(v)for v in r]for r in pri['computational_gram']]
    normals={};gates={r['wall']:r for r in pri['native_gates']}
    for f in pri['facet_witnesses']:
        hp=dec(pri['center_products'][str(f['wall'])]);assert hp[1]<0
        inverse=(1/hp[1],1/hp[0]);v=[]
        for a,b in zip(h,map(dec,f['point'])):
            v.append(times((2*(a[0]-b[1]),2*(a[1]-b[0])),inverse))
        normals[f['wall']]=v
    def pair(i,j):
        if i==j:return Q(2),Q(2)
        if 66 in(i,j):
            r=gates[j if i==66 else i]
            return(Q(0),Q(0))if r['label']==2 else dec(r['product'])
        out=(Q(0),Q(0))
        for a,x in enumerate(normals[i]):
            for b,y in enumerate(normals[j]):out=plus(out,times(times(x,G[a][b]),y))
        return out
    profiles=[[0.,0.]]
    for k in(0,1):
        for m in(3,4,5):
            a=[0.,0.];a[k]=-2*np.cos(np.pi/m);profiles.append(a)
    J=np.array([[2.,1.],[1.,2.]])/3;seen=set()
    for row in saved['actual_side_profiles']:
        i=row['wall'];assert i not in seen;seen.add(i)
        p=np.array(profiles[row['profile']]);assert max(abs(p-np.array(row['stored_A2_contacts'])))<2e-5
        assert abs(row['intrinsic_square']-(2-p@J@p))<2e-5
        assert [r['white']for r in row['original_native_checks']]==[0,1,2,3,5,6]
        for r in row['original_native_checks']:
            q=pair(i,r['white']);assert tuple(map(str,q))==tuple(r['exact_native_enclosure'])
            assert float(q[0])-2e-5<=r['numerical_product']<=float(q[1])+2e-5
    assert {r['wall']:r['profile']for r in saved['actual_side_profiles']}[10]==1
    assert [r['walls']for r in saved['all_old_frame_pairs']]==[list(x)for x in combinations([7,10,55,66],2)]
    for r in saved['all_old_frame_pairs']:
        q=pair(*r['walls']);assert tuple(map(str,q))==tuple(r['exact_native_enclosure'])
        assert float(q[0])-2e-5<=r['native_product']<=float(q[1])+2e-5
        assert relation(r['native_product'])==r['native_relation']
    return profiles,seen,pair


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    profiles,sideids,native_interval=exact_interface(pri,d['old_frame_interface'])
    assert len(sideids)==18 and [r['walls']for r in d['added_old44_frame_checks']]==[[7,44],[10,44],[44,55],[44,66]]
    for r in d['added_old44_frame_checks']:
        q=native_interval(*r['walls']);assert tuple(map(str,q))==tuple(r['exact_native_enclosure'])
        assert float(q[0])-2e-5<=r['native_product']<=float(q[1])+2e-5 and relation(r['native_product'])==r['native_relation']
    vs=[tuple(v)for v in pri['retained_old_finite_supports']]+[tuple(r['support'])for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost=[i for i,v in enumerate(vs)if set(v)&{30,33}];assert len(lost)==15 and all(8 not in vs[i]for i in lost)
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,ii in faces.items():
        if not(set(ii)&set(lost))or len(ii)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i]if j!=k)]==2 for i in ii for k in vs[i]if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i])for i in ii]))
    assert d['exact_deletion_debit']['lost_finite_supports']==[list(vs[i])for i in lost]
    assert touched==d['exact_deletion_debit']['compact_faces_touched']and Counter(r['dimension']for r in touched)=={2:10,3:1}
    assert d['deleted_walls']==[30,33]and d['retained_old_walls']==61 and d['deletion_score']is None
    assert len(d['unresolved'])==13
    for r in d['unresolved']:
        assert r['stage']in('delete30_33','new33')and r['boundary']['warnings']
        assert any(set(a['active'])>={10,44,55,63,64}and a['kind']=='exterior'for a in r['boundary']['rays'])
    for l in d['boundary_ledgers']:
        assert not l['warnings']and l['PD_status']=='numerically_separated_acute_incidence'
        assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
        assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
        assert {tuple(r['active'])for r in l['rays']if r['kind']=='finite'}=={tuple(r)for r in l['PD_finite_supports']}
    assert len(d['boundary_ledgers'])==4 and d['boundary_ledgers'][d['control_score']]['open_ends']==10
    counts=Counter();oldkeys=None
    expected=[dict(frame=[10,44,7],duties=[[0,1]],extra_corners=[[7,55,67]],new_wall=67,replaced_old_wall=30),
              dict(frame=[7,10,66],duties=[[0,1],[0,2]],extra_corners=[],new_wall=68,replaced_old_wall=33)]
    assert len(d['banks'])==2 and d['low_labels']==[2,3,4,5]
    for bi,(b,spec)in enumerate(zip(d['banks'],expected)):
        assert all(b[k]==v for k,v in spec.items())
        assert b['complete_low_recipe_count']==len(b['recipes'])==[64,21][bi]
        assert len({(r['profile'],*r['labels'])for r in b['recipes']})==len(b['recipes'])
        assert len(b['caps'])==[1,12][bi]
        assert [r['walls']for r in b['frame_pairs']]==[list(v)for v in combinations(b['frame'],2)]
        for r in b['frame_pairs']:
            lo,hi=native_interval(*r['walls']);assert float(lo)-2e-5<=r['native_product']<=float(hi)+2e-5
            assert relation(r['native_product'])==r['relation']
        refs=[]
        for ri,r in enumerate(b['recipes']):
            counts['norm_equations']+=1
            assert r['index']==ri and 0<=r['profile']<7 and all(v in(2,3,4,5)for v in r['labels'])
            assert r['norm']['status']=='two_real_norm_branches'and sorted(v['branch']for v in r['branches'])==[-1,1]
            for br in r['branches']:
                counts['norm_branches']+=1
                good=not br['blockers']and br['center_product']<-1e-8 and br['norm_residual']<2e-5 and not br.get('extra_corner_failures')
                assert good==('cap'in br or br.get('duplicate',False))
                assert all(relation(v)in(None,'parallel_or_unresolved')for _,v in br['blockers'])
                if 'cap'in br:
                    c=b['caps'][br['cap']];refs.append(br['cap'])
                    assert(c['recipe'],c['branch'],c['profile'],c['labels'])==(ri,br['branch'],r['profile'],r['labels'])
        assert sorted(refs)==list(range(len(b['caps'])))
        for ci,c in enumerate(b['caps']):
            assert c['index']==ci and c['center_product']<-1e-8 and c['norm_residual']<2e-5
            vals=dict(c['native_products']);assert len(vals)==len(c['native_products'])==55
            if oldkeys is None:oldkeys=set(vals)
            assert set(vals)==oldkeys and not set(vals)&{30,33}
            assert all(relation(v)not in(None,'parallel_or_unresolved')for v in vals.values())
            assert all(relation(vals[i])==m for i,m in zip(b['frame'],c['labels']))
            a=np.array(profiles[c['profile']]);J=np.array([[2.,1.],[1.,2.]])/3;v=np.array(c['normal'])
            assert abs(v@M@v-(2-a@J@a))<2e-5
            assert c['cone_validated']==(bi==0)and(c['score']is not None)==(bi==0)
            if bi==1:assert len(c['cone_diagnostic'])==1
    assert counts==Counter(d['counts'])==dict(norm_equations=85,norm_branches=170)
    assert [c['index']for c in d['banks'][0]['caps']if c['restores_old_numerically']]==[0]
    assert [c['index']for c in d['banks'][1]['caps']if c['restores_old_numerically']]==[11]
    assert d['banks'][0]['caps'][0]['labels']==[2,2,5]
    r=d['banks'][0]['recipes'][3];assert r['labels']==[2,2,5]and r['profile']==0
    br=next(b for b in r['branches']if b['branch']==1);assert len(br['blockers'])==1 and br['blockers'][0][0]==8 and abs(br['blockers'][0][1]-2)<2e-7
    for br in r['branches']:
        match=br['retained_control_comparison'];assert match['old_wall']==(8 if br['branch']==1 else 30)and match['profile']==0
        difference=float(np.max(abs(np.array(match['normal'])-np.array(match['old_normal']))))
        assert difference==match['maximum_coordinate_difference']and difference<2e-6
        assert abs(br['center_product']-match['old_center_product'])<2e-6 and match['old_center_product']<0
    digest=sha256();hist=Counter();coherent={r['pair']:r for r in d['coherent_pairs']}
    assert d['complete_pair_count']==len(d['pairs'])==12 and len(coherent)==len(d['coherent_pairs'])==3
    assert d['required_corners']==[[10,44,67],[7,55,67],[7,10,68],[7,66,68],[10,67,68]]
    for ix,x in enumerate(d['banks'][0]['caps']):
        for iy,y in enumerate(d['banks'][1]['caps']):
            i=ix*12+iy;r=d['pairs'][i]
            J=np.array([[2.,1.],[1.,2.]])/3
            val=float(np.array(x['normal'])@M@np.array(y['normal']))+float(np.array(profiles[x['profile']])@J@np.array(profiles[y['profile']]))
            rel=relation(val);hist[str(rel)]+=1
            assert(r['index'],r['new30'],r['new33'],r['native_product'],r['relation'])==(i,ix,iy,val,rel)
            digest.update(json.dumps([ix,iy,val,rel],separators=(',',':')).encode())
            if rel is None:assert i not in coherent and r['score']is None;continue
            assert i in coherent;a=coherent[i];assert a['score']==r['score']
            l=d['boundary_ledgers'][a['score']]
            assert a['old30_restored']and a['old33_restored']==(iy==11)
            assert a['exterior_supports']==[z['active']for z in l['rays']if z['kind']=='exterior']
            assert a['open_edges']==[e for e in l['edges']if e['kind']in('one_open','two_open')]
            assert len(a['whole_face_records'])==4
            for f,w in zip(a['whole_face_records'],[10,55,67,68]):
                assert f['wall']==w and f['rays']==[z for z in l['rays']if w in z['active']]
                assert f['all_rays_finite']==all(z['kind']=='finite'for z in f['rays'])
            assert a['compact10_preserved']==a['whole_face_records'][0]['all_rays_finite']==(iy==11)
            assert a['compact55_preserved']==a['whole_face_records'][1]['all_rays_finite']==True
            assert a['finite_rays']==l['counts']['finite']==(26 if iy==11 else 25)
            assert a['open_ends']==l['open_ends']==(10 if iy==11 else 11)
            assert a['failed_corners']==([]if iy==11 else[[10,67,68]])
            if iy!=11:assert[10,67,68]in a['exterior_supports']
    assert hist==Counter(d['pair_histogram'])=={'None':9,'dotted':2,'2':1}and digest.hexdigest()==d['complete_pair_digest']
    return dict(status='PASS_COMPACT10_FIRST_SAVED',exact_debit=[15,10,1],original_side_profiles=18,old_frame_pairs=10,
        recipes=85,branches=170,banks=[1,12],complete_pairs=12,pair_hashes=1,coherent_pairs=3,new11_end_pairs=2,old_control=1,
        scope='No project imports. Independent saved original-native interval identities,199-support/compact debit, domain/branch/native values,all12 mutual products/hash, all3 final boundary and whole-face records, stored cone/PD agreement. No classifier, norm-root, cone or ambient-volume replay.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))

