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
TAG='fresh_joint30_33_hinge_native_pilot_20261018'
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
    return profiles,seen


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    profiles,sideids=exact_interface(pri,d['old_frame_interface'])
    vs=[tuple(v)for v in pri['retained_old_finite_supports']]+[tuple(r['support'])for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost=[i for i,v in enumerate(vs)if set(v)&{30,33}];assert len(lost)==15
    assert all(8 not in vs[i]for i in lost)
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,ii in faces.items():
        if not(set(ii)&set(lost))or len(ii)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i]if j!=k)]==2 for i in ii for k in vs[i]if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i])for i in ii]))
    debit=d['exact_deletion_debit']
    assert debit['lost_finite_supports']==[list(vs[i])for i in lost]
    assert touched==debit['compact_faces_touched']and Counter(r['dimension']for r in touched)=={2:10,3:1}
    assert d['deleted_walls']==[30,33]and d['retained_old_walls']==61
    assert d['deletion_score']is None and len(d['deletion_diagnostic'])==1 and len(d['unresolved'])==13
    for u in d['unresolved']:
        assert u['stage']in('delete30_33','replace33')and u['boundary']['warnings']
        assert u['PD_incidence']['open_ends']==13
        assert any(set(r['active'])>={10,44,55,63,64}and r['kind']=='exterior'for r in u['boundary']['rays'])
    for l in d['boundary_ledgers']:
        assert not l['warnings']and l['PD_status']=='numerically_separated_acute_incidence'
        assert Counter(r['kind']for r in l['rays'])==Counter(l['counts'])
        assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in l['edges'])
        assert {tuple(r['active'])for r in l['rays']if r['kind']=='finite'}=={tuple(r)for r in l['PD_finite_supports']}
    assert len(d['boundary_ledgers'])==5
    assert d['boundary_ledgers'][d['control_score']]['counts']==dict(finite=26,exterior=6)
    counts=Counter();nativekeys=None
    def bank(b,first):
        nonlocal nativekeys
        assert b['complete_low_recipe_count']==len(b['recipes'])
        keys={(r['profile'],*r['labels'])for r in b['recipes']};assert len(keys)==len(b['recipes'])
        refs=[]
        for ri,r in enumerate(b['recipes']):
            counts['norm_equations']+=1
            assert r['index']==ri and 0<=r['profile']<7 and all(m in(2,3,4,5)for m in r['labels'])
            assert r['norm']['status']=='two_real_norm_branches'and sorted(v['branch']for v in r['branches'])==[-1,1]
            for br in r['branches']:
                counts['norm_branches']+=1
                good=not br['blockers']and br['center_product']<-1e-8 and br['norm_residual']<2e-5
                assert all(relation(v)in(None,'parallel_or_unresolved')for _,v in br['blockers'])
                assert good==('cap'in br or br.get('duplicate',False))
                if 'cap'in br:
                    c=b['caps'][br['cap']];refs.append(br['cap'])
                    assert(c['recipe'],c['branch'],c['profile'],c['labels'])==(ri,br['branch'],r['profile'],r['labels'])
        assert sorted(refs)==list(range(len(b['caps'])))
        for ci,c in enumerate(b['caps']):
            assert c['index']==ci and c['center_product']<-1e-8 and c['norm_residual']<2e-5
            vals=dict(c['native_products']);assert len(vals)==len(c['native_products'])==55+(not first)
            if nativekeys is None:nativekeys=set(vals)
            assert set(vals)==nativekeys|(set()if first else{67})and not set(vals)&{30,33}
            assert all(relation(v)not in(None,'parallel_or_unresolved')for v in vals.values())
            assert all(relation(vals[i])==m for i,m in zip(b['frame'],c['labels']))
            v=np.array(c['normal']);a=np.array(profiles[c['profile']]);J=np.array([[2.,1.],[1.,2.]])/3
            assert abs(v@M@v-(2-a@J@a))<2e-5
            if first:assert not c['cone_validated']and c['score']is None and len(c['cone_diagnostic'])==1
            else:assert c['score']is not None
    b=d['first_bank'];assert b['frame']==[7,10,66]and b['duties']==[[0,1],[0,2]]
    bank(b,True);assert len(b['caps'])==12 and b['complete_low_recipe_count']==21
    assert [c['index']for c in b['caps']if c['restores_old33_numerically']]==[11]
    ranked=sorted([c for c in b['caps']if not c['restores_old33_numerically']],key=lambda c:(c['cone_diagnostic'][0]['PD_incidence']['open_ends'],-c['cone_diagnostic'][0]['PD_incidence']['finite'],c['profile'],c['labels'],c['branch']))
    assert d['second_selection']==[c['index']for c in ranked[:8]]+[11]==[r['first_cap']for r in d['records']]
    allcaps=set();rejected={}
    for ri,r in enumerate(d['records']):
        assert r['bank']['frame']==[10,55,67]and r['bank']['duties']==[[0,2]]
        bank(r['bank'],False)
        assert [x['walls']for x in r['frame_native_checks']]==[[10,55],[10,67],[55,67]]
        for x in r['frame_native_checks']:assert abs(x['native_product']-x['reference_product'])<2e-5 and relation(x['native_product'])==relation(x['reference_product'])
        for c in r['bank']['caps']:allcaps.add((ri,c['index']))
    assert counts==Counter(d['counts'])==dict(norm_equations=397,norm_branches=794)
    for r in d['required_corner_failures']:
        key=(r['record'],r['cap']);assert key not in rejected;rejected[key]=r
        assert r['failed_corners']==[[7,55,68]]
        c=d['records'][r['record']]['bank']['caps'][r['cap']]
        assert relation(dict(c['native_products'])[7])=='dotted'
        l=d['boundary_ledgers'][c['score']];assert l['counts']==dict(finite=20,exterior=12)and l['open_ends']==16
        assert any(10 in x['active']and x['kind']=='exterior'for x in l['rays'])
    assert len(rejected)==14 and len(d['accepted'])==1
    a=d['accepted'][0];assert allcaps==set(rejected)|{(a['record'],a['cap'])}and len(allcaps)==15
    assert a['old33_restored']and a['old30_restored']and a['open_ends']==10 and a['finite_rays']==26
    l=d['boundary_ledgers'][a['score']]
    assert a['exterior_supports']==[r['active']for r in l['rays']if r['kind']=='exterior']
    assert a['open_edges']==[e for e in l['edges']if e['kind']in('one_open','two_open')]
    for key,w in [('whole10',10),('whole55',55),('new33_face',67),('new30_face',68)]:
        f=a[key];assert f['wall']==w and f['rays']==[r for r in l['rays']if w in r['active']]
        assert f['all_rays_finite']==all(r['kind']=='finite'for r in f['rays'])
    assert a['whole10']['all_rays_finite']and a['whole55']['all_rays_finite']
    return dict(status='PASS_JOINT30_33_SAVED',actual_side_profiles=len(sideids),old_frame_pairs=6,second_frame_pairs=27,
        exact_debit=[15,10,1],recipes=397,branches=794,first_caps=12,selected=9,final_native_pairs=15,
        closing_corner_failures=14,old_geometry_restorations=1,provisional_exterior_incidence_issues=13,
        scope='No project imports. Independent original-native rational interval identities and saved-profile comparisons;199-support compact debit; saved recipe/branch/cap domain and all native values; selection, all15 final pair ledgers and stored cone/PD equality. No norm-root, classifier, numerical cone or ambient finite-volume reconstruction.')


if __name__=='__main__':print(json.dumps(run(),sort_keys=True))
