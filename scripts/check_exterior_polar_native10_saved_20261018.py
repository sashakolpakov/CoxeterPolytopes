#!/usr/bin/env python3
"""Independent saved numerical interfaces and ledgers; no project imports."""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_exterior_polar_native10_frame_20261018'
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


def native_interval(p):
    dec=lambda x:tuple(map(F,x))
    def mul(a,b):
        z=[x*y for x in a for y in b];return min(z),max(z)
    h=list(map(dec,p['common_center']));g=[[dec(x) for x in r] for r in p['computational_gram']];ns={}
    for row in p['facet_witnesses']:
        hp=dec(p['center_products'][str(row['wall'])]);assert hp[1]<0
        ns[row['wall']]=[mul((2*(a[0]-b[1]),2*(a[1]-b[0])),(1/hp[1],1/hp[0])) for a,b in zip(h,map(dec,row['point']))]
    gates={r['wall']:r for r in p['native_gates']}
    def pair(i,j):
        if i==j:return F(2),F(2)
        if 66 in (i,j):
            gate=gates[j if i==66 else i]
            return (F(),F()) if gate['label']==2 else dec(gate['product'])
        q=(F(),F())
        for a,x in enumerate(ns[i]):
            for b,y in enumerate(ns[j]):
                z=mul(mul(x,g[a][b]),y);q=q[0]+z[0],q[1]+z[1]
        return q
    return pair


def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text())
    assert d['status']=='BOUNDED_EXTERIOR_POLAR_AND_NATIVE10_FRAME'
    for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    pair=native_interval(pri);snap=d['numerical_snapshot'];old={int(k):np.array(v) for k,v in snap['old_normals'].items()}
    center=np.array(snap['center']);pp=snap['profiles'];assert len(pp)==7 and len(old)==57
    profiles=[np.zeros(2)]
    for node in (0,1):
        for label in (3,4,5):
            v=np.zeros(2);v[node]=-2*np.cos(np.pi/label);profiles.append(v)
    J=np.array([[2.,1.],[1.,2.]])/3
    off=lambda a,b:float(profiles[a]@J@profiles[b])
    for pi,p in enumerate(pp):
        assert max(abs(np.array(p['contacts'])-profiles[pi]))<2e-7
        assert abs(p['square']-(2-off(pi,pi)))<2e-7
        assert len(p['native_old_offsets'])==len(old) and set(dict(p['native_old_offsets']))==set(old)
    types={r['wall']:r['profile'] for r in d['old_frame_interface']['actual_side_profiles']}
    assert len(types)==18 and types[7]==4 and types[10]==1
    for r in d['old_frame_interface']['actual_side_profiles']:
        assert np.max(abs(profiles[r['profile']]-r['stored_A2_contacts']))<2e-5
        assert abs(r['intrinsic_square']-pp[r['profile']]['square'])<2e-5
        for z in r['original_native_checks']:
            q=pair(r['wall'],z['white']);assert list(map(str,q))==z['exact_native_enclosure']
            assert float(q[0])-2e-5<=z['numerical_product']<=float(q[1])+2e-5
    expected=set(tuple(sorted(e)) for e in list(combinations([7,55,10],2))+list(combinations([8,11,33,44],2))+[(10,i) for i in [7,8,11,33,44]])
    assert [r['walls'] for r in d['original_native_checks']]==list(map(list,sorted(expected)))
    for r in d['old_frame_interface']['all_old_frame_pairs']+d['original_native_checks']:
        i,j=r['walls'];q=pair(i,j);assert list(map(str,q))==r['exact_native_enclosure']
        value=float(old[i]@M@old[j])+off(types[i],types[j])
        assert abs(value-r['native_product'])<2e-7 and float(q[0])-2e-5<=value<=float(q[1])+2e-5
        assert relation(value)==r.get('relation',r.get('native_relation'))
    # Independent old certificate support/compact-face debit.
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
    assert debit['lost_finite_supports']==[list(vs[i]) for i in sorted(lost)] and debit['compact_faces_touched']==touched
    assert Counter(r['dimension'] for r in touched)=={2:10,3:1}
    assert d['deleted_walls']==[30] and d['retained_old_walls']==62
    def native(v,pi,keys):
        offsets=dict(pp[pi]['native_old_offsets'])
        return {i:float(v@M@old[i])+offsets[i] for i in keys}
    def check_values(v,pi,row,keys):
        values=native(v,pi,keys);assert len(row['native_products'])==len(keys) and set(dict(row['native_products']))==set(keys)
        assert max(abs(values[i]-x) for i,x in row['native_products'])<2e-7
        assert row['blockers']==[[i,values[i]] for i in keys if relation(values[i]) in (None,'parallel_or_unresolved')]
        assert abs(row['center_product']-(float(v@M@center)+pp[pi]['center_offset']))<2e-7
        return values
    pol=d['polar'];assert pol['anchors']==[8,11,33,44] and pol['spanning_triple']==[8,11,44]
    gram=np.array([[old[i]@M@old[j] for j in pol['spanning_triple']] for i in pol['spanning_triple']])
    assert np.max(abs(gram-pol['spanning_gram']))<2e-7 and np.linalg.det(gram)<-1
    assert max(abs(float(old[10]@M@old[i])) for i in pol['anchors'])<2e-7
    assert abs(float(old[10]@M@old[10])-4/3)<2e-7 and abs(float(old[10]@M@old[7])+1/3)<2e-7
    assert len(pol['rows'])==14 and [(r['profile'],r['sign']) for r in pol['rows']]==list(product(range(7),(-1,1)))
    polar_keys=[i for i in old if i not in (10,30,65)];accepted=[]
    for row in pol['rows']:
        pi,sign=row['profile'],row['sign'];scale=np.sqrt(3*pp[pi]['square']/4);v=sign*scale*old[10]
        assert np.max(abs(v-row['normal']))<2e-7 and abs(scale-row['scale'])<2e-7
        values=check_values(v,pi,row,polar_keys);form=off(pi,4)-sign*scale/3
        assert abs(form-row['native7_formula'])<2e-7 and abs(values[7]-row['native7'])<2e-7
        if sign==-1:assert row['center_product']>0
        elif pi==0:assert -1<values[7]<0 and relation(values[7]) is None
        elif pi!=1:assert values[7]>0
        if not row['blockers'] and row['center_product']<-1e-8:accepted.append((pi,sign))
    assert accepted==[(1,1)]
    seed=json.loads((ROOT/'results/fresh_joint10_30_65_compact_pair_20261018.json').read_text())
    assert len(pol['changed_R_conflicts'])==2
    for r in pol['changed_R_conflicts']:
        rec=next(x for x in seed['records'] if x['selection']['R_prefix']==r['R_prefix'])['native'][0]
        val=float(old[10]@M@np.array(rec['normal']))+off(1,rec['profile'])
        assert abs(val-r['native10_product'])<2e-7 and -1<val<0
    # Independently exhaust the declared low cube using its finite5x5 Gram.
    b=d['bank'];assert b['frame']==[7,55,10] and b['duties']==[[0,1]] and b['low_labels']==[2,3,4,5]
    domain=[]
    for pi in range(7):
        for ms in product(range(2,6),repeat=3):
            gram=2*np.eye(5);gram[0,1]=gram[1,0]=-1
            for k,p in enumerate((4,0,pi)):
                gram[k+2,:2]=profiles[p];gram[:2,k+2]=profiles[p]
            gram[2,4]=gram[4,2]=-2*np.cos(np.pi/ms[0]);gram[3,4]=gram[4,3]=-2*np.cos(np.pi/ms[1])
            if np.linalg.eigvalsh(gram).min()>1e-8:domain.append((pi,list(ms)))
    assert [(r['profile'],r['labels']) for r in b['recipes']]==domain and len(domain)==64
    assert b['symbolic_high_families_not_tested']==21
    keys=[i for i in old if i!=30];refs=[];digest=sha256();hist=Counter();branches=0
    for ri,r in enumerate(b['recipes']):
        assert r['index']==ri and r['norm']['status']=='two_real_norm_branches'
        assert sorted(x['branch'] for x in r['branches'])==[-1,1]
        for br in r['branches']:
            branches+=1;pi=r['profile'];v=np.array(br['normal']);vals=check_values(v,pi,br,keys)
            assert abs(v@M@v-pp[pi]['square'])<2e-5 and br['norm_residual']<2e-5
            assert all(relation(vals[i])==m for i,m in zip(b['frame'],r['labels']))
            digest.update(json.dumps([ri,br['branch'],pi,r['labels'],br['native_products'],br['center_product']],separators=(',',':')).encode())
            hist[str(tuple(i for i,_ in br['blockers']))]+=1
            good=not br['blockers'] and br['center_product']<-1e-8
            assert good==('cap' in br or br.get('duplicate',False))
            if 'cap' in br:refs.append(br['cap'])
    assert branches==128 and refs==list(range(35)) and digest.hexdigest()==b['complete_branch_digest'] and dict(hist)==b['branch_histogram']
    assert len(b['caps'])==35;clean=Counter();provisional=[]
    for ci,c in enumerate(b['caps']):
        assert c['index']==ci
        r=b['recipes'][c['recipe']];br=next(x for x in r['branches'] if x['branch']==c['branch'])
        assert br['cap']==ci and c['profile']==r['profile'] and c['labels']==r['labels'] and c['normal']==br['normal']
        if c['score'] is None:
            assert len(c['cone_diagnostic'])==1;u=c['cone_diagnostic'][0];raw=u['boundary'];provisional.append(ci)
            assert u['PD_incidence']['status']=='unresolved_possible_affine_supports' and u['PD_incidence']['ambiguous_supports']==[[10,67]]
        else:raw=d['boundary_ledgers'][c['score']];clean[c['open_ends']]+=1
        assert len(c['whole_face_records'])==11
        for f,wall in zip(c['whole_face_records'],(8,10,11,33,44,55,63,64,65,66,67)):
            rays=[r for r in raw['rays'] if wall in r['active']]
            assert f==dict(wall=wall,rays=rays,all_rays_finite=bool(rays) and all(r['kind']=='finite' for r in rays))
        assert c['finite_rays']==raw['counts']['finite'] and c['open_ends']==raw['open_ends']
        assert c['compact8'] and c['compact10']==c['restores_old30']
        assert c['exterior_supports']==[r['active'] for r in raw['rays'] if r['kind']=='exterior']
        assert c['open_edges']==[e for e in raw['edges'] if e['kind'] in ('one_open','two_open')]
    assert clean=={10:1,14:19,16:12} and provisional==[22,29,32]
    assert len(d['boundary_ledgers'])==5 and len(d['unresolved'])==4 and d['deletion_score'] is None
    for l in d['boundary_ledgers']:
        assert not l['warnings'] and l['PD_status']=='numerically_separated_acute_incidence'
        assert l['open_ends']==l['PD_open_ends']==sum(1 if e['kind']=='one_open' else 2 if e['kind']=='two_open' else 0 for e in l['edges'])
        assert {tuple(r['active']) for r in l['rays'] if r['kind']=='finite'}=={tuple(r) for r in l['PD_finite_supports']}
    print('PASS_EXTERIOR_POLAR_NATIVE10_SAVED',dict(old_supports=199,debit=[13,10,1],polar_profiles=7,recipes=64,branches=128,caps=35,clean=32,provisional=3,ledgers=5))


if __name__=='__main__':run()
