#!/usr/bin/env python3
"""Independent numerical saved-domain checks; NumPy, no project imports.

Checks all new physical branches, native gates, finite supports and saved
cone accounting. Does not reconstruct cones or the historical scaffold.
"""
import json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_x_s_diagonal_native_pilot_20261018'
PREV='fresh_oblique_high_n_right_r_join_20261018'
M=np.diag([-1.,1.,1.,1.])
TRI=[[8,10,67],[8,68,69],[8,67,69],[8,52,68],[52,68,72],
     [10,67,73],[67,69,73],[69,70,73],[68,69,72],[69,70,72]]
BOUNDARY=[8,10,73,70,72,52]

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

def topology(t):
    assert t['triangles']==TRI and t['boundary']==BOUNDARY
    vs=sorted(set(sum(TRI,[])));ee=Counter(tuple(sorted(e))for f in TRI for e in combinations(f,2))
    assert vs==t['vertices'] and(len(vs),len(ee),len(TRI))==(t['V'],t['E'],t['F'])==(9,18,10)
    assert t['euler']==1 and 71 not in vs and(68,70)not in ee and(69,72)in ee
    assert {e for e,n in ee.items()if n==1}=={tuple(sorted(e))for e in zip(BOUNDARY,BOUNDARY[1:]+BOUNDARY[:1])}
    assert set(ee.values())=={1,2} and set(t['links'])==set(map(str,vs))
    for v in vs:
        edges=[sorted(set(f)-{v})for f in TRI if v in f];assert t['links'][str(v)]==edges
        degrees=Counter(x for e in edges for x in e);seen={next(iter(degrees))}
        while True:
            new=seen|{x for e in edges if seen.intersection(e)for x in e}
            if new==seen:break
            seen=new
        assert seen==set(degrees)
        assert sorted(degrees.values())==([1,1]+[2]*(len(degrees)-2)if v in BOUNDARY else[2]*len(degrees))

def debit(d):
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    vs=[tuple(v)for v in pri['retained_old_finite_supports']]+[tuple(r['support'])for r in pri['new_finite_vertices']]
    assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs)if set(v)&set(d['deleted_original_walls'])}
    assert len(lost)==19 and sum(8 in vs[i]for i in lost)==3
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,ii in faces.items():
        if not(lost&set(ii))or len(ii)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i]if j!=k)]==2 for i in ii for k in vs[i]if k not in f):
            touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i])for i in ii]))
    z=d['exact_deletion_debit']
    assert z['lost_finite_supports']==[list(vs[i])for i in sorted(lost)] and z['compact_faces_touched']==touched
    assert Counter(r['dimension']for r in touched)=={2:12,3:1}

class Check:
    def __init__(self,d):
        self.d=d;s=d['numerical_snapshot'];assert s['metric']==M.tolist()
        self.old={int(k):np.array(v)for k,v in s['old_normals'].items()};self.center=np.array(s['center'])
        self.pp=s['profiles'];self.labels=s['labels'];assert len(self.old)==57 and set(self.old)==set(self.labels)
        self.base={i:v for i,v in self.old.items()if i not in d['deleted_original_walls']}
        self.ids=[i for i in s['actual_sides']if i not in d['deleted_original_walls']]
        assert len(self.base)==52 and len(self.ids)==13
        self.types={r['wall']:r['profile']for r in d['old_frame_interface']['actual_side_profiles']}
        self.profiles=[np.zeros(2)]
        for i in(0,1):
            for m in(3,4,5):
                v=np.zeros(2);v[i]=-2*np.cos(np.pi/m);self.profiles.append(v)
        self.J=np.array([[2.,1.],[1.,2.]])/3
        for pi,p in enumerate(self.pp):
            assert np.max(abs(self.profiles[pi]-p['contacts']))<2e-7
            assert abs(p['square']-(2-self.off(pi,pi)))<2e-7
        for r in d['old_frame_interface']['actual_side_profiles']:
            assert np.max(abs(self.profiles[r['profile']]-r['stored_A2_contacts']))<2e-5
        self.branches=0;self.recipes=0
    def off(self,a,b):return float(self.profiles[a]@self.J@self.profiles[b])
    def pt(self,i,tt):return tt.get(i,self.types.get(i,0))
    def native(self,v,pi,nd,tt):
        ofs=dict(zip(self.labels,self.pp[pi]['offsets']))
        vals={i:float(v@M@w)+(ofs[i]if i in ofs else self.off(pi,tt[i]))for i,w in nd.items()}
        bad=[i for i in nd if relation(vals[i])in(None,'parallel_or_unresolved')]
        return vals,bad,float(v@M@self.center)+self.pp[pi]['center_offset']
    def gates(self,rows,nd,tt):
        assert [r['wall']for r in rows]==list(tt)
        for r in rows:
            i=r['wall'];vals,bad,cp=self.native(nd[i],tt[i],{j:v for j,v in nd.items()if j!=i},tt)
            assert [j for j,_ in r['blockers']]==bad and all(abs(vals[j]-q)<2e-7 for j,q in r['blockers'])
            assert abs(r['center_product']-cp)<2e-7
        return all(not r['blockers']and r['center_product']<-1e-8 for r in rows)
    def bank(self,r,nd,tt,w,frame,ms):
        assert(r['wall'],r['profile'],r['frame'],r['labels'])==(w,0,frame,ms)
        assert r['norm']['status']=='two_real_norm_branches'
        assert sorted(b['branch']for b in r['branches'])==[-1,1]and len(r['branches'])==2
        g=np.array([[float(nd[i]@M@nd[j])for j in frame]for i in frame])
        assert abs(np.linalg.det(g))>1e-7 and np.max(abs(g-np.array(r['norm']['gram'])))<2e-7
        roots=[np.array(b['normal'])for b in r['branches']];assert np.linalg.norm(roots[0]-roots[1])>1e-6
        for b,v in zip(r['branches'],roots):
            vals,bad,cp=self.native(v,0,nd,tt)
            assert len(b['native_products'])==len(nd)and set(dict(b['native_products']))==set(nd)
            assert all(abs(vals[i]-q)<2e-7 for i,q in b['native_products'])
            assert [i for i,_ in b['blockers']]==bad and all(abs(vals[i]-q)<2e-7 for i,q in b['blockers'])
            err=max(abs(float(v@M@v)-2),max(abs(float(v@M@nd[i])+2*np.cos(np.pi/m))for i,m in zip(frame,ms)))
            assert err<2e-5 and abs(err-b['norm_residual'])<2e-7 and abs(cp-b['center_product'])<2e-7
            assert bool(not bad and cp<-1e-8 and err<2e-5)==b['admissible']
        self.branches+=2;self.recipes+=1
    def finite(self,face,nd,tt):
        g=2*np.eye(5);g[0,1]=g[1,0]=-1
        for k,w in enumerate(face):g[k+2,:2]=self.profiles[self.pt(w,tt)];g[:2,k+2]=g[k+2,:2]
        for i,j in combinations(range(3),2):
            a,b=face[i],face[j];q=float(nd[a]@M@nd[b])+self.off(self.pt(a,tt),self.pt(b,tt))
            assert isinstance(relation(q),int);g[i+2,j+2]=g[j+2,i+2]=q
        return np.linalg.eigvalsh(g).min()>1e-8
    def ledger(self,a,nd,tt,ids):
        assert a['score']is not None and a['cone_diagnostic']is None
        bd=self.d['boundary_ledgers'][a['score']]
        assert not bd['warnings']and bd['PD_status']=='numerically_separated_acute_incidence'
        assert a['finite_rays']==bd['counts']['finite']and a['open_ends']==bd['open_ends']==bd['PD_open_ends']
        finite={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
        assert finite=={tuple(sorted(r))for r in bd['PD_finite_supports']}
        g=np.array([[float(nd[i]@M@nd[j])for j in ids]for i in ids])
        found={tuple(sorted(ids[i]for i in ix))for ix in combinations(range(len(ids)),3)if np.linalg.eigvalsh(g[np.ix_(ix,ix)]).min()>1e-8}
        assert finite==found
        assert a['exterior_supports']==[r['active']for r in bd['rays']if r['kind']=='exterior']
        assert a['open_edges']==[r for r in bd['edges']if r['kind']in('one_open','two_open')]
        assert len(a['whole_faces'])==len(ids)and{r['wall']for r in a['whole_faces']}==set(ids)
        for f in a['whole_faces']:
            rr=[r for r in bd['rays']if f['wall']in r['active']]
            assert f==dict(wall=f['wall'],rays=rr,all_rays_finite=bool(rr)and all(r['kind']=='finite'for r in rr))
        return bd

def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text());old=json.loads((ROOT/f'results/{PREV}.json').read_text())
    for doc in(d,old):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert d['status']=='BOUNDED_X_S_DIAGONAL_HIGH12_NATIVE_PILOT'
    assert d['deleted_original_walls']==[30,44,63,64,65]and d['removed_new_wall']==71
    assert d['N_labels']==[3,2,12]and d['S_labels']==[12,2,2]and d['RS_labels']==[2,3,4,5]
    assert d['old_frame_interface']==old['old_frame_interface']and d['exact_deletion_debit']==old['exact_deletion_debit']
    assert d['numerical_snapshot']==old['numerical_snapshot']
    assert d['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,retained_old_actual_boundary_sides=13)
    topology(d['topology']);debit(d);c=Check(d)
    assert len(d['omitted_timelike_constraints'])==39 and float(c.center@M@c.center)<0
    assert [r['wall']for r in d['omitted_timelike_constraints']]==[i for i in c.base if i not in c.ids]
    for r in d['omitted_timelike_constraints']:
        v=c.old[r['wall']];assert abs(r['projected_square']-v@M@v)<2e-7
        assert abs(r['projected_center_product']-v@M@c.center)<2e-7
        assert r['projected_square']<-1e-5 and r['projected_center_product']<-1e-5
    c.ledger(d['control'],c.old,{},d['numerical_snapshot']['actual_sides'])
    assert(d['control']['finite_rays'],d['control']['open_ends'])==(26,10)
    assert len(d['records'])==len(old['records'])==2;accepted=[]
    for row,prior in zip(d['records'],old['records']):
        assert row['index']==prior['index']and row['origins']==prior['origins']
        assert len(row['reconstruction'])==3 and set(row['prefix_normals'])==set(row['prefix_profiles'])=={'67','68','69'}
        nd=dict(c.base);tt={}
        for j,(step,prev)in enumerate(zip(row['reconstruction'],prior['reconstruction'][:3])):
            assert all(step[k]==prev[k]for k in('wall','frame','profile','labels','branch'))
            w=step['wall'];v=np.array(row['prefix_normals'][str(w)]);pi=row['prefix_profiles'][str(w)]
            assert np.max(abs(v-prior['prefix_normals'][str(w)]))<2e-6 and pi==step['profile']
            assert abs(v@M@v-c.pp[pi]['square'])<2e-5
            vals,_,_=c.native(v,pi,c.old if j<2 else nd,{}if j<2 else tt)
            assert all(relation(vals[i])==m for i,m in zip(step['frame'],step['labels']))
            nd[w]=v;tt[w]=pi
        c.bank(row['N'],nd,tt,73,[10,67,69],[3,2,12]);c.bank(row['S'],nd,tt,72,[52,68,69],[12,2,2])
        nbs=[b for b in row['N']['branches']if b['admissible']];sbs=[b for b in row['S']['branches']if b['admissible']]
        assert len(nbs)==len(sbs)==len(row['pairs'])==1
        pair=row['pairs'][0];n,s=nbs[0],sbs[0]
        assert(pair['N_branch'],pair['S_branch'])==(n['branch'],s['branch'])
        base={**nd,73:np.array(n['normal']),72:np.array(s['normal'])};bt={**tt,73:0,72:0}
        q=float(base[72]@M@base[73]);assert abs(q-pair['SN_product'])<2e-7 and relation(q)==pair['SN_relation']=='dotted'
        assert c.gates(pair['prefix_all_new_gates'],base,bt)and len(pair['R'])==4
        for m,rec in zip((2,3,4,5),pair['R']):
            c.bank(rec,base,bt,70,[69,73,72],[2,2,m])
            for b in rec['branches']:
                if b['admissible']:accepted.append((row['index'],n['branch'],s['branch'],m,b['branch'],base,bt,b))
    assert len(accepted)==len(d['accepted'])==8
    for (parent,nb,sb,m,rb,base,bt,b),a in zip(accepted,d['accepted']):
        assert(a['parent'],a['N_branch'],a['S_branch'],a['RS_label'],a['R_branch'])==(parent,nb,sb,m,rb)
        full={**base,70:np.array(b['normal'])};tt={**bt,70:0};ids=c.ids+list(tt)
        assert c.gates(a['all_new_gates'],full,tt)and len(a['corner_types'])==10 and all(a['corner_types'])
        assert all(c.finite(face,full,tt)for face in TRI)and a['all_disk_corners']
        bd=c.ledger(a,full,tt,ids);assert(a['finite_rays'],a['open_ends'])==(23,13)
        assert a['whole_compact_faces']=={str(f['wall']):f['all_rays_finite']for f in a['whole_faces']}
        for w in(67,68,69):
            f=next(f for f in a['whole_faces']if f['wall']==w)
            assert f['all_rays_finite']and sorted(sorted(r['active'])for r in f['rays'])==sorted(sorted(t)for t in TRI if w in t)
        assert a['exact_whole_interior_caps']=={'67':True,'68':True,'69':True}
        raw=a['full_projected_boundary'];assert not raw['warnings']and raw['counts']=={'finite':23,'exterior':15}and raw['open_ends']==13
        assert len(raw['rays'])==38 and [r['index']for r in raw['rays']]==list(range(38))
        assert len({tuple(r['active'])for r in raw['rays']})==38
        assert {tuple(sorted(r['active']))for r in raw['rays']if r['kind']=='finite'}=={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
        for r in raw['rays']:
            v=np.array(r['coordinates']);products={i:float(w@M@v)/np.linalg.norm(w@M)for i,w in full.items()}
            assert max(products.values())<2e-8
            assert r['active']==[i for i in full if abs(products[i])<2e-10]
            assert abs(r['square']-v@M@v)<2e-9
            assert r['kind']==('finite'if r['square']<-1e-12 else'exterior'if r['square']>1e-12 else'unresolved')
        pairs=defaultdict(list)
        for r in raw['rays']:
            for e in combinations(r['active'],2):pairs[e].append(r['index'])
        assert all(len(ii)==2 for ii in pairs.values()) and len(raw['edges'])==len(pairs)
        rebuilt=[]
        for e,ii in pairs.items():
            u,v=[raw['rays'][i]for i in ii];uv=float(np.array(u['coordinates'])@M@np.array(v['coordinates']))
            nfinite=sum(r['kind']=='finite'for r in(u,v))
            kind='compact'if nfinite==2 else'one_open'if nfinite==1 else'two_open'if uv<-1e-12 and u['square']*v['square']-uv*uv<-1e-15 else'outside'
            rebuilt.append(dict(sides=list(e),rays=ii,kind=kind))
        assert rebuilt==raw['edges']
        assert sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in rebuilt)==13
    assert c.recipes==12 and c.branches==24
    assert d['summary']==dict(parents=2,origins=4,N_recipes=2,S_recipes=2,N_branches=4,S_branches=4,joint_prefixes=2,
        R_recipes=8,R_branches=16,native_R=8,accepted_scored=8,complete_disks=8,unresolved_norm_recipes=0)
    print('PASS_SAVED_X_S_DIAGONAL',dict(new_recipes=c.recipes,new_branches=c.branches,native_configurations=8,whole_interior_caps=3,finite=23,open_ends=13),flush=True)

if __name__=='__main__':run()
