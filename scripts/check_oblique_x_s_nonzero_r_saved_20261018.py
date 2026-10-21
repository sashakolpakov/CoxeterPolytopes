#!/usr/bin/env python3
"""Independent nonzero-R saved-domain checks; NumPy, no project imports.

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
TAG='fresh_oblique_x_s_nonzero_r_profiles_20261018'
PREV='fresh_oblique_x_s_diagonal_native_pilot_20261018'
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

def finite_duty(c,pi,m,high):
    g=2*np.eye(5);g[0,1]=g[1,0]=-1
    g[4,:2]=c.profiles[pi];g[:2,4]=c.profiles[pi]
    for i,j,label in((2,3,12 if high else 2),(2,4,2),(3,4,2 if high else m)):
        g[i,j]=g[j,i]=-2*np.cos(np.pi/label)
    return bool(np.linalg.eigvalsh(g).min()>1e-8)

def branch(c,b,pi,nd,tt,frame,ms):
    v=np.array(b['normal']);vals,bad,cp=c.native(v,pi,nd,tt)
    assert len(b['native_products'])==len(nd)and set(dict(b['native_products']))==set(nd)
    assert all(abs(vals[i]-q)<2e-7 for i,q in b['native_products'])
    assert [i for i,_ in b['blockers']]==bad and all(abs(vals[i]-q)<2e-7 for i,q in b['blockers'])
    rhs=[-2*np.cos(np.pi/m)-c.off(pi,tt[i])for i,m in zip(frame,ms)]
    err=max(abs(float(v@M@v)-c.pp[pi]['square']),max(abs(float(v@M@nd[i])-q)for i,q in zip(frame,rhs)))
    assert err<2e-5 and abs(err-b['norm_residual'])<2e-7 and abs(cp-b['center_product'])<2e-7
    assert bool(not bad and cp<-1e-8 and err<2e-5)==b['native_admissible']
    return v

def raw_boundary(raw,full,finite):
    assert not raw['warnings']and raw['counts']=={'finite':23,'exterior':15}and raw['open_ends']==13
    assert len(raw['rays'])==38 and[r['index']for r in raw['rays']]==list(range(38))
    assert len({tuple(r['active'])for r in raw['rays']})==38
    assert {tuple(sorted(r['active']))for r in raw['rays']if r['kind']=='finite'}==finite
    pairs=defaultdict(list)
    for r in raw['rays']:
        v=np.array(r['coordinates']);ps={i:float(w@M@v)/np.linalg.norm(w@M)for i,w in full.items()}
        assert max(ps.values())<2e-8 and r['active']==[i for i in full if abs(ps[i])<2e-10]
        assert abs(r['square']-v@M@v)<2e-9
        assert r['kind']==('finite'if r['square']<-1e-12 else'exterior'if r['square']>1e-12 else'unresolved')
        for e in combinations(r['active'],2):pairs[e].append(r['index'])
    assert all(len(ii)==2 for ii in pairs.values())and len(raw['edges'])==len(pairs)
    rebuilt=[]
    for e,ii in pairs.items():
        u,v=[raw['rays'][i]for i in ii];uv=float(np.array(u['coordinates'])@M@np.array(v['coordinates']))
        nf=sum(r['kind']=='finite'for r in(u,v))
        kind='compact'if nf==2 else'one_open'if nf==1 else'two_open'if uv<-1e-12 and u['square']*v['square']-uv*uv<-1e-15 else'outside'
        rebuilt.append(dict(sides=list(e),rays=ii,kind=kind))
    assert rebuilt==raw['edges']
    assert sum(1 if e['kind']=='one_open'else 2 if e['kind']=='two_open'else 0 for e in rebuilt)==13

def run():
    d=json.loads((ROOT/f'results/{TAG}.json').read_text());old=json.loads((ROOT/f'results/{PREV}.json').read_text())
    for doc in(d,old):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert d['status']=='BOUNDED_X_S_NONZERO_R_PROFILE_COMPLETION'and d['predecessor']==PREV
    assert d['deleted_original_walls']==old['deleted_original_walls']==[30,44,63,64,65]
    for k in('exact_deletion_debit','topology','old_frame_interface','numerical_snapshot','omitted_timelike_constraints','wall_count_breakdown'):
        assert d[k]==old[k]
    topology(d['topology']);debit(d);c=Check(d)
    assert len(d['records'])==len(old['records'])==2
    statuses=Counter();blockers=Counter();expected=[];nbranches=0
    for row,prior in zip(d['records'],old['records']):
        assert row['parent']==prior['index']and row['origins']==prior['origins']
        assert len(row['reconstruction'])==5 and set(row['normals'])==set(row['profiles'])=={'67','68','69','73','72'}
        nd=dict(c.base);tt={}
        instructions=prior['reconstruction'][:]
        saved={int(k):v for k,v in prior['prefix_normals'].items()}
        for key,w,bkey in(('N',73,'N_branch'),('S',72,'S_branch')):
            r=prior[key];b=next(b for b in r['branches']if b['branch']==prior['pairs'][0][bkey]);assert b['admissible']
            instructions.append(dict(wall=w,frame=r['frame'],profile=r['profile'],labels=r['labels'],branch=b['branch']))
            saved[w]=b['normal']
        for j,(step,prev)in enumerate(zip(row['reconstruction'],instructions)):
            assert all(step[k]==prev[k]for k in('wall','frame','profile','labels','branch'))
            w=step['wall'];v=np.array(row['normals'][str(w)]);pi=row['profiles'][str(w)]
            assert np.max(abs(v-saved[w]))<2e-6 and pi==step['profile']
            assert abs(v@M@v-c.pp[pi]['square'])<2e-5
            vals,_,_=c.native(v,pi,c.old if j<2 else nd,{}if j<2 else tt)
            assert all(relation(vals[i])==m for i,m in zip(step['frame'],step['labels']))
            nd[w]=v;tt[w]=pi
        assert 70 not in nd and 71 not in nd
        assert len(row['recipes'])==24
        domain=[(pi,m)for pi in range(1,7)for m in range(2,6)]
        for ri,((pi,m),r)in enumerate(zip(domain,row['recipes'])):
            assert(r['index'],r['profile'],r['frame'],r['labels'])==(ri,pi,[69,73,72],[2,2,m])
            assert r['required_corners']==[[69,73,70],[69,72,70]]
            finite=[finite_duty(c,pi,m,True),finite_duty(c,pi,m,False)]
            assert [x is not None for x in r['finite_duty_types']]==finite
            # The explicit finite chains are A/B/H3 at XNR and
            # A4/B4/F4/H4 or their disconnected subdiagrams at XRS.
            allowed=m in({2,3,4,5}if pi in(1,4)else{2,3}if pi in(2,5)else{2})
            assert finite==[True,allowed]
            statuses[r['status']]+=1
            if not allowed:
                assert r['status']=='finite_duty_rejected'and r['norm']is None and not r['branches'];continue
            assert r['status']==r['norm']['status']=='two_real_norm_branches'
            assert len(r['branches'])==2 and sorted(b['branch']for b in r['branches'])==[-1,1]
            g=np.array([[float(nd[i]@M@nd[j])for j in r['frame']]for i in r['frame']])
            assert abs(np.linalg.det(g))>1e-7 and np.max(abs(g-r['norm']['gram']))<2e-7
            assert np.linalg.norm(np.array(r['branches'][0]['normal'])-r['branches'][1]['normal'])>1e-6
            for b in r['branches']:
                v=branch(c,b,pi,nd,tt,r['frame'],r['labels']);nbranches+=1;blockers[str(tuple(i for i,_ in b['blockers']))]+=1
                if not b['native_admissible']:
                    assert b['accepted_index']is None and b['retained_corner_types']is None and b['all_new_gates']is None;continue
                full={**nd,70:v};ft={**tt,70:pi}
                assert len(b['retained_corner_types'])==2 and all(b['retained_corner_types'])
                assert all(c.finite(f,full,ft)for f in r['required_corners'])and c.gates(b['all_new_gates'],full,ft)
                assert b['accepted_index']==len(expected);expected.append((row['parent'],ri,pi,m,b['branch'],full,ft))
    assert len(expected)==len(d['accepted'])==4 and nbranches==56
    for (parent,ri,pi,m,br,full,tt),a in zip(expected,d['accepted']):
        assert(a['parent'],a['recipe'],a['profile'],a['RS_label'],a['branch'])==(parent,ri,pi,m,br)
        assert pi in(1,4)and m==5 and br==1
        assert len(a['corner_types'])==10 and all(a['corner_types'])and a['all_disk_corners']
        assert all(c.finite(f,full,tt)for f in TRI)
        bd=c.ledger(a,full,tt,c.ids+list(tt));assert(a['finite_rays'],a['open_ends'])==(23,13)
        assert a['whole_compact_faces']=={str(f['wall']):f['all_rays_finite']for f in a['whole_faces']}
        for w in(67,68,69):
            f=next(f for f in a['whole_faces']if f['wall']==w)
            assert f['all_rays_finite']and sorted(sorted(r['active'])for r in f['rays'])==sorted(sorted(f)for f in TRI if w in f)
        assert a['exact_whole_interior_caps']=={'67':True,'68':True,'69':True}
        finite={tuple(sorted(r['active']))for r in bd['rays']if r['kind']=='finite'}
        raw_boundary(a['full_projected_boundary'],full,finite)
        before=next(v for v in old['accepted']if v['parent']==parent and v['RS_label']==m)
        assert a['whole_compact_faces']==before['whole_compact_faces']
        assert a['exterior_supports']==before['exterior_supports']and a['open_edges']==before['open_edges']
        assert [(r['active'],r['kind'])for r in a['full_projected_boundary']['rays']]==[(r['active'],r['kind'])for r in before['full_projected_boundary']['rays']]
    assert not d['unresolved_cones']
    assert d['summary']==dict(parents=2,origins=4,profiles=list(range(1,7)),RS_labels=list(range(2,6)),recipes=48,
        status_counts=dict(statuses),branches=56,native_branches=4,accepted=4,complete_disks=4,
        boundary_patterns={'(True, 23, 13)':4},unresolved_norm_recipes=0,blocker_histogram=dict(blockers))
    print('PASS_SAVED_X_S_NONZERO_R',dict(recipes=48,finite_rejections=20,norm_branches=56,survivors=4,unchanged_open_ends=13),flush=True)

if __name__=='__main__':run()
