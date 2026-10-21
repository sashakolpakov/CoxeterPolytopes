#!/usr/bin/env python3
"""Independent remote60 bridge saved-domain checks; NumPy, no project imports.

Checks all new physical branches, native gates, finite supports and saved
cone accounting. Does not reconstruct cones or the historical scaffold.
"""
import json
from collections import Counter,defaultdict
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_remote60_compact_face_bridge_20261018'
CAP='fresh_oblique_n_r_s_polar_cap_20261018'
NONZERO='fresh_oblique_x_s_nonzero_r_profiles_20261018'
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

def finite_duty(c,profiles,m,ma,mb):
    g=2*np.eye(5);g[0,1]=g[1,0]=-1
    for k,pi in enumerate(profiles,2):g[k,:2]=c.profiles[pi];g[:2,k]=c.profiles[pi]
    for i,j,z in((2,3,m),(2,4,ma),(3,4,mb)):g[i,j]=g[j,i]=-2*np.cos(np.pi/z)
    return bool(np.linalg.eigvalsh(g).min()>1e-8)

def saved_bank(c,bk,nd,tt,frame):
    assert bk['frame']==frame and len(bk['recipes'])==448
    ps=[c.pt(i,tt)for i in frame]
    bonds=[relation(float(nd[frame[0]]@M@nd[frame[k]])+c.off(ps[0],ps[k]))for k in(1,2)]
    assert ps==bk['frame_profiles']and bonds==bk['frame_bonds']
    counts=Counter();native=[];nroots=0
    for ix,((pi,ms),r)in enumerate(zip(product(range(7),product(range(2,6),repeat=3)),bk['recipes'],strict=True)):
        assert(r['index'],r['profile'],r['frame'],r['labels'])==(ix,pi,frame,list(ms))
        assert r['required_corners']==[[frame[0],frame[k],75]for k in(1,2)]
        duties=[finite_duty(c,[ps[0],ps[k],pi],bonds[k-1],ms[0],ms[k])for k in(1,2)]
        assert duties==[x is not None for x in r['finite_duty_types']]
        counts[r['status']]+=1
        if not all(duties):assert r['status']=='finite_duty_rejected'and r['norm']is None and not r['branches'];continue
        fs=np.array([nd[i]for i in frame]);g=fs@M@fs.T
        if 'gram' in r['norm']:assert np.max(abs(g-np.array(r['norm']['gram'])))<2e-7
        eig=np.linalg.eigvalsh(g);assert min(abs(eig))>1e-8
        ofs=dict(zip(c.labels,c.pp[pi]['offsets']))
        rhs=np.array([-2*np.cos(np.pi/m)-(ofs[i]if i in ofs else c.off(pi,tt[i]))for i,m in zip(frame,ms)])
        part=rhs@np.linalg.inv(g)@fs;remaining=float(c.pp[pi]['square']-part@M@part)
        perp_sign=-1 if all(eig>0)else 1
        assert sum(eig<0)in(0,1)and abs(remaining)>1e-8
        if perp_sign*remaining<0:
            assert r['status']=='no_real_norm_root'and not r['branches'];continue
        assert r['status']=='two_real_norm_branches'and len(r['branches'])==2
        assert sorted(x['branch']for x in r['branches'])==[-1,1]
        assert np.linalg.norm(np.array(r['branches'][0]['normal'])-r['branches'][1]['normal'])>1e-6
        for b in r['branches']:
            nroots+=1;v=np.array(b['normal']);vv,bb,cp=c.native(v,pi,nd,tt)
            assert len(b['native_products'])==len(nd)and set(dict(b['native_products']))==set(nd)
            assert all(abs(vv[i]-z)<2e-7 for i,z in b['native_products'])
            assert [i for i,_ in b['blockers']]==bb and all(abs(vv[i]-z)<2e-7 for i,z in b['blockers'])
            err=max(abs(v@M@v-c.pp[pi]['square']),max(abs(fs@M@v-rhs)))
            assert err<2e-5 and abs(err-b['norm_residual'])<2e-7 and abs(cp-b['center_product'])<2e-7
            assert b['native_admissible']==bool(not bb and cp<-1e-8 and err<2e-5)
            dup=[i for i in c.d['numerical_snapshot']['actual_sides']if c.pt(i,{})==pi and np.max(abs(v-c.old[i]))<2e-5]
            assert dup==b['original_duplicates']
            if b['native_admissible']:native.append((r,b))
    assert bk['summary']==dict(recipes=448,status_counts=dict(counts),branches=nroots,native_branches=len(native))
    return native

def raw_check(raw,full):
    assert raw['counts']=={'finite':27,'exterior':15}and raw['open_ends']==13 and not raw['warnings']
    assert len(raw['rays'])==42 and[r['index']for r in raw['rays']]==list(range(42))
    assert len({tuple(r['active'])for r in raw['rays']})==42
    pairs=defaultdict(list)
    for r in raw['rays']:
        v=np.array(r['coordinates']);ps={i:float(w@M@v)/np.linalg.norm(w@M)for i,w in full.items()}
        assert abs(np.linalg.norm(v)-1)<2e-8 and max(ps.values())<2e-8
        assert r['active']==[i for i in full if abs(ps[i])<2e-10]and len(r['active'])==3
        assert abs(r['square']-v@M@v)<2e-9 and abs(r['square'])>1e-12
        assert r['kind']==('finite'if r['square']<0 else'exterior')
        for e in combinations(r['active'],2):pairs[e].append(r['index'])
    edges=[]
    for e,ii in pairs.items():
        assert len(ii)==2
        a,b=[raw['rays'][i]for i in ii];ab=float(np.array(a['coordinates'])@M@np.array(b['coordinates']));nf=sum(x['kind']=='finite'for x in(a,b))
        kind='compact'if nf==2 else'one_open'if nf==1 else'two_open'if ab<-1e-12 and a['square']*b['square']-ab*ab<-1e-15 else'outside'
        edges.append(dict(sides=list(e),rays=ii,kind=kind))
    assert edges==raw['edges']and sum(1 if x['kind']=='one_open'else 2 if x['kind']=='two_open'else 0 for x in edges)==13

def run():
    read=lambda stem:json.loads((ROOT/f'results/{stem}.json').read_text())
    d,prior,original=read(TAG),read(CAP),read(PREV)
    for x in(d,prior,original):
        for p,h in x['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    assert d['status']=='BOUNDED_REMOTE60_WHOLE_FACE_BRIDGE'and d['predecessor']==CAP and d['new_wall']==75
    assert d['deleted_original_walls']==[30,44,63,64,65]
    for key in('numerical_snapshot','old_frame_interface','exact_deletion_debit'):
        assert d[key]==prior[key]==original[key]
    assert d['main_topology']==prior['topology']and 60 not in d['main_topology']['vertices']
    tri=[[58,59,60],[58,60,75],[59,60,75]]
    assert d['remote_topology']==dict(triangles=tri,boundary=[58,59,75],interior_wall=60,interior_link=[58,59,75],V=4,E=6,F=3,euler=1,
        separate_from_main_T_disk=True,forbidden_finite_duty=[60,66,75])
    assert not set(sum(tri,[]))&set(d['main_topology']['vertices'])
    assert Counter(tuple(sorted(e))for f in tri for e in combinations(f,2))=={(58,59):1,(58,60):2,(59,60):2,(58,75):1,(60,75):2,(59,75):1}
    assert d['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,retained_new_before_U=7,all_facets_after_U=66)
    debit(d);c=Check(d)
    assert len(d['preflight_parents'])==len(original['records'])==2 and len(d['failed_preflights'])==6
    for ix,(row,p)in enumerate(zip(d['preflight_parents'],original['records'],strict=True)):
        assert(row['parent'],row['origins'])==(p['index'],p['origins'])
        steps=p['reconstruction'][:];normals=dict(p['prefix_normals'])
        for key,w,bk in(('N',73,'N_branch'),('S',72,'S_branch')):
            r=p[key];b=next(b for b in r['branches']if b['branch']==p['pairs'][0][bk])
            steps.append(dict(wall=w,frame=r['frame'],profile=r['profile'],labels=r['labels'],branch=b['branch']));normals[str(w)]=b['normal']
        assert len(row['reconstruction'])==len(steps)==5
        nd=dict(c.base);tt={}
        for j,(st,ex)in enumerate(zip(row['reconstruction'],steps,strict=True)):
            assert all(st[k]==ex[k]for k in('wall','frame','profile','labels','branch'))
            w=st['wall'];v=np.array(row['normals'][str(w)]);pi=st['profile']
            assert np.max(abs(v-normals[str(w)]))<2e-6 and row['profiles'][str(w)]==pi
            vv,_,_=c.native(v,pi,c.old if j<2 else nd,{}if j<2 else tt)
            assert abs(v@M@v-c.pp[pi]['square'])<2e-5 and all(relation(vv[i])==m for i,m in zip(st['frame'],st['labels']))
            nd[w]=v;tt[w]=pi
        assert set(row['normals'])==set(row['profiles'])==set(map(str,tt))
        for j,(name,fr)in enumerate([('old10',[10,33,73]),('old52',[52,47,72]),('old55',[55,7,59])]):
            bk=d['failed_preflights'][3*ix+j];assert(bk['parent'],bk['target'])==(p['index'],name)
            assert saved_bank(c,bk,nd,tt,fr)==[]
    roots=saved_bank(c,d['remote_bank'],c.base,{},[60,58,59]);assert len(roots)==15
    assert all(r['profile']==0 and r['labels'][0]==2 and b['branch']==-1 for r,b in roots)
    assert sum(len(r['branches'])for r in d['remote_bank']['recipes'])==30
    assert len(d['configurations'])==len(prior['accepted'])==28 and len(d['accepted'])==420
    expected=0;direct=0
    for ix,(row,oldcap)in enumerate(zip(d['configurations'],prior['accepted'],strict=True)):
        p=prior['records'][oldcap['configuration']];tr=p['recipes'][oldcap['profile']];tb=next(x for x in tr['branches']if x['branch']==oldcap['branch'])
        assert(row['index'],row['source_configuration'],row['T_profile'],row['T_branch'])==(ix,oldcap['configuration'],oldcap['profile'],oldcap['branch'])
        exnorm={**p['normals'],'74':tb['normal']};extt={**p['profiles'],'74':oldcap['profile']}
        assert set(row['normals'])==set(row['profiles'])==set(exnorm)=={'67','68','69','70','72','73','74'}
        assert row['profiles']==extt and all(np.max(abs(np.array(row['normals'][i])-v))<2e-6 for i,v in exnorm.items())
        nd={**c.base,**{int(i):np.array(v)for i,v in row['normals'].items()}};tt={int(i):v for i,v in row['profiles'].items()}
        assert len(row['pairs'])==15
        for (r,b),z in zip(roots,row['pairs'],strict=True):
            assert(z['recipe'],z['branch'])==(r['index'],b['branch']);v=np.array(b['normal']);vv,bb,cp=c.native(v,0,nd,tt)
            assert not bb and cp<-1e-8 and z['native_admissible']and not z['blockers']
            assert len(z['native_products'])==len(nd)and set(dict(z['native_products']))==set(nd)
            assert all(abs(vv[i]-q)<2e-7 for i,q in z['native_products'])and abs(cp-z['center_product'])<2e-7
            assert z['accepted_index']==expected;a=d['accepted'][expected];expected+=1
            assert(a['configuration'],a['recipe'],a['branch'])==(ix,r['index'],b['branch'])
            full={**nd,75:v};ft={**tt,75:0};raw=a['full_projected_boundary'];proof=a['clipping'];old=oldcap['full_projected_boundary']
            cov=v@M/np.linalg.norm(v@M);products=[float(cov@np.array(q['coordinates']))for q in old['rays']]
            assert min(map(abs,products))>1e-10 and np.max(abs(np.array(products)-proof['normalized_ray_products']))<2e-12
            keep=[i for i,q in enumerate(products)if q<0];removed=[i for i,q in enumerate(products)if q>0]
            assert keep==proof['retained_old_rays']and removed==proof['removed_old_rays']and proof['complete_input_cone_inherited']
            assert {tuple(old['rays'][i]['active'])for i in removed}=={(58,60,66),(59,60,66)}
            vectors=[old['rays'][i]['coordinates']for i in keep];cross=[]
            for ei,e in enumerate(old['edges']):
                i,j=e['rays'];x,y=products[i],products[j]
                if x*y>=0:continue
                if x<0:i,j=j,i;x,y=y,x
                vv=x*np.array(old['rays'][j]['coordinates'])-y*np.array(old['rays'][i]['coordinates'])
                cross.append(dict(old_edge=ei,positive_ray=i,negative_ray=j,new_ray=len(vectors)));vectors.append(vv)
            assert len(cross)==4 and cross==proof['crossing_edges']and len(vectors)==42
            assert all(np.max(abs(np.array(vec)/np.linalg.norm(vec)-q['coordinates']))<2e-7 for vec,q in zip(vectors,raw['rays'],strict=True))
            raw_check(raw,full)
            ids=c.ids+list(ft);g=np.array([[full[i]@M@full[j]for j in ids]for i in ids])
            pd=sorted(sorted(ids[i]for i in zz)for zz in combinations(range(len(ids)),3)if np.linalg.eigvalsh(g[np.ix_(zz,zz)]).min()>1e-8)
            finite=sorted(sorted(q['active'])for q in raw['rays']if q['kind']=='finite');assert finite==pd==a['PD_finite_supports']
            assert len(a['whole_faces'])==len(ids)and[f['wall']for f in a['whole_faces']]==ids
            for f in a['whole_faces']:
                rr=[q for q in raw['rays']if f['wall']in q['active']]
                assert f['rays']==rr and f['all_rays_finite']==all(q['kind']=='finite'for q in rr)
            sixty=next(f for f in a['whole_faces']if f['wall']==60)
            assert sixty['all_rays_finite']and sorted(sorted(q['active'])for q in sixty['rays'])==tri
            assert all(c.finite(f,full,ft)for f in tri)and len(a['finite_corner_types'])==3 and all(a['finite_corner_types'])
            oldf={tuple(sorted(q['active']))for q in old['rays']if q['kind']=='finite'};olde={tuple(sorted(q['active']))for q in old['rays']if q['kind']=='exterior'}
            newf={tuple(q)for q in finite};newe={tuple(sorted(q['active']))for q in raw['rays']if q['kind']=='exterior'}
            assert newf-oldf=={(58,60,75),(59,60,75)}and not oldf-newf
            assert olde-newe=={(58,60,66),(59,60,66)}and newe-olde=={(58,66,75),(59,66,75)}
            assert a['finite_supports_added']==[[58,60,75],[59,60,75]]and a['finite_supports_removed']==[]
            assert a['exterior_supports_removed']==[[58,60,66],[59,60,66]]and a['exterior_supports_added']==[[58,66,75],[59,66,75]]
            if ix==0:
                direct+=1;raw_check(a['direct_cone'],full)
                assert sorted((sorted(q['active']),q['kind'])for q in raw['rays'])==sorted((sorted(q['active']),q['kind'])for q in a['direct_cone']['rays'])
            else:assert a['direct_cone']is None
    assert expected==420 and direct==15
    assert d['summary']==dict(preflight_families=6,preflight_recipes=2688,preflight_native=0,remote_recipes=448,remote_eligible=17,
        remote_real_roots=30,remote_native=15,T_configurations=28,complete_pairs=420,accepted=420,direct_cones=15,
        finite_rays=27,exterior_rays=15,open_ends=13,new_whole_old_face=60,unresolved=0)
    print('PASS_STANDALONE_REMOTE60_BRIDGE',d['summary'],flush=True)

if __name__=='__main__':run()

