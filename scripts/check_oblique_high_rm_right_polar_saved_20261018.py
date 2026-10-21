#!/usr/bin/env python3
"""Independent right-polar saved-data and normal reconstruction checks; NumPy, no project imports.

Reconstructs the positive polar line by cofactors from saved N/S frames.
It does not solve the preceding N/S norm equations or reconstruct cones.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_high_rm_right_polar_20261018'
PREV='fresh_oblique_joint_n_r_pilot_20261018'
A='fresh_joint44_65_compact_chain_20261018'
B='fresh_joint44_63_65_link_flip_20261018'
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



class Check:
    def __init__(self,tag):
        self.d=json.loads((ROOT/f'results/{tag}.json').read_text());d=self.d
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
        snap=d['numerical_snapshot'];self.old={int(k):np.array(v) for k,v in snap['old_normals'].items()};self.center=np.array(snap['center']);self.pp=snap['profiles'];self.labels=snap['labels']
        assert len(self.old)==len(self.labels)==57 and set(self.old)==set(self.labels) and snap['metric']==M.tolist()
        self.profiles=[np.zeros(2)]
        for node in (0,1):
            for m in (3,4,5):
                v=np.zeros(2);v[node]=-2*np.cos(np.pi/m);self.profiles.append(v)
        self.J=np.array([[2.,1.],[1.,2.]])/3
        self.types={r['wall']:r['profile'] for r in d['old_frame_interface']['actual_side_profiles']}
        for r in d['old_frame_interface']['actual_side_profiles']:assert max(abs(self.profiles[r['profile']]-r['stored_A2_contacts']))<2e-5
        for pi,p in enumerate(self.pp):
            assert max(abs(self.profiles[pi]-p['contacts']))<2e-7 and abs(p['square']-(2-self.off(pi,pi)))<2e-7 and len(p['offsets'])==57
        self.base={i:v for i,v in self.old.items() if i not in d['deleted_walls']};assert len(self.base)+6==58
        self.totals=Counter()
    def off(self,a,b):return float(self.profiles[a]@self.J@self.profiles[b])
    def pt(self,i,tt):return tt.get(i,self.types.get(i,0))
    def native(self,v,pi,nd,tt):
        ofs=dict(zip(self.labels,self.pp[pi]['offsets']))
        vals={i:float(v@M@w)+(ofs[i] if i in ofs else self.off(pi,tt[i])) for i,w in nd.items()}
        bad=[i for i in nd if relation(vals[i]) in (None,'parallel_or_unresolved')]
        cp=float(v@M@self.center)+self.pp[pi]['center_offset']
        return vals,bad,cp
    def ledger(self,rec):
        d=self.d
        if rec['score'] is None:
            assert len(rec['cone_diagnostic'])==1;bd=rec['cone_diagnostic'][0]['boundary'];assert rec['cone_diagnostic'][0] in d['unresolved'] and bd['warnings']
        else:
            bd=d['boundary_ledgers'][rec['score']];assert not bd['warnings'] and bd['PD_status']=='numerically_separated_acute_incidence'
            assert bd['open_ends']==bd['PD_open_ends']
            assert {tuple(r['active']) for r in bd['rays'] if r['kind']=='finite'}=={tuple(r) for r in bd['PD_finite_supports']}
        assert rec['finite_rays']==bd['counts']['finite'] and rec['open_ends']==bd['open_ends']
        assert rec['exterior_supports']==[r['active'] for r in bd['rays'] if r['kind']=='exterior']
        assert rec['open_edges']==[e for e in bd['edges'] if e['kind'] in ('one_open','two_open')]
        assert len({f['wall'] for f in rec['whole_faces']})==len(rec['whole_faces'])
        for f in rec['whole_faces']:
            rr=[r for r in bd['rays'] if f['wall'] in r['active']]
            assert f==dict(wall=f['wall'],rays=rr,all_rays_finite=bool(rr) and all(r['kind']=='finite' for r in rr))
    def finite(self,fs,nd,tt):
        g=2*np.eye(5);g[0,1]=g[1,0]=-1
        for j,w in enumerate(fs):g[j+2,:2]=self.profiles[self.pt(w,tt)];g[:2,j+2]=g[j+2,:2]
        for i,j in combinations(range(3),2):g[i+2,j+2]=g[j+2,i+2]=float(nd[fs[i]]@M@nd[fs[j]])+self.off(self.pt(fs[i],tt),self.pt(fs[j],tt))
        if any(not isinstance(relation(g[i+2,j+2]),int) for i,j in combinations(range(3),2)):return False
        return bool(np.linalg.eigvalsh(g).min()>1e-8)
    def branch(self,b,pi,nd,tt,frame,ms):
        v=np.array(b['normal']);vals,bad,cp=self.native(v,pi,nd,tt)
        assert len(b['native_products'])==len(nd) and set(dict(b['native_products']))==set(nd)
        assert max(abs(vals[i]-x) for i,x in b['native_products'])<2e-7
        assert [i for i,_ in b['blockers']]==bad and all(abs(vals[i]-x)<2e-7 for i,x in b['blockers'])
        assert abs(v@M@v-self.pp[pi]['square'])<2e-5 and b['norm_residual']<2e-5 and abs(b['center_product']-cp)<2e-7
        assert all(relation(vals[i])==m for i,m in zip(frame,ms))
        self.totals['branches']+=1
        return v,bad,cp


def debit(d,pri):
    vs=[tuple(v) for v in pri['retained_old_finite_supports']]+[tuple(r['support']) for r in pri['new_finite_vertices']];assert len(vs)==len(set(vs))==199
    lost={i for i,v in enumerate(vs) if set(v)&set(d['deleted_walls'])};assert len(lost)==19 and sum(8 in vs[i] for i in lost)==3
    edges=Counter(e for v in vs for e in combinations(v,8));faces=defaultdict(list)
    for i,v in enumerate(vs):
        for k in range(8):
            for f in combinations(v,k):faces[f].append(i)
    touched=[]
    for f,ii in faces.items():
        if not(lost&set(ii)) or len(ii)<10-len(f):continue
        if all(edges[tuple(j for j in vs[i] if j!=k)]==2 for i in ii for k in vs[i] if k not in f):touched.append(dict(support=list(f),dimension=9-len(f),finite_supports=[list(vs[i]) for i in ii]))
    z=d['exact_deletion_debit'];assert z['lost_finite_supports']==[list(vs[i]) for i in sorted(lost)] and z['compact_faces_touched']==touched
    assert Counter(r['dimension'] for r in touched)=={2:12,3:1}


def topology(t):
    fs=[[8,10,67],[70,71,69],[71,68,69],[8,68,69],[8,67,69],[8,52,68],
        [52,68,72],[68,71,72],[71,70,72],[10,67,73],[67,69,73],[69,70,73]]
    bc=[8,10,73,70,72,52];assert t['triangles']==fs and t['boundary']==bc
    es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2));vs=sorted(set(sum(fs,[])))
    assert (len(vs),len(es),len(fs))==(t['vertices'],t['edges'],t['faces'])==(10,21,12)
    assert set(es.values())=={1,2} and all(e not in es for e in [(10,70),(67,70),(10,69),(67,68),(52,71)])
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])}
    assert [r['vertex'] for r in t['links']]==vs
    for r in t['links']:
        v=r['vertex'];ee=[list(i for i in f if i!=v) for f in fs if v in f];assert r['edges']==ee
        ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            nxt=seen|{i for e in ee if seen.intersection(e) for i in e}
            if nxt==seen:break
            seen=nxt
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))


def run():
    c=Check(TAG);d=c.d
    old=json.loads((ROOT/'results/fresh_oblique_high_n_right_r_join_20261018.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_HIGH_RM_RIGHT_POLAR_DIAGNOSTIC'
    assert d['deleted_walls']==[30,44,63,64,65]and d['exact_deletion_debit']==old['exact_deletion_debit']
    assert d['old_frame_interface']==old['old_frame_interface']and d['numerical_snapshot']==old['numerical_snapshot']
    assert d['necessary_conditions']==dict(N_profile=0,R_profile=0,RX=2,RN=2,SM=2,SR=2,RM_label_lower_bound=6,RM_squared_lower_bound=3)
    assert d['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_R=6,explicit_R_gates=58,total_fixed_facets=64)
    debit(d,pri);topology(d['topology'])
    assert len(d['parent_reconstruction'])==2
    for p,pr in zip(d['parent_reconstruction'],old['records']):
        for k in('index','origins','instruction_key','reconstruction','removed'):assert p[k]==pr[k]
    expected=[(p,n)for p in old['records']for n in p['caps']]
    assert len(d['records'])==len(expected)==20
    hist=Counter();futurehist=Counter();npol=nbranches=nfuture=0;absprod=[]
    for ri,(r,(pr,nc))in enumerate(zip(d['records'],expected)):
        assert r['index']==ri and r['parent']==pr['index']and r['parent_origins']==pr['origins']and r['N_index']==nc['index']and r['N_labels']==nc['labels']
        base={**c.base,**{int(i):np.array(v)for i,v in pr['prefix_normals'].items()}}
        bt={int(i):v for i,v in pr['prefix_profiles'].items()}
        nd={int(i):np.array(v)for i,v in r['prefix_normals'].items()};tt={int(i):v for i,v in r['prefix_profiles'].items()}
        assert set(nd)==set(tt)=={67,68,69,71,73}and all(np.max(abs(nd[i]-base[i]))<2e-6 and tt[i]==bt[i]for i in bt)
        assert tt[73]==0 and np.max(abs(nd[73]-nc['normal']))<2e-6
        ch=r['N_reconstruction'];assert[ch[k]for k in('frame','profile','labels','branch')]==[[10,67,69],0,nc['labels'],nc['branch']]
        assert ch['norm']['status']=='two_real_norm_branches'and ch['max_saved_difference']<2e-6
        vals,bad,cp=c.native(nd[73],0,base,bt);assert not bad and cp<-1e-8 and abs(nd[73]@M@nd[73]-2)<2e-5
        assert all(relation(vals[i])==m for i,m in zip([10,67,69],nc['labels']))
        nd={**c.base,**nd}
        selected=[z for z in nc['S_bank']['caps']if z['labels'][2]==2]
        assert r['complete_S_bank_size']==len(nc['S_bank']['caps'])and r['selected_S_indices']==[z['index']for z in selected]
        assert len(r['polars'])==len(selected)
        assert sha256(json.dumps(r['polars'],separators=(',',':')).encode()).hexdigest()==r['complete_polar_digest']
        for p,sc in zip(r['polars'],selected):
            npol+=1
            assert[p[k]for k in('S_index','S_profile','S_labels','S_branch')]==[sc['index'],sc['profile'],sc['labels'],sc['branch']]
            sv=np.array(p['S_normal']);assert np.max(abs(sv-sc['normal']))<2e-6
            sch=p['S_reconstruction'];assert[sch[k]for k in('frame','profile','labels','branch')]==[[52,68,71],sc['profile'],sc['labels'],sc['branch']]
            assert sch['norm']['status']=='two_real_norm_branches'and sch['max_saved_difference']<2e-6
            vals,bad,cp=c.native(sv,sc['profile'],nd,tt);assert not bad and cp<-1e-8
            assert all(relation(vals[i])==m for i,m in zip([52,68,71],sc['labels']))
            assert abs(sv@M@sv-c.pp[sc['profile']]['square'])<2e-5
            full={**nd,72:sv};ft={**tt,72:sc['profile']}
            assert p['frame']==[69,73,72]and p['status']==p['norm']['status']=='two_real_norm_branches'
            vv=np.array([full[i]for i in p['frame']]);G=vv@M@vv.T;inv=np.linalg.inv(G);eig=np.linalg.eigvalsh(G)
            assert np.max(abs(G-p['frame_Gram']))<2e-7 and np.max(abs(inv-p['frame_inverse']))<2e-7
            assert np.max(abs(eig-p['frame_eigenvalues']))<2e-7 and(sum(eig>1e-8),sum(eig< -1e-8))==(2,1)
            cov=vv@M@full[71];q=float(full[71]@M@full[71]-cov@inv@cov)
            assert np.max(abs(cov-p['M_covector']))<2e-7 and abs(q-p['M_complement_square'])<2e-7
            assert q>0 and abs(2*q-p['forced_RM_square'])<2e-7 and 2*q<3
            # Reconstruct the entire positive polar line independently by
            # signed maximal minors, not the producer's norm-root routine.
            eq=vv@M;line=np.array([(-1)**j*np.linalg.det(np.delete(eq,j,axis=1))for j in range(4)])
            assert np.max(abs(eq@line))<2e-7 and line@M@line>0
            pole=line*np.sqrt(2/(line@M@line))
            assert len(p['branches'])==2 and sorted(b['branch']for b in p['branches'])==[-1,1]
            seen=[]
            for b in p['branches']:
                v,bad,cp=c.branch(b,0,full,ft,p['frame'],[2,2,2]);nbranches+=1
                side=1 if np.max(abs(v-pole))<2e-6 else -1
                assert np.max(abs(v-side*pole))<2e-6;seen.append(side)
                vals,_,_=c.native(v,0,full,ft);rm=vals[71]
                assert abs(rm-b['RM_product'])<2e-7 and relation(rm)==b['RM_relation']and not b['finite_high_RM']
                assert abs(rm*rm-p['forced_RM_square'])<2e-6 and abs(rm)<np.sqrt(3)
                hist[str(tuple(bad))]+=1;absprod.append(abs(rm))
                if cp<-1e-8:
                    assert bad==[68,71]and rm>0;nfuture+=1;futurehist[str(tuple(bad))]+=1
                else:assert cp>0 and rm<0
            assert sorted(seen)==[-1,1]
    s=d['summary'];assert npol==192 and nbranches==384 and nfuture==192
    assert s==dict(parent_origins=4,distinct_parents=2,N_geometries=20,complete_S_caps=712,right_S_polars=192,branches=384,future_branches=192,native_future_branches=0,finite_high_RM_branches=0,accepted=0,new_cones=0,unresolved=0,minimum_abs_RM=min(absprod),maximum_abs_RM=max(absprod),future_blocker_histogram=dict(futurehist),complete_blocker_histogram=dict(hist))
    assert not d['accepted']and not d['boundary_ledgers']and not d['unresolved']
    print('PASS_SAVED_HIGH_RM_RIGHT_POLAR','192cofactorlines','384branches','all_abs_RM_below_sqrt3','19/12/1debit',flush=True)


if __name__=='__main__':run()
