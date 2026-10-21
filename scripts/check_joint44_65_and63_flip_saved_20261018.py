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
        self.base={i:v for i,v in self.old.items() if i not in d['deleted_walls']};assert d['retained_old_walls']==len(self.base)+6
        self.totals=Counter();self.ledger(d['control']);self.ledger(d['deletion'])
        assert d['control']['finite_rays']==26 and d['control']['open_ends']==10
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
    def bank(self,bank,nd,tt):
        frame=bank['frame'];duties=bank['duties'];bonds=bank['duty_bonds'];domain=[]
        assert bank['frame_profiles']==[self.pt(i,tt) for i in frame]
        for a,b in duties:assert relation(float(nd[frame[a]]@M@nd[frame[b]])+self.off(self.pt(frame[a],tt),self.pt(frame[b],tt)))==bonds[3*a+b]
        for pi in range(7):
            for ms in product(range(2,6),repeat=3):
                good=True
                for a,b in duties:
                    g=2*np.eye(5);g[0,1]=g[1,0]=-1
                    for k,p in enumerate((self.pt(frame[a],tt),self.pt(frame[b],tt),pi)):g[k+2,:2]=self.profiles[p];g[:2,k+2]=self.profiles[p]
                    for i,j,m in ((2,3,bonds[3*a+b]),(2,4,ms[a]),(3,4,ms[b])):g[i,j]=g[j,i]=-2*np.cos(np.pi/m)
                    good &= np.linalg.eigvalsh(g).min()>1e-8
                if good:domain.append((pi,list(ms)))
        assert [(r['profile'],r['labels']) for r in bank['recipes']]==domain and len(domain)==bank['complete_low_recipe_count']
        digest=sha256();hist=Counter();seen=set();ci=0
        for ri,r in enumerate(bank['recipes']):
            assert r['index']==ri and r['norm']['status']=='two_real_norm_branches' and sorted(b['branch'] for b in r['branches'])==[-1,1]
            pi=r['profile']
            for b in r['branches']:
                v,bad,cp=self.branch(b,pi,nd,tt,frame,r['labels'])
                hist[str(tuple(bad))]+=1;digest.update(json.dumps([ri,b['branch'],pi,r['labels'],b['native_products'],b['center_product']],separators=(',',':')).encode())
                if bad or cp>=-1e-8:assert 'cap' not in b;continue
                key=(pi,)+tuple(np.round(v,8))
                if key in seen:assert b['duplicate'] and 'cap' not in b;continue
                seen.add(key);assert b['cap']==ci;c=bank['caps'][ci];ci+=1
                assert c['index']==b['cap'] and c['recipe']==ri and c['branch']==b['branch'] and c['profile']==pi and c['labels']==r['labels'] and c['normal']==b['normal']
                assert c['native_products']==b['native_products'] and c['center_product']==b['center_product'];self.ledger(c)
        assert ci==len(bank['caps']) and digest.hexdigest()==bank['complete_branch_digest'] and dict(hist)==bank['branch_histogram']
        self.totals['recipes']+=len(domain);self.totals['caps']+=ci


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


def run():
    ca=Check(A);a=ca.d;cb=Check(B);b=cb.d
    assert a['status']=='BOUNDED_JOINT44_65_COMPACT_CHAIN' and b['status']=='BOUNDED_JOINT44_63_65_LINK_FLIP'
    assert a['deleted_walls']==[30,44,65] and b['deleted_walls']==[30,44,63,65]
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text());debit(a,pri);debit(b,pri)
    assert a['additional63_deletion_debit']==b['exact_deletion_debit']
    pairs=list(combinations([8,10,33,44,52,63,64,65],2));iv=original_intervals(pri,pairs)
    assert [r['walls'] for r in a['original_native_intervals']]==list(map(list,pairs))
    for r in a['original_native_intervals']:
        i,j=r['walls'];x=float(ca.old[i]@M@ca.old[j])+ca.off(ca.types[i],ca.types[j]);q=iv[i,j]
        assert r['exact_native_enclosure']==q and abs(x-r['native_product'])<2e-7 and relation(x)==r['relation']
        assert float(F(q[0]))-2e-5<=x<=float(F(q[1]))+2e-5
    assert len(a['J_recipes'])==4 and len(a['J_caps'])==4
    for ri,r in enumerate(a['J_recipes']):
        assert r['index']==ri and r['profile']==0 and r['labels']==[2,ri+2,2] and sorted(x['branch'] for x in r['branches'])==[-1,1]
        for br in r['branches']:
            v,bad,cp=ca.branch(br,0,ca.base,{},[8,10,63],r['labels'])
            if bad or cp>=-1e-8:assert 'cap' not in br;continue
            c=a['J_caps'][br['cap']];assert c['normal']==br['normal'] and c['recipe']==ri and c['branch']==br['branch'];ca.ledger(c)
            assert c['restores_old44']==bool(np.max(abs(v-ca.old[44]))<2e-6)
    ca.bank(a['L_bank'],ca.base,{})
    assert len(a['joins'])==4*87 and [(r['J'],r['L']) for r in a['joins']]==list(product(range(4),range(87)))
    pc=0
    for r in a['joins']:
        j=a['J_caps'][r['J']];l=a['L_bank']['caps'][r['L']];nd={**ca.base,67:np.array(j['normal']),68:np.array(l['normal'])};tt={67:0,68:l['profile']}
        vals,bad,cp=ca.native(nd[68],tt[68],{i:v for i,v in nd.items() if i!=68},tt)
        assert abs(r['native_product']-vals[67])<2e-7 and relation(vals[67])==r['relation'] and [i for i,_ in r['blockers']]==bad and abs(cp-r['center_product'])<2e-7
        assert [x is not None for x in r['corner_types']]==[ca.finite(fs,nd,tt) for fs in a['L_five_corner_duties']]
        if bad or cp>=-1e-8:assert 'pair' not in r;continue
        p=a['pairs'][pc];assert r['pair']==pc and p['index']==pc and p['J']==r['J'] and p['L']==r['L'];pc+=1;ca.ledger(p)
        assert p['all_five_corners']==all(x is not None for x in r['corner_types'])
    assert pc==345 and len(a['R_banks'])==len(a['final'])==1
    rb=a['R_banks'][0];p=a['pairs'][rb['pair']];j=a['J_caps'][p['J']];l=a['L_bank']['caps'][p['L']]
    assert p['both_controls'] and p['all_five_corners']
    ca.bank(rb['bank'],{**ca.base,67:np.array(j['normal']),68:np.array(l['normal'])},{67:0,68:l['profile']})
    assert len(rb['bank']['caps'])==1 and a['final'][0]['restores_old30'] and a['final'][0]['clean'] and all(a['final'][0]['whole_compact_faces'].values()) and a['final'][0]['open_ends']==10
    ca.bank(a['necessary63_64_bank'],ca.base,{})
    nc=a['necessary63_64_bank']['caps'];assert len(nc)==20
    for c in nc:
        finite=ca.finite([8,52,68],{**ca.base,68:np.array(c['normal'])},{68:c['profile']})
        assert finite==(c['remaining8_52_corner'] is not None)==c['restores_old65']
    assert sum(c['restores_old65'] for c in nc)==1 and len(a['unresolved'])==457
    top=b['topology'];fs=top['triangles'];bc=top['boundary'];es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2))
    assert len(set(sum(fs,[])))==8 and len(es)==16 and len(fs)==9 and (67,68) not in es
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])} and set(es.values())=={1,2}
    assert len(top['links'])==8
    for rec in top['links']:
        v=rec['vertex'];ee=[list(i for i in f if i!=v) for f in fs if v in f];assert ee==rec['edges'];ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            ss=seen|{i for e in ee if seen.intersection(e) for i in e}
            if ss==seen:break
            seen=ss
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))
    selected=[]
    for j in a['J_caps']:
        if j['restores_old44']:continue
        for pi in range(7):selected.append(min((x for x in a['joins'] if x['J']==j['index'] and not x['blockers'] and a['L_bank']['caps'][x['L']]['profile']==pi),key=lambda x:(abs(x['native_product']+2),x['L'])))
    assert b['selected']==selected and len(b['records'])==len(selected)==21
    x64=[]
    for ri,(r,sel) in enumerate(zip(b['records'],selected)):
        assert r['index']==ri and r['selected_join']==sel and r['J']==sel['J'] and r['L']==sel['L'];j=a['J_caps'][r['J']];l=a['L_bank']['caps'][r['L']]
        nd={**cb.base,67:np.array(j['normal']),68:np.array(l['normal'])};tt={67:0,68:l['profile']}
        for i in (67,68):
            _,bad,cp=cb.native(nd[i],tt[i],{z:v for z,v in nd.items() if z!=i},tt);assert not bad and cp<-1e-8
        cb.ledger(r['prefix']);assert r['bank']['frame']==[8,67,68] and r['bank']['duties']==[[0,1],[0,2]]
        cb.bank(r['bank'],nd,tt)
        for c in r['bank']['caps']:
            x=dict(c['native_products'])[64];x64.append(x);assert relation(x)=='dotted' and c['corner_L64'] is None
            assert c['score'] is None and c['finite_rays']==18 and c['open_ends']==14
            assert [f['wall'] for f in c['whole_faces'] if f['all_rays_finite']]==[8,11]
            f8=next(f for f in c['whole_faces'] if f['wall']==8);assert sorted(r['active'] for r in f8['rays'])==sorted([[7,8,11],[7,8,12],[8,10,11],[8,10,67],[8,12,52],[8,52,68],[8,67,69],[8,68,69]])
            fx=next(f for f in c['whole_faces'] if f['wall']==69);assert sorted(r['active'] for r in fx['rays'] if r['kind']=='exterior')==[[64,67,69],[64,68,69]]
    assert len(x64)==270 and min(x64)<-6.88 and -2.009<max(x64)<-2.008
    assert not b['R_banks'] and not b['accepted'] and len(b['unresolved'])==292
    assert ca.totals==dict(branches=340,recipes=166,caps=108) and cb.totals==dict(branches=1584,recipes=792,caps=270)
    assert a['summary']==dict(J_recipes=4,J_branches=8,J_native=4,J_changed=3,L_recipes=96,L_branches=192,L_native=87,joins=348,native_pairs=345,five_corner_pairs=1,changed_five_corner_pairs=0,R_banks=1,R_recipes=42,R_branches=84,native_triples=1,changed_triples=0,clean_triples=1,restored_compact_faces=1,fewer_than_ten_ends=0,necessary_recipes=28,necessary_branches=56,necessary_native=20,necessary_completed=1,unresolved=457)
    assert b['summary']==dict(selected_prefixes=21,X_recipes=792,X_branches=1584,X_native=270,X_L64_finite=0,R_banks=0,R_recipes=0,R_branches=0,native_quadruples=0,clean_quadruples=0,completed_disks=0,fewer_than_ten_ends=0,unresolved=292)
    print('PASS_JOINT44_65_AND63_FLIP_SAVED',dict(old_supports=199,debit=[19,12,1],on8=3,first=dict(ca.totals),joins=348,selected=21,flip=dict(cb.totals),dotted_X64=270,provisional=[457,292]))


if __name__=='__main__':run()
