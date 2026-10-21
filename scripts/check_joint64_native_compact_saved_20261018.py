#!/usr/bin/env python3
"""Independent moving64 saved-domain/native/incidence checks; numpy, no project imports.

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
TAG='fresh_joint64_native_compact_closure_20261018'
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


def topology(t):
    fs=t['triangles'];bc=t['boundary']
    expected=[[8,10,67],[10,67,70],[67,70,69],[70,71,69],[71,68,69],[8,68,69],[8,67,69],[8,52,68],[52,71,68],[52,70,71]]
    assert fs==expected and bc==[8,10,70,52]
    es=Counter(tuple(sorted(e)) for f in fs for e in combinations(f,2))
    assert len(set(sum(fs,[])))==t['vertices']==8 and len(es)==t['edges']==17 and len(fs)==t['faces']==10
    assert (67,68) not in es and (52,69) not in es
    assert {e for e,n in es.items() if n==1}=={tuple(sorted(e)) for e in zip(bc,bc[1:]+bc[:1])} and set(es.values())=={1,2}
    assert [r['vertex'] for r in t['links']]==sorted(set(sum(fs,[])))
    for r in t['links']:
        v=r['vertex'];ee=[list(i for i in f if i!=v) for f in fs if v in f];assert ee==r['edges']
        ds=Counter(i for e in ee for i in e);seen={next(iter(ds))}
        while True:
            ss=seen|{i for e in ee if seen.intersection(e) for i in e}
            if ss==seen:break
            seen=ss
        assert seen==set(ds) and sorted(ds.values())==([1,1]+[2]*(len(ds)-2) if v in bc else [2]*len(ds))


def seed(c,check,frame,pi,ms,branch,v,nd,tt):
    assert check['frame']==frame and check['profile']==pi and check['labels']==ms and check['branch']==branch
    assert check['norm']['status']=='two_real_norm_branches' and check['max_saved_difference']<2e-6
    vals,_,_=c.native(v,pi,nd,tt)
    assert abs(v@M@v-c.pp[pi]['square'])<2e-5
    assert all(relation(vals[i])==m for i,m in zip(frame,ms))


def run():
    c=Check(TAG);d=c.d
    a=json.loads((ROOT/f'results/{A}.json').read_text());b=json.loads((ROOT/f'results/{B}.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_JOINT64_NATIVE_COMPACT_CLOSURE' and d['deleted_walls']==[30,44,63,64,65]
    assert d['replacement64_id']==71 and d['replacement30_id']==70 and d['retained_old_walls']==58
    assert all(i not in c.old for i in (67,68,69,70,71))
    debit(d,pri)
    assert {k:v for k,v in d['exact_deletion_debit'].items() if k!='deleted_walls'}=={k:v for k,v in b['exact_deletion_debit'].items() if k!='deleted_walls'}
    assert d['old_frame_interface']==b['old_frame_interface'];topology(d['topology'])
    assert d['deletion']['score'] is not None and d['deletion']['finite_rays']==13 and d['deletion']['open_ends']==11
    cc=d['original64_recipe'];assert cc['frame']==[52,65,63] and cc['profile']==0 and cc['labels']==[2,2,27]
    assert cc['norm']['status']=='two_real_norm_branches' and sorted(r['branch'] for r in cc['branches'])==[-1,1]
    restored=0
    for r in cc['branches']:
        v=np.array(r['normal']);nd={i:w for i,w in c.old.items() if i!=64};vals,bad,cp=c.native(v,0,nd,{})
        assert set(dict(r['native_products']))==set(nd) and len(r['native_products'])==len(nd)
        assert max(abs(vals[i]-x) for i,x in r['native_products'])<2e-7
        assert [i for i,_ in r['blockers']]==bad and all(abs(vals[i]-x)<2e-7 for i,x in r['blockers'])
        assert abs(cp-r['center_product'])<2e-7 and abs(v@M@v-c.pp[0]['square'])<2e-5
        assert all(relation(vals[i])==m for i,m in zip(cc['frame'],cc['labels']))
        assert r['restores_old64']==bool(np.max(abs(v-c.old[64]))<2e-6)
        restored+=r['restores_old64'] and not bad and cp<-1e-8
    assert restored==1
    classes={}
    for rec in b['records']:
        if a['J_caps'][rec['J']]['label10']!=3:continue
        for x in rec['bank']['caps']:
            key=(rec['L_profile'],x['profile']);value=dict(x['native_products'])[64]
            cand=dict(prefix=rec['index'],X=x['index'],L_profile=rec['L_profile'],X_profile=x['profile'],old_X64=value)
            order=lambda z:(abs(z['old_X64']+2),z['prefix'],z['X'])
            if key not in classes or order(cand)<order(classes[key]):classes[key]=cand
    selection=[classes[k] for k in sorted(classes)]
    assert d['selected']==selection and len(d['records'])==len(selection)==9
    contexts={};patterns=Counter()
    for ri,(r,sel) in enumerate(zip(d['records'],selection)):
        assert r['index']==ri and r['selection']==sel
        old=b['records'][sel['prefix']];j=a['J_caps'][old['J']];l=a['L_bank']['caps'][old['L']];x=old['bank']['caps'][sel['X']]
        assert r['J']==old['J'] and r['L']==old['L'] and len(r['reconstruction'])==3
        jv=np.array(j['normal']);lv=np.array(l['normal']);xv=np.array(x['normal'])
        seed(c,r['reconstruction'][0],[8,10,63],0,[2,3,2],j['branch'],jv,c.old,{})
        seed(c,r['reconstruction'][1],[8,52,64],l['profile'],l['labels'],l['branch'],lv,c.old,{})
        nd={**c.base,67:jv,68:lv};tt={67:0,68:l['profile']}
        seed(c,r['reconstruction'][2],[8,67,68],x['profile'],x['labels'],x['branch'],xv,nd,tt)
        nd[69]=xv;tt[69]=x['profile']
        for i in (67,68,69):
            _,bad,cp=c.native(nd[i],tt[i],{k:v for k,v in nd.items() if k!=i},tt);assert not bad and cp<-1e-8
        c.ledger(r['prefix']);assert r['prefix']['score'] is not None and r['prefix']['finite_rays']==17 and r['prefix']['open_ends']==13
        bank=r['bank'];assert bank['frame']==[52,68,69] and bank['duties']==[[0,1],[1,2]]
        c.bank(bank,nd,tt)
        for m in bank['caps']:
            nn={**nd,71:np.array(m['normal'])};types={**tt,71:m['profile']};contexts[ri,m['index']]=(nn,types)
            assert m['score'] is not None and m['finite_rays']==19 and m['open_ends']==13
            assert [f['wall'] for f in m['whole_faces'] if f['all_rays_finite']]==[8,11,68]
            assert c.finite([52,68,71],nn,types) and c.finite([68,69,71],nn,types)
            for i in types:
                _,bad,cp=c.native(nn[i],types[i],{k:v for k,v in nn.items() if k!=i},types);assert not bad and cp<-1e-8
            faces={f['wall']:f for f in m['whole_faces']}
            assert sorted(r['active'] for r in faces[8]['rays'])==sorted([[7,8,11],[7,8,12],[8,10,11],[8,10,67],[8,12,52],[8,52,68],[8,67,69],[8,68,69]])
            assert sorted(r['active'] for r in faces[68]['rays'])==sorted([[8,52,68],[8,68,69],[52,68,71],[68,69,71]])
            assert not faces[71]['all_rays_finite'] and not faces[69]['all_rays_finite']
            patterns[tuple(tuple(s) for s in m['exterior_supports'])]+=1
    assert len(contexts)==len(d['R_banks'])==63
    assert [(r['prefix'],r['M']) for r in d['R_banks']]==list(contexts)
    common=set(c.base)
    for r in d['R_banks']:
        nn,tt=contexts[r['prefix'],r['M']];bank=r['bank']
        assert bank['frame']==[10,33,67] and bank['duties']==[[0,1],[0,2]] and len(bank['recipes'])==21
        c.bank(bank,nn,tt);assert not bank['caps']
        for row in bank['recipes']:
            for br in row['branches']:
                bad={i for i,_ in br['blockers']};assert bad&set(c.base);common&=bad
    assert common=={7,29,31,32,52,55} and sorted(patterns.values())==[9,22,32]
    assert not d['accepted'] and not d['unresolved'] and len(d['boundary_ledgers'])==7
    assert c.totals==dict(branches=2950,recipes=1475,caps=63)
    assert d['summary']==dict(selected_prefixes=9,M_recipes=152,M_branches=304,M_native=63,R_banks=63,R_recipes=1323,R_branches=2646,native_five_caps=0,clean_five_caps=0,completed_disks=0,fewer_than_ten_ends=0,unresolved=0)
    print('PASS_JOINT64_NATIVE_COMPACT_SAVED',dict(old_supports=199,debit=[19,12,1],on8=3,selected=9,**dict(c.totals),clean_M=63,exterior_patterns=sorted(patterns.values()),common_R_blockers=sorted(common),unresolved=0))


if __name__=='__main__':run()

