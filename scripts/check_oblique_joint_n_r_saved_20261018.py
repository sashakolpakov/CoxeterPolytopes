#!/usr/bin/env python3
"""Independent joint N/R pilot saved-data checks; NumPy, no project imports.

Does not solve norm equations or reconstruct geometric cones. The frozen
S packet supplies its previously verified full recipe-chain provenance.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_joint_n_r_pilot_20261018'
PREV='fresh_oblique_left_boundary_connector_20261018'
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
        self.totals=Counter();self.ledger(d['control'])
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
        assert [(r['profile'],r['labels']) for r in bank['recipes']]==domain
        digest=sha256();hist=Counter();seen=set();ci=0
        for ri,r in enumerate(bank['recipes']):
            assert r['index']==ri
            if r['norm']['status']=='no_real_norm_root':
                assert not r['branches'];self.totals['no_real_recipes']+=1;continue
            assert r['norm']['status']=='two_real_norm_branches' and sorted(b['branch'] for b in r['branches'])==[-1,1]
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
    c=Check(TAG);d=c.d;old=json.loads((ROOT/f'results/{PREV}.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_OBLIQUE_JOINT_N_R_PILOT' and d['new_id']==70
    assert d['deleted_walls']==[30,44,63,64,65] and d['retained_old_walls']==58
    assert d['old_frame_interface']==old['old_frame_interface'] and d['exact_deletion_debit']==old['exact_deletion_debit']
    assert d['wall_count_breakdown']==old['wall_count_breakdown']
    debit(d,pri);topology(d['topology']);candidates=[]
    assert len(old['records'])==376
    for row in old['records']:
        for q in row['recipes']:
            for b in q['branches']:
                if [i for i,_ in b['blockers']]==[70] and b['center_product']<-1e-8:
                    candidates.append(dict(prefix=row['index'],recipe=q['index'],branch=b['branch'],profile=q['profile'],
                        labels=q['labels'],normal=b['normal'],former_NR=b['blockers'][0][1]))
    assert len(candidates)==1100;key=lambda p:(abs(p['former_NR']),p['prefix'],p['recipe'],p['branch'])
    selected=[];used=set()
    for label in (3,4,5):
        p=min((p for p in candidates if p['labels'][2]==label and p['prefix'] not in used),key=key)
        selected.append(p);used.add(p['prefix'])
    selected.append(min((p for p in candidates if p['prefix'] not in used),key=key))
    assert selected==d['selections'] and [r['prefix']for r in selected]==[140,133,166,348]
    assert len(d['records'])==4;expected=[];patterns=Counter()
    for ri,(sel,r)in enumerate(zip(selected,d['records'])):
        assert r['index']==ri and r['selection']==sel
        prior=old['records'][sel['prefix']]
        instructions=[[q['wall'],q['frame'],q['profile'],q['labels'],q['branch']] for q in prior['reconstruction']]
        instructions.append([73,[10,67,69],sel['profile'],sel['labels'],sel['branch']])
        normals=[prior['old_R_normal']if j==3 else prior['prefix_normals'][str(st[0])]for j,st in enumerate(instructions[:-1])]+[sel['normal']]
        assert len(r['reconstruction'])==len(instructions)==8
        assert set(r['prefix_normals'])==set(r['prefix_profiles'])==set(map(str,[67,68,69,71,72,73]))
        nd=dict(c.base);tt={}
        for j,(step,saved,ch)in enumerate(zip(instructions,normals,r['reconstruction'])):
            if j==6:nd.pop(70);tt.pop(70)
            i,frame,pi,ms,branch=step
            v=np.array(prior['old_R_normal']if j==3 else r['former_R']if j==6 else r['prefix_normals'][str(i)])
            assert (prior['old_R_profile']if j==3 else r['former_R_profile']if j==6 else r['prefix_profiles'][str(i)])==pi
            assert ch['wall']==i and[ch[k]for k in('frame','profile','labels','branch')]==[frame,pi,ms,branch]
            assert ch['norm']['status']=='two_real_norm_branches'and ch['max_saved_difference']<2e-6
            assert max(abs(v-np.array(saved)))<2e-6 and abs(v@M@v-c.pp[pi]['square'])<2e-5
            vals,_,_=c.native(v,pi,c.old if j<2 else nd,{}if j<2 else tt)
            assert all(relation(vals[w])==m for w,m in zip(frame,ms))
            nd[i]=v;tt[i]=pi
        nr=float(nd[73]@M@nd[70])+c.off(tt[73],tt[70]);assert abs(nr-sel['former_NR'])<2e-7
        nd.pop(70);tt.pop(70)
        for i in tt:
            _,bad,cp=c.native(nd[i],tt[i],{k:v for k,v in nd.items()if k!=i},tt);assert not bad and cp<-1e-8
        c.ledger(r['before']);assert(r['before']['finite_rays'],r['before']['open_ends'])==((23,13)if ri==0 else(22,14))
        assert r['frame']==[69,71,73]and r['duties']==[[0,1],[0,2]]
        c.bank(r,nd,tt)
        for cap in r['caps']:
            full={**nd,70:np.array(cap['normal'])};ft={**tt,70:cap['profile']}
            for i in ft:
                _,bad,cp=c.native(full[i],ft[i],{k:v for k,v in full.items()if k!=i},ft);assert not bad and cp<-1e-8
            cs=[c.finite(f,full,ft)for f in d['topology']['triangles']]
            assert [z is not None for z in cap['corner_types']]==cs and cap['all_disk_corners']==all(cs)
            assert(cap['deferred_M_R_S']is not None)==c.finite([71,70,72],full,ft)
            faces={f['wall']:f for f in cap['whole_faces']};assert set(faces)=={7,8,10,11,33,52,55,66,67,68,69,70,71,72,73}
            cc={str(w):faces[w]['all_rays_finite']for w in faces};assert cap['whole_compact_faces']==cc
            caps={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in d['topology']['triangles']if w in f)for w in(67,68,69,71)}
            assert cap['exact_whole_interior_caps']==caps
            assert cap['score']is not None and(cap['finite_rays'],cap['open_ends'])==(24,14)
            assert[f for f,z in zip(d['topology']['triangles'],cs)if not z]==[[71,70,72]]
            assert{k for k,v in caps.items()if v}=={'67','68','69'}and{k for k,v in cc.items()if v}=={'8','11','67','68','69'}
            assert relation(dict(cap['native_products'])[72])=='dotted'
            patterns[tuple(tuple(z)for z in cap['exterior_supports'])]+=1;expected.append(dict(prefix=ri,**cap))
    assert expected==d['accepted']and len(expected)==80 and sorted(patterns.values())==[34,46]
    assert d['summary']==dict(forward_sole_R_branches=1100,tested_prefixes=4,recipes=140,branches=264,native_R=80,clean_R=80,completed_disks=0,fewer_than_ten_ends=0,unresolved=0)
    assert dict(c.totals)==dict(branches=264,no_real_recipes=8,recipes=140,caps=80)and not d['unresolved']
    print('PASS_SAVED_JOINT_N_R',dict(c.totals),'4prefixes','80clean24finite14end','19/12/1debit','two_exterior_patterns',flush=True)


if __name__=='__main__':run()
