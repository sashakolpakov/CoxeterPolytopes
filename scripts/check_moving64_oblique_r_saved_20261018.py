#!/usr/bin/env python3
"""Independent changed-R/reversed-M saved checks; NumPy, no project imports.

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
TAG='fresh_moving64_oblique_r_closure_20261018'
PREV='fresh_joint64_native_compact_closure_20261018'
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
    old=json.loads((ROOT/f'results/{PREV}.json').read_text())
    a=json.loads((ROOT/f'results/{A}.json').read_text());b=json.loads((ROOT/f'results/{B}.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_MOVING64_OBLIQUE_R_CLOSURE'
    assert d['deleted_walls']==[30,44,63,64,65] and d['retained_old_walls']==58
    assert d['replacement30_id']==70 and d['replacement64_id']==71
    assert d['exact_deletion_debit']==old['exact_deletion_debit'] and d['old_frame_interface']==old['old_frame_interface']
    assert d['prior_selected']==old['selected'] and len(old['records'])==9
    debit(d,pri);topology(d['topology'])
    assert d['deletion']['score'] is not None and (d['deletion']['finite_rays'],d['deletion']['open_ends'])==(13,11)
    contexts={};expected=[]
    for pr in old['records']:
        sel=pr['selection'];fr=b['records'][sel['prefix']]
        j=a['J_caps'][fr['J']];l=a['L_bank']['caps'][fr['L']];x=fr['bank']['caps'][sel['X']]
        nd={**c.base,67:np.array(j['normal']),68:np.array(l['normal']),69:np.array(x['normal'])}
        tt={67:0,68:l['profile'],69:x['profile']}
        for m in pr['bank']['caps']:
            expected.append((pr['index'],m['index']))
            contexts[pr['index'],m['index']]=(nd,tt,j,l,x,m)
    assert len(expected)==len(d['records'])==63 and [(r['prefix'],r['M'])for r in d['records']]==expected
    selected={};duplicates=[];common=set(c.base);only_m=[]
    for ri,r in enumerate(d['records']):
        assert r['index']==ri and r['selection']==old['records'][r['prefix']]['selection']
        nd,tt,j,l,x,m=contexts[r['prefix'],r['M']];checks=r['reconstruction'];assert len(checks)==4
        seed(c,checks[0],[8,10,63],0,[2,3,2],j['branch'],nd[67],c.old,{})
        seed(c,checks[1],[8,52,64],l['profile'],l['labels'],l['branch'],nd[68],c.old,{})
        seed(c,checks[2],[8,67,68],x['profile'],x['labels'],x['branch'],nd[69],{k:v for k,v in nd.items()if k!=69},{k:v for k,v in tt.items()if k!=69})
        mv=np.array(m['normal']);seed(c,checks[3],[52,68,69],m['profile'],m['labels'],m['branch'],mv,nd,tt)
        nn={**nd,71:mv};types={**tt,71:m['profile']}
        for i in types:
            _,bad,cp=c.native(nn[i],types[i],{k:v for k,v in nn.items()if k!=i},types);assert not bad and cp<-1e-8
        c.ledger(r['before']);assert r['before']['score'] is not None and (r['before']['finite_rays'],r['before']['open_ends'])==(19,13)
        bank=r['bank'];assert bank['frame']==[10,67,69] and bank['duties']==[[0,1],[1,2]]
        c.bank(bank,nn,types);assert not bank['caps']
        for q in bank['recipes']:
            for br in q['branches']:
                bad=[i for i,_ in br['blockers']];common&=set(bad)
                if bad==[8]:
                    diff=float(np.max(abs(np.array(br['normal'])-c.old[8])))
                    duplicates.append(dict(record=ri,recipe=q['index'],branch=br['branch'],profile=q['profile'],
                        old8_profile=c.pt(8,{}),max_projected_difference=diff,same_full_profile=q['profile']==c.pt(8,{}),
                        duplicates_old8=q['profile']==c.pt(8,{})and diff<2e-6))
                if bad==[71]and br['center_product']<-1e-8:
                    only_m.append((ri,q['index'],br['branch']))
                    key=(r['prefix'],q['profile'],*np.round(br['normal'],8))
                    if key not in selected:
                        selected[key]=dict(record=ri,recipe=q['index'],branch=br['branch'],prefix=r['prefix'],
                            profile=q['profile'],labels=q['labels'],normal=br['normal'],removed_M_product=dict(br['native_products'])[71])
    assert d['old8_duplicate_checks']==duplicates and len(duplicates)==123 and sum(r['duplicates_old8']for r in duplicates)==63
    assert not common and len(only_m)==51 and len(selected)==11
    assert len(d['reverse_records'])==11 and [r['selection']for r in d['reverse_records']]==list(selected.values())
    completions=[];patterns=Counter();r52=[];m52=[]
    for ri,r in enumerate(d['reverse_records']):
        assert r['index']==ri;sel=r['selection'];source=d['records'][sel['record']]
        nd,tt,*_=contexts[source['prefix'],source['M']]
        nn={**nd,70:np.array(sel['normal'])};types={**tt,70:sel['profile']}
        for i in types:
            _,bad,cp=c.native(nn[i],types[i],{k:v for k,v in nn.items()if k!=i},types);assert not bad and cp<-1e-8
        c.ledger(r['before']);assert r['before']['score'] is not None and (r['before']['finite_rays'],r['before']['open_ends'])==(19,13)
        assert [f['wall']for f in r['before']['whole_faces']if f['all_rays_finite']]==[8,11,67]
        bank=r['bank'];assert bank['frame']==[68,69,70]and bank['duties']==[[0,1],[1,2]]
        c.bank(bank,nn,types)
        for m in bank['caps']:
            full={**nn,71:np.array(m['normal'])};ft={**types,71:m['profile']}
            faces={f['wall']:f for f in m['whole_faces']}
            comps={str(w):faces[w]['all_rays_finite']for w in (8,10,11,33,52,55,66,67,68,69,70,71)}
            capsets={str(w):faces[w]['all_rays_finite']and sorted(sorted(z['active'])for z in faces[w]['rays'])==sorted(sorted(f)for f in d['topology']['triangles']if w in f)for w in (67,68,69,71)}
            corners=[c.finite(f,full,ft)for f in d['topology']['triangles']]
            extra=[[7,55,70],[55,59,70],[10,33,70]];ec=[c.finite(f,full,ft)for f in extra]
            vals,_,_=c.native(full[70],ft[70],{k:v for k,v in full.items()if k!=70},ft)
            completions.append((ri,m['index'],m,corners,extra,ec,comps,capsets,vals[33]))
            assert m['score'] is not None and (m['finite_rays'],m['open_ends'])==(21,13)
            assert [f['wall']for f in m['whole_faces']if f['all_rays_finite']]==[8,11,67,69]
            assert [f for f,v in zip(d['topology']['triangles'],corners)if not v]==[[52,71,68],[52,70,71]]
            assert capsets['67']and capsets['69']and not capsets['68']and not capsets['71']
            patterns[tuple(tuple(s)for s in m['exterior_supports'])]+=1
            r52.append(vals[52]);m52.append(dict(m['native_products'])[52])
    assert len(completions)==len(d['reverse_accepted'])==52
    for z,a in zip(completions,d['reverse_accepted']):
        ri,mi,m,cs,extra,ec,compact,capsets,p33=z
        assert a['reverse']==ri and a['M']==mi and [v is not None for v in a['corner_types']]==cs
        assert a['diagnostic_external_duties']==extra and [v is not None for v in a['diagnostic_external_corner_types']]==ec
        assert a['whole_compact_faces']==compact and a['exact_whole_new_caps']==capsets and a['whole_J_closed']
        assert not a['all_disk_corners']and a['clean']and abs(a['old33_product']-p33)<2e-7
        for k in ('score','finite_rays','open_ends','exterior_supports','open_edges','whole_faces','cone_diagnostic'):assert a[k]==m[k]
    assert max(r52)<-2-2e-7 and max(m52)<-2-2e-7 and len(patterns)==1
    assert not d['accepted']and not d['unresolved']and len(d['boundary_ledgers'])==7
    assert c.totals==dict(branches=1442,recipes=721,caps=52)
    assert d['summary']==dict(frozen_selected_JLX_prefixes=9,reconstructed_M_prefixes=63,R_recipes=647,R_branches=1294,
        native_R=0,clean_R=0,whole_J_closed=0,completed_disks=0,fewer_than_ten_ends=0,common_retained_old_blockers=[],
        sole_old8_branches=123,old8_duplicates=63,reverse_prefixes=11,reverse_M_recipes=74,reverse_M_branches=148,
        reverse_native_M=52,reverse_clean_M=52,reverse_completed_disks=0,reverse_fewer_than_ten_ends=0,unresolved=0)
    print('PASS_MOVING64_OBLIQUE_R_SAVED',dict(old_supports=199,debit=[19,12,1],on8=3,original_M=63,
        **dict(c.totals),only_M_failures=51,reverse_prefixes=11,clean_joint=52,finite=21,ends=13,
        R52_range=[min(r52),max(r52)],M52_range=[min(m52),max(m52)],old8_duplicates=63),flush=True)


if __name__=='__main__':run()

