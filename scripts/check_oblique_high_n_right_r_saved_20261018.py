#!/usr/bin/env python3
"""Independent high-N/right-R saved-data checks; NumPy, no project imports.

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
TAG='fresh_oblique_high_n_right_r_join_20261018'
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
                assert c['native_products']==b['native_products'] and c['center_product']==b['center_product']
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
    c=Check(TAG);d=c.d
    prev=json.loads((ROOT/f'results/{PREV}.json').read_text())
    forward=json.loads((ROOT/'results/fresh_oblique_left_boundary_connector_20261018.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_HIGH_N_RIGHT_R_JOIN'and len(d['records'])==2
    assert d['deleted_walls']==[30,44,63,64,65]and d['high_labels']==[6,7,8,10,12]
    assert d['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_N=4,fixed_new_before_R_S=5,explicit_N_gates=56,explicit_R_or_S_gates=57)
    assert d['necessary_R_rights']==dict(NX_squared_lower_bound=3,other_nonright_squared_lower_bound=1,determinant_upper_bound=0,RX=2,RN=2)
    assert 8-2*(3+1)==0
    debit(d,pri);topology(d['topology'])
    assert d['old_frame_interface']==prev['old_frame_interface']and d['exact_deletion_debit']==prev['exact_deletion_debit']
    origins=[o['index']for r in d['records']for o in r['origins']];assert origins==[0,1,2,3]
    relations=Counter();paircount=0;diagnostics=0;bankcounts=Counter();allcaps=0
    for ri,r in enumerate(d['records']):
        assert r['index']==ri and len(r['origins'])==[3,1][ri]
        for o in r['origins']:
            src=prev['records'][o['index']]
            assert o['forward_prefix']==src['selection']['prefix']
            key=json.dumps([[z['wall'],z['frame'],z['profile'],z['labels'],z['branch']]for z in src['reconstruction'][:5]],separators=(',',':'))
            assert key==r['instruction_key']
            for i in(67,68,69,71):assert max(abs(np.array(r['prefix_normals'][str(i)])-np.array(src['prefix_normals'][str(i)])))<2e-6
        src=prev['records'][r['origins'][0]['index']];prior=forward['records'][src['selection']['prefix']]
        instructions=[[q['wall'],q['frame'],q['profile'],q['labels'],q['branch']]for q in src['reconstruction']]
        assert len(r['reconstruction'])==len(instructions)==8 and set(r['prefix_normals'])==set(r['prefix_profiles'])==set(map(str,[67,68,69,71]))
        assert set(r['removed'])=={'70','72','73'};base=dict(c.base);bt={}
        for j,(st,ch)in enumerate(zip(instructions,r['reconstruction'])):
            if j==6:base.pop(70);bt.pop(70)
            i,frame,pi,ms,branch=st
            saved=prior['old_R_normal']if j==3 else src['former_R']if j==6 else src['prefix_normals'][str(i)]
            v=np.array(prior['old_R_normal']if j==3 else r['removed'][str(i)]['normal']if i in(70,72,73)else r['prefix_normals'][str(i)])
            assert[ch[k]for k in('wall','frame','profile','labels','branch')]==[i,frame,pi,ms,branch]
            assert ch['norm']['status']=='two_real_norm_branches'and ch['max_saved_difference']<2e-6
            assert max(abs(v-np.array(saved)))<2e-6 and abs(v@M@v-c.pp[pi]['square'])<2e-5
            vals,_,_=c.native(v,pi,c.old if j<2 else base,{}if j<2 else bt)
            assert all(relation(vals[w])==m for w,m in zip(frame,ms))
            base[i]=v;bt[i]=pi
        for i in(70,72,73):assert r['removed'][str(i)]['profile']==bt.pop(i);base.pop(i)
        assert{str(i):p for i,p in bt.items()}==r['prefix_profiles']
        for i in bt:
            _,bad,cp=c.native(base[i],bt[i],{j:v for j,v in base.items()if j!=i},bt);assert not bad and cp<-1e-8
        c.ledger(r['before']);assert r['frame']==[10,67,69]and r['duties']==[[0,1],[1,2]]
        assert r['frame_profiles']==[c.pt(i,bt)for i in r['frame']]
        for a,b in r['duties']:assert relation(float(base[r['frame'][a]]@M@base[r['frame'][b]])+c.off(c.pt(r['frame'][a],bt),c.pt(r['frame'][b],bt)))==r['duty_bonds'][3*a+b]
        assert sorted(r['symbolic_high_families'])==[[0,[2,2,6],[2]],[0,[3,2,6],[2]]]
        G=np.array([[float(base[i]@M@base[j])for j in r['frame']]for i in r['frame']]);inv=np.linalg.inv(G)
        assert np.max(abs(G-r['projected_Gram']))<2e-7 and np.max(abs(inv-r['projected_inverse']))<2e-7
        eig=np.linalg.eigvalsh(G);assert np.max(abs(eig-r['eigenvalues']))<2e-7 and sum(eig>1e-8)==2 and sum(eig< -1e-8)==1
        assert len(r['recipes'])==len(r['caps'])==10
        assert[(q['profile'],q['labels'])for q in r['recipes']]==[(0,[a,2,m])for a,m in product((3,2),d['high_labels'])]
        digest=sha256();capidx=0
        for ni,nr in enumerate(r['recipes']):
            assert nr['index']==ni and nr['norm']['status']=='two_real_norm_branches'and sorted(b['branch']for b in nr['branches'])==[-1,1]
            theta=2*np.cos(np.pi/nr['labels'][2]);b=np.array([-2*np.cos(np.pi/nr['labels'][0]),0.,0.]);v=np.array([0.,0.,-1.])
            coeff=[float(2-b@inv@b),float(-2*b@inv@v),float(-v@inv@v)]
            assert np.max(abs(np.array(coeff)-nr['residual_square_polynomial']))<2e-7 and abs(nr['theta']-theta)<2e-7
            assert abs(nr['residual_square']-sum(a*theta**k for k,a in enumerate(coeff)))<2e-7 and nr['residual_square']>0
            for br in nr['branches']:
                nv,bad,cp=c.branch(br,0,base,bt,r['frame'],nr['labels'])
                digest.update(json.dumps([ni,br['branch'],0,nr['labels'],br['native_products'],br['center_product']],separators=(',',':')).encode())
                if bad or cp>=-1e-8:
                    assert bad==[8]and br['branch']==1 and 'cap'not in br;continue
                assert br['branch']==-1 and br['cap']==capidx
                cap=r['caps'][capidx];assert cap['index']==capidx and cap['recipe']==ni and cap['normal']==br['normal']and cap['labels']==nr['labels']and cap['profile']==0;capidx+=1;allcaps+=1
                nd={**base,73:nv};tt={**bt,73:0};c.ledger(cap['before'])
                assert all(c.finite([r['frame'][a],r['frame'][b],73],nd,tt)for a,b in r['duties'])
                assert cap['R_bank']['frame']==[69,71,73]and cap['R_bank']['duties']==[[0,1],[0,2]]
                assert cap['S_bank']['frame']==[52,68,71]and cap['S_bank']['duties']==[[0,1],[1,2]]
                for k in('R_bank','S_bank'):
                    bk=cap[k];c.bank(bk,nd,tt);bankcounts[k+'_recipes']+=len(bk['recipes']);bankcounts[k+'_branches']+=sum(len(q['branches'])for q in bk['recipes']);bankcounts[k+'_caps']+=len(bk['caps'])
                assert all(q['labels'][0]==q['labels'][2]==2 for q in cap['R_bank']['recipes'])
                aa=cap['R_bank']['caps'];bb=cap['S_bank']['caps']
                assert len(cap['pairs'])==len(aa)*len(bb)and[(q['R'],q['S'])for q in cap['pairs']]==list(product(range(len(aa)),range(len(bb))))
                pdigest=sha256()
                for q in cap['pairs']:
                    a,b=aa[q['R']],bb[q['S']];value=float(np.array(a['normal'])@M@np.array(b['normal']))+c.off(a['profile'],b['profile'])
                    assert abs(value-q['product'])<2e-7 and relation(value)==q['relation']
                    full={**nd,70:np.array(a['normal']),72:np.array(b['normal'])};ft={**tt,70:a['profile'],72:b['profile']}
                    assert q['finite_MRS']is None and not c.finite([71,70,72],full,ft)
                    relations[str(q['relation'])]+=1;paircount+=1;pdigest.update(json.dumps(q,separators=(',',':')).encode())
                assert pdigest.hexdigest()==cap['complete_pair_digest']and not cap['accepted']
                native=[q for q in cap['pairs']if q['relation']=='dotted'or isinstance(q['relation'],int)]
                best=max(native,key=lambda p:(p['product'],-p['R'],-p['S']));z=cap['closest_native_diagnostic']
                assert(z['R'],z['S'])==(best['R'],best['S'])and z['diagnostic_only']
                a,b=aa[z['R']],bb[z['S']];full={**nd,70:np.array(a['normal']),72:np.array(b['normal'])};ft={**tt,70:a['profile'],72:b['profile']}
                for i in ft:
                    _,bad,cp=c.native(full[i],ft[i],{j:v for j,v in full.items()if j!=i},ft);assert not bad and cp<-1e-8
                c.ledger(z);assert(z['finite_rays'],z['open_ends'])==(24,14)
                cs=[c.finite(f,full,ft)for f in d['topology']['triangles']]
                assert[v is not None for v in z['corner_types']]==cs and not z['all_disk_corners']
                assert[f for f,yes in zip(d['topology']['triangles'],cs)if not yes]==[[71,70,72]]
                faces={f['wall']:f for f in z['whole_faces']};assert set(faces)=={7,8,10,11,33,52,55,66,67,68,69,70,71,72,73}
                cc={str(w):faces[w]['all_rays_finite']for w in faces};assert z['whole_compact_faces']==cc
                exact={str(w):faces[w]['all_rays_finite']and sorted(sorted(q['active'])for q in faces[w]['rays'])==sorted(sorted(f)for f in d['topology']['triangles']if w in f)for w in(67,68,69,71)}
                assert z['exact_whole_interior_caps']==exact and{k for k,v in exact.items()if v}=={'67','68','69'};diagnostics+=1
        assert capidx==10 and digest.hexdigest()==r['complete_branch_digest']
    assert allcaps==20 and paircount==6104 and relations=={'dotted':6088,'None':16}and diagnostics==20
    assert d['summary']==dict(parent_origins=4,distinct_parents=2,N_recipes=20,N_branches=40,native_N=20,R_recipes=360,S_recipes=840,R_branches=720,S_branches=1680,R_caps=172,S_caps=712,pairs=6104,pair_relations=dict(relations),finite_MRS=0,completed_disks=0,diagnostic_cones=20,unresolved=0)
    assert dict(c.totals)==dict(branches=2440,recipes=1200,caps=884)and not d['accepted']and not d['unresolved']
    assert bankcounts=={'R_bank_recipes':360,'S_bank_recipes':840,'R_bank_branches':720,'S_bank_branches':1680,'R_bank_caps':172,'S_bank_caps':712}
    print('PASS_SAVED_HIGH_N_RIGHT_R',dict(c.totals),'20N','6104pairs','20diagnostics','19/12/1debit',flush=True)


if __name__=='__main__':run()
