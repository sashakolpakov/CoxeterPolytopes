#!/usr/bin/env python3
"""Independent continuous-N event saved checks; NumPy, no project imports.

Checks complete saved S/R low cubes, numerical native products and all event
polynomials/roots. Historical physical scaffold and interval claims are inherited.
"""
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_oblique_low_s_continuous_n_events_20261018'
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


def triangle(c,profiles,bonds):
    g=2*np.eye(5);g[0,1]=g[1,0]=-1
    for i,p in enumerate(profiles):g[i+2,:2]=c.profiles[p];g[:2,i+2]=c.profiles[p]
    for (i,j),m in zip(((0,1),(0,2),(1,2)),bonds):g[i+2,j+2]=g[j+2,i+2]=-2*np.cos(np.pi/m)
    return np.linalg.eigvalsh(g).min()>1e-8


def recipes(c,nd,tt,kind):
    frame=[52,68,71]if kind=='S'else[69,71,72]
    profiles=[c.pt(i,tt)for i in frame]
    pairs=[(0,1),(1,2)]
    bonds={pair:relation(float(nd[frame[pair[0]]]@M@nd[frame[pair[1]]])+c.off(profiles[pair[0]],profiles[pair[1]]))for pair in pairs}
    assert all(isinstance(x,int)for x in bonds.values())
    return [(p,list(ms))for p in range(7)for ms in product((2,3,4,5),repeat=3)
        if(kind=='S'or ms[0]==2)and all(triangle(c,[profiles[a],profiles[b],p],[bonds[a,b],ms[a],ms[b]])for a,b in pairs)]


def complete_branches(c,q,nd,tt,frame):
    assert q['norm']['status']=='two_real_norm_branches'
    assert len(q['branches'])==2 and [b['branch']for b in q['branches']]==[-1,1]
    ff=np.array([nd[i]for i in frame]);gram=ff@M@ff.T;assert abs(np.linalg.det(gram))>1e-7
    vs=[]
    for b in q['branches']:
        v,bad,cp=c.branch(b,q['profile'],nd,tt,frame,q['labels']);vs.append(v)
    # A nonsingular three-frame leaves one affine line. Distinct saved roots
    # meeting both its equations and norm constraint exhaust the norm branches.
    assert np.max(abs(vs[0]-vs[1]))>1e-5


def run():
    c=Check(TAG);d=c.d
    prior=json.loads((ROOT/'results/fresh_oblique_high_n_right_r_join_20261018.json').read_text())
    pri=json.loads((ROOT/'results/fresh_nonarithmetic_release49_62_replacement_20261017.json').read_text())
    assert d['status']=='BOUNDED_LOW_S_CONTINUOUS_HIGH_N_EVENT_SEARCH'
    assert d['deleted_walls']==[30,44,63,64,65]and d['exact_deletion_debit']==prior['exact_deletion_debit']
    assert d['numerical_snapshot']==prior['numerical_snapshot']and d['old_frame_interface']==prior['old_frame_interface']
    assert d['wall_count_breakdown']==dict(projected_originals=52,fixed_original_whites=6,fixed_new_before_S=4,fixed_new_before_R=5,fixed_new_before_N=6)
    debit(d,pri);topology(d['topology'])
    assert len(d['records'])==2 and d['N_families']==[[2,2,'m>=6'],[3,2,'m>=6']]
    assert d['parameter_domain']==[float(np.sqrt(3)),2.]
    nSrecipes=nSbranches=nScaps=nRrecipes=nRbranches=nRc=neq=nroots=0;status=Counter();hist=Counter();inrange=[]
    for row,old in zip(d['records'],prior['records']):
        for k in('index','origins','instruction_key','reconstruction','prefix_normals','prefix_profiles'):assert row[k]==old[k]
        nd={**c.base,**{int(i):np.array(v)for i,v in row['prefix_normals'].items()}}
        tt={int(i):p for i,p in row['prefix_profiles'].items()};assert set(tt)=={67,68,69,71}and not({70,72,73}&set(nd))
        for i in tt:
            _,bad,cp=c.native(nd[i],tt[i],{j:v for j,v in nd.items()if j!=i},tt);assert not bad and cp<-1e-8
        assert len(row['N_families'])==2
        for a,nf in enumerate(row['N_families']):
            assert nf['a']==a and nf['frame']==[10,67,69]
            ff=np.array([nd[i]for i in nf['frame']]);g=ff@M@ff.T;inv=np.linalg.inv(g)
            assert np.max(abs(g-nf['Gram']))<2e-7 and np.max(abs(inv-nf['inverse']))<2e-7
            u=np.array([-a,0,0])@inv@ff;v=np.array([0,0,-1])@inv@ff
            cv=ff@M@c.old[64];b=c.old[64]-cv@inv@ff;sig=float(b@M@b);b/=np.sqrt(sig)
            assert sig>0 and abs(sig-nf['old64_residual_square'])<2e-7
            for key,x in(('u',u),('v',v),('b',b)):assert np.max(abs(x-nf[key]))<2e-7
            rho=[float(2-u@M@u),float(-2*u@M@v),float(-v@M@v)]
            assert np.max(abs(np.array(rho)-nf['rho']))<2e-7
            lo,hi=d['parameter_domain'];points=[lo,hi]
            if rho[2]>0 and lo<-rho[1]/(2*rho[2])<hi:points.append(-rho[1]/(2*rho[2]))
            assert min(sum(q*t**k for k,q in enumerate(rho))for t in points)>0
        sb=row['S_bank'];assert sb['frame']==[52,68,71]
        want=recipes(c,nd,tt,'S');assert [(q['profile'],q['labels'])for q in sb['recipes']]==want
        assert [q['index']for q in sb['recipes']]==list(range(len(want)))
        nSrecipes+=len(want);digest=sha256();native=[]
        for q in sb['recipes']:
            complete_branches(c,q,nd,tt,sb['frame']);nSbranches+=2
            for b in q['branches']:
                digest.update(json.dumps([q['index'],b['branch'],q['profile'],q['labels'],b['native_products'],b['center_product']],separators=(',',':')).encode())
                if not b['blockers']and b['center_product']<-1e-8:
                    assert not b.get('duplicate',False);native.append((q,b))
        assert digest.hexdigest()==sb['complete_branch_digest']and len(native)==len(sb['caps'])==len(row['R_banks'])
        for si,((sq,sbr),sc,bank)in enumerate(zip(native,sb['caps'],row['R_banks'])):
            assert sc['index']==bank['S']==si and sc['recipe']==sq['index']and sc['branch']==sbr['branch']and sc['normal']==sbr['normal']
            nScaps+=1;full={**nd,72:np.array(sc['normal'])};ft={**tt,72:sc['profile']}
            want=recipes(c,full,ft,'R');assert bank['frame']==[69,71,72]and[(q['profile'],q['labels'])for q in bank['recipes']]==want
            assert [q['index']for q in bank['recipes']]==list(range(len(want)))
            nRrecipes+=len(want);nativeR=[]
            for q in bank['recipes']:
                complete_branches(c,q,full,ft,bank['frame']);nRbranches+=2
                for b in q['branches']:
                    hist[str(tuple(i for i,v in b['blockers']))]+=1
                    if not b['blockers']and b['center_product']<-1e-8:nativeR.append((q,b))
            assert len(nativeR)==len(bank['native_caps'])
            for(q,b),cap in zip(nativeR,bank['native_caps']):
                nRc+=1
                assert [cap[k]for k in('recipe','branch','profile','labels')]==[q['index'],b['branch'],q['profile'],q['labels']]
                rv=np.array(b['normal']);chamber={**full,70:rv};ct={**ft,70:q['profile']}
                assert c.finite([69,71,70],chamber,ct)and c.finite([71,70,72],chamber,ct)
                assert len(cap['events'])==2
                for nf,ev in zip(row['N_families'],cap['events']):
                    neq+=1;assert nf['a']==ev['a']
                    u,v,bv=(np.array(nf[k])for k in('u','v','b'));ab=[float(rv@M@w)for w in(u,v,bv)];aa,bb,cc=ab
                    assert np.max(abs(np.array(ab)-ev['coefficients']))<2e-7
                    rho=nf['rho'];pol=[aa*aa-cc*cc*rho[0],2*aa*bb-cc*cc*rho[1],bb*bb-cc*cc*rho[2]]
                    assert np.max(abs(np.array(pol)-ev['polynomial']))<2e-7
                    assert ev['root_analysis']['status']=='two_real_events'and len(ev['roots'])==2
                    expected=sorted(float(x)for x in np.roots(pol[::-1]));assert abs(expected[1]-expected[0])>1e-5
                    for er,rt in zip(expected,ev['roots']):
                        nroots+=1;theta=rt['theta'];assert abs(er-theta)<2e-7
                        inside=np.sqrt(3)-1e-9<=theta<=2+1e-9;assert inside==rt['in_closed_high_domain']
                        status[rt['status']]+=1
                        if not inside:assert rt['status']=='outside_closed_high_domain';continue
                        rem=sum(z*theta**i for i,z in enumerate(rho));nv=u+theta*v+np.sqrt(rem)*bv
                        assert rem>0 and abs(rem-rt['norm_remainder'])<2e-7 and np.max(abs(nv-rt['N_normal']))<2e-7
                        assert abs(aa+bb*theta+cc*np.sqrt(rem))<1e-6 and abs(rv@M@nv)<1e-6
                        assert rt['status']=='unsquared_match'and rt['NX_relation']is None
                        vals,bad,cp=c.native(nv,0,chamber,ct)
                        assert len(rt['native_products'])==len(chamber)and set(dict(rt['native_products']))==set(chamber)
                        assert max(abs(vals[i]-x)for i,x in rt['native_products'])<2e-7
                        assert bad==[i for i,x in rt['blockers']]and 69 in bad and abs(cp-rt['center_product'])<2e-7
                        lo,hi=rt['adjacent_integer_labels'];assert hi==lo+1 and(lo,hi)in((8,9),(10,11))
                        gaps=[theta-2*np.cos(np.pi/lo),2*np.cos(np.pi/hi)-theta]
                        assert min(gaps)>0.002 and np.max(abs(np.array(gaps)-rt['adjacent_theta_gaps']))<2e-7
                        inrange.append((row['index'],nf['a'],theta))
    s=d['summary'];assert(nSrecipes,nSbranches,nScaps,nRrecipes,nRbranches,nRc,neq,nroots)==(84,168,76,1238,2476,60,120,240)
    assert s['parents']==2 and s['parent_origins']==4
    for k,v in zip(('S_recipes','S_branches','native_S','R_recipes','R_branches','native_R','event_equations','event_roots'),(nSrecipes,nSbranches,nScaps,nRrecipes,nRbranches,nRc,neq,nroots)):assert s[k]==v
    assert s['event_status']==dict(status)=={'outside_closed_high_domain':236,'unsquared_match':4}and s['R_blocker_histogram']==dict(hist)
    assert s['accepted']==s['new_cones']==s['unresolved']==0 and not d['accepted']and not d['boundary_ledgers']and not d['unresolved']
    assert len(inrange)==4 and{(p,a)for p,a,t in inrange}==set(product(range(2),repeat=2))
    print('PASS_SAVED_CONTINUOUS_N_EVENTS','84/1238complete_recipes','2644norm_branches','120quadratics','240roots','4nonnative_NX_events','19/12/1debit')


if __name__=='__main__':run()

