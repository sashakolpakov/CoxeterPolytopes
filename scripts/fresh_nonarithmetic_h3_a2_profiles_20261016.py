#!/usr/bin/env python3
"""Complete finite A2 contact profiles, with a bounded H3 endpoint scout."""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import sympy as sp
from scipy.linalg import null_space
import fresh_nonarithmetic_h3_two_contact_20261016 as C
import fresh_nonarithmetic_h3_two_contact_20261016_pairs as B
import fresh_nonarithmetic_h3_edge_completion_20261016 as L
ROOT=C.ROOT;ETA=C.ETA;TAG='fresh_nonarithmetic_h3_a2_profiles_20261016'

def setup():
    data=json.loads((ROOT/f'results/{C.TAG}.json').read_text());NC=np.array(data['old_projected_normals'])
    Q,N,H,E=C.P.geometry();old=C.P.P.B.read(C.H3.CAP);new=C.P.P.B.read(C.CAPS)
    N.append(mp.matrix([C.P.mid(x)for x in old['normal_interval']]))
    N.extend(mp.matrix([C.P.mid(x)for x in r['normal_interval']])for r in new['caps']);assert len(N)==55
    basis=[mp.matrix([mp.mpf(x)for x in v])for v in data['basis']];tau=-(H.T*Q*basis[0])[0]
    A=mp.matrix([[N[j][i]for j in(5,6)]for i in range(10)]);J=(A.T*Q*A)**-1
    assert abs(J[0,0]-mp.mpf(2)/3)<mp.mpf('1e-40')and abs(J[0,1]-mp.mpf(1)/3)<mp.mpf('1e-40')
    profiles=[]
    for node in(5,6):
        for label in(3,4,5):
            c=mp.matrix([-2*mp.cos(mp.pi/label)if k==node else 0 for k in(5,6)])
            p=A*J*c;q=2-(p.T*Q*p)[0];offset=np.array([float((p.T*Q*n)[0])for n in N])
            profiles.append(dict(node=node,label=label,residual_square=float(q),native_offsets=offset,
                center_lower_bound=float((p.T*Q*H)[0]/tau),physical_projection=[mp.nstr(x,65)for x in p]))
    rays=[dict(mask=sum(1<<i for i in r['active']),v=np.array(r['point']),kind=r['kind'])for r in data['actual_face_rays']]
    return data,NC,profiles,rays

def algebra():
    phi=(1+sp.sqrt(5))/2;J=sp.Matrix([[2,1],[1,2]])/3;rows=[]
    H4=sp.Matrix([[2,-phi,0,0],[-phi,2,-1,0],[0,-1,2,-1],[0,0,-1,2]])
    H4i=sp.simplify(H4.inv());assert all(z>0 for z in H4i)
    assert all(sp.simplify(H4i[i,i]-2)>0 for i in range(4))
    for label,c in[(3,sp.Integer(1)),(4,sp.sqrt(2)),(5,phi)]:
        for node in(5,6):
            v=sp.Matrix([c if node==5 else 0,c if node==6 else 0]);q=sp.simplify(2-(v.T*J*v)[0]);assert q>0
            rows.append(dict(node=node,label=label,residual_square=str(q)))
    assert (sp.Matrix([1,1]).T*J*sp.Matrix([1,1]))[0]==2
    assert sp.simplify(2-(sp.Matrix([sp.sqrt(3),0]).T*J*sp.Matrix([sp.sqrt(3),0]))[0])==0
    inv=json.loads((ROOT/'results/fresh_h9_persistent_cusp_screen_20261016.json').read_text())['pair_matrix']
    assert [(i,j,inv[i][j])for i,j in itertools.combinations(C.S,2)if inv[i][j]!=2]==[(0,1,5),(1,2,3),(2,3,3),(5,6,3)]
    return dict(support_components='H4+A2',H4_nodes=[0,1,2,3],H4_inverse=[[str(sp.simplify(H4i[i,j]))for j in range(4)]for i in range(4)],
        H4_inverse_diagonal=[str(sp.simplify(H4i[i,i]))for i in range(4)],A2_nodes=[5,6],A2_inverse=[['2/3','1/3'],['1/3','2/3']],
        nonzero_profiles=rows,zero_profile_residual_square='2',two_positive_contacts_minimum_cost='2',single_label6_minimum_cost='2',
        proof='An actual additional facet of the S6 face has positive projected norm. The H4 inverse is positive and each diagonal exceeds2, so any nonzero Coxeter contact of magnitude at least1 already has projection cost>2. Thus every facet is right to H4, without a classification import. Its two A2 products are zero or -2cos(pi/m); parallel/dotted products also have magnitude>=2. The positive inverse gives cost >=2 if both are nonzero, and >=2 for a singleton label>=6. Thus precisely the zero profile and the six singleton labels3,4,5 remain. All products here use norm-two outward normals.')

def run(maxlabel):
    data,N,profiles,rays=setup();edges,ends,sides=L.ledger(rays,N)
    pairs=sorted({tuple(e['sides'])for e in ends});points=np.array([e['point']for e in ends]);allrows=[];totals=Counter()
    for pr in profiles:
        q=pr['residual_square'];off=pr['native_offsets'];rows=[];counts=Counter();pool={}
        def admit(f,t,recipe):
            if t<=0 or not math.isfinite(t):return
            v=C.evaluate(f,t);g=N@ETA@v+off
            if v[0]<=pr['center_lower_bound']+1e-7 or abs(v@ETA@v-q)>3e-5:return
            if not all(C.compatible(x,4e-5)for x in g):return
            cover=np.flatnonzero(points@ETA@v>2e-7).tolist()
            if not cover:return
            key=tuple(np.round(v,7));pool.setdefault(key,dict(vector=v,recipe=recipe,open_end_coverage=cover,native_products=g))
        for pair in pairs:
            A=N[list(pair)];G=A@ETA@A.T;Ji=np.linalg.inv(G);K=null_space(A@ETA);ee,U=np.linalg.eigh(K.T@ETA@K);K=K@U@np.diag(1/np.sqrt(abs(ee)))
            assert ee[0]<0<ee[1]
            for labels in itertools.product(range(2,maxlabel+1),repeat=2):
                counts['endpoint_recipes']+=1;c=np.array([L.product(m)-off[i]for i,m in zip(pair,labels)]);p=c@Ji@A;rho=q-p@ETA@p
                if rho<=1e-7:counts['nonpositive_endpoint_schur']+=1;continue
                for sign in(-1,1):
                    counts['finite_endpoint_branches']+=1
                    f=dict(kind='hyperbola',p=p,a=sign*math.sqrt(rho)*(K[:,1]+K[:,0])/2,b=sign*math.sqrt(rho)*(K[:,1]-K[:,0])/2)
                    polys=[C.poly(f,np.array([-1.,0,0,0]),pr['center_lower_bound'])];fixed=[];bad=False
                    for j,n in enumerate(N):
                        if j in C.S or j in pair:continue
                        cov=ETA@n;vary=max(abs(f['a']@cov),abs(f['b']@cov));constant=p@cov+off[j]
                        if vary<1e-8:
                            if not C.compatible(constant,4e-5):bad=True;break
                            fixed.append(j)
                        else:polys.append(C.poly(f,cov,off[j]+2))
                    if bad:counts['fixed_contact_gap']+=1;continue
                    good=[]
                    for interval in C.intervals(polys,True):
                        lo,hi=interval;cuts=sorted({t for point in points for t in C.roots(C.poly(f,ETA@point))if t>0 and(lo is None or t>lo+1e-8)and(hi is None or t<hi-1e-8)})
                        cuts=[lo]+cuts+[hi];ns=0
                        for a,b in zip(cuts,cuts[1:]):
                            for t in C.parameter_samples([a,b],True):
                                before=len(pool);admit(f,t,dict(kind='continuous_two_contact',anchors=pair,labels=labels,branch=sign,parameter=t));ns+=len(pool)>before
                        if ns:good.append(dict(interval=interval,new_samples=ns))
                    if good:rows.append(dict(**f,anchors=pair,labels=labels,branch=sign,particular_square=float(p@ETA@p),intervals=good))
                    # Isolated third finite contacts lie outside the open
                    # all-other-dotted intervals, and are solved separately.
                    for third in sides:
                        if third in pair:continue
                        for m in range(2,maxlabel+1):
                            counts['third_contact_equations']+=1
                            poly=C.poly(f,ETA@N[third],off[third]-L.product(m))
                            if max(abs(poly))<1e-10:continue
                            for t in C.roots(poly):
                                admit(f,t,dict(kind='isolated_three_contact',anchors=list(pair)+[third],labels=list(labels)+[m],branch=sign,parameter=t))
        counts['admissible_endpoint_families']=len(rows);counts['candidate_count']=len(pool);totals.update(counts)
        best=[]
        for cand in pool.values():
            try:rr=B.clip(rays,cand['vector'],55)
            except ArithmeticError:continue
            cc=dict(Counter(r['kind']for r in rr));score=(cc.get('hyperideal_rays',0)+cc.get('unresolved_null',0),-len(cand['open_end_coverage']))
            best.append((score,dict(candidate_recipe=cand['recipe'],vector=cand['vector'],counts=cc,open_end_coverage=cand['open_end_coverage'])))
        best.sort(key=lambda z:z[0])
        allrows.append(dict(**pr,counts=dict(counts),families=rows,candidates=list(pool.values()),best_single_clips=[v for _,v in best[:12]]))
        print('A2_PROFILE',pr['node'],pr['label'],dict(counts),flush=True)
    paths=[f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'scripts/{L.TAG}.py',f'results/{C.TAG}.json',C.H3.CAP,C.CAPS,
           'results/fresh_h9_persistent_cusp_screen_20261016.json']
    return C.clean(dict(status='COMPLETE_FINITE_A2_PROFILE_LEMMA_WITH_BOUNDED_ENDPOINT_SCOUT',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
        exact_profile_lemma=algebra(),max_endpoint_label=maxlabel,retained_walls=55,initial_open_edge_ends=len(ends),counts=dict(totals),profiles=allrows,
        scope='The seven finite S6 facet profiles are classified for all labels by exact positive H4/A2 inverse arguments. The endpoint search is floating point, bounded to old open-edge pairs and labels2..max_endpoint_label, retains the saved common center, and requires other varying native contacts strictly dotted. Additional isolated third finite contacts are separately solved. No continuous-family exclusion, new exact cap, or compactness claim.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--maxlabel',type=int,default=10);p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run(a.maxlabel);out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved A2 profile scout')
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
