#!/usr/bin/env python3
"""Bounded nine-contact proposals over rank-seven supports, without an old-edge endpoint requirement."""
import argparse,gzip,json,math
from collections import Counter
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import mpmath as mp
import numpy as np
import fresh_h9_polar_endpoint_profiles_20261015 as P

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_joint_pair_20261015_scout'
BASE='results/fresh_h9_nine_contact_extension_20261015.json'
ANCHORS=[2,3,4,5,6,8,10,12,16,24,32,64]
TAILS=[6,7,8,10,12,16,24,32,64]

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def scalar(x):return mp.mpf(x[0][0])+mp.mpf(x[0][1])*mp.sqrt(5)+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*mp.sqrt(5))*mp.sqrt(6+16*mp.sqrt(5)/5)
def qval(x):return float(mp.mpf(str(x.a))+mp.mpf(str(x.b))*mp.sqrt(5))
def finite(S,labels):
    if any(tuple(p)not in labels for p in combinations(S,2)):return False
    es={p:labels[p]for p in combinations(S,2)}
    cs,ad=P.P.components(list(S),es)
    try:
        for c in cs:
            if max(len(ad[i])for i in c)>3:return False
            P.P.component_type(c,es)
    except AssertionError:return False
    return True

def run(all_faces=False,right_only=False):
    mp.mp.dps=100
    d=json.loads((ROOT/BASE).read_text())
    for p,h in d['sources'].items():assert sha(p)==h,p
    Q=np.array([[float(scalar(x))for x in row]for row in d['form']])
    T=np.zeros((10,10));T[0,0]=math.sqrt(-Q[0,0]);T[1:,1:]=np.linalg.cholesky(Q[1:,1:])
    metric=np.diag([-1.]+[1.]*9)
    assert np.max(abs(T@metric@T.T-Q))<1e-12
    N=np.array([[float(scalar(x))for x in row]for row in d['normals']])
    scales=np.array([math.sqrt(2/float(scalar(x)))for x in d['norms']])
    N=(N@T)*scales[:,None]
    H=np.array([float(scalar(x))for x in d['center']])@T
    E=np.array([float(scalar(x))for x in d['strict_escape']['point']])@T
    G=np.array([[float(scalar(x))*scales[i]*scales[j]for j,x in enumerate(row)]for i,row in enumerate(d['gram'])])
    assert np.max(abs(N@metric@N.T-G))<1e-8
    labels={};fvalues=-2*np.cos(np.pi/np.arange(2,257))
    for i,j in combinations(range(len(N)),2):
        k=int(np.argmin(abs(G[i,j]-fvalues)))
        if abs(G[i,j]-fvalues[k])<1e-8:labels[i,j]=k+2
    active=d['strict_escape']['active'];pool=list(range(len(N)))if all_faces else active
    supports=[S for S in combinations(pool,7)if finite(S,labels)]
    counts=Counter();rows=[];found={};anchor_values=np.array(list(product([-2*math.cos(math.pi/m)for m in ANCHORS],repeat=2)))
    for sid,S in enumerate(supports):
        S=list(S);outside=[i for i in range(len(N))if i not in S]
        GS=G[np.ix_(S,S)];J=np.linalg.inv(GS);assert np.linalg.eigvalsh(GS)[0]>1e-7
        if right_only:profiles=[[]]
        else:
            es={(i,j):labels[tuple(sorted((S[i],S[j])))] if j<7 else 2 for i,j in combinations(range(8),2)}
            data=P.profiles(list(range(8)),es)
            profiles=[r['contacts']for r in data['finite_profiles']if all(i!=7 for i,m in r['contacts'])]
            profiles += [[[r['wall'],m]]for r in data['dihedral_tail_profiles']if r['wall']!=7 for m in TAILS]
        cv=np.zeros((len(profiles),7))
        for k,cs in enumerate(profiles):
            for i,m in cs:cv[k,i]=-2*math.cos(math.pi/m)
        coeff=cv@J;rho=2-np.sum(coeff*cv,axis=1);assert min(rho)>0
        parts=coeff@N[S];old_pairs=coeff@G[S]
        _,ss,vh=np.linalg.svd(N[S]@metric,full_matrices=True)
        if ss[-1]<1e-9:counts['numerical_support_skip']+=1;continue
        C=vh[7:];g=C@metric@C.T;vals,U=np.linalg.eigh(g)
        if min(abs(vals))<1e-12:counts['numerical_frame_skip']+=1;continue
        assert list(np.sign(vals))==[-1,1,1]
        C=(U/np.sqrt(abs(vals))).T@C;j3=np.diag([-1.,1.,1.])
        assert np.max(abs(C@metric@C.T-j3))<1e-6
        CN=C@metric@N.T;CH=C@metric@H;CE=C@metric@E
        HP=parts@metric@H;EP=parts@metric@E
        before=len(found);escape_events=compatible=0
        for a,b in combinations(outside,2):
            D=CN[:,[a,b]].T;dd=D@D.T
            if np.linalg.det(dd)<1e-12*max(1.,np.linalg.norm(dd)**2):counts['dependent_anchor_skips']+=1;continue
            R=D.T@np.linalg.inv(dd);w=np.cross(D[0],D[1]);w/=np.linalg.norm(w)
            wa=float(w@j3@w)
            if abs(wa)<1e-10:counts['null_kernel_skips']+=1;continue
            rhs=anchor_values[None,:,:]-old_pairs[:,None,[a,b]]
            z=rhs@R.T
            wb=2*z@j3@w;wc=np.sum((z@j3)*z,axis=2)-rho[:,None]
            disc=wb*wb-4*wa*wc;good=disc>1e-9
            counts['nine_contact_parameter_tests']+=int(good.size)
            for branch in(-1,1):
                lam=(-wb+branch*np.sqrt(np.maximum(disc,0)))/(2*wa)
                zz=z+lam[:,:,None]*w
                ep=EP[:,None]+zz@CE;hp=HP[:,None]+zz@CH
                candidate=good&(ep>1e-7)&(hp<-1e-7)
                pi,mi=np.where(candidate);escape_events+=len(pi)
                if not len(pi):continue
                vv=old_pairs[pi]+zz[pi,mi]@CN
                errors=np.min(abs(vv[:,:,None]-fvalues[None,None,:]),axis=2)
                accepted=np.all((vv<-2-1e-7)|(errors<1e-7),axis=1)
                for ix in np.flatnonzero(accepted):
                    ip,im=int(pi[ix]),int(mi[ix]);vec=vv[ix];key=tuple(np.round(vec,7))
                    compatible+=1
                    if key in found:continue
                    ia,ib=divmod(im,len(ANCHORS));contacts=[[S[i],m]for i,m in profiles[ip]]
                    normal=parts[ip]+zz[ip,im]@C
                    found[key]=dict(support=S,support_contacts=contacts,anchors=[[a,ANCHORS[ia]],[b,ANCHORS[ib]]],
                        quadratic_branch=branch,normalized_wall_pairings=vec.tolist(),
                        escape_pairing=float(ep[ip,im]),center_pairing=float(hp[ip,im]),
                        physical_normal=normal.tolist(),normal_square=float(normal@metric@normal),
                        finite_neighbors=[i for i,t in enumerate(vec)if t>-2+1e-7])
        row=dict(id=sid,support=S,sampled_profiles=len(profiles),escape_cutting_events=escape_events,
                 coxeter_compatible_events=compatible,new_distinct_candidates=len(found)-before)
        rows.append(row);print(row,flush=True)
    bank=sorted(found.values(),key=lambda r:-r['escape_pairing'])
    paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h9_polar_endpoint_profiles_20261015.py',
           'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py']
    return dict(status='BOUNDED_RANK7_TWO_ANCHOR_ESCAPE_SCOUT',sources={p:sha(p)for p in paths},
        all_rank7_faces=all_faces,right_profiles_only=right_only,support_pool=pool,
        support_count=len(supports),anchors=ANCHORS,tail_samples=TAILS,counts=dict(counts),supports=rows,
        candidate_count=len(bank),candidates=bank,
        scope='Bounded numerical proposals with a spherical rank7 support plus a new wall, two old finite anchors and a norm equation. No old-edge finite endpoint is required. All finite pair gates use labels2..256 and1e-7 tolerance after100-digit scalar evaluation. Numerical singular cases are skipped and counted; symbolic I2 tails sampled. No exact wall, pair completion or all-parameter exclusion is asserted.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--all-faces',action='store_true');p.add_argument('--right-only',action='store_true');a=p.parse_args()
    d=run(a.all_faces,a.right_only);bank=d.pop('candidates');suffix=('_all'if a.all_faces else'_escape')+('_right'if a.right_only else'')
    bp=ROOT/f'data/enumerations/{TAG}{suffix}_bank.json.gz';blob=gzip.compress(json.dumps(bank,separators=(',',':')).encode(),mtime=0)
    d['candidate_bank']=dict(path=str(bp.relative_to(ROOT)),sha256=sha256(blob).hexdigest(),count=len(bank));d['samples']=bank[:20]
    out=ROOT/f'results/{TAG}{suffix}.json'
    if a.verify:assert blob==bp.read_bytes()and d==json.loads(out.read_text())
    else:bp.write_bytes(blob);out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['support_count'],d['candidate_count'])
