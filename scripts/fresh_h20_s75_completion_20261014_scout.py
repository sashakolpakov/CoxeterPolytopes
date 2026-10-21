#!/usr/bin/env python3
"""Bounded numerical nonright-profile scout on the untouched S75 face.

All acceptance in this file is numerical. Frozen predecessors are read-only.
"""
import argparse,json,math,itertools
from pathlib import Path
from fractions import Fraction as F
from decimal import Decimal,localcontext
from hashlib import sha256
from collections import Counter
import numpy as np
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_s75_completion_20261014_scout'
FILES=['results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json','results/fresh_h20_e7_bridge_face_20261014_exact.json','results/fresh_h20_bridge_edge_closure_20261014.json','results/fresh_h20_bridge_faces_20261014_review.json']
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def load():
    old,b,edge,face=[json.loads((ROOT/p).read_text())for p in FILES]
    for doc in[old,b,edge,face]:
        for p,h in doc['sources'].items():assert sha(p)==h,p
    f=next(x for x in face['faces']if 75 in x['support']);S=f['support']
    with localcontext()as ctx:
        ctx.prec=120;rr=[Decimal(d).sqrt()for d in old['field']['base_square_root_basis']]
        def v0(x):return sum((Decimal(F(a).numerator)/Decimal(F(a).denominator))*r for a,r in zip(x,rr))
        D=v0(old['field']['discriminant']).sqrt()
        def v1(x):return v0(x[0])+D*v0(x[1])
        K=v1(b['field']['new_positive_sqrt_radical']).sqrt()
        def v2(x):return v1(x[0])+K*v1(x[1])
        Delta=v2(edge['field']['new_positive_sqrt_radical']).sqrt();A=v2(edge['recipe']['A'])
        def v3(x):return v2(x[0])+Delta*v2(x[1])
        ND=[[v1(x)for x in v]for v in old['normals']]+[[v2(x)for x in b['new_normal']]]
        ND +=[[v3(x)/A for x in v]for v in edge['new_normals']]
        HD=[v2(x)for x in edge['strict_center']];ED=[v2(x)for x in edge['strict_escape']['point']]
        def dec(x):return Decimal(F(x).numerator)/Decimal(F(x).denominator)
        QD=[[dec(x)for x in row]for row in old['form']];JD=[[dec(x)for x in row]for row in f['inverse']]
        def qv(v):return[sum(a*b for a,b in zip(row,v))for row in QD]
        def dp(v,w):return sum(a*b for a,b in zip(v,qv(w)))
        nq=[qv(v)for v in ND]
        GD=[[sum(a*b for a,b in zip(v,w))for w in nq]for v in ND]
        EN=[sum(a*b for a,b in zip(ED,w))for w in nq];HN=[sum(a*b for a,b in zip(HD,w))for w in nq]
        ids=[0,11,28]
        basis=[]
        for j in ids:
            c=[sum(JD[k][ell]*GD[S[ell]][j]for ell in range(18))for k in range(18)]
            basis.append([ND[j][a]-sum(c[k]*ND[S[k]][a]for k in range(18))for a in range(21)])
        rawpairs=[[dp(v,w)for w in basis]for v in ND]
        # Exact support projections vanish; never normalize floating noise.
        for i in S:rawpairs[i]=[Decimal(0)]*3
        f['_highprecision']=dict(gram=np.array(GD,dtype=float),EN=np.array(EN,dtype=float),HN=np.array(HN,dtype=float),basis=np.array(basis,dtype=float).T,basis_gram=np.array([[dp(v,w)for w in basis]for v in basis],dtype=float),pairs=np.array(rawpairs,dtype=float),E=np.array([dp(ED,w)for w in basis],dtype=float),H=np.array([dp(HD,w)for w in basis],dtype=float))
        N=np.array(ND,dtype=float);H=np.array(HD,dtype=float);E=np.array(ED,dtype=float)
    Q=np.array([[float(F(x))for x in row]for row in old['form']]);norms=np.array(list(map(float,old['norms']))+[2.,2.,2.])
    GI=np.array([[float(F(x))for x in row]for row in f['inverse']])
    return np.array(N),Q,norms,H,E,S,GI,f
def run():
    N,Q,norms,H,E,S,GI,exactface=load();hpdata=exactface['_highprecision'];G=hpdata['gram']
    vals,V=np.linalg.eigh(hpdata['basis_gram']);assert vals[0]<0<vals[1]
    T=V@np.diag(1/np.sqrt(abs(vals)));U=hpdata['basis']@T;J=np.diag([-1.,1.,1.]);assert np.max(abs(T.T@hpdata['basis_gram']@T-J))<1e-10
    AA=hpdata['pairs']@T;ee=hpdata['E']@T;EE=ee/np.linalg.norm(ee);HH=hpdata['H']@T
    profilelabels={1:[3,4],2:[3],9:[3],10:[3,4],12:[3],17:[3,4,5,6,7,8,10,12,16,24,30,60],31:[3],87:[3]}
    profiles=[dict(contacts=[],coefficients=[0.]*18,cost=0.)]
    for size in[1,2]:
        for ids in itertools.combinations(profilelabels,size):
            for ls in itertools.product(*(profilelabels[i]for i in ids)):
                b=np.zeros(18)
                for i,m in zip(ids,ls):b[S.index(i)]=2*math.cos(math.pi/m)
                cost=b@GI@b
                if cost<2-1e-10:profiles.append(dict(contacts=list(zip(ids,ls)),coefficients=b.tolist(),cost=float(cost)))
    # Fixed, compact anchor pool from old S75 polar rays plus the three new walls.
    anchors=[0,61,65,76,79,80,95,96,97]
    labels=[2,3,4,5,6,8,10,12,'parallel']
    pairs=list(itertools.combinations(anchors,2));hitkeys=set();hits=[];stats=Counter();profile_stats=[]
    for pi,pr in enumerate(profiles):
        b=np.array(pr['coefficients']);coeff=GI@b;p=-coeff@N[S];q=2-pr['cost'];pn=-G[:,S]@coeff;pn[S]=-b;ep=-hpdata['EN'][S]@coeff;hp=-hpdata['HN'][S]@coeff
        # Linear feasibility is only a relaxed diagnostic; it omits norm=q.
        scales=np.maximum(np.linalg.norm(AA,axis=1),1)
        lp=linprog(-EE,A_ub=AA/scales[:,None],b_ub=-pn/scales,bounds=[(None,None)]*3,method='highs')
        report=dict(profile=pi,contacts=pr['contacts'],cost=pr['cost'],linear_relaxation_status=int(lp.status),linear_max_escape=None if lp.fun is None else float(-lp.fun+ep/np.linalg.norm(ee)))
        local=Counter()
        for i,j in pairs:
            B=AA[[i,j]];cr=np.cross(*B);cn=np.linalg.norm(cr)
            if cn<1e-10:continue
            cr/=cn;bb=B@B.T
            if np.linalg.cond(bb)>1e15:local['ill_conditioned_anchor_pairs']+=1;continue
            pinv=B.T@np.linalg.inv(bb);qa=cr@J@cr
            for mi,mj in itertools.product(labels,repeat=2):
                rhs=-np.array([math.sqrt(2*norms[k])*(1 if m=='parallel'else math.cos(math.pi/m))for k,m in[(i,mi),(j,mj)]])-pn[[i,j]]
                part=pinv@rhs;qb=2*part@J@cr;qc=part@J@part-q;disc=qb*qb-4*qa*qc
                if disc< -1e-8:continue
                roots=[-qc/qb]if abs(qa)<1e-13 and abs(qb)>1e-13 else[]if abs(qa)<1e-13 else[(-qb+sg*math.sqrt(max(0,disc)))/(2*qa)for sg in[-1,1]]
                for t in roots:
                    u=part+t*cr;local['norm_candidates']+=1
                    cut=ee@u+ep
                    if cut<=1e-8*np.linalg.norm(ee):continue
                    local['escape_cutters']+=1;g=AA@u+pn
                    if np.max(g/np.maximum(np.linalg.norm(N,axis=1),1))>1e-8:continue
                    if HH@u+hp>=0:continue
                    local['acute_cutters']+=1;bad=[];relations=[]
                    for k,c in enumerate(-g/np.sqrt(2*norms)):
                        if abs(c)<2e-6:relations.append(2)
                        elif c>=1-2e-6:relations.append('parallel_or_dotted')
                        elif c<.5-2e-6:bad.append((k,float(c)));break
                        else:
                            m=math.pi/math.acos(min(1,c));mi0=round(m)
                            if abs(m-mi0)>2e-5:bad.append((k,float(c)));break
                            relations.append(mi0)
                    if bad:continue
                    r=p+U@u;key=tuple(np.round(r,6))
                    if key in hitkeys:continue
                    hitkeys.add(key);local['accepted_Coxeter_cutters']+=1
                    rec=dict(profile=pi,contacts=pr['contacts'],anchors=[i,j],labels=[mi,mj],normal=r.tolist(),residual_coordinates=u.tolist(),escape_pairing=float(cut),normalization_error=float(u@J@u-q),all_relations=relations)
                    hits.append(rec);print('NUMERIC_HIT',pi,pr['contacts'],[i,j],[mi,mj],cut,flush=True)
        stats.update(local);report['counts']=dict(local);profile_stats.append(report)
    paths=FILES+[f'scripts/{TAG}.py']
    return dict(status='NUMERICAL_S75_LOW_CONTACT_BOUNDED_ANCHOR_SCOUT',sources={p:sha(p)for p in paths},support=S,inverse_diagonal=list(zip(S,[exactface['inverse'][i][i]for i in range(18)])),profiles=profiles,profile_count=len(profiles),anchors=anchors,anchor_labels=labels,counts=dict(stats),profile_diagnostics=profile_stats,candidates=hits,arithmetic='Physical Gram, residual basis and all projected wall/escape pairings evaluated at120 decimal digits before float eigensolver,LP and quadratic scout; support projection zeros imposed by exact identity; knownnormal norms retained.',scope='Numerical only. Necessary low-contact profiles with cost<2; isolated17 I2(m) tail sampled at7,8,10,12,16,24,30,60 only. Two finite anchors drawn from specified9-wall pool and specified labels. Linear relaxation omits norm constraint. No exhaustive real-parameter or all-label exclusion; no exact accepted construction.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert json.loads(path.read_text())==out
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['profile_count'],out['counts'],flush=True)
