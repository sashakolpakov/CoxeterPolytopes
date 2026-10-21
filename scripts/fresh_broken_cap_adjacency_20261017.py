#!/usr/bin/env python3
"""Numerical five-cap disk scout after subdividing the forbidden A--R edge.

All old geometry is inherited unchanged. Compact summary output records the
tested domains, exact deletion ledger, blockers, and any coherent local bridge.
"""
import argparse,json
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from pathlib import Path
import numpy as np
import fresh_new_cavity_construction_20261017 as N

S,H,B,G,M,ROOT=N.S,N.H,N.B,N.G,N.M,N.ROOT
TAG='fresh_broken_cap_adjacency_20261017'
TGRID=['6/5','5/4','4/3','3/2','2','3','4','6','8','12','20']
HIGH_LABELS=[6,10,24]

@lru_cache(None)
def corner(pi,m8,mo,otherpi,oldbond):
    g=[[2]*5 for _ in range(5)]
    for i in range(5):g[i][i]=1
    g[0][1]=g[1][0]=3
    for i,m in enumerate(G.profile_labels(otherpi)):g[i][3]=g[3][i]=m
    for i,m in enumerate(G.profile_labels(pi)):g[i][4]=g[4][i]=m
    g[2][3]=g[3][2]=oldbond;g[2][4]=g[4][2]=m8;g[3][4]=g[4][3]=mo
    return N.A.components(tuple(v for row in g for v in row))

@lru_cache(None)
def domain(rpi,r8):
    rows=[]
    for pi,m8,ma,mr in product(range(7),range(2,7),range(2,7),range(2,7)):
        ca=corner(pi,m8,ma,0,2);cr=corner(pi,m8,mr,rpi,r8)
        if ca is None or cr is None:continue
        rows.append(dict(profile=pi,labels=[m8,ma,mr],high_positions=[j for j,m in enumerate([m8,ma,mr])if m==6],corner_types=[ca,cr]))
    return rows

def a_curve(ns,h):
    f=np.array([ns[7],ns[8]]);inv=np.linalg.inv(f@M@f.T)
    def proj(v):return v-(v@M@f.T)@inv@f
    k0=proj(h);assert k0@M@k0<0;k0=k0/np.sqrt(-k0@M@k0)
    k1=proj(ns[30]);k1=k1+float(k1@M@k0)*k0;assert k1@M@k1>0;k1=k1/np.sqrt(k1@M@k1)
    assert abs(k0@M@k0+1)<1e-8 and abs(k1@M@k1-1)<1e-8 and abs(k0@M@k1)<1e-8
    return k0,k1

def native_values(y,pi,labels,ns,profiles):
    p=profiles[pi];return {i:float(y@M@ns[i]+p['offsets'][labels.index(i)])for i in labels if i not in S.DELETE}

def gates(vals,assigned=()):
    assigned=set(assigned);bad=[]
    for i,v in vals.items():
        if i in assigned:continue
        rel=B.relation(v)
        if rel is None or rel=='parallel_or_unresolved':bad.append((i,v))
    return bad

def sample_recipe(r):
    high=r['high_positions']
    for hs in product(HIGH_LABELS,repeat=len(high)):
        ms=list(r['labels'])
        for j,m in zip(high,hs):ms[j]=m
        yield ms

def run():
    frozen=json.loads((ROOT/'results'/f'{S.TAG}.json').read_text());rb=frozen['accepted_low_caps'];assert len(rb)==94
    z,labels,ns,h,profiles,primary=H.setup();actual=[i for i in z['actual_sides']if i not in S.DELETE]
    k0,k1=a_curve(ns,h);acandidates=[];atrials=[]
    for sign,ts in product([-1,1],TGRID):
        t=float(Fraction(ts));y=np.sqrt(.5)*(t*(k0+sign*k1)+(-k0+sign*k1)/t)
        assert abs(y@M@y-2)<1e-7 and max(abs(y@M@ns[i])for i in[7,8])<1e-7
        vals=native_values(y,0,labels,ns,profiles);bad=gates(vals,[7,8]);cp=float(y@M@h)
        rec=dict(sign=sign,t=ts,projection=y.tolist(),center_product=cp,blockers=bad)
        atrials.append(rec)
        if not bad and cp<-1e-8:acandidates.append(rec)
    counts=Counter();blockers=Counter();perR=[];accepted=[];best=[];unresolved=[]
    for rr in rb:
        ry=np.array(rr['projection']);rpi=rr['recipe']['profile'];r8=rr['recipe']['labels'][0]
        rq=profiles[rpi]['square'];local=Counter()
        for ar in acandidates:
            ay=np.array(ar['projection']);av=float(ay@M@ry)
            rel=B.relation(av)
            if rel is None or rel=='parallel_or_unresolved':counts['AR_native_incompatible']+=1;continue
            if av>=-np.sqrt(2*rq)-1e-7:counts['AR_not_intrinsically_dotted']+=1;continue
            counts['AR_admissible_sample_pairs']+=1
            frame=[ns[8],ay,ry]
            for rec in domain(rpi,r8):
                xpi=rec['profile'];p=profiles[xpi]
                for ms in sample_recipe(rec):
                    off=H.A.native(np.zeros(4),xpi,np.zeros(4),rpi,profiles)
                    rhs=[-B.theta(ms[0]),-B.theta(ms[1]),-B.theta(ms[2])-off]
                    roots,norm=N.norm_roots(frame,rhs,p['square']);counts['norm_equations']+=1;local['norm_equations']+=1
                    if norm['status'].startswith('unresolved'):
                        unresolved.append(dict(R=rr['index'],A_sign=ar['sign'],A_t=ar['t'],X=rec,sampled_labels=ms,norm=norm));continue
                    if not roots:counts['no_real_norm']+=1;continue
                    for branch,xy in roots:
                        counts['real_branches']+=1;local['real_branches']+=1
                        assert abs(xy@M@xy-p['square'])<2e-6
                        assert max(abs(float(xy@M@v)-w)for v,w in zip(frame,rhs))<2e-6
                        vals=native_values(xy,xpi,labels,ns,profiles);bad=gates(vals,[8]);center=float(xy@M@h+p['center_offset'])
                        if center>=-1e-8:counts['wrong_center']+=1;continue
                        counts['future_branches']+=1;blockers.update(i for i,v in bad)
                        item=dict(R_index=rr['index'],R_recipe=rr['recipe'],A_sign=ar['sign'],A_t=ar['t'],AR_product=av,
                                  X_profile=xpi,X_labels=ms,X_symbolic_recipe=rec,X_branch=branch,X_projection=xy.tolist(),
                                  center_product=center,blockers=bad)
                        if len(bad)<=2:
                            best.append(item)
                        if not bad:
                            nd=dict(ns);nd.update({67:ay,68:ry,69:xy})
                            item['boundary']=S.C.detailed_cone(nd,actual+[67,68,69])
                            squares=[profiles[0]['square'],rq,p['square']]
                            oldsq={r['wall']:B.mid(r['square'])for r in z['projections']}
                            item['PD_incidence']=H.incidence_score([nd[i]for i in actual+[67,68,69]],actual+[67,68,69],[oldsq[i]for i in actual]+squares)
                            accepted.append(item);counts['coherent_X_bridges']+=1;local['coherent_X_bridges']+=1
        perR.append(dict(index=rr['index'],recipe=rr['recipe'],counts=dict(local)))
        if len(perR)%10==0:print('BROKEN_ADJACENCY',len(perR),'of94',dict(counts),'best',len(best),flush=True)
    paths=[Path(__file__),Path(N.__file__),Path(S.__file__),Path(H.__file__),Path(N.A.__file__),ROOT/'results'/f'{S.TAG}.json',H.PARENT,H.BOUNDARY]
    return dict(status='BOUNDED_NUMERICAL_STELLAR_SUBDIVISION_FIVE_CAP_DISK_DOMAIN',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        deleted_walls=S.DELETE,boundary_cycle=S.CYCLE,exact_old_deletion_cost=frozen['exact_old_deletion_cost'],
        deletion_boundary_summary={k:frozen['deletion_boundary'][k]for k in['counts','open_ends','warnings']},
        topology=dict(old_triangles_replaced=[[8,'A','R'],['A','R','B']],
                      new_triangles=[[8,'A','X'],[8,'X','R'],['A','X','B'],['X','R','B']],
                      other_central_triangle=['A','B','C'],removed_edge=['A','R'],
                      X_neighbors=[8,'A','R','B'],vertices=10,edges=22,triangles=13,x=5,y=5,z=3),
        A_curve=dict(frame=[7,8],timelike_unit=k0.tolist(),spacelike_unit=k1.tolist(),samples=atrials,
                     equation='sqrt(2)/2*[t(k0+s*k1)+(-k0+s*k1)/t], t>0, s=+/-1'),
        X_domains=[dict(R_profile=pi,R8_label=m,recipes=domain(pi,m))for pi,m in sorted({(r['recipe']['profile'],r['recipe']['labels'][0])for r in rb})],
        counts=dict(counts,A_samples=len(atrials),native_coherent_A_samples=len(acandidates),R_low_seeds=94),
        per_R=perR,retained_wall_blocker_counts=sorted(blockers.items()),near_candidates=best,accepted_bridges=accepted,unresolved_events=unresolved,
        sampling=dict(A_parameters=TGRID,symbolic_high_labels_sampled=HIGH_LABELS),
        scope='Bounded numerical necessary-bridge search for the five-cap disk obtained by subdividing edge A-R. A-R is required intrinsically dotted, so the frozen finite A-R gate does not apply. All retained old native contacts and the common center are checked for A,R,X. Both norm branches are tested. The complete finite endpoint recipe shapes are listed, but A parameters, R low roots and high labels are sampled; no all-parameter obstruction follows. B,C completion, permanent-cusp screening, exact root promotion and ambient finite-volume proof are required for any eventual full patch. No new chamber is certified here.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'accepted',len(out['accepted_bridges']),flush=True)
