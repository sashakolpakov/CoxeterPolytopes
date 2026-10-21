#!/usr/bin/env python3
"""Discovery only: H9 seeds with three/four sum blocks, finite-angle caps.

All geometry and norm checks here use floats. A positive candidate requires a
separate exact certificate. Parameter and selection restrictions are recorded.
"""
import argparse,itertools,json,math
from pathlib import Path
import numpy as np
from scipy.spatial import HalfspaceIntersection,QhullError

ROOT=Path(__file__).resolve().parents[1]
PHI=(1+5**.5)/2
VALS=np.array([0.,1.,PHI])
LABELS=np.array([0.]+[2*np.cos(np.pi/m) for m in range(3,121)])

def blocktypes(n):
    out=[]
    for lo,hi in itertools.combinations_with_replacement(range(len(VALS)),2):
        for count in ([n] if lo==hi else range(n+1)):
            q=VALS[lo]-VALS[hi]
            out.append((lo,hi,count,VALS[lo]+VALS[hi],(2*count-n)*q,n*q*q,math.comb(n,count)))
    return out

def profiles(ns):
    types=[blocktypes(n) for n in ns]
    inds=np.array(list(itertools.product(*[range(len(t)) for t in types])),dtype=int)
    tv=np.array([[types[g][i][3] for g,i in enumerate(row)] for row in inds])
    dv=np.array([[types[g][i][4] for g,i in enumerate(row)] for row in inds])
    sq=np.array([sum(types[g][i][5] for g,i in enumerate(row)) for row in inds])
    mult=np.array([math.prod(types[g][i][6] for g,i in enumerate(row)) for row in inds])
    return types,inds,tv,dv,sq,mult

def parameters(ns,b,P):
    types,inds,tv,dv,sq,mult=P;n=np.array(ns);k=len(ns)
    denominators=4-b*n
    C=np.sum(tv*tv,axis=1)/(2*b)+sq/8+b*np.sum(dv*dv/denominators,axis=1)/8-2
    A=np.sum(tv,axis=1)**2/(2*b)
    B=np.sum(dv/denominators,axis=1)**2/2
    L=np.sum(n/denominators)
    qa=C*k*L-A*L-B*k;qb=C*(k-b*L)-A+B*b;qc=-C*b
    points={}
    for i in range(len(inds)):
        if abs(qa[i])<1e-10:
            roots=[] if abs(qb[i])<1e-10 else [-qc[i]/qb[i]]
        else:
            disc=qb[i]**2-4*qa[i]*qc[i]
            if disc< -1e-8:continue
            roots=[(-qb[i]+s*max(0,disc)**.5)/(2*qa[i]) for s in [-1,1]]
        for a in roots:
            if not (b+2+1e-7<a<20):continue
            key=round(float(a),8)
            if key not in points: points[key]={'a':float(a),'profiles':set()}
            points[key]['profiles'].add(i)
    # Re-evaluate complete orbit equations to avoid rounded-key omissions.
    for point in points.values():
        a=point['a'];res=C-a*A/(k*a-b)-a*B/(1+a*L)
        ii=np.flatnonzero(abs(res)<1e-7)
        point['profiles']=ii.tolist();point['cap_count']=int(sum(mult[ii]))
    return sorted(points.values(),key=lambda p:(-p['cap_count'],p['a']))

def seed(ns,a,b):
    Z=np.zeros((sum(ns),len(ns)));start=0
    for g,n in enumerate(ns):Z[start:start+n,g]=1;start+=n
    K=Z@(a*np.ones((len(ns),len(ns)))-b*np.eye(len(ns)))@Z.T
    C=-2*np.eye(sum(ns))-K
    return np.block([[2*np.eye(sum(ns)),C],[C,2*np.eye(sum(ns))]]),Z,K

def variants(n,typ):
    lo,hi,count,*_=typ
    for I in itertools.combinations(range(n),count):
        x=np.full(n,VALS[hi]);y=np.full(n,VALS[lo]);x[list(I)]=VALS[lo];y[list(I)]=VALS[hi]
        yield x,y

def bank(ns,a,b,P,ii):
    types,inds,*_=P;G,Z,K=seed(ns,a,b);inv=np.linalg.pinv(G,rcond=1e-10)
    cov=[];refs=[]
    for i in ii:
        ts=[types[g][j] for g,j in enumerate(inds[i])]
        for parts in itertools.product(*[list(variants(n,t)) for n,t in zip(ns,ts)]):
            x=np.concatenate([p[0] for p in parts]);y=np.concatenate([p[1] for p in parts]);cov.append(-np.r_[x,y]);refs.append([x.tolist(),y.tolist()])
    cov=np.array(cov);R=cov@inv
    assert np.max(abs(R@G-cov))<1e-6
    assert np.max(abs(np.sum((R@G)*R,axis=1)-2))<1e-6
    return G,R,refs

def compat(H):
    out=H<=-2+1e-7
    for q in LABELS:out|=abs(H+q)<1e-7
    np.fill_diagonal(out,False)
    return out

def geometry(G,R,selection):
    val,Q=np.linalg.eigh(G);idx=list(np.where(val< -1e-7)[0])+list(np.where(val>1e-7)[0])
    assert len(idx)==10 and sum(val < -1e-7)==1
    N=Q[:,idx]*np.sqrt(abs(val[idx]));N[:,0]*=np.sign(sum(N[:,0]));W=np.vstack([N,R[selection]@N])
    assert min(W[:,0])>1e-7
    try:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9)).intersections
    except QhullError:V=HalfspaceIntersection(np.c_[W[:,1:],-W[:,0]],np.zeros(9),qhull_options='Qs').intersections
    q=np.sum(V*V,axis=1)
    return {'vertices':len(V),'outside':int(sum(q>1+1e-7)),'ideal':int(sum(abs(q-1)<1e-7)),'max_klein_square':float(max(q)),'worst_vertex':V[np.argmax(q)].tolist()}

def run(out,per_slice,attempts):
    records=[];census=[];rng=np.random.default_rng(610132)
    bs=[.5,PHI-1,1.,PHI/2,(PHI+1)/2,PHI]
    for ns in [[2,2,1,1],[3,2,2]]:
        P=profiles(ns)
        for b in bs:
            if b>=4/max(ns):continue
            pts=parameters(ns,b,P);census.append({'groups':ns,'b':b,'parameter_count':len(pts),'top_cap_count':pts[0]['cap_count'] if pts else 0})
            print(census[-1],flush=True)
            for point in pts[:per_slice]:
                a=point['a'];G,R,refs=bank(ns,a,b,P,point['profiles']);C=compat(R@G@R.T)
                selections=set()
                for trial in range(attempts):
                    remaining=np.ones(len(R),dtype=bool);sel=[]
                    while np.any(remaining):
                        avail=np.flatnonzero(remaining)
                        if trial==0:j=avail[np.argmax(np.sum(C[np.ix_(avail,avail)],axis=1))]
                        else:j=rng.choice(avail)
                        sel.append(int(j));remaining&=C[j]
                    selections.add(tuple(sorted(sel)))
                best=None
                for sel in selections:
                    geo=geometry(G,R,list(sel));key=(geo['max_klein_square'],geo['outside'])
                    if best is None or key<best[0]:best=(key,sel,geo)
                rec={'groups':ns,'b':b,'a':a,'cap_count':len(R),'edges':int(sum(sum(C))//2),'selections':len(selections),'max_selected':max(map(len,selections)),'selected_indices':list(best[1]),'selected_covectors':[refs[i] for i in best[1]],'geometry':best[2]}
                records.append(rec)
                print({k:rec[k] for k in ['groups','a','b','cap_count','max_selected','geometry']},flush=True)
                out.write_text(json.dumps({'scope':'Numerical only; K=Z(aJ-bI)Z^T; groups2,2,1,1 and3,2,2; six b values subject to signature; a roots in(b+2,20) of all three-label orbit norm equations; top parameter points by cap count per slice; random greedy compatible selections; labels2..120 for cap pairs. No exact acceptance or exclusion claim.','parameters_per_slice':per_slice,'selection_attempts':attempts,'census':census,'records':records},indent=2)+'\n')
    print({'done':len(records),'best':min(r['geometry']['max_klein_square'] for r in records)},flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',default='data/enumerations/fresh_h9_multiblock_seed_20261013_scout.json');p.add_argument('--per-slice',type=int,default=12);p.add_argument('--attempts',type=int,default=24);args=p.parse_args();run(ROOT/args.output,args.per_slice,args.attempts)
