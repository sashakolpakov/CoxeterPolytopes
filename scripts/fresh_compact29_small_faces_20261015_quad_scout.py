#!/usr/bin/env python3
"""Bounded quadrilateral-star existence scout; every failure remains open."""
import json
from itertools import combinations
import numpy as np
from fresh_compact29_small_faces_20261015_scout import ROOT,diagram,extensions,gram

TAG='fresh_compact29_small_faces_20261015_quad_scout'


def quad_ok(R):
    A=R[:3,:3]
    if np.linalg.det(A)>=-1e-8: return None
    d=np.linalg.solve(A,R[:3,3]); I=np.linalg.inv(A)
    vv=[]
    for i,j in [(0,1),(1,2),(2,3),(3,0)]:
        rows=[]
        for k in [i,j]: rows.append(A[k] if k<3 else R[:3,3])
        k=next(k for k in range(3) if k not in [i,j])
        rows.append(A[k]); x=np.linalg.solve(np.array(rows),[0,0,-1])
        ps=np.append(A@x,np.dot(R[:3,3],x))
        if max(ps)>1e-8 or x@A@x>=-1e-8: return None
        vv.append(x)
    dots=np.array(vv)@A@np.array(vv).T
    if np.max(dots)>=-1e-8:return None
    return dots.tolist()


def find_quad(kind,n,a,b):
    e=diagram(kind,n); W=[i for i in range(n) if i not in [a,b]]
    g0=gram(n,e); J=np.linalg.inv(g0[np.ix_(W,W)])
    op_c=extensions(W+[b],e,n); op_d=extensions(W+[a],e,n+1)
    for c in op_c:
        if any(m>5 for m in c.values()):continue
        for d in op_d:
            if any(m>5 for m in d.values()):continue
            G=np.zeros((n+2,n+2)); G[:n,:n]=g0; G[n,n]=G[n+1,n+1]=2
            for i,m in c.items():G[i,n]=G[n,i]=-2*np.cos(np.pi/m)
            for i,m in d.items():G[i,n+1]=G[n+1,i]=-2*np.cos(np.pi/m)
            order=[a,b,n,n+1]
            P=G[np.ix_(order,W)]@J@G[np.ix_(W,order)]
            for cd in [2,3,4,5]:
                G[n,n+1]=G[n+1,n]=-2*np.cos(np.pi/cd) if cd>2 else 0
                R=G[np.ix_(order,order)]-P
                if np.linalg.det(R[2:,2:])<=1e-8:continue
                for ac in [0,1,np.sqrt(2),(1+np.sqrt(5))/2,np.sqrt(3),2.1,2.5,3,4,6,10]:
                    G[a,n]=G[n,a]=-ac
                    R=G[np.ix_(order,order)]-P
                    def det(y):
                        Z=R.copy();Z[1,3]=Z[3,1]=-y-P[1,3]
                        return np.linalg.det(Z)
                    f0,f1,fm=det(0),det(1),det(-1)
                    A=(f1+fm)/2-f0;B=(f1-fm)/2;C=f0
                    if abs(A)<1e-10:continue
                    disc=B*B-4*A*C
                    if disc<-1e-9:continue
                    for bd in [(-B+np.sqrt(max(0,disc)))/(2*A),(-B-np.sqrt(max(0,disc)))/(2*A)]:
                        if bd<=2+1e-6:continue
                        G[b,n+1]=G[n+1,b]=-bd
                        RR=G[np.ix_(order,order)]-P
                        v=quad_ok(RR)
                        if v is None:continue
                        bad=False
                        for size in range(3,n+2):
                            for inds in combinations(range(n+2),size):
                                ev=np.linalg.eigvalsh(G[np.ix_(inds,inds)])
                                if abs(ev[0])<1e-8 and ev[1]>1e-8:bad=True;break
                            if bad:break
                        if bad:continue
                        return dict(c=sorted(c.items()),d=sorted(d.items()),cd=cd,ac=float(ac),bd=float(bd),
                                    bd_quadratic=[float(A),float(B),float(C)],vertex_pairings=v)
    return None


def run():
    old=json.loads((ROOT/'results/fresh_compact29_small_faces_20261015_scout.json').read_text())
    out={}
    for kind,n in [('E',7),('E',8)]:
        rows=[]
        for r in old['diagrams'][kind+str(n)]:
            if r['triangle'] is not None:continue
            a,b=r['marked'];w=find_quad(kind,n,a,b)
            rows.append(dict(marked=[a,b],quad=w));print(kind+str(n),a,b,bool(w),flush=True)
        out[kind+str(n)]=rows
    (ROOT/f'results/{TAG}.json').write_text(json.dumps(dict(scope='Bounded numerical existence scout; unfilled cases remain open.',diagrams=out),indent=2)+'\n')


if __name__=='__main__':run()
