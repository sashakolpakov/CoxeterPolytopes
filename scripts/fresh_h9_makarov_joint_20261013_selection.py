#!/usr/bin/env python3
"""Numerical compatibility-conditioned selection, never a discovery certificate."""
from pathlib import Path
from fractions import Fraction
import argparse,json,time
import numpy as np
from scipy.optimize import linprog
from fresh_h9_makarov_joint_20261013 import read

ROOT=Path(__file__).resolve().parents[1];PREFIX='fresh_h9_makarov_joint_20261013'
def setup():
    data,Q,W,c=read();caps=[];refs=[];keys=set()
    for suffix in ['', '_profiles']:
        D=json.loads((ROOT/f'results/{PREFIX}{suffix}.json').read_text())
        for row in D['caps']:
            r=np.array(row['normal']);key=tuple(np.rint(r*1e6).astype('int64'))
            if key in keys:continue
            keys.add(key);caps.append(r);refs.append([suffix,row['id']])
    R=np.array(caps);F=json.loads((ROOT/'results/fresh_h9_makarov_frontier_20261013.json').read_text())
    def f(v):return sum(float(Fraction(a))*np.sqrt(int(k))for k,a in v.items())
    null=[];edges=[]
    for ei,e in enumerate(F['outgoing_edges']):
        if e['kind']!='open':continue
        p=np.array([f(x)for x in F['finite_vertices'][e['vertex']]['point']])
        u=np.array([f(x)for x in e['dual']]);ray=p-np.sqrt(-(p@Q@p)/(u@Q@u))*u
        assert abs(ray@Q@ray)<1e-5
        null.append(ray/(-(ray@Q@c)));edges.append(ei)
    return Q,W,c,R,refs,np.array(null),edges

def run(trials=1000):
    start=time.time();Q,W,c,R,refs,null,edges=setup();RQ=R@Q
    covers=RQ@null.T>=-1e-10;N=len(R);rng=np.random.default_rng(611309)
    finite=np.array([-2*np.cos(np.pi/m)for m in [2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]])
    counts=covers.sum(axis=1);print('unique',N,'covers',dict(zip(*np.unique(counts,return_counts=True))),flush=True)
    good=[];best=0;seen=set();maxq=-1e30
    for trial in range(trials):
        available=np.ones(N,dtype=bool);uncovered=np.ones(40,dtype=bool);chosen=[]
        weight=np.exp(rng.uniform(-1,1,size=40))
        while np.any(uncovered):
            score=covers[:,uncovered]@weight[uncovered]
            score[~available]=-1
            if max(score)<=0:break
            top=np.flatnonzero(score>=max(score)*(.75 if trial else 1))
            k=int(rng.choice(top));chosen.append(k);uncovered[covers[k]]=False
            p=RQ@R[k];ok=(p<-2+3e-7)|(np.min(abs(p[:,None]-finite),axis=1)<3e-7)
            available &=ok;available[k]=False
        n=int(40-uncovered.sum())
        if n>best:best=n;print('covered',best,'trial',trial,'size',len(chosen),flush=True)
        if n<40:continue
        key=tuple(sorted(chosen))
        if key in seen:continue
        seen.add(key);rows=np.vstack([W,R[chosen]]);A=rows@Q;bestv=(-1e30,None)
        for step in range(24):
            d=Q@bestv[1]if step%3==0 and bestv[1]is not None else rng.normal(size=10)
            sol=linprog(-d,A_ub=A,b_ub=np.zeros(len(A)),A_eq=[c@Q],b_eq=[-1],bounds=[(None,None)]*10,method='highs')
            if not sol.success:continue
            q=sol.x@Q@sol.x
            if q>bestv[0]:bestv=(float(q),sol.x)
        good.append({'selected':chosen,'references':[refs[i]for i in chosen], 'best_norm':bestv[0],
                     'best_point':None if bestv[1]is None else bestv[1].tolist()})
        print('40-cover',len(good),'size',len(chosen),'volume',bestv[0],flush=True)
        if bestv[0]<-1e-8:break
    out={'status':'NUMERICAL_COMPATIBILITY_CONDITIONED_SCOUT','entries':N,'trials':trials,
         'best_old_edge_coverage':best,'cover_histogram':{str(k):int(v)for k,v in zip(*np.unique(counts,return_counts=True))},
         'all40_edge_covers':good,'old_edges':edges,'elapsed_seconds':time.time()-start,
         'scope':'Bounded heuristic and numerical angle gates only; no exhaustive failure or exact geometry claim.'}
    (ROOT/f'results/{PREFIX}_selection.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',time.time()-start,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--trials',type=int,default=1000);a=p.parse_args();run(a.trials)
