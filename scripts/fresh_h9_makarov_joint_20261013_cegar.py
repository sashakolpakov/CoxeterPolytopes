#!/usr/bin/env python3
"""Bounded numerical joint selection with actual gap constraints."""
import argparse,json,time
import numpy as np
from scipy.optimize import linprog
from fresh_h9_makarov_joint_20261013_selection import setup,ROOT,PREFIX
from fresh_h9_makarov_joint_20261013_cover_dfs import bits

def run(rounds=30,seconds=120):
    start=time.time();Q,W,c,R,refs,null,edges=setup();N=len(R);RQ=R@Q
    cv=RQ@null.T>=-1e-10;masks=[bits(cv[:,j])for j in range(40)]
    capm=[sum(1<<int(j)for j in np.flatnonzero(row))for row in cv]
    cache={};trace=[];rng=np.random.default_rng(816);
    def compat(k):
        if k not in cache:
            ps=RQ@R[k];ok=ps<-2+3e-7
            for m in [2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]:ok|=abs(ps+2*np.cos(np.pi/m))<3e-7
            ok[k]=False;cache[k]=bits(ok)
        return cache[k]
    for it in range(rounds):
        nodes=0;tic=time.time()
        def dfs(avail,missing,chosen):
            nonlocal nodes
            nodes+=1
            if nodes%100==0 and time.time()-tic>seconds:raise TimeoutError
            if not missing:return chosen
            opts=None
            for j in range(len(masks)):
                if missing>>j&1:
                    a=avail&masks[j]
                    if not a:return None
                    if opts is None or a.bit_count()<opts.bit_count():opts=a
            ks=[]
            while opts:
                one=opts&-opts;k=one.bit_length()-1;opts-=one;ks.append(k)
            ks.sort(key=lambda k:-(capm[k]&missing).bit_count())
            for k in ks:
                out=dfs(avail&compat(k),missing&~capm[k],chosen+[k])
                if out is not None:return out
            return None
        try:sel=dfs((1<<N)-1,(1<<len(masks))-1,[])
        except TimeoutError:trace.append({'round':it,'status':'TIME_LIMIT','nodes':nodes});break
        if sel is None:trace.append({'round':it,'status':'NUMERICAL_UNSAT','nodes':nodes});break
        A=np.vstack([W,R[sel]])@Q;best=(-1e100,None)
        for j in range(40):
            direction=Q@best[1]if j%4==0 and best[1]is not None else rng.normal(size=10)
            lp=linprog(-direction,A_ub=A,b_ub=np.zeros(len(A)),A_eq=[c@Q],b_eq=[-1],bounds=[(None,None)]*10,method='highs')
            if not lp.success:continue
            q=lp.x@Q@lp.x
            if q>best[0]:best=(float(q),lp.x)
        trace.append({'round':it,'status':'NUMERICAL_SELECTION','nodes':nodes,'selected':sel,'references':[refs[k]for k in sel],
                      'best_norm':best[0],'best_point':best[1].tolist()})
        print('round',it,'nodes',nodes,'size',len(sel),'norm',best[0],flush=True)
        if best[0]<1e-10:break
        # Make the spacelike proposal strict on its chosen walls.
        x=best[1]*.999+c/(-c@Q@c)*.001
        assert x@Q@x>0 and max(A@x)<0
        cuts=RQ@x>1e-10;masks.append(bits(cuts));j=len(masks)-1
        for k in np.flatnonzero(cuts):capm[int(k)] |=1<<j
        trace[-1]['strict_escape']=x.tolist();trace[-1]['bank_cutters']=int(cuts.sum())
        if not cuts.any():trace[-1]['whole_bank_escape']=True;break
    out={'status':'BOUNDED_NUMERICAL_JOINT_COUNTEREXAMPLE_SEARCH','entries':N,'round_budget':rounds,'node_time_budget':seconds,
         'constraints':len(masks),'trace':trace,'seconds':time.time()-start,'cached_pair_rows':len(cache),
         'scope':'No exact compatibility, graph exclusion or volume claim; all selected recipes and remaining gap proposals are saved.'}
    (ROOT/f'results/{PREFIX}_cegar.json').write_text(json.dumps(out,indent=2)+'\n');print('DONE',out['seconds'],flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=30);p.add_argument('--seconds',type=int,default=120);a=p.parse_args();run(a.rounds,a.seconds)
