#!/usr/bin/env python3
"""Numerical lazy conflict search for compatible old-edge endpoint coverage."""
import json,time,argparse
import numpy as np
import z3
from fresh_h9_makarov_joint_20261013_selection import setup,ROOT,PREFIX

def run(rounds=1000,budget=30):
    start=time.time();Q,W,c,R,refs,null,edges=setup();cov=R@Q@null.T>=-1e-10
    ids=np.flatnonzero(cov.any(axis=1));R=R[ids];cov=cov[ids];refs=[refs[i]for i in ids]
    v=[z3.Bool('r'+str(i))for i in range(len(ids))];s=z3.Solver();s.set(timeout=30000)
    for j in range(40):s.add(z3.Or(*[v[i]for i in np.flatnonzero(cov[:,j])]))
    if budget:s.add(z3.PbLe([(x,1)for x in v],budget))
    finite=np.array([-2*np.cos(np.pi/m)for m in [2,3,4,5,6,7,8,9,10,12,15,16,17,20,24,30]])
    conflicts=set();trace=[];answer=None
    for step in range(rounds):
        status=s.check()
        if status!=z3.sat:answer=str(status);break
        model=s.model();sel=[i for i,x in enumerate(v)if z3.is_true(model.eval(x,model_completion=True))]
        P=R[sel]@Q@R[sel].T;ok=(P<-2+3e-7)
        for f in finite:ok|=abs(P-f)<3e-7
        bad=[(sel[i],sel[j])for i in range(len(sel))for j in range(i)if not ok[i,j]]
        trace.append({'round':step,'selected':len(sel),'new_conflicts':len(bad)})
        if not bad:
            answer='NUMERIC_COMPATIBLE_COVER';break
        for i,j in bad:
            s.add(z3.Or(z3.Not(v[i]),z3.Not(v[j])));conflicts.add((min(i,j),max(i,j)))
        if step%20==0:print(step,'sel',len(sel),'conflicts',len(conflicts),'seconds',time.time()-start,flush=True)
    out={'status':answer or 'ROUND_LIMIT','numerical_only':True,'round_budget':rounds,'cap_budget':budget,
         'candidate_count':len(ids),'trace':trace,'conflicts_learned':len(conflicts),'seconds':time.time()-start,
         'selected':[int(ids[i])for i in sel]if answer=='NUMERIC_COMPATIBLE_COVER' else None,
         'references':[refs[i]for i in sel]if answer=='NUMERIC_COMPATIBLE_COVER' else None}
    (ROOT/f'results/{PREFIX}_cover.json').write_text(json.dumps(out,indent=2)+'\n');print(out['status'],len(conflicts),time.time()-start,flush=True)
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=1000);p.add_argument('--budget',type=int,default=30);a=p.parse_args();run(a.rounds,a.budget)
