#!/usr/bin/env python3
"""Numerical mixed old-vertex/18-short bank and compatible completion tests."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
from fresh_h20_targeted_escape_20261013_second_volume import extended_setup

ROOT=Path(__file__).resolve().parents[1]
PREFIX='fresh_h20_edge_family_completion_20261013'
OLD='data/enumerations/fresh_h20_targeted_escape_20261013_vertex_bank.json'
EDGE='data/enumerations/fresh_h20_targeted_escape_20261013_edge_directions.json'
OUT='data/enumerations/'+PREFIX+'_bank.json'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def gate(rat):
    finite=np.cos(np.pi/np.arange(2,121))
    ix=np.clip(np.searchsorted(finite,rat),1,len(finite)-1)
    err=np.minimum(abs(rat-finite[ix]),abs(rat-finite[ix-1]))
    return (rat>=1-1e-7)|(err<1e-7)

def geometry(Q,R,h,objectives):
    ns=np.einsum('ij,jk,ik->i',R,Q,R);assert min(ns)>0
    A=R@Q/np.sqrt(ns[:,None]);c=h/np.sqrt(-h@Q@h);best=None
    for i,obj in enumerate(objectives):
        p=linprog(obj,A_ub=A,b_ub=np.zeros(len(A)),A_eq=(c@Q)[None,:],b_eq=[-1],bounds=[(None,None)]*21,method='highs')
        if not p.success:return {'status':'LP_UNRESOLVED','objective':i,'lp_status':int(p.status)}
        q=float(p.x@Q@p.x)
        if best is None or q>best['norm']:
            best={'norm':q,'point':p.x.tolist(),'objective':i,'active_indices':np.flatnonzero(abs(A@p.x)<1e-7).tolist()}
        if q>1e-7:best['status']='NUMERICAL_SPACELIKE_ESCAPE';return best
    best['status']='NO_ESCAPE_IN_BOUNDED_OBJECTIVES';return best

def run(rounds=24,extra_sources=(),out_path=OUT):
    Q,R,h,*_=extended_setup();assert len(R)==87
    sources=[OLD,EDGE]+list(extra_sources);caps=[]
    for p in sources:
        for i,c in enumerate(read(p)['caps']):caps.append({'source':p,'source_index':i,'normal':c['normal']})
    C=np.array([c['normal']for c in caps]);norms=np.einsum('ij,jk,ik->i',C,Q,C)
    assert max(abs(norms-2))<1e-5
    assert np.all(gate(-(C@Q@R.T)/np.sqrt(norms[:,None]*np.einsum('ij,jk,ik->i',R,Q,R))))
    assert np.max(C@Q@h)<0
    compat=gate(-(C@Q@C.T)/np.sqrt(norms[:,None]*norms[None,:]));np.fill_diagonal(compat,False)
    assert np.all(compat==compat.T)
    rng=np.random.default_rng(881693);objectives=rng.normal(size=(160,21))
    relaxed=geometry(Q,np.vstack([R,C]),h,objectives);print({'relaxed':relaxed['status'],'norm':relaxed.get('norm')},flush=True)
    checks=[];edge_offset=len(read(OLD)['caps'])
    seeds=[edge_offset+1]+[edge_offset+i for i in(2,8,9,17,18,25,26,44)]+[None]*max(0,rounds-9)
    for iteration,seed in enumerate(seeds[:rounds]):
        allowed=np.ones(len(C),dtype=bool);selected=[]
        if seed is not None:selected=[seed];allowed&=compat[seed]
        while np.any(allowed):
            ids=np.flatnonzero(allowed);deg=np.sum(compat[np.ix_(ids,ids)],axis=1)
            j=int(ids[np.argmax(deg+rng.uniform(0,3,len(ids)))]);selected.append(j);allowed&=compat[j]
        assert all(compat[i,j]for a,i in enumerate(selected)for j in selected[:a])
        g=geometry(Q,np.vstack([R,C[selected]]),h,objectives[:80])
        checks.append({'seed':seed,'cap_indices':selected,'walls':87+len(selected),'geometry':g})
        print({'selection':iteration,'walls':87+len(selected),'status':g['status'],'norm':g.get('norm')},flush=True)
    out={'status':'NUMERICAL_MIXED_DIRECTION_BANK_COMPLETION_SCOUT','sources':{p:sha(p)for p in sources},
         'retained_walls':87,'caps':caps,'compatible_pairs':int(np.sum(compat)//2),'whole_bank_geometry':relaxed,'selections':checks,
         'scope':'Combine1593old-vertex height proposals,76 initial changed18-short cutters, and explicitly listed later targeted families; test full relaxation and bounded greedy compatible selections, including alternative choices to the88th cap. Numerical angles/dedup and LP geometry only; no exact volume or all-angle completeness claim.'}
    (ROOT/out_path).write_text(json.dumps(out,indent=2)+'\n');return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--rounds',type=int,default=24);p.add_argument('--extra',action='append',default=[]);p.add_argument('--output',default=OUT);a=p.parse_args();run(a.rounds,a.extra,a.output)
