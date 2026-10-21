#!/usr/bin/env python3
"""Numerical larger-endpoint caps cutting the inter-layer annular gap."""
import argparse,json,time
import numpy as np
from fresh_h20_intrinsic_caps_20261014 import ROOT,TAG,sha,data
from fresh_h20_intrinsic_caps_20261014_quarter_face_scout import build
from fresh_h20_intrinsic_caps_20261014_central_contact_scout import allowed

def run(allprobes):
 t=time.time();D,V,X,Q,_,_,_,_=build('contiguous');quarter=json.loads((ROOT/f'results/{TAG}_quarter_face_scout.json').read_text());Q=Q[quarter['selected']]
 bank=f'data/enumerations/{TAG}_endpoint_bank.npz';meta=f'results/{TAG}_completion_scout.json';a=np.load(ROOT/bank);C=a['centers'];rr=a['radii'];old=json.loads((ROOT/meta).read_text());recipes=old['recipes']
 target=V*.75 if allprobes else V[[0]]*.75;indices=[];bestpower=[]
 for start in range(0,len(C),2000):
  cc=C[start:start+2000];rs=rr[start:start+2000];p=np.sum(cc*cc,axis=1)[:,None]+np.sum(target*target,axis=1)[None,:]-2*cc@target.T-rs[:,None]**2;good=np.flatnonzero(p.min(axis=1)<-1e-7);indices.extend((start+good).tolist());bestpower.extend(p[good].min(axis=1).tolist())
 labels=list(range(2,65));cs=list(np.cos(np.pi/np.array(labels)))+[1.];fixedC=np.concatenate([np.zeros((1,19)),V/2,Q]);fixedR=np.concatenate([[2],np.ones(252),np.full(len(Q),1.5)]);fixed2=np.sum(fixedC*fixedC,axis=1)
 rows=[];survive=[];best=None;counts=dict(endpoint_bank=len(C),target_cutters=len(indices),allright_target_cutters=0,central_compatible=0,half_layer_compatible=0,quarter_layer_compatible=0,all_interior_layers_compatible=0)
 for start in range(0,len(indices),1000):
  ii=np.array(indices[start:start+1000]);cc=C[ii];rs=rr[ii];dd=np.maximum(np.sum(cc*cc,axis=1)[:,None]+fixed2[None,:]-2*cc@fixedC.T,0);co=(dd-rs[:,None]**2-fixedR[None,:]**2)/(2*rs[:,None]*fixedR[None,:]);ok=allowed(co,cs)
  counts['central_compatible']+=int(ok[:,0].sum());counts['half_layer_compatible']+=int(ok[:,1:253].all(axis=1).sum());counts['quarter_layer_compatible']+=int(ok[:,253:].all(axis=1).sum());good=ok.all(axis=1);counts['all_interior_layers_compatible']+=int(good.sum())
  for j,k0 in enumerate(ii):
   k=int(k0);counts['allright_target_cutters']+=int(recipes[k][1]==0)
   conflict=np.flatnonzero(~ok[j]);row=dict(endpoint_bank_index=k,recipe=recipes[k],radius=float(rs[j]),center=cc[j].tolist(),fixed_conflict_count=len(conflict),fixed_conflicts=conflict.tolist(),fixed_conflict_cosines=co[j,conflict].tolist(),target_minimum_power=float(bestpower[start+j]))
   if best is None or(len(conflict),recipes[k][1]!=0)<(best['fixed_conflict_count'],best['recipe'][1]!=0):best=row
   if good[j]:
    vv=np.maximum(np.sum((V-cc[j])**2,axis=1),0);vp=(vv-rs[j]**2-1)/(2*rs[j]);drop=np.flatnonzero(~allowed(vp,cs));row.update(old_vertex_caps_to_replace=drop.tolist(),old_vertex_conflict_cosines=vp[drop].tolist());survive.append(row)
 out=dict(status='NUMERICAL_ANNULUS_ENDPOINT_CAP_SCOUT',sources={p:sha(p)for p in[bank,meta,f'scripts/{TAG}_annulus_endpoint_scout.py',f'scripts/{TAG}_quarter_face_scout.py',f'scripts/{TAG}_central_contact_scout.py',f'results/{TAG}_quarter_face_scout.json']},target='all252 three-quarter vertices'if allprobes else'first three-quarter vertex',cross_labels_tested=labels,counts=counts,surviving_candidates=survive,best_replacement_candidate=best,elapsed_seconds=time.time()-t,scope='Numerical finite endpoint bank. Seed,central,half andquarter layers tested; conflicting original radius-one vertex caps may be replaced. No exact root/pair/finite-volume claim.')
 path=ROOT/('results/'+TAG+'_annulus_endpoint'+('_all'if allprobes else'')+'_scout.json');path.write_text(json.dumps(out,indent=2)+'\n')
 print(counts,'bestconflicts',best['fixed_conflict_count']if best else None,'elapsed',time.time()-t,flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--all-probes',action='store_true');a=p.parse_args();run(a.all_probes)
