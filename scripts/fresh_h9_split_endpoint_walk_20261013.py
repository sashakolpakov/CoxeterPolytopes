#!/usr/bin/env python3
"""Numerical spherical-edge walk from the selected split-seed chamber."""
import argparse,hashlib,itertools,json,os
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_endpoint_selected_geometry_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_endpoint_walk_20261013_scout.json'

def run(all_anchors=False):
 Path('/tmp/fresh_h9_split_endpoint_walk_20261013.pid').write_text(str(os.getpid()))
 old=json.loads(SOURCE.read_text());W=np.array(old['physical_normals']);N=np.array(old['N']);eta=np.diag([-1]+[1]*9);G=W@eta@W.T;p=(1+5**.5)/2
 finite={i for row in old['finite_vertices']for i in row['active']if i>=16};mandatory=list(range(16))+sorted(finite)
 edge_sets=set()
 for v in old['finite_vertices']:
  for removed in v['active']:
   block=tuple(sorted(set(v['active'])-{removed}))
   if len(block)==8 and sum(i>=16 for i in block)==1:edge_sets.add(block)
 allprofiles=[];alllabels=[]
 for count in range(4):
  for inds in itertools.combinations(range(8),count):
   for labels in itertools.product([3,5],repeat=count):
    row=np.zeros(8)
    for i,m in zip(inds,labels):row[i]=-1 if m==3 else-p
    allprofiles.append(row);alllabels.append(list(zip(inds,labels)))
 profiles=np.array(allprofiles);fvals=-2*np.cos(np.pi/np.array([2,3,4,5,6,8,10,12,15,16,17,20,24,30,32,34,40,48,60]));seen=set();rows=[];tested=0;seeds_pass=0;mandatory_pass=0
 escaped=np.r_[1,old['leading_escape']['point']]
 def allowed(v):return(v< -2-1e-7)|np.any(abs(v[...,None]-fvals)<1e-7,axis=-1)
 for bi,block in enumerate(sorted(edge_sets)):
  B=G[np.ix_(block,block)];assert min(np.linalg.eigvalsh(B))>1e-7;ib=np.linalg.inv(B);qs=2-np.einsum('ij,jk,ik->i',profiles,ib,profiles);valid=np.flatnonzero(qs>1e-7)
  for anchor in range(len(W)if all_anchors else 16):
   if anchor in block:continue
   A=W[list(block)+[anchor]]@eta;uu,ss,vv=np.linalg.svd(A,full_matrices=True)
   if min(ss)<1e-8:continue
   normal=vv[-1];inv=np.linalg.pinv(A);baseprofiles=np.repeat(profiles[valid],3,axis=0);targets=np.column_stack([baseprofiles,np.tile([0.,-1.,-p],len(valid))]);r0=targets@inv.T
   aa=normal@eta@normal;bb=2*r0@eta@normal;cc=np.einsum('ij,jk,ik->i',r0,eta,r0)-2;dd=bb*bb-4*aa*cc
   for sign in [-1,1]:
    inds=np.flatnonzero(dd>=-1e-8)
    if abs(aa)<1e-10:
     if sign==1:continue
     inds=inds[abs(bb[inds])>1e-8];t=-cc[inds]/bb[inds]
    else:t=(-bb[inds]+sign*np.sqrt(np.maximum(0,dd[inds])))/(2*aa)
    roots=r0[inds]+t[:,None]*normal;tested+=len(roots);vals=roots@eta@W.T
    good=np.all(allowed(vals[:,:16]),axis=1)&(roots[:,0]>1e-9);seeds_pass+=sum(good).item()
    good&=np.all(allowed(vals[:,mandatory]),axis=1);mandatory_pass+=sum(good).item()
    for local in np.flatnonzero(good):
     v=vals[local];key=tuple(np.round(v[:16],7))
     if key in seen:continue
     seen.add(key);r=roots[local]
     # An affine subset among a seed orthogonal block prevents compactness.
     affine=False
     for start in [0,8]:
      finite_values=[z*z for z in v[start:start+8]if z>=-2-1e-7]
      sums={0.}
      for value in finite_values:sums|={s+value for s in list(sums)}
      if any(abs(s-4)<1e-7 for s in sums):affine=True;break
     if affine:continue
     conflicts=np.flatnonzero(~allowed(v)).tolist()
     if len(conflicts)>4:continue
     point=W[list(block)].sum(axis=0);point=np.r_[sum(W[:,0]),np.zeros(9)]
     cov=W[list(block)]@eta@point;projected=point-(ib@cov)@W[list(block)]
     full=W[list(block)+[]];pr=projected@eta@r
     vertex=projected-pr/qs[valid[inds[local]//3]]*(r-(ib@(W[list(block)]@eta@r))@W[list(block)])
     assert vertex@eta@vertex<0
     outside=float(r@eta@escaped);idx=int(inds[local]);profile_index=int(valid[idx//3]);angle=[2,3,5][idx%3]
     rows.append({'edge':list(block),'anchor':anchor,'finite_profile':alllabels[profile_index],'anchor_label':angle,'quadratic_branch':sign,'normal_physical':r.tolist(),'normal_seed_coefficients':(np.linalg.pinv(N.T)@r).tolist(),'seed_covector':v[:16].tolist(),'conflicts':conflicts,'cuts_leading_escape':outside>1e-7,'leading_escape_pairing':outside,'finite_vertex':(vertex/vertex[0]).tolist()})
  print({'edges_done':bi+1,'candidates':len(rows),'leading_cutters':sum(r['cuts_leading_escape']for r in rows)},flush=True)
 out={'status':'NUMERICAL_TARGETED_EDGE_WALK','scope':'Finite-angle equations at32selected spherical rank8 edges, seed anchors, up to3finite attachments usingpi/3,pi/5; numerical geometry and compatibility only. All16seed and4seed-edge endpoint caps retained; up to4other optional-wall conflicts allowed.','source_snapshot':old,'source_snapshot_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'edge_count':len(edge_sets),'tested_roots':tested,'seed_compatible_roots':seeds_pass,'mandatory_compatible_roots':mandatory_pass,'candidates':rows}
 out['anchor_scope']='all selected walls'if all_anchors else'sixteen seed walls'
 out['scope']=out['scope'].replace('seed anchors',out['anchor_scope'])
 target=OUT.with_name('fresh_h9_split_endpoint_walk_optional_anchor_20261013_scout.json')if all_anchors else OUT
 target.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items()if k not in['source_snapshot','candidates','scope']},flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--all-anchors',action='store_true');a=p.parse_args();run(a.all_anchors)
