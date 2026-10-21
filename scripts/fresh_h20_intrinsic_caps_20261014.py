#!/usr/bin/env python3
"""Intrinsic upper-half-space cap construction scout, with no ambient rationality gate.

This file uses floating geometry for proposals and escape tests. Exact verification
is required before claiming a finite-volume Coxeter polytope.
"""
from pathlib import Path
from itertools import product
import argparse,json,hashlib,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_caps_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
SEED='results/fresh_h20_quartic_orbit_cap_gate_20261014.json'

def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def data(s):
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64)
 seed=json.loads((ROOT/SEED).read_text())['seed_cusp'];factors=[seed['unchanged_affine_factors'][0]['walls'],[23,50,200],seed['unchanged_affine_factors'][1]['walls']]
 walls=sum(factors,[]);G=W[walls]@Q@W[walls].T;norm=np.sqrt(np.diag(G));U=G/np.outer(norm,norm);U[12:15,12:15]=1.5*np.eye(3)-.5*np.ones((3,3));radii=[1.,s,1.]
 vectors=np.zeros((22,19));offset=0;row=0
 for f in factors:
  n=len(f);g=U[row:row+n,row:row+n];vals,vec=np.linalg.eigh(g);assert abs(vals[0])<1e-10 and min(vals[1:])>0
  vectors[row:row+n,offset:offset+n-1]=vec[:,1:]*np.sqrt(vals[1:]);offset+=n-1;row+=n
 assert np.max(abs(vectors@vectors.T-U))<1e-10
 marks=[np.ones(12),np.ones(3),np.array([3,1,1,2,1,2,2],dtype=float)]
 assert all(np.linalg.norm(m@vectors[sum(map(len,factors[:i])):sum(map(len,factors[:i+1]))])<1e-10 for i,m in enumerate(marks))
 b=np.concatenate([np.full(len(f),r)for f,r in zip(factors,radii)])
 vertices=[]
 for omitted in product(*[range(sum(map(len,factors[:i])),sum(map(len,factors[:i+1])))for i in range(3)]):
  keep=[i for i in range(22)if i not in omitted];z=np.linalg.solve(vectors[keep],b[keep]);d=b-vectors@z
  assert min(d)>-1e-9 and max(abs(d[keep]))<1e-9
  vertices.append({'omitted':list(omitted),'omitted_walls':[walls[i]for i in omitted],'retained':keep,'point':z,'heights':d[list(omitted)]})
 return {'factors':factors,'walls':walls,'unit_Gram':U,'vectors':vectors,'offsets':b,'marks':marks,'vertices':vertices}

def compatible_ratios(x,cosines,tol=1e-8):
 if min(x)<-tol:return False
 return all(v>=1-tol or min(abs(v-c)for c in cosines)<tol for v in x)
def allright_caps(D,labels):
 C=np.array([math.cos(math.pi/m)for m in labels]+[1.]);caps=[]
 for vi,v in enumerate(D['vertices']):
  A=v['heights'];radii=[]
  for a,c in product(A,C):
   if c<1e-10:continue
   R=a/c
   if compatible_ratios(A/R,C) and not any(abs(R-z)<1e-8 for z in radii):radii.append(R)
  for R in sorted(radii):
   caps.append({'vertex':vi,'omitted':v['omitted_walls'],'center':v['point'],'radius':R,'labels':{},'profile':'all_right','center_escape_margin':float(v['point']@v['point']-R*R)})
 return caps

def pair_compat(caps,labels):
 C=np.array([math.cos(math.pi/m)for m in labels]+[1.]);W=np.array([v['center']for v in caps]);rr=np.array([v['radius']for v in caps]);dist=np.maximum(np.sum(W*W,axis=1)[:,None]+np.sum(W*W,axis=1)[None,:]-2*W@W.T,0)
 cos=(dist-rr[:,None]**2-rr[None,:]**2)/(2*rr[:,None]*rr[None,:]);ok=cos>=1-1e-8
 for c in C:ok|=abs(cos-c)<1e-8
 np.fill_diagonal(ok,False);return ok,cos

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--triangle-radius',type=float,default=1.9378680925704);ap.add_argument('--max-label',type=int,default=12);args=ap.parse_args()
 D=data(args.triangle_radius);labels=list(range(2,args.max_label+1));caps=allright_caps(D,labels);ok,cos=pair_compat(caps,labels)
 # Independent small caps always exist: take one common radius below half every
 # vertex separation and every opposite-facet height. Then all cap pairs and
 # omitted-wall pairs are strictly divergent, and all19 retained contacts right.
 V=np.array([v['point']for v in D['vertices']]);ds=np.linalg.norm(V[:,None,:]-V[None,:,:],axis=2);ds[ds<1e-8]=np.inf
 small=.4*min(float(ds.min()),min(float(v['heights'].min())for v in D['vertices']))
 assert min(float(np.sum(z*z)-small*small)for z in V)>0
 out={'status':'NUMERICAL_INTRINSIC_CUSP_VERTEX_CAP_CONSTRUCTION','sources':{p:sha(p)for p in[RAW,SEED,f'scripts/{TAG}.py']},'triangle_inradius':args.triangle_radius,'other_factor_inradii':[1,1],'finite_labels_tested':labels,'cusp_vertices':len(V),'allright_candidate_caps':len(caps),'candidate_pairs':len(caps)*(len(caps)-1)//2,'compatible_candidate_pairs':int(ok.sum()//2),'minimum_center_escape_margin':min(v['center_escape_margin']for v in caps),'small_cap_radius':small,'small_cap_count':len(V),'small_caps_close_all_cusp_edges_but_leave_incenter_escape':True,'vertex_norm_range':[float(np.sum(V*V,axis=1).min()),float(np.sum(V*V,axis=1).max())],'caps':[{**v,'center':v['center'].tolist()}for v in caps],'scope':'Floating candidate construction, labels bounded. Full finite volume is not established. Small allright vertex caps close all252 old cusp edges but leave a strict ideal-boundary escape at the Euclidean incenter.'}
 p=ROOT/'results'/f'{TAG}_allright_scout.json';p.write_text(json.dumps(out,indent=2)+'\n');print({k:v for k,v in out.items()if k not in['caps','sources']})
if __name__=='__main__':main()
