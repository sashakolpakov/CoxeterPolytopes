#!/usr/bin/env python3
"""Exact complete symmetry-orbit repair of the H22 k=6 cusp links.

Native integer operations have checked overflow bounds. Fraction/SymPy derive
local cusp geometry. No numerical scout or producer is imported.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from math import isqrt
from pathlib import Path
import numpy as np
import sympy as sp
P=Path(__file__).resolve().parents[1];PRE='fresh_h22_new_cusp_completion_20261013'
OUT='results/'+PRE+'_verify.json';PRIOR='results/fresh_h22_gap_face_family_20261013.json'
FRAMES='data/enumerations/'+PRE+'_frames.npz';HIGHEST='data/enumerations/'+PRE+'_highest_roots.npz'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';SYMS='data/enumerations/fresh_h20_global_interior_sections_20261013.json'
DEN=5208;SHIFT=10

def read(n):return json.loads((P/n).read_text())
def sha(n):return hashlib.sha256((P/n).read_bytes()).hexdigest()
def ff(x):return F(int(x.p),int(x.q))
def sqrt_upper(q):
 k=isqrt(q.numerator//q.denominator)
 return k if k*k==q else k+1

def setup():
 prior=read(PRIOR);assert all(sha(n)==h for n,h in prior['sources'].items())
 old=read('results/fresh_h22_cusp_orbits_20261013_verify.json');assert all(sha(n)==h for n,h in old['sources'].items())
 olddata=read('data/enumerations/fresh_h22_cusp_orbits_20261013.json');assert all(sha(n)==h for n,h in olddata['sources'].items())
 oc=read('results/fresh_h22_symmetry_orbit_verify_20261013.json');assert all(sha(n)==h for n,h in oc['source_sha256'].items())
 raw=read(RAW)['coxeter_input'];Q=-np.array(raw['basis_gram'],dtype=np.int64);W=np.array(raw['wall_coordinates_as_columns_transposed'],dtype=np.int64);h=W.sum(0);h//=np.gcd.reduce(h)
 meta=read(SYMS)['ambient_symmetry_generators'];gens=[np.array(g['matrix'],dtype=np.int64)for g in meta]
 z=np.load(P/FRAMES);D=z['six_cusp_old'];V=z['first_tangent_numerator'];B=z['second_tangent_numerator']
 z=np.load(P/HIGHEST);A=z['old_bases'];ids=z['cusp_indices'];sizes=z['component_sizes']
 return prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes

def families(D,V,B,A,ids):
 N=len(D)
 yield 'highest_roots',A*DEN,ids
 for name,rows in [('first+',V*(DEN//28)),('first-',-V*(DEN//28)),('second+',B),('second-',-B)]:
  yield name,rows,np.arange(N)

def run():
 prior,Q,W,h,meta,gens,D,V,B,A,ids,sizes=setup();N=len(D)
 assert D.dtype==V.dtype==B.dtype==A.dtype==np.int64 and D.shape==V.shape==B.shape==(120960,22)
 assert A.shape==(7*N,22) and np.array_equal(ids,np.repeat(np.arange(N),7))
 maxval=max(int(abs(x).max())for x in(Q,W,D,V,B,A,*gens));assert 22**2*maxval**3<2**63
 assert h@Q@h==-14
 seed=np.array(prior['six_old_projection'],dtype=np.int64);partner=np.array(prior['new_cusp']['six_partner'],dtype=np.int64)
 assert np.array_equal(D[0],seed+partner)
 lookup={v.tobytes():i for i,v in enumerate(D)};assert len(lookup)==N
 reached=np.zeros(N,dtype=bool);reached[0]=True
 for g,m in zip(gens,meta):
  assert np.array_equal(g.T@Q@g,Q)and np.array_equal(g@h,h)
  assert np.array_equal(W@g.T,W[m['wall_permutation']])
  perm=np.array([lookup[x.tobytes()]for x in D@g.T]);assert len(np.unique(perm))==N
  assert np.array_equal(V@g.T,V[perm])and np.array_equal(B@g.T,B[perm])
  mask=perm>np.arange(N);reached[perm[mask]]=True
 assert np.all(reached)
 # Reconstruct both residual tangent directions independently from the old face.
 common=prior['new_cusp']['common_old_walls'];SM=sp.Matrix(Q.tolist())
 ns=sp.Matrix(np.vstack([W[common]@Q,seed@Q,partner@Q]).tolist()).nullspace();assert len(ns)==2
 v=ns[0];w=ns[1]-v*((v.T*SM*ns[1])[0]/(v.T*SM*v)[0])
 assert [ff(x)for x in v]==[F(int(x),28)for x in V[0]]
 assert [ff(x)for x in w]==[F(int(x),DEN)for x in B[0]]
 assert V[0]@Q@V[0]==15624 and B[0]@Q@B[0]==23870 and V[0]@Q@B[0]==0
 assert D[0]@Q@D[0]==-1584 and h@Q@D[0]==-432
 assert seed@Q@seed==partner@Q@partner==-372 and seed@Q@partner==-420
 assert (seed-partner)@Q@(seed-partner)==96
 assert all(x@Q@y==0 for x in(V[0],B[0])for y in(seed,partner))
 assert sp.Matrix((W[common]@Q@W[common].T).tolist()).rank()==18
 # Reconstruct every highest root from every cusp's complete incident old set.
 G=W@Q@W.T;adj=[sum(1<<int(j)for j in np.flatnonzero(row))for row in G]
 component_records=[]
 for start in range(0,N,2048):
  pairs=D[start:start+2048]@Q@W.T
  for off,row in enumerate(pairs):
   js=np.flatnonzero(row==0).tolist();assert len(js)==18 and max(row)==0
   unseen=sum(1<<j for j in js);local=[]
   while unseen:
    bit=unseen&-unseen;unseen^=bit;todo=bit;component=0
    while todo:
     bit=todo&-todo;todo^=bit;j=bit.bit_length()-1;component|=bit
     added=unseen&adj[j];unseen^=added;todo|=added
    cs=[]
    while component:
     bit=component&-component;component^=bit;cs.append(bit.bit_length()-1)
    local.append(cs);i=7*(start+off)+len(local)-1
    assert np.array_equal(A[i],-W[cs].sum(0)) and sizes[i]==len(cs)
    # A finite chain, or a singleton of norm two/four: highest root is the sum.
    sub=G[np.ix_(cs,cs)]
    if len(cs)==1:assert sub[0,0]in(2,4)
    else:
     assert np.all(np.diag(sub)==2)
     assert set(sub.ravel())<={-1,0,2}
     degree=np.sum(sub==-1,axis=0)
     assert sorted(degree.tolist())==[1,1]+[2]*(len(cs)-2)
   assert sorted(map(len,local))==[1,1,1,3,3,4,5]
   if start+off==0:component_records=local
 # Representatives of all seven highest-root families and four tangent families.
 bases=np.vstack([A[:7]*DEN,V[0]*(DEN//28),-V[0]*(DEN//28),B[0],-B[0]])
 C=bases+SHIFT*(DEN//6)*D[0]
 cov=C@Q
 # Bound native dot products before using large arrays.
 normnums=np.einsum('ij,ij->i',C,cov)+44*SHIFT**2*DEN**2
 norms=[F(int(x),DEN**2)for x in normnums]
 assert norms[:7]==[F(2)]*6+[F(4)] and norms[7:]==[F(279,14)]*2+[F(55,62496)]*2
 assert all(0<x<20 for x in norms)and max(C@(Q@h))<0
 local=bases@Q@bases.T
 for i in range(11):
  for j in range(i):
   assert local[i,j]==(-int(normnums[i])if(i,j)in((8,7),(10,9))else 0)
 extrema={}
 for name,base,ix in families(D,V,B,A,ids):
  roots=base+SHIFT*(DEN//6)*D[ix]
  bound=22*int(abs(roots).max())*int(abs(cov).max())+44*SHIFT**2*DEN**2
  assert bound<2**63
  pair=roots@cov.T+44*SHIFT**2*DEN**2
  maximum=pair[ix!=0].max(axis=0)
  assert np.all(maximum < -20*DEN**2)
  for root,row in zip(roots[ix==0],pair[ix==0]):
   match=np.flatnonzero(np.all(C==root,axis=1));assert len(match)==1
   assert np.array_equal(row,local[match[0]])
  extrema[name]=[str(F(int(t),DEN**2))for t in maximum]
 # By full orbit transfer, these representative gates cover all new/old pairs.
 oldnorm=np.diag(G);pair=W@cov.T;oldkinds=[]
 for j,q in enumerate(norms):
  count=Counter()
  for value,n in zip(pair[:,j],oldnorm):
   p=F(int(value),DEN);ratio=p*p/(q*int(n))
   if p==0:kind='right'
   elif p<0 and ratio==F(1,4):kind='pi/3'
   elif p<0 and ratio==F(1,2):kind='pi/4'
   elif p<0 and ratio==F(3,4):kind='pi/6'
   elif p<0 and ratio==1:kind='parallel'
   else:assert p<0 and ratio>1;kind='dotted'
   count[kind]+=1
  oldkinds.append(dict(count))
 U=np.load(P/'data/enumerations/fresh_h22_symmetry_orbit_20261013.npz')['twice_old_projections']
 z=np.load(P/'data/enumerations/fresh_h22_cusp_orbits_20261013.npz');K=z['old_null_projections'];four=z['four_heights']
 R=np.load(P/'data/enumerations/fresh_h22_gap_face_orbit_20261013.npz')['six_old_projections']
 maxima=[]
 for j,q in enumerate(norms):
  cc=cov[j];upper=sqrt_upper(q)
  assert 22*max(int(abs(x).max())for x in(U,K,R))*int(abs(cc).max())<2**63
  mu=int(max(U@cc));assert mu+21*SHIFT*DEN < -2*DEN*upper # sqrt110<21/2
  kp=K@cc;nullmax={}
  for f in(1,4):
   value=int(max(kp[four==f]));assert 12*value+20*SHIFT*DEN*f < -3*DEN*f*upper # 2sqrt11<20/3
   nullmax[str(f)]=value
  hp=int(h@cc);assert 3*hp+80*SHIFT*DEN < -3*DEN*sqrt_upper(2*q)
  rp=R@cc+132*SHIFT*DEN;zero=np.flatnonzero(rp==0).tolist()
  assert zero==[0,prior['new_cusp']['partner_index']]
  non=rp[rp!=0];assert max(non)<-6*DEN*sqrt_upper(F(2,3)*q)
  maxima.append({'twice_old_timelike_pair_numerator':mu,'old_null_pair_numerators_by_four_height':nullmax,'center_pair_numerator':hp,'new_k6_nonincident_max_numerator':int(max(non))})
 # Each old A component plus its negative highest root has a positive all-one
 # affine relation M*delta. The existing parallel pair sums to delta. The two
 # residual pairs sum to 2M*delta. Orthogonal positive ranks add to 18+1+2=21.
 old_escape=np.array([F(x)for x in prior['uniform_strict_escape']['old_coordinates']],dtype=object)
 cuts=[sum(F(int(a),DEN)*b for a,b in zip(cc,old_escape))+44*SHIFT for cc in cov]
 assert max(cuts)>0
 return {'status':'PASS_EXACT_COMPLETE_H22_NEW_CUSP_LINK_ORBIT_REPAIR',
  'sources':{n:sha(n)for n in(PRIOR,RAW,SYMS,FRAMES,HIGHEST)},'shift':SHIFT,'common_denominator':DEN,
  'positive_cusp_orbit':N,'completed_signed_cusps':2*N,'new_positive_walls':11*N,'new_signed_walls':22*N,'total_facets':prior['total_facets']+22*N,
  'formula':'For each cusp delta=(D/6,+/-2sqrt11), add base+10delta for seven negative highest roots and four signed residual tangent directions.',
  'seed_cusp_six_old':D[0].tolist(),'old_components':component_records,'base_old_numerators':bases.tolist(),'new_cap_old_numerators':C.tolist(),'new_cap_norms':list(map(str,norms)),
  'old210_pair_counts':oldkinds,'new_pair_nonlocal_maxima':extrema,'prior_family_bounds':maxima,
  'affine_link_type':'affine A5 + affine A4 + 2 affine A3 + 6 affine A1','component_ranks':[5,4,3,3,1,1,1,1,1,1],'horosphere_rank':21,'incident_walls_per_new_cusp':31,
  'old_escape_new_pairings':list(map(str,cuts)),'old_escape_cut':True,
  'scope':'Exact compatible complete orbit at new label k=6 and every old center label m>=3; all241920 new cusp links bounded. Actual facets follow from strict center and acute projection. No finite-volume conclusion without a separate whole-chamber certificate.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();result=run()
 if a.verify:assert result==read(OUT)
 else:(P/OUT).write_text(json.dumps(result,indent=2)+'\n')
 print(result['status'],result['total_facets'],result['completed_signed_cusps'])
