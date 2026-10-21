#!/usr/bin/env python3
"""Exact strict escape and certified finite vertices of the835-wall branch."""
import argparse,contextlib,io,itertools,json
from fractions import Fraction as F
import numpy as np
from scipy.optimize import linprog
import fresh_h22_incidence_prism_20261013_verify as V
from fresh_h22_borcherds_extension_20261013 import prepare
OUT='results/'+V.PRE+'_geometry.json'
def div(a,b):
 assert not b[0]<=0<=b[1]
 if b[1]<0:return div(V.neg(a),V.neg(b))
 return V.mul(a,(V.S*V.S//b[1],-((-V.S*V.S)//b[0])))
def pd(A):
 A=[[x for x in row]for row in A];piv=[]
 for i in range(len(A)):
  z=A[i][i];assert z[0]>0;piv.append(list(map(str,z)))
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]=V.add(A[j][k],V.neg(div(V.mul(A[j][i],A[i][k]),z)))
 return piv
def escape_check(nums,den,M,R,h,caps):
 x=[F(a,den)for a in nums];q=V.pair(M,x,x);assert q>0 and V.pair(M,x,h)<0;old=[V.pair(M,x,w)for w in R];assert max(old)<0
 cov=V.mv(M,x);pairs=[V.isum(V.scale(a,b)for a,b in zip(r,cov))for r in caps];assert max(p[1]for p in pairs)<0
 return{'numerator':nums,'denominator':den,'norm':str(q),'center_pairing':str(V.pair(M,x,h)),'maximum_base_pairing':str(max(old)),'new_cap_pairings':[list(map(str,p))for p in pairs]}
def run(verify=False):
 with contextlib.redirect_stdout(io.StringIO()):cert,state=V.run()
 assert cert==V.read(V.OUT);M,R,h,caps,qn=state;rn=np.array(R,dtype=float);mn=np.array(M,dtype=float);cn=np.array([[(a+b)/(2*V.S)for a,b in r]for r in caps]);normals=np.vstack([rn,cn]);roots=V.read(V.REC)['recipes'];vertices=[];ri=np.array(R,dtype=np.int64);mi=np.array(M,dtype=np.int64);assert 23**2*int(abs(mi).max())*int(abs(ri).max())**2<2**63;G=ri@mi@ri.T
 for k,rec in enumerate(roots):
  ff=[x['wall']for x in rec['finite_contacts']];count=0;cp=[tuple(map(int,x))for x in cert['new_caps'][k]['old_pair_intervals']]
  for ss in itertools.combinations(ff,21):
   ix=list(ss);gram=G[np.ix_(ix,ix)];nn=np.diag(gram)
   if np.any(np.triu(gram*gram>=nn[:,None]*nn[None,:]-1e-8,1)):continue
   pair=rn[ix]@mn@cn[k];a=np.block([[gram,pair[:,None]],[pair[None,:],np.array([[float(qn[k])]])]])
   if np.linalg.eigvalsh(a)[0]<1e-7:continue
   exact=[[V.rat(int(G[i,j]))for j in ix]+[cp[i]]for i in ix]+[[cp[i]for i in ix]+[V.rat(qn[k])]];piv=pd(exact);vertices.append({'new_cap_index':k,'candidate':rec['candidate_index'],'base_walls':ix,'positive_Gram_pivots':piv});count+=1
  print('finite vertices',rec['candidate_index'],count,flush=True)
 if verify:
  saved=V.read(OUT);escape=escape_check(saved['strict_escape']['numerator'],saved['strict_escape']['denominator'],M,R,h,caps)
 else:
  Q,W,GG,hh,X,A,b,Qn,_=prepare();B=np.zeros((23,22));B[:22,:21]=X;B[-1,-1]=1/np.sqrt(2);center=np.array(h)/np.sqrt(14);lim=-(normals@mn@center);H=(normals@mn@B)/lim[:,None];rng=np.random.default_rng(8357)
  for trial in range(100):
   lp=linprog(rng.normal(size=22),A_ub=H,b_ub=np.ones(len(H)),bounds=[(None,None)]*22,method='highs');assert lp.success
   if lp.x@lp.x>1+1e-5:
    x=center+.999*(B@lp.x);den=10**7;nums=np.rint(x*den).astype(np.int64).tolist();escape=escape_check(nums,den,M,R,h,caps);break
  else:raise RuntimeError('no escape proposal; requires full volume analysis')
 out={'status':'PASS_EXACT_835_FACET_GEOMETRY_AND_STRICT_ESCAPE','sources':{V.OUT:V.sha(V.OUT),V.REC:V.sha(V.REC)},'facets':835,'new_finite_vertices':vertices,'strict_escape':escape,'scope':'Each listed22-wall Gram is rigorously positive definite. Pairwise acuteness and the strict timelike center make its orthogonal center projection an actual finite vertex, strict on every other wall. The vertex list is not claimed exhaustive. A separate rational point is strictly feasible and spacelike for all835walls, proving infinite volume.'}
 if verify:assert out==V.read(OUT)
 else:(V.P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'vertices',len(vertices),'escape',escape['norm'],flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');run(p.parse_args().verify)
