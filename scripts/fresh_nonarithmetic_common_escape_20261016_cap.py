#!/usr/bin/env python3
"""Exact nonpolar wall cutting the common escape of the48-wall chamber."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_joint_regions_20261016_batch as B
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_common_escape_20261016_cap'
BATCH=f'results/{B.TAG}.json';CONE='results/fresh_nonarithmetic_joint_regions_20261016_cone.json'
IDS=[0,1,2,4,5,6,8,10,17,27]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def decode(a):return tuple(A.B.db(x)for x in a)
def run():
 j=B.read(BATCH);cone=B.read(CONE)
 assert cone['counts']==dict(finite_vertices=172,cusps=0,hyperideal_rays=2109,opposite_sheet_rays=0)
 Q=[[A.B.db(x)for x in row]for row in j['form']]
 N=[tuple(decode(a)for a in row)for row in j['batch_normals']]
 H=tuple(decode(a)for a in j['center']);E=tuple(decode(a)for a in j['common_strict_escape']['point'])
 cov=[tuple(D.mv(Q,n))for n in N];norms=[D.dot(n,c)for n,c in zip(N,cov)]
 rows=[]
 for i in IDS:
  assert norms[i]==D.lift(2*O) and all(not any(x[1:])for x in cov[i]);rows.append([x[0]for x in cov[i]])
 assert len(C.eliminate(rows)[1])==10
 target=[-O if i in[8,10]else -(O+R.BASIS[1])/2 if i==17 else Z for i in IDS]
 inverse=C.inverse(rows);v=C.mv(inverse,target);r=tuple(D.lift(x)for x in v)
 rcov=tuple(D.mv(Q,r));q=D.dot(r,rcov);assert q==D.lift(2*O)
 assert D.sign(D.dot(H,rcov))<0 and D.sign(D.dot(E,rcov))>0
 vals=[D.dot(r,c)for c in cov];relations=[dict(wall=i,pairing=g,**B.relation(g,D.mul(q,nn)))for i,(g,nn)in enumerate(zip(vals,norms))]
 assert all(a['kind']in['right','finite','dotted']for a in relations)
 right=[a['wall']for a in relations if a['kind']=='right'];finite=[[a['wall'],a['label']]for a in relations if a['kind']=='finite']
 assert right==[0,1,2,4,5,6,27,40] and finite==[[8,3],[10,3],[17,5]]
 N=N+[r];cov=cov+[rcov];norms=norms+[q];hp=[D.dot(H,c)for c in cov];hh=D.pair(H,H,Q)
 assert D.sign(hh)<0 and all(D.sign(x)<0 for x in hp)
 facets=[]
 for i,(n,nn,hn)in enumerate(zip(N,norms,hp)):
  V=B.vecsub(B.vecscale(H,nn),B.vecscale(n,hn));vv=D.sub(D.mul(D.mul(nn,nn),hh),D.mul(nn,D.mul(hn,hn)))
  assert D.sign(vv)<0 and D.sign(D.dot(V,cov[i]))==0
  # Every other facet pairing is nn*hp[j]-hn*Gram[i,j]<0:
  # nn>0, hp[j]<0, hn<0, and Gram[i,j]<=0.
  for k,c in enumerate(cov):
   if k!=i:assert B.dot_sign(n,c)<=0
  facets.append(dict(wall=i,point=V,norm=vv))
 paths=[f'scripts/{TAG}.py',f'scripts/{B.TAG}.py',BATCH,CONE]
 return B.encode(dict(status='PASS_EXACT_NONPOLAR_49TH_WALL_CUTS_COMMON_ESCAPE',sources={p:sha(p)for p in paths},form=Q,center=H,normals=N,norms=norms,actual_facets=49,facet_witnesses=facets,
  cap=dict(index=48,linear_walls=IDS,linear_targets=target,linear_rank=10,linear_inverse=inverse,normal=r,norm=q,relations=relations,right_neighbors=right,finite_neighbors=finite,counts=dict(Counter(a['kind']for a in relations)),center_pairing=hp[48],old_escape_pairing=D.dot(E,rcov)),
  cut_escape=dict(point=E,norm=D.pair(E,E,Q),new_pairing=D.dot(E,rcov)),
  discovery=dict(bank='fresh_h9_octic_germ_scout_20261015',index=1928,recipe='8-wall endpoint profile plus right anchor4; literal10-contact linear equations above are the certificate'),
  scope='Exact compatible nonpolar cap with eight right, three finite and37 dotted contacts to the48-wall chamber. It cuts the common escape missed by every admissible current polar. All49 facets are actual. This certificate alone asserts neither complete new incidence nor finite volume; complete cone clipping is separate.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact49wall cap replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
