#!/usr/bin/env python3
"""Bounded numerical ideal-boundary covering with lazy pair constraints.

This is discovery only. A finite cone without exterior/ideal rays would still
require exact reconstruction of all selected normals and the entire face.
"""
import argparse,json,itertools,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
from scipy.optimize import milp,Bounds,LinearConstraint
from scipy.sparse import csr_matrix
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;S=H3.S;TAG='fresh_nonarithmetic_h3_joint_cover_20261016'
TRIAL='results/fresh_nonarithmetic_h3_closure_20261016_trial.json'
REPAIR='results/fresh_nonarithmetic_h3_closure_20261016_repair_scout.json'
NEW='results/fresh_nonarithmetic_nonright_s6_20261016_caps.json'
ETA=np.diag([-1.,1.,1.,1.]);MAXLABEL=10
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return[i for i in range(m.bit_length())if m>>i&1]
def fit(g,tol=2e-5):
 if g<-2-1e-5:return True
 if abs(g)<tol:return True
 if not(-2+1e-5<g<-1+tol):return False
 m=math.pi/math.acos(-g/2)
 return round(m)>=3 and abs(m-round(m))<tol
def setup():
 Q,N,H,E=P.geometry();old=P.P.B.read(H3.CAP);tr=P.P.B.read(TRIAL);repair=P.P.B.read(REPAIR);new=P.P.B.read(NEW)
 N.append(mp.matrix([P.mid(x)for x in old['normal_interval']]))
 N.extend(mp.matrix([mp.mpf(x)for x in row])for row in tr['normal_intervals_or_proposals'])
 extra=[mp.matrix([P.mid(x)for x in r['normal_interval']])for r in new['caps']]
 W=mp.matrix([[N[j][i]for j in S]for i in range(10)]);J=(W.T*Q*W)**-1
 T=H-W*J*W.T*Q*H;T/=mp.sqrt(-(T.T*Q*T)[0]);basis=[T]
 for unused in range(3):
  choices=[]
  for n in N:
   v=n-W*J*W.T*Q*n;v+=((v.T*Q*T)[0])*T
   for z in basis[1:]:v-=((v.T*Q*z)[0])*z
   choices.append(((v.T*Q*v)[0],v))
  q,v=max(choices,key=lambda p:p[0]);assert q>mp.mpf('1e-30');basis.append(v/mp.sqrt(q))
 def coord(v):return np.array([float(-(v.T*Q*basis[0])[0])]+[float((v.T*Q*z)[0])for z in basis[1:]])
 NC=np.array([coord(n)for n in N]);XC=np.array([coord(n)for n in extra]);
 assert min(NC[i,0]for i in range(len(N))if i not in S)>0 and min(XC[:,0])>0
 base_rel=[]
 for i,n in enumerate(N[53:],55):
  for j,e in enumerate(extra,53):
   g=(n.T*Q*e)[0];base_rel.append(dict(trial_wall=i,new_wall=j,pairing=mp.nstr(g,45),compatible=fit(float(g))))
 assert all(r['compatible']for r in base_rel),base_rel
 cache={};pool=[];vectors=[];rejected=Counter()
 for idx,r in enumerate(repair['candidates']):
  if max(r['labels'])>MAXLABEL:continue
  f=tuple(S+r['anchors'])
  if f not in cache:
   F=mp.matrix([[N[j][i]for j in f]for i in range(10)]);Jf=(F.T*Q*F)**-1
   k=H-F*Jf*F.T*Q*H;qk=(k.T*Q*k)[0]
   cache[f]=(np.array(Jf.tolist(),float),np.array([coord(N[i])for i in f]),coord(k),float(qk))
  j,ns,k,qk=cache[f];c=np.array([0.]*6+[-2*math.cos(math.pi/m)if m!=2 else 0 for m in r['labels']]);pc=c@j;p=pc@ns;rho=(2-pc@c)/qk
  assert rho>0
  v=p+r['branch']*math.sqrt(rho)*k;g=XC@ETA@v
  if not all(fit(x)for x in g):rejected['new55_contact_gap']+=1;continue
  oldg=NC@ETA@v
  if np.max(oldg)>2e-5 or not all(fit(x)for x in oldg):rejected['reconstruction_gate']+=1;continue
  pool.append(dict(repair_index=idx,target=r['target'],anchors=r['anchors'],labels=r['labels'],branch=r['branch'],new53_54_products=g.tolist()));vectors.append(v)
 rays=[]
 for r in tr['rays']:
  v=coord(mp.matrix([mp.mpf(x)for x in r['point']]));v/=v[0]
  # Trial indices53..58 become55..60 after inserting the two exact cutters.
  active=[i+2 if i>=53 else i for i in r['active']]
  rays.append(dict(mask=sum(1<<i for i in active),v=v,kind=r['kind']))
 return np.array(vectors),pool,rays,base_rel,rejected
def null_point(ray):
 u=ray['v'][1:];return np.r_[1.,u/np.linalg.norm(u)]
def clip(current,n,wall):
 # Scale cap products before sign tests; exactification must settle all zeros.
 cov=ETA@n;scale=max(abs(n));g=np.array([r['v']@cov/scale for r in current]);tol=2e-11
 neg=np.flatnonzero(g<-tol);pos=np.flatnonzero(g>tol);zero=np.flatnonzero(abs(g)<=tol)
 out={current[i]['mask']:current[i]for i in neg}
 for i in zero:
  row=dict(current[i]);row['mask']|=1<<wall;out[row['mask']]=row
 crossings=[]
 for i,j in itertools.product(neg,pos):
  m=current[i]['mask']&current[j]['mask']
  if m.bit_count()<8:continue
  face=[k for k,r in enumerate(current)if r['mask']&m==m]
  if face!=sorted([int(i),int(j)]):continue
  v=g[j]*current[i]['v']-g[i]*current[j]['v'];v/=v[0];q=v@ETA@v
  kind='finite_vertices'if q<-2e-10 else'hyperideal_rays'if q>2e-10 else'unresolved_null'
  nm=m|(1<<wall);assert nm not in out
  out[nm]=dict(mask=nm,v=v,kind=kind);crossings.append(bits(nm))
 return sorted(out.values(),key=lambda r:bits(r['mask'])),dict(negative=len(neg),positive=len(pos),zero=len(zero),crossings=len(crossings))
def run(iterations):
 V,pool,initial,base_rel,rejected=setup();assert len(pool)>0
 print('pool',len(pool),'rejected',dict(rejected),'base55compatible',flush=True)
 def fits(G):
  labs=np.pi/np.arccos(np.clip(-G/2,-1+1e-15,1-1e-15))
  return (G<-2-1e-5)|(abs(G)<2e-5)|((G<-1+2e-5)&(G>-2+1e-5)&(abs(labs-np.rint(labs))<2e-5))
 groups=[[i for i,r in enumerate(pool)if r['target']==t]for t in range(6)]
 fixedrows=[];compatstats=[]
 for t,ids in enumerate(groups):
  G=V[ids]@ETA@V[ids].T;ok=fits(G);np.fill_diagonal(ok,False)
  assert not ok.any(),('same-region compatible pair needs a different encoding',t)
  fixedrows.append({i:1 for i in ids})
 for a,b in itertools.combinations(range(6),2):
  aa,bb=groups[a],groups[b];G=V[aa]@ETA@V[bb].T;bad=~fits(G)
  compatstats.append(dict(targets=[a,b],tested=int(G.size),incompatible=int(bad.sum())))
  # With at most one candidate from target b, xi+sum(bad neighbors)<=1
  # is precisely the conjunction of the corresponding pair exclusions.
  for i,row in zip(aa,bad):
   js=np.flatnonzero(row)
   if len(js):fixedrows.append({i:1,**{bb[int(j)]:1 for j in js}})
 print('complete numerical compatibility rows',len(fixedrows),flush=True)
 samples=[null_point(r)for r in initial if r['kind']=='hyperideal_rays'];conflicts=set();history=[];proposal=None
 for it in range(iterations):
  rows=list(fixedrows);lo=[-np.inf]*len(rows);hi=[1]*len(rows);sample_counts=[]
  for p in samples:
   g=V@ETA@p;cover=np.flatnonzero(g>2e-7)
   rows.append({int(i):1 for i in cover});lo.append(1);hi.append(np.inf);sample_counts.append(len(cover))
  for i,j in sorted(conflicts):rows.append({i:1,j:1});lo.append(-np.inf);hi.append(1)
  rr=[];cc=[];vv=[]
  for i,row in enumerate(rows):
   for j,x in row.items():rr.append(i);cc.append(j);vv.append(x)
  matrix=csr_matrix((vv,(rr,cc)),shape=(len(rows),len(pool)))
  ans=milp(c=np.ones(len(pool)),integrality=np.ones(len(pool)),bounds=Bounds(0,1),constraints=LinearConstraint(matrix,lo,hi),options=dict(time_limit=45,mip_rel_gap=0))
  rec=dict(iteration=it,samples=len(samples),conflicts=len(conflicts),sample_cover_counts=sample_counts,solver_status=int(ans.status),solver_message=ans.message)
  if ans.x is None:history.append(rec);print(rec,flush=True);break
  chosen=np.flatnonzero(ans.x>.5).tolist();rec['selected_pool_indices']=chosen;rec['selected_recipes']=[pool[i]for i in chosen]
  newbad=[]
  for i,j in itertools.combinations(chosen,2):
   g=V[i]@ETA@V[j]
   if not fit(g):conflicts.add((i,j));newbad.append(dict(pair=[i,j],product=float(g)))
  rec['new_pair_conflicts']=newbad
  if newbad:history.append(rec);print('iteration',it,'chosen',len(chosen),'newconflicts',len(newbad),flush=True);continue
  current=[dict(r)for r in initial];stages=[]
  for k,i in enumerate(chosen,61):current,stage=clip(current,V[i],k);stages.append(stage)
  counts=dict(Counter(r['kind']for r in current));rec['face_counts']=counts;rec['clipping_stages']=stages
  rec['exterior_supports']=[bits(r['mask'])for r in current if r['kind']!='finite_vertices'];history.append(rec)
  print('iteration',it,'chosen',len(chosen),'counts',counts,flush=True)
  if not counts.get('hyperideal_rays',0)and not counts.get('unresolved_null',0):proposal=rec;break
  added=0
  for r in current:
   if r['kind']!='hyperideal_rays':continue
   p=null_point(r)
   if all(np.linalg.norm(p-q)>1e-11 for q in samples):samples.append(p);added+=1
  if not added:break
 paths=[f'scripts/{TAG}.py',f'scripts/{H3.TAG}.py',H3.CAP,H3.CONE,TRIAL,REPAIR,NEW]
 return dict(status='BOUNDED_NUMERICAL_IDEAL_BOUNDARY_COVER_SEARCH',sources={p:sha(p)for p in paths},max_anchor_label=MAXLABEL,pool_count=len(pool),rejected=dict(rejected),initial_trial_extra55_contacts=base_rel,pool=pool,complete_pair_compatibility=compatstats,compatibility_rows=len(fixedrows),initial_face_counts=dict(Counter(r['kind']for r in initial)),history=history,compact_face_proposal=proposal,
 scope='Discovery only. The six prior trial sides are retained and checked against both newexact55walls. Candidate labels2..10 on the six trial exterior-ray anchor triples are screened against all61retained normals; selected mutual incompatibilities become binary exclusions. Every coherent proposed cover is checked by complete numerical4Dface clipping, and any remaining exterior ray supplies a new ideal sample. Solver infeasibility excludes only this finite numerical model; positive compactness requires independent exact reconstruction, affine filtering, and cone certification.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--iterations',type=int,default=60);ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run(a.iterations);p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved numerical cover replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'proposal',d['compact_face_proposal']is not None)
