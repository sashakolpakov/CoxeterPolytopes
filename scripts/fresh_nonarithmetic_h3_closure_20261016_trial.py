#!/usr/bin/env python3
"""High-precision complete three-face clipping for a concrete six-cap trial."""
import argparse,itertools,json
from collections import Counter
from hashlib import sha256
import mpmath as mp
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_trial';S=H3.S
RECIPES=[([7,30,49],[2,3,3],-1),([8,44,52],[2,10,5],-1),([30,47,55],[2,2,3],-1)]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def bits(m):return[i for i in range(m.bit_length())if m>>i&1]
def run():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);sc=P.P.B.read(f'results/{H3.TAG}.json');N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));hc=Q*H
 for i in[1,3,6]:
  ids=sc['proposals'][i]['active'];v=mp.matrix(P.P.null([list(Q*N[j])for j in ids]));v/=-(v.T*hc)[0];v*=mp.sqrt(2/(v.T*Q*v)[0]);N.append(v)
 for anchors,labels,branch in RECIPES:
  frame=S+anchors;W=mp.matrix([[N[i][j]for i in frame]for j in range(10)]);J=(W.T*Q*W)**-1;c=mp.matrix([0]*6+[-2*mp.cos(mp.pi/m)if m!=2 else 0 for m in labels]);p=W*J*c;k=H-W*J*W.T*Q*H;qk=(k.T*Q*k)[0];rho=(2-(p.T*Q*p)[0])/qk;assert rho>0;r=p+branch*mp.sqrt(rho)*k;assert abs((r.T*Q*r)[0]-2)<mp.mpf('1e-60');N.append(r)
 cov=[Q*n for n in N];relations=[]
 for i in range(53,len(N)):
  row=[]
  for j in range(i):
   g=(N[i].T*cov[j])[0];rel=P.relation(g);assert rel['kind']in['finite','dotted'],(i,j,rel);row.append(dict(wall=j,pairing=mp.nstr(g,45),**rel))
  relations.append(dict(wall=i,relations=row))
 current=[]
 for kind,sets in sc['face_census'].items():
  for ids in sets:
   v=mp.matrix(P.P.null([list(cov[j])for j in ids]));v/=-(v.T*hc)[0];current.append(dict(mask=sum(1<<j for j in ids),v=v,kind=kind))
 stages=[]
 def sg(x):return 1 if x>mp.mpf('1e-28')else-1 if x<-mp.mpf('1e-28')else 0
 def kind(v):
  q=(v.T*Q*v)[0];assert abs(q)>mp.mpf('1e-25'),('unresolved norm',mp.nstr(q,40));return'finite_vertices'if q<0 else'hyperideal_rays'
 for wall in range(53,len(N)):
  gs=[(r['v'].T*cov[wall])[0]for r in current];neg=[i for i,g in enumerate(gs)if sg(g)<0];pos=[i for i,g in enumerate(gs)if sg(g)>0];zero=[i for i,g in enumerate(gs)if sg(g)==0];new={current[i]['mask']:current[i]for i in neg}
  for i in zero:
   rr=dict(current[i]);rr['mask']|=1<<wall;assert rr['mask']not in new;new[rr['mask']]=rr
  crosses=[]
  for i,j in itertools.product(neg,pos):
   mask=current[i]['mask']&current[j]['mask']
   if mask.bit_count()<8:continue
   face=[k for k,r in enumerate(current)if r['mask']&mask==mask]
   if face!=sorted([i,j]):continue
   v=gs[j]*current[i]['v']-gs[i]*current[j]['v'];v/=-(v.T*hc)[0];nm=mask|(1<<wall);assert nm not in new;new[nm]=dict(mask=nm,v=v,kind=kind(v));crosses.append(bits(nm))
  current=sorted(new.values(),key=lambda r:bits(r['mask']));stages.append(dict(wall=wall,negative=len(neg),positive=len(pos),zero=len(zero),crossings=crosses,counts=dict(Counter(r['kind']for r in current))));print(stages[-1],flush=True)
 rays=[dict(active=bits(r['mask']),kind=r['kind'],norm=mp.nstr((r['v'].T*Q*r['v'])[0],50),point=[mp.nstr(x,65)for x in r['v']])for r in current]
 polars=[]
 for i,r in enumerate(current):
  if r['kind']!='hyperideal_rays':continue
  v=r['v']*mp.sqrt(2/(r['v'].T*Q*r['v'])[0]);rel=[dict(wall=j,**P.relation((v.T*c)[0]))for j,c in enumerate(cov)];bad=[a for a in rel if a['kind']not in['finite','dotted']];polars.append(dict(ray_index=i,active=bits(r['mask']),bad=bad,relations=rel,point=[mp.nstr(x,65)for x in v]));print('next polar',bits(r['mask']),'bad',bad,flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{H3.TAG}.py',H3.CAP,H3.CONE,f'results/{H3.TAG}.json',P.CL,P.CO,P.P.BATCH]
 return dict(status='NUMERICAL_COMPLETE_H3_FACE_SIX_CAP_TRIAL',sources={p:sha(p)for p in paths},white_support=S,polar_seed_indices=[1,3,6],nonpolar_recipes=RECIPES,normal_intervals_or_proposals=[[mp.nstr(x,65)for x in n]for n in N[53:]],relations=relations,clipping_stages=stages,counts=dict(Counter(r['kind']for r in current)),rays=rays,next_polars=polars,scope='Numerical high-precision complete three-face clipping only. All53 ambient old-wall gates and all new pair gates are checked numerically. Exact positivity, affine-circuit filtering, actual facets and complete ambient cone remain separate requirements; no compactness or finite-volume claim is made.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved H3 six-cap trial')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
