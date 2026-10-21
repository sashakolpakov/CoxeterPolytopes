#!/usr/bin/env python3
"""Complete49-wall clipping, with rational interval DAG arithmetic.

Every node denotes an exact positive combination. Outward rational intervals
certify strict signs; ambiguous signs use the frozen formal field arithmetic.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_joint_regions_20261016_batch as B
import fresh_nonarithmetic_escape49_20261016_cluster as U
import sympy as sy
D,A,R,C,K,Z,O=B.D,B.A,B.R,B.C,B.K,B.Z,B.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape49_20261016_cone'
BATCH='results/fresh_nonarithmetic_common_escape_20261016_cap.json';PRIOR48='results/fresh_nonarithmetic_joint_regions_20261016_cone.json';PRIOR49='results/fresh_nonarithmetic_common_escape_20261016_cone.json';CLUSTER='results/fresh_nonarithmetic_escape49_20261016_cluster.json';I=R.L.I
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def decode(a):return tuple(A.B.db(x)for x in a)
def ii(a):
 z=I()
 for x,b in zip(a,D.BIV):
  if x:z=z+D.biv(x)*b
 return B.rnd(z,256)
def idot(v,w):return B.rnd(sum((a*b for a,b in zip(v,w)),I()),256)
def imv(Q,v):return[idot(row,v)for row in Q]
def power_scale(v):
 b=max(max(abs(x.l),abs(x.h))for x in v);assert b>0
 p=b.numerator.bit_length()-b.denominator.bit_length()
 return F(1,1<<p)if p>=0 else F(1<<(-p))
def bits(mask):return[i for i in range(mask.bit_length())if mask>>i&1]
def run():
 batch=B.read(BATCH);prior=B.read(PRIOR48);prior49=B.read(PRIOR49);cluster=B.read(CLUSTER);local,start,Q,H,N36,C36,reconstruct=B.geometry()
 N=[tuple(decode(x)for x in n)for n in batch['normals']];assert len(N)==49 and N[:36]==N36
 cov=[tuple(D.mv(Q,n))for n in N];hcov=tuple(D.mv(Q,H));QI=[[I(str(x))for x in row]for row in [[0]*10]*10]
 QI=[[ii(D.lift(x))for x in row]for row in Q];hI=[ii(x)for x in hcov];covI=[[ii(x)for x in c]for c in cov]
 newnormals=[[I(*x)for x in r['normal_interval']]for r in cluster['caps']];covI.extend(imv(QI,n)for n in newnormals)
 nodes=[];stats=Counter()
 def create(mask,kind,vec=None,parents=None,wall=None):
  if vec is not None:vI=[ii(x)for x in vec]
  else:
   a,b=parents;ga=idot(nodes[a]['iv'],covI[wall]);gb=idot(nodes[b]['iv'],covI[wall])
   # Exact opposite signs have already been certified. Tightening their
   # intervals is needed only if cancellation hid an endpoint from interval.
   if not ga.h<0:ga=ii(value(a,wall))
   if not gb.l>0:gb=ii(value(b,wall))
   assert ga.h<0<gb.l
   vI=[B.rnd(gb*x-ga*y,256)for x,y in zip(nodes[a]['iv'],nodes[b]['iv'])]
  s=power_scale(vI);vI=[B.rnd(s*x,256)for x in vI]
  row=dict(mask=mask,kind=kind,iv=vI,parents=parents,wall=wall,scale=s)
  row['exact']=tuple(D.scale(x,K.convert(s))for x in vec)if vec is not None else None
  nodes.append(row);return len(nodes)-1
 @lru_cache(maxsize=None)
 def value(i,wall):
  row=nodes[i]
  if row['mask']>>wall&1:return D.ZERO
  if row['exact']is not None:return D.dot(row['exact'],cov[wall])
  a,b=row['parents'];g1=value(a,row['wall']);g2=value(b,row['wall'])
  return D.scale(D.sub(D.mul(g2,value(a,wall)),D.mul(g1,value(b,wall))),K.convert(row['scale']))
 def sign_value(i,wall):
  if nodes[i]['mask']>>wall&1:return 0
  z=idot(nodes[i]['iv'],covI[wall])
  if z.l>0:return 1
  if z.h<0:return-1
  assert wall<49,('UNRESOLVED new-wall sign',wall,bits(nodes[i]['mask']),z.enc())
  stats['exact_pairing_fallbacks']+=1;return D.sign(value(i,wall))
 def exact(i):
  row=nodes[i]
  if row['exact']is None:
   a,b=row['parents'];g1=value(a,row['wall']);g2=value(b,row['wall'])
   out=B.vecsub(B.vecscale(exact(a),g2),B.vecscale(exact(b),g1))
   row['exact']=tuple(D.scale(x,K.convert(row['scale']))for x in out)
  return row['exact']
 def classify(i):
  row=nodes[i];vI=row['iv'];q=idot(vI,imv(QI,vI))
  if q.l>0:return'hyperideal_rays'
  if q.h<0:sg=-1
  else:
   stats['exact_norm_fallbacks']+=1
   if row['wall'] is not None and row['wall']>=49:
    assert affine_rank(tuple(bits(row['mask'])))==8;sg=0;stats['exact_affine_gram_cusps']+=1
   else:sg=D.sign(D.pair(exact(i),exact(i),Q))
  if sg>0:return'hyperideal_rays'
  hp=idot(vI,hI)
  if not hp.h<0:
   assert D.sign(D.dot(exact(i),hcov))<0,'opposite-sheet or zero direction'
  return'finite_vertices'if sg<0 else'cusps'
 @lru_cache(maxsize=None)
 def relation(i,j):
  if i==j:return 2
  if i<j:i,j=j,i
  if i>=49:
   if j>=49:return None
   if j in U.S+[U.ANCHOR]:return U.PROFILES[i-49].get(j,2)
   return None
  g=D.pair(N[i],N[j],Q);nn=D.mul(D.pair(N[i],N[i],Q),D.pair(N[j],N[j],Q));rel=B.relation(g,nn)
  if rel['kind']=='dotted':return None
  assert rel['kind']in['right','finite'];return rel['label']
 @lru_cache(maxsize=None)
 def affine_rank(ids):
  def theta(m):return {2:sy.Integer(0),3:sy.Integer(1),4:sy.sqrt(2),5:(1+sy.sqrt(5))/2,6:sy.sqrt(3),10:sy.sqrt(10+2*sy.sqrt(5))/2}[m]
  M=[]
  for i in ids:
   row=[]
   for j in ids:
    if i==j:row.append(sy.Integer(2))
    else:
     m=relation(i,j);assert m is not None,('dotted pair in unresolved norm',ids,i,j);row.append(-theta(m))
   M.append(row)
  rank=0
  for k in range(len(M)):
   pivot=sy.simplify(M[k][k])
   if pivot==0:assert all(sy.simplify(M[k][j])==0 for j in range(k,len(M)));continue
   assert pivot.is_positive is True,(ids,k,pivot);rank+=1
   for i in range(k+1,len(M)):
    for j in range(i,len(M)):
     M[i][j]=sy.simplify(M[i][j]-M[i][k]*M[k][j]/pivot);M[j][i]=M[i][j]
  assert rank==8,('not affine rank8',ids,rank)
  return rank
 current=[]
 for kind,sets in start['complete_census'].items():
  for ids in sets:
   v=reconstruct(35,tuple(ids));i=create(sum(1<<j for j in ids),kind,vec=v)
   assert classify(i)==kind;current.append(i)
 assert len(current)==1168
 stages=[];cluster_coverage=[]
 for wall in range(36,52):
  if wall==49:
   for nw in range(49,52):
    ss=[sign_value(i,nw)for i in current];cluster_coverage.append(dict(wall=nw,positive_active_sets=[bits(nodes[i]['mask'])for i,sg in zip(current,ss)if sg>0],zero_active_sets=[bits(nodes[i]['mask'])for i,sg in zip(current,ss)if sg==0],counts={str(k):v for k,v in sorted(Counter(ss).items())}))
  signs=[sign_value(i,wall)for i in current];neg=[k for k,s in enumerate(signs)if s<0];pos=[k for k,s in enumerate(signs)if s>0];zero=[k for k,s in enumerate(signs)if s==0]
  incidence=[0]*wall
  for k,i in enumerate(current):
   for j in bits(nodes[i]['mask']):incidence[j]|=1<<k
  new={nodes[current[k]]['mask']:current[k]for k in neg}
  for k in zero:
   i=current[k];old=nodes[i];j=len(nodes);row=dict(old);row['mask']|=1<<wall;nodes.append(row)
   assert row['mask']not in new;new[row['mask']]=j
  crossings=[];candidates=0
  for ni in neg:
   a=current[ni]
   for pi in pos:
    b=current[pi];mask=nodes[a]['mask']&nodes[b]['mask']
    if mask.bit_count()<8:continue
    candidates+=1;face=(1<<len(current))-1
    for j in bits(mask):face&=incidence[j]
    if face!=((1<<ni)|(1<<pi)):continue
    newmask=mask|(1<<wall);assert newmask not in new
    i=create(newmask,None,parents=(a,b),wall=wall);nodes[i]['kind']=classify(i);new[newmask]=i
    crossings.append(dict(negative_support=bits(nodes[a]['mask']),positive_support=bits(nodes[b]['mask']),new_support=bits(newmask),kind=nodes[i]['kind']))
  current=sorted(new.values(),key=lambda i:bits(nodes[i]['mask']));counts=Counter(nodes[i]['kind']for i in current)
  stage=dict(wall=wall,candidate_index=wall-36,negative=len(neg),positive=len(pos),zero=len(zero),adjacency_candidates=candidates,crossing_edges=len(crossings),counts=dict(counts),total_rays=len(current),crossings=crossings)
  if wall in[47,48]:
   check=prior if wall==47 else prior49
   for kind in check['complete_census']:assert sorted(bits(nodes[i]['mask'])for i in current if nodes[i]['kind']==kind)==check['complete_census'][kind]
  stages.append(stage);print('clip',wall,{k:v for k,v in stage.items()if k!='crossings'},'fallbacks',dict(stats),flush=True)
 census={kind:sorted(bits(nodes[i]['mask'])for i in current if nodes[i]['kind']==kind)for kind in['finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays']}
 oldvertices={tuple(v)for v in prior49['complete_census']['finite_vertices']};assert oldvertices<={tuple(v)for v in census['finite_vertices']}
 raycert=[]
 for i in current:
  row=nodes[i];q=idot(row['iv'],imv(QI,row['iv']));raycert.append(dict(active=bits(row['mask']),norm_interval=q.enc(),kind=row['kind']))
 # Choose an exact surviving hyperideal ray and perturb toward the strict center.
 escape_node=next(i for i in current if nodes[i]['kind']=='hyperideal_rays' and(nodes[i]['wall']is None or nodes[i]['wall']<49))
 V=exact(escape_node) # Feasibility follows from the complete clipping census; final escape signs are certified below.
 for power in range(13,81):
  eps=K.convert(F(1,1<<power));E=B.vecadd(V,tuple(D.scale(h,eps)for h in H));en=D.pair(E,E,Q)
  if D.sign(en)>0:break
 else:raise AssertionError('No strict spacelike perturbation')
 ep=[D.dot(E,c)for c in cov];assert all(D.sign(x)<0 for x in ep) and D.sign(D.dot(E,hcov))<0
 eiv=[ii(x)for x in E];new_ep=[idot(eiv,c)for c in covI[49:]];assert all(x.h<0 for x in new_ep)
 escape=B.encode(dict(point=E,norm=en,pairings=ep,ray=V,active=bits(nodes[escape_node]['mask']),epsilon=f'1/{1<<power}',new_cap_pairing_intervals=[x.enc()for x in new_ep]))
 paths=[CLUSTER,PRIOR49,f'scripts/{U.TAG}.py',f'scripts/{TAG}.py',f'scripts/{B.TAG}.py','scripts/fresh_nonarithmetic_common_escape_20261016_cap.py',BATCH,PRIOR48,B.CONE,B.PRIMARY,D.START_REVIEW]
 return dict(status='PASS_COMPLETE_EXACT_52_WALL_CONE_AFTER_THREE_CAP_CLUSTER',sources={p:sha(p)for p in paths},initial_counts={k:len(v)for k,v in start['complete_census'].items()},clipping_stages=stages,counts={k:len(v)for k,v in census.items()},complete_census=census,final_ray_norm_intervals=raycert,arithmetic_counts=dict(stats),dag_nodes=len(nodes),strict_escape=escape,old49_cluster_coverage=cluster_coverage,old49_cluster_cut_union_count=len({tuple(ids)for row in cluster_coverage for ids in row['positive_active_sets']}),
  proof='Each exact node is a positive crossing combination divided by a positive rational power of two. Rational outward interval propagation certifies strict signs; uncertain values fall back to the exact formal field. Full common-zero-face adjacency makes the clipping exhaustive. Zero-ray masks are updated explicitly; mask uniqueness is asserted. All178finite vertices of the retained49-wall chamber remain.',
  scope='Complete52-wall cone. The three mutually dotted caps cut the previous common strict escape; the explicit renewed strict spacelike escape still proves infinite volume. No approximate sign is used as an acceptance gate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved complete52 cone replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
