#!/usr/bin/env python3
"""Complete exact32-to36 cone clipping with active-face adjacency.

Field coefficients are sparse arrays over the frozen octic field. No division
in any new radical extension and no enumeration of all9-wall subsets is used.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
import fresh_nonarithmetic_escape35_20261016 as A
R,C,K,Z,O=A.R,A.C,A.K,A.Z,A.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_escape35_20261016_cone'
PRIMARY='results/fresh_nonarithmetic_escape35_20261016.json'
START=A.BASE32;START_REVIEW='results/fresh_h9_joint_pair_20261015_two_review.json'
D=[A.B.P.D,A.B.D,A.D]
ZERO=(Z,)*8;ONE=(O,)+(Z,)*7
FACT=[O]
for mask in range(1,8):
 q=O
 for b in range(3):
  if mask>>b&1:q*=D[b]
 FACT.append(q)
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def lift(a):return(a,)+(Z,)*7
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def sub(a,b):return tuple(x-y for x,y in zip(a,b))
def scale(a,q):return tuple(x*q for x in a)
def mul(a,b):
 if a==ZERO or b==ZERO:return ZERO
 aa=[(i,x)for i,x in enumerate(a)if x];bb=[(i,x)for i,x in enumerate(b)if x];r=[Z]*8
 for i,x in aa:
  for j,y in bb:r[i^j]+=x*y*FACT[i&j]
 return tuple(r)
@lru_cache(maxsize=100000)
def biv(a):
 z=R.L.I()
 for x in a.to_list():z=z*R.THETA+R.L.I(F(str(x)))
 return z
RIV=[biv(x).sqrt()for x in D];BIV=[]
for mask in range(8):
 z=R.L.I(1)
 for b in range(3):
  if mask>>b&1:z=z*RIV[b]
 BIV.append(z)
def recursive_sign(a,level=3):
 if not level:return R.sg(a[0])
 n=1<<(level-1);lo=a[:n];hi=a[n:2*n]
 sa=recursive_sign(lo,level-1);sb=recursive_sign(hi,level-1)
 if not sa:return sb
 if not sb or sa==sb:return sa
 al=tuple(lo)+(Z,)*(8-n);ah=tuple(hi)+(Z,)*(8-n)
 diff=sub(mul(al,al),scale(mul(ah,ah),D[level-1]))
 return sa*recursive_sign(diff[:n],level-1)
def sign(a):
 if a==ZERO:return 0
 z=R.L.I()
 for x,b in zip(a,BIV):
  if x:z=z+biv(x)*b
 if z.l>0:return 1
 if z.h<0:return-1
 return recursive_sign(a)
def flatten_raw(x):
 # Each supplied scalar is the complete three-level formal tower.
 return tuple(A.B.db(x[b2][b1][b0])for b2 in range(2)for b1 in range(2)for b0 in range(2))
def dot(a,b):
 out=ZERO
 for x,y in zip(a,b):out=add(out,mul(x,y))
 return out
def mv(Q,v):
 return[sumfield(scale(v[j],q)for j,q in enumerate(row)if q)for row in Q]
def sumfield(xs):
 out=ZERO
 for x in xs:out=add(out,x)
 return out
def pair(v,w,Q):return dot(v,mv(Q,w))
def encfield(a):return [R.enc(x)for x in a]
def bit_ids(mask):return[i for i in range(mask.bit_length())if mask>>i&1]
def classify(v,Q,hcov):
 q=sign(pair(v,v,Q))
 if q>0:return'hyperideal_rays'
 if sign(dot(v,hcov))>=0:return'opposite_sheet_rays'
 return'finite_vertices'if q<0 else'cusps'
def run():
 d=A.read(PRIMARY);start=A.read(START);review=A.read(START_REVIEW)
 assert review['status']=='PASS_INDEPENDENT_TWO_CAPS_AND_COMPLETE_32_WALL_CONE'
 assert review['sources'][START]==sha(START) and review['complete_census']==start['complete_census']
 Q=[[A.B.db(x)for x in row]for row in d['form']];H=list(map(A.B.db,d['center']));hcov=tuple(lift(x)for x in C.mv(Q,H))
 N=[tuple(flatten_raw(x)for x in row)for row in d['normals']];assert len(N)==36
 cov=[tuple(mv(Q,n))for n in N]
 oldN=R.m4(start['normals']);oldcov=[C.mv(Q,n)for n in oldN]
 assert [[x[0]for x in n]for n in N[:32]]==oldN
 assert all(all(not any(x[1:])for x in n)for n in N[:32])
 rays=[];oldcounts=Counter()
 for kind,sets in start['complete_census'].items():
  for ids in sets:
   v=C.null_one([oldcov[i]for i in ids]);vals=[C.dot(v,c)for c in oldcov]
   g=next(x for x in vals if x)
   if R.sg(g)>0:v=[-x for x in v];vals=[-x for x in vals];g=-g
   assert all(R.sg(x)<=0 for x in vals) and ids==[i for i,x in enumerate(vals)if not x]
   v=[x/-g for x in v];fv=tuple(lift(x)for x in v)
   assert classify(fv,Q,hcov)==kind
   rays.append(dict(v=fv,mask=sum(1<<i for i in ids),kind=kind));oldcounts[kind]+=1
   if len(rays)%100==0:print('reconstructed32 rays',len(rays),flush=True)
 assert len(rays)==972 and len({x['mask']for x in rays})==972
 stages=[]
 for wall in range(32,36):
  values=[dot(r['v'],cov[wall])for r in rays];signs=[sign(x)for x in values]
  neg=[i for i,s in enumerate(signs)if s<0];pos=[i for i,s in enumerate(signs)if s>0];zero=[i for i,s in enumerate(signs)if s==0]
  new={rays[i]['mask']:rays[i]for i in neg}
  for i in zero:
   row=dict(rays[i]);row['mask']|=1<<wall;assert row['mask']not in new;new[row['mask']]=row
  incidence=[0]*wall
  for i,r in enumerate(rays):
   for j in bit_ids(r['mask']):incidence[j]|=1<<i
  crossing=[];tested=0
  for i in neg:
   for j in pos:
    common=rays[i]['mask']&rays[j]['mask']
    if common.bit_count()<8:continue
    tested+=1;face=(1<<len(rays))-1
    for k in bit_ids(common):face&=incidence[k]
    if face!=((1<<i)|(1<<j)):continue
    v=tuple(sub(mul(values[j],a),mul(values[i],b))for a,b in zip(rays[i]['v'],rays[j]['v']))
    assert any(x!=ZERO for x in v) and sign(dot(v,cov[wall]))==0
    mask=common|(1<<wall);assert mask not in new
    kind=classify(v,Q,hcov);new[mask]=dict(v=v,mask=mask,kind=kind)
    crossing.append(dict(negative_support=bit_ids(rays[i]['mask']),positive_support=bit_ids(rays[j]['mask']),new_support=bit_ids(mask),kind=kind))
  rays=sorted(new.values(),key=lambda r:bit_ids(r['mask']))
  counts=Counter(r['kind']for r in rays)
  assert counts['finite_vertices']=={32:92,33:95,34:98,35:100}[wall]
  stage=dict(wall=wall,negative=len(neg),positive=len(pos),zero=len(zero),adjacency_candidates=tested,crossing_edges=len(crossing),counts=dict(counts),total_rays=len(rays),crossings=crossing)
  stages.append(stage);print('clip',wall,{k:v for k,v in stage.items()if k!='crossings'},flush=True)
 census={kind:sorted(bit_ids(r['mask'])for r in rays if r['kind']==kind)for kind in('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
 assert census['finite_vertices']==sorted(d['finite_vertex_incidence'])
 counts={k:len(v)for k,v in census.items()}
 digest=sha256()
 for r in rays:digest.update(json.dumps([bit_ids(r['mask']),[encfield(x)for x in r['v']]],separators=(',',':')).encode())
 paths=[f'scripts/{TAG}.py',PRIMARY,f'scripts/{A.TAG}.py',START,START_REVIEW,'scripts/fresh_h9_joint_corner_20261015_bridge_review.py']
 return dict(status='PASS_COMPLETE_EXACT_36_WALL_CONE_BY_FOUR_CLIPS',sources={p:sha(p)for p in paths},inherited_complete_cone=START_REVIEW,initial_counts=dict(oldcounts),clipping_stages=stages,counts=counts,complete_census=census,coordinate_digest=digest.hexdigest(),
   arithmetic='Eight sparse coefficients over the frozen octic field for three independently formal positive radicals; exact interval signs with recursive norm fallback; no extension division.',
   adjacency_proof='For each crossing pair the complete common-zero face contains exactly its two rays, hence is a two-dimensional face of the pointed cone. New rays are positive combinations, so all old signs and the exact common active set follow. The new-wall pairing is checked exactly. The inherited cone is full-dimensional and pointed; strict H preserves dimension after each clip.',
   scope='Complete36-wall polyhedral cone, inherited from the independently complete32-wall cone and four exhaustive clipping steps. Finite-corner list agrees with the exact local update. The separately certified strict spacelike escape proves infinite volume. No finite-volume or construction novelty claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved complete36 cone replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
