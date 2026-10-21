#!/usr/bin/env python3
"""Exact bounded transplant of an I_12,1 A12 corner into the E8 face.

Both A12 chain orientations are tested. This is not a search over arbitrary
H12 chambers or a claim that either singly compatible cap closes the face.
"""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from hashlib import sha256
from collections import Counter
import fresh_h20_contact_closure_20261014_replacement as b
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_contact_closure_20261014_h12_transfer'
FACE='results/fresh_h20_contact_closure_20261014_e8_face.json'
BASE='results/fresh_h20_contact_closure_20261014_replacement.json'
Z=b.ZZ
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 f=json.loads((ROOT/FACE).read_text());d=json.loads((ROOT/BASE).read_text())
 for data in [f,d]:
  for p,h in data['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']]
 N=[[tuple(map(F,x))for x in row]for row in d['normals']]
 H=list(map(lambda x:b.c(F(x)),d['strict_center']))
 T=[tuple(map(F,x))for x in f['normalized_timelike_basis']]
 A=f['orthogonal_face_walls'];E=f['E8_walls'];C=[[F(x)for x in row]for row in f['A12_gram']];CI=b.rat.inv(C)
 cache={}
 def dot(x,y):
  key=tuple(y)
  if key not in cache:cache[key]=[tuple(sum(Q[i][j]*y[j][k]for j in range(21))for k in range(4))for i in range(21)]
  out=Z
  for a,v in zip(x,cache[key]):out=b.add(out,b.mul(a,v))
  return out
 assert dot(T,T)==b.c(F(-1,13))
 norms=list(map(F,d['norms']));rows=[];allroots=[]
 for orientation,indices in [('forward',(3,11)),('reverse',(10,2))]:
  roots=[];entries=[]
  for j,k,q in [(indices[0],2,F(2)),(indices[1],3,F(1))]:
   r=[b.scale(x,k)for x in T]
   for i,a in enumerate(A):r=[b.add(x,b.scale(y,-CI[j-1][i]))for x,y in zip(r,N[a])]
   assert dot(r,r)==b.c(q)
   assert all(dot(r,N[a])==b.c(-int(i==j-1))for i,a in enumerate(A))
   assert all(dot(r,N[e])==Z for e in E)
   assert b.sign(dot(r,H))<0
   counts=Counter();pairs=[];bad=[]
   for i,w in enumerate(N):
    g=dot(r,w);assert b.sign(g)<=0
    cosine_squared=b.scale(b.mul(g,g),1/(q*norms[i]))
    if g==Z:kind='right'
    elif b.sign(b.add(cosine_squared,b.c(-1)))>0:kind='dotted'
    elif cosine_squared==b.c(1):kind='parallel'
    elif cosine_squared==b.c(F(1,4)):kind='pi3'
    elif cosine_squared==b.c(F(1,2)):kind='pi4'
    elif cosine_squared==b.c(F(3,4)):kind='pi6'
    else:
     assert b.sign(cosine_squared)>0 and b.sign(b.add(cosine_squared,b.c(-1)))<0
     assert cosine_squared[1]==cosine_squared[3]==0 and cosine_squared[2]!=0
     # A nonrational squared Coxeter cosine in Q(sqrt79) is impossible:
     # phi(m)<=4 leaves only m=5,8,10,12, with fields Qsqrt5,Qsqrt2,Qsqrt3.
     kind='forbidden_all_labels_Qsqrt79';bad.append(i)
    counts[kind]+=1;pairs.append(dict(wall=i,pairing=g,cosine_squared=cosine_squared,relation=kind,original60=i<60))
   assert sum(counts.values())==76
   entries.append(dict(fundamental_weight_index=j,timelike_coefficient=k,norm=q,normal=r,pair_counts=dict(counts),forbidden_walls=bad,pairs=pairs,strict_center_pairing=dot(r,H)))
   roots.append(r);allroots.append(r)
  assert dot(*roots)==Z
  assert entries[0]['forbidden_walls']==([0,18,20,31]if orientation=='forward'else[0])
  assert not entries[1]['forbidden_walls'] and entries[1]['pair_counts']==dict(right=19,pi4=1,dotted=56)
  assert all(i<60 for e in entries for i in e['forbidden_walls'])
  rows.append(dict(orientation=orientation,caps=entries,mutual_pairing=dot(*roots),whole_pair_compatible=False))
 # Small all-label separation witnesses, independent of the field-degree argument.
 forward=rows[0]['caps'][0]['pairs'][0]['cosine_squared']
 reverse=rows[1]['caps'][0]['pairs'][0]['cosine_squared']
 assert b.sign(forward)>0 and b.sign(b.add(forward,b.c(F(-1,4))))<0
 assert b.sign(b.add(reverse,b.c(F(-3,4))))>0 and b.sign(b.add(reverse,b.c(F(-4,5))))<0
 paths=[FACE,BASE,f'scripts/{TAG}.py','scripts/fresh_h20_contact_closure_20261014_e8_face.py','scripts/fresh_h20_contact_closure_20261014_replacement.py','scripts/fresh_h20_coupled_lift_20261014.py']
 return b.enc(dict(status='PASS_BOTH_A12_SOURCE_TRANSPLANTS_FAIL_RETAINED_OLD60_ANGLE_GATES',sources={p:sha(p)for p in paths},field_basis=d['field_basis'],A12_chain=A,E8_support=E,orientations=rows,standalone_norm1_caps_mutual_pairing=dot(allroots[1],allroots[3]),all_label_argument='A nonrational squared cosine in Q(sqrt79) has degree2. For m>=3 its degree is phi(m)/2; the degree2 labels5,8,10,12 give quadratic fields5,2,5,3, never79. Labels2,3,4,6 are rational and directly checked.',scope='Exactly two chain orientations of the two fixed I12,1 source caps. Each norm1 cap passes all76 pair gates individually, but each norm2 cap fails retained original60 walls. No arbitrary H12-source or changed-core exclusion; no completed H20 construction.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
