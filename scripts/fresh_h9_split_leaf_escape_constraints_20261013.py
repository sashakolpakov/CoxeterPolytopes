#!/usr/bin/env python3
"""Exact necessary escape clauses for compatible finite-volume completion.

Every listed rational-coefficient vector is spacelike and strictly inside the
sixteen seed halfspaces. A finite-volume completion must choose some cap with
positive pairing. Integer interval signs give a conservative cutter superset;
uncertain signs are included, never discarded. This does not itself prove
that the clauses are jointly unsatisfiable.
"""
import argparse,hashlib,json,math
from fractions import Fraction as F
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
NORMALS=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
SCOUT=ROOT/'data/enumerations/fresh_h9_split_leaf_z3_20261013_exact_source_scout.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_escape_constraints_20261013.json'
S=10**24;DEN=10**12
sqrt5=math.isqrt(5*S*S);PHI=((S+sqrt5)//2,(S+sqrt5+2)//2+1)

def floor(x):return x.numerator//x.denominator
def ceil(x):return -((-x.numerator)//x.denominator)
def gsign(a,b):
 r=2*a+b
 if not b:return (r>0)-(r<0)
 if not r:return (b>0)-(b<0)
 if (r>0)==(b>0):return 1 if r>0 else -1
 delta=r*r-5*b*b
 return ((r>0)-(r<0))*((delta>0)-(delta<0))
def gint(v):
 a=F(v['a']);b=F(v['b']);ends=[a*S+b*t for t in PHI];return floor(min(ends)),ceil(max(ends))
def mul(x,y):
 vals=[a*b for a in x for b in y];return min(vals)//S,-((-max(vals))//S)
def add(x,y):return x[0]+y[0],x[1]+y[1]
def root_int(v):
 assert gsign(F(v['a']),F(v['b']))>=0,'radicand must be exactly nonnegative'
 lo,hi=gint(v);assert hi>=0
 return math.isqrt(max(0,lo)*S),math.isqrt(hi*S)+1

def setup():
 normal=json.loads(NORMALS.read_text());source=ROOT/normal['source'];old=json.loads(source.read_text());assert hashlib.sha256(source.read_bytes()).hexdigest()==normal['source_sha256']
 A=[[0]*16 for _ in range(16)];B=[[0]*16 for _ in range(16)]
 for i in range(16):
  for j in range(16):
   if i==j:A[i][j]=22
   elif i//8!=j//8:
    same=(i%8)//4==(j%8)//4;A[i][j]=-(3 if same else 8)-(22 if i%8==j%8 else 0);B[i][j]=-(12 if same else 10)
   assert F(normal['gram'][i][j]['a'])==F(A[i][j],11) and F(normal['gram'][i][j]['b'])==F(B[i][j],11)
 pats=[];dints=[]
 def pair(v):return F(v['a']),F(v['b'])
 for cap in normal['caps']:
  masks=[[],[],[]];d0,d1=cap['dotted_parameter'];negative_d=tuple(tuple(-x for x in pair(v)) for v in (d0,d1))
  for j,v in enumerate(cap['covector']):
   v0,v1=map(pair,v)
   if v1==(0,0) and v0==(0,0):continue
   if v1==(0,0) and v0==(-1,0):masks[0].append(j)
   elif v1==(0,0) and v0==(0,-1):masks[1].append(j)
   else:
    assert (v0,v1)==negative_d
    masks[2].append(j)
  di=add(gint(d0),mul(gint(d1),root_int(cap['radicand'])))
  pats.append(masks);dints.append(di)
 return A,B,pats,dints

def certify(w,A,B,pats,dints):
 assert isinstance(w,(list,tuple)) and len(w)==16 and all(type(x) is int for x in w),'escape vector must contain exactly sixteen integer coefficients'
 aw=[sum(a*x for a,x in zip(row,w)) for row in A];bw=[sum(a*x for a,x in zip(row,w)) for row in B];qa=sum(x*y for x,y in zip(w,aw));qb=sum(x*y for x,y in zip(w,bw));ca=sum(aw);cb=sum(bw)
 assert gsign(qa,qb)>0 and gsign(ca,cb)<0 and all(gsign(a,b)<0 for a,b in zip(aw,bw))
 cutters=0;ambiguous=0
 for i,(masks,di) in enumerate(zip(pats,dints)):
  u,v,t=[sum(w[j] for j in inds) for inds in masks]
  # cap pairing is -(u+v*phi+t*d)/DEN. Its upper bound is
  # minus the lower bound of the parenthesis; include if possibly positive.
  ph=min(v*PHI[0],v*PHI[1]);dot=min(t*di[0],t*di[1]);lo=u*S+ph+dot
  hi=u*S+max(v*PHI[0],v*PHI[1])+max(t*di[0],t*di[1])
  if lo<0:cutters|=1<<i
  if lo<=0<=hi:ambiguous+=1
 return {'integer_coefficients':w,'denominator':DEN,'norm_numerator':[str(qa),str(qb)],'norm_denominator':11*DEN*DEN,'center_pairing_numerator':[str(ca),str(cb)],'center_pairing_denominator':11*DEN,'possible_cutters_hex':hex(cutters),'interval_ambiguous_cap_pairs':ambiguous}

def input_regressions(valid,A,B,pats,dints):
 message='escape vector must contain exactly sixteen integer coefficients'
 malformed={'too_short':valid[:-1],'too_long':valid+[0],
            'boolean':[False]+valid[1:],
            'integral_fraction':[F(valid[0])]+valid[1:],
            'floating_coefficient':[float(valid[0])]+valid[1:]}
 for name,w in malformed.items():
  try:certify(w,A,B,pats,dints)
  except AssertionError as error:assert str(error)==message,(name,str(error))
  else:raise AssertionError('malformed vector accepted: '+name)
 # This negative number rounds to an interval with upper endpoint zero.
 # The old upper-bound-only test would incorrectly accept its square root.
 negative={'a':str(F(-1,S*S)),'b':'0'}
 assert gint(negative)[1]>=0 and gsign(F(negative['a']),F(negative['b']))<0
 try:root_int(negative)
 except AssertionError as error:assert str(error)=='radicand must be exactly nonnegative'
 else:raise AssertionError('negative radicand accepted')
 for value,root in((0,0),(1,S)):
  lo,hi=root_int({'a':str(value),'b':'0'});assert lo<=root<=hi
 return {'malformed_vectors_rejected':list(malformed),'negative_radicand_with_nonnegative_interval_upper_rejected':True,'zero_and_unit_radicands_retained':True}

def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');p.add_argument('--source',type=Path,default=OUT);args=p.parse_args()
 if not args.source.is_absolute():args.source=ROOT/args.source
 if not args.verify and args.source!=OUT:p.error('--source is only supported with --verify')
 A,B,pats,dints=setup()
 if args.verify:
  old=json.loads(args.source.read_text());assert old['normal_source']==str(NORMALS.relative_to(ROOT));assert old['normal_sha256']==hashlib.sha256(NORMALS.read_bytes()).hexdigest();records=[]
  for record in old['constraints']:
   check=certify(record['integer_coefficients'],A,B,pats,dints);assert check==record;records.append(check)
  try:certify([0]*16,A,B,pats,dints)
  except AssertionError:pass
  else:raise AssertionError('zero vector accepted')
  assert records
  controls=input_regressions(records[0]['integer_coefficients'],A,B,pats,dints)
  print({'status':'PASS_EXACT_NECESSARY_ESCAPE_CLAUSES','constraints':len(records),'forged_zero_rejected':True,**controls});return
 import numpy as np
 from scripts.fresh_h9_split_leaf_caps_20261013 import seed
 source=json.loads(SCOUT.read_text());phi=(1+5**.5)/2;a=phi+.5;b=-a/(12*a+1);G=seed(a,b);ev,Q=np.linalg.eigh(G);ii=list(np.where(ev< -1e-7)[0])+list(np.where(ev>1e-7)[0]);N=Q[:,ii]*np.sqrt(abs(ev[ii]));N[:,0]*=np.sign(sum(N[:,0]));I=np.linalg.pinv(N.T);records=[];seen=set()
 for old in source['witnesses']:
  c=I@np.r_[1,old['klein_point']];w=[int(round(float(x)*DEN)) for x in c];record=certify(w,A,B,pats,dints)
  if record['possible_cutters_hex'] in seen:continue
  seen.add(record['possible_cutters_hex']);records.append(record)
 out={'status':'PASS_EXACT_NECESSARY_ESCAPE_CLAUSES','scope':'Necessary Boolean cutter clauses for finite-volume completion from this archived bank; no satisfiability/nonexistence conclusion by itself. Each noncutter has certified nonpositive pairing; uncertain interval signs remain allowed cutters.','normal_source':str(NORMALS.relative_to(ROOT)),'normal_sha256':hashlib.sha256(NORMALS.read_bytes()).hexdigest(),'scout_source':str(SCOUT.relative_to(ROOT)),'scout_sha256_at_extraction':hashlib.sha256(SCOUT.read_bytes()).hexdigest(),'interval_scale':S,'phi_interval':list(PHI),'constraints':records}
 OUT.write_text(json.dumps(out,indent=2)+'\n');print({'status':out['status'],'constraints':len(records),'ambiguous_cap_pairings':sum(r['interval_ambiguous_cap_pairs'] for r in records)},flush=True)
if __name__=='__main__':main()
