#!/usr/bin/env python3
"""Exact reconstruction of the3040 algebraic cap normals used by the SAT scout.
Finite-angle compatibility and finite volume are not inferred from this census.
"""
import argparse,hashlib,json,math
from fractions import Fraction as Q
from pathlib import Path
from scripts.fresh_h9_quartic_completion_20261013_arithmetic import K,Quadratic,dot,matvec,total
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_leaf_caps_20261013_dotted_scout.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_leaf_sat_20261013_exact_normals.json'
Z=K.zero;O=K.one;P=(Q(0),Q(1));a=K.add(P,K.rational(Q(1,2)));b=K.neg(K.div(a,K.add(K.mul(K.rational(12),a),O)));eps=[1]*4+[-1]*4

def scale(x,q):return K.mul(x,K.rational(q))
def numeric(x):return float(x[0])+float(x[1])*(1+5**.5)/2

def gram():
 G=[]
 for i in range(16):
  G.append([K.rational(2) if i==j else (Z if i//8==j//8 else K.neg(K.add(K.add(a,scale(b,eps[i%8]*eps[j%8])),K.rational(2 if i%8==j%8 else 0)))) for j in range(16)])
 return G

def coefficients(x,y):
 delta=[K.sub(u,v) for u,v in zip(x,y)];t=scale(K.add(K.add(x[0],y[0]),K.add(x[4],y[4])),Q(1,2));h=scale(K.sub(K.add(x[0],y[0]),K.add(x[4],y[4])),Q(1,2))
 ds=total(K,delta);de=total(K,[scale(v,e) for v,e in zip(delta,eps)])
 r=[];s=[]
 for v,e in zip(delta,eps):
  w=K.add(scale(v,Q(-1,8)),K.add(K.div(K.mul(a,ds),scale(K.add(O,scale(a,2)),32)),scale(K.div(K.mul(b,de),scale(K.add(O,scale(b,2)),32)),e)))
  c=K.add(K.div(t,scale(a,16)),scale(K.div(h,scale(b,16)),e));r.append(K.add(c,w));s.append(K.sub(c,w))
 return tuple(r+s)

def run():
 source=json.loads(SOURCE.read_text());G=gram();caps=[];seen=set();fields=set();sqrt_branches={-1:0,0:0,1:0}
 for index,xy in enumerate(source['all_profiles']):
  values=sum(xy,[]);fixed=[];linear=[]
  for value in values:
   close=[abs(value-z)<1e-8 for z in [0.,1.,(1+5**.5)/2]]
   if any(close):fixed.append([Z,O,P][close.index(True)]);linear.append(Z)
   else:assert value>=2-1e-8;fixed.append(Z);linear.append(O)
  r0=coefficients(fixed[:8],fixed[8:]);r1=coefficients(linear[:8],linear[8:]);g0=matvec(K,G,r0);g1=matvec(K,G,r1)
  assert g0==tuple(K.neg(x) for x in fixed) and g1==tuple(K.neg(x) for x in linear)
  A=dot(K,r1,g1);B=scale(dot(K,r0,g1),2);C=K.sub(dot(K,r0,g0),K.rational(2))
  if not any(v!=Z for v in linear):
   assert C==Z;D=Z;d0=Z;d1=Z;branch=0
  elif A==Z:
   assert B!=Z;D=Z;d0=K.neg(K.div(C,B));d1=Z;branch=0
  else:
   D=K.sub(K.mul(B,B),scale(K.mul(A,C),4));assert K.sign(D)>=0
   target=next(v for v,l in zip(values,linear) if l!=Z);d0=K.neg(K.div(B,scale(A,2)));den=scale(A,2)
   plus=numeric(d0)+math.sqrt(max(0,numeric(D)))/numeric(den);minus=numeric(d0)-math.sqrt(max(0,numeric(D)))/numeric(den);branch=1 if abs(plus-target)<abs(minus-target) else -1
   assert min(abs(plus-target),abs(minus-target))<1e-6
   d1=K.div(K.rational(branch),den)
  F=Quadratic(K,D);d=(d0,d1);r=tuple((K.add(u,K.mul(v,d0)),K.mul(v,d1)) for u,v in zip(r0,r1));g=tuple((K.add(u,K.mul(v,d0)),K.mul(v,d1)) for u,v in zip(g0,g1))
  assert F.sign(F.sub(dot(F,r,g),F.rational(2)))==0
  assert all(F.sign(x)<=0 for x in g) and F.sign(total(F,g))<0
  assert F.sign(F.sub(d,F.rational(2)))>=0 or not any(v!=Z for v in linear)
  # Selected branch independently matches every saved seed pairing.
  assert max(abs(numeric(v[0])+numeric(v[1])*math.sqrt(max(0,numeric(D)))+float(old)) for v,old in zip(g,values))<1e-6
  encoded={'index':index,'radicand':K.encode(D),'dotted_parameter':F.encode(d),'normal':[[K.encode(v[0]),K.encode(v[1])] for v in r],'covector':[[K.encode(v[0]),K.encode(v[1])] for v in g],'norm_polynomial':[K.encode(x) for x in [A,B,C]],'branch':branch}
  caps.append(encoded);fields.add(tuple(D));sqrt_branches[branch]+=1
  key=tuple(round(float(v),9) for v in values);assert key not in seen;seen.add(key)
  if index%500==0:print({'verified_normals':index},flush=True)
 assert len(caps)==3040
 return {'status':'PASS_EXACT_3040_NORMALS','scope':'Every archived3040 normal has norm2, Coxeter seed angles/dotted pairs, and a common strict timelike interior. Cap-cap compatibility and finite volume require separate checks. Numeric deduplication is not a completeness proof for arbitrary caps.','source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'a':K.encode(a),'b':K.encode(b),'gram':[[K.encode(x) for x in row] for row in G],'center_coefficients':['1']*16,'center_square':K.encode(dot(K,(O,)*16,matvec(K,G,(O,)*16))),'caps':caps,'distinct_literal_radicands':len(fields),'branches':{str(k):v for k,v in sqrt_branches.items()}}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run()
 if args.verify:assert json.loads(OUT.read_text())==out
 else:OUT.write_text(json.dumps(out,indent=2)+'\n')
 print({k:out[k] for k in ['status','distinct_literal_radicands','branches']},flush=True)
