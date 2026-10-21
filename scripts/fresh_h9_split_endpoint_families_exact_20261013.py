#!/usr/bin/env python3
"""Exact reconstruction of the small-label endpoint-cap discovery subset."""
import argparse,hashlib,itertools,json,math
from fractions import Fraction as F
from pathlib import Path
from scripts.fresh_h9_split_leaf_sat_20261013_exact import K,Quadratic,coefficients,gram,scale,numeric,dot,matvec,total,Z,O,P
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'data/enumerations/fresh_h9_split_endpoint_families_20261013_scout.json'
OUT=ROOT/'data/enumerations/fresh_h9_split_endpoint_families_20261013_exact.json'

def recognize(t):
 for val,k in [(0,Z),(1,O),((1+5**.5)/2,P)]:
  if abs(t-val)<1e-7:return k
 return None

def build():
 data=json.loads(SOURCE.read_text());G=gram();caps=[]
 for index,row in enumerate(data['caps']):
  mags=[-v for v in row['seed_covector']]
  if any(recognize(v)is None for v in mags if v<2-1e-6):continue
  if any(abs(v-2)<1e-6 for v in mags):continue
  if any(sum(abs(v-1)<1e-6 for v in mags[s:s+8])>=4 for s in [0,8]):continue
  side=row['seed_side'];x=[recognize(v) for v in mags[side*8:side*8+8]];assert all(v is not None for v in x)
  y=mags[(1-side)*8:(1-side)*8+8];ts=[]
  for group in range(2):
   known=[K.add(x[i],recognize(y[i])) for i in range(group*4,group*4+4)if recognize(y[i])is not None]
   assert not known or all(v==known[0]for v in known)
   ts.append(known[0]if known else None)
  assert any(t is not None for t in ts)
  fixed=[K.sub(ts[i//4],x[i])if ts[i//4]is not None else K.neg(x[i])for i in range(8)]
  linear=[Z if ts[i//4]is not None else O for i in range(8)]
  if side==0:r0=coefficients(x,fixed);r1=coefficients([Z]*8,linear)
  else:r0=coefficients(fixed,x);r1=coefficients(linear,[Z]*8)
  g0=matvec(K,G,r0);g1=matvec(K,G,r1)
  A=dot(K,r1,g1);B=scale(dot(K,r0,g1),2);C=K.sub(dot(K,r0,g0),scale(O,2))
  if all(t is not None for t in ts):assert C==Z;D=d0=d1=Z
  elif A==Z:D=Z;d0=K.neg(K.div(C,B));d1=Z
  else:
   D=K.sub(K.mul(B,B),scale(K.mul(A,C),4));assert K.sign(D)>=0
   d0=K.neg(K.div(B,scale(A,2)));den=scale(A,2);unknown=ts.index(None);target=y[unknown*4]+numeric(x[unknown*4])
   choices=[K.div(K.rational(s),den)for s in [-1,1]];d1=min(choices,key=lambda z:abs(numeric(d0)+numeric(z)*math.sqrt(max(0,numeric(D)))-target))
  E=Quadratic(K,D);normal=tuple((K.add(v,K.mul(u,d0)),K.mul(u,d1))for v,u in zip(r0,r1));cov=tuple((K.add(v,K.mul(u,d0)),K.mul(u,d1))for v,u in zip(g0,g1))
  assert E.sign(E.sub(dot(E,normal,cov),E.rational(2)))==0
  assert E.sign(total(E,cov))<0
  finite=[]
  for v in cov:
   match=next((label for label,k in [(2,Z),(3,K.neg(O)),(5,K.neg(P))]if E.sign(E.sub(v,(k,Z)))==0),None)
   if match is None:assert E.sign(E.add(v,E.rational(2)))<0
   finite.append(match)
  assert all(sum(m==3 for m in finite[s:s+8])<4 for s in [0,8])
  q=K.sub(scale(O,2),scale(total(K,[K.mul(v,v)for v in x]),F(1,2)));assert K.sign(q)>0
  witness=[8,8,8,1,-3,-3,-3,-3,-3,-3,-3,5,5,5,5,5]
  wp=total(E,[E.mul(E.rational(v),z)for v,z in zip(witness,cov)])
  assert max(abs(numeric(v[0])+numeric(v[1])*math.sqrt(max(0,numeric(D)))+m)for v,m in zip(cov,mags))<1e-6
  caps.append({'index':len(caps),'scout_index':index,'seed_side':side,'finite_profile':row['finite_profile'],'radicand':K.encode(D),'normal':[[K.encode(v[0]),K.encode(v[1])]for v in normal],'covector':[[K.encode(v[0]),K.encode(v[1])]for v in cov],'finite_seed_labels':finite,'endpoint_projection_norm':K.encode(q),'old_universal_escape_pairing':E.encode(wp),'old_universal_escape_sign':E.sign(wp)})
 assert len(caps)==728
 return {'status':'PASS_EXACT_ENDPOINT_CAPS_WITHOUT_SINGLE_CAP_AFFINE_SUBDIAGRAMS','scope':'728 archived entries with finite seed labels2,3,5, strictly dotted other seed relations, norm2 and a common strict center. Each extends a seed rank8 orthogonal block elliptically. Cap-cap compatibility and compactness require separate checks. Not a completeness claim for all endpoint caps.','source':str(SOURCE.relative_to(ROOT)),'source_sha256':hashlib.sha256(SOURCE.read_bytes()).hexdigest(),'gram':[[K.encode(v)for v in row]for row in G],'caps':caps,'positive_old_escape_cutters':sum(r['old_universal_escape_sign']>0 for r in caps)}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=build()
 if a.verify:assert r==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(r,indent=2)+'\n')
 print({'status':r['status'],'entries':len(r['caps']),'old_escape_cutters':r['positive_old_escape_cutters']},flush=True)
