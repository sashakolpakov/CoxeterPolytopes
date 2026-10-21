#!/usr/bin/env python3
"""Rational reconstruction of a strict spacelike witness for the final roof."""
import json,sys
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_higher_face_20261013_final_verify as F
V=F.V;K=F.K

def run():
 _,(G,I,h,g,roots,names,center,pairs)=F.build();scout=V.read(F.F.PREFIX+'_adaptive_next_scout.json');x=scout['remaining_point'];gf=lambda z:float(z[0])+float(z[1])*(1+5**.5)/2;raw=[-sum(gf(G[i][j])*x[j]for j in range(10))for i in range(10)];base=tuple((Q(round(a*10**12),10**12),Q(0))for a in raw)
 for eps in[Q(1,10**n)for n in range(12,3,-1)]:
  y=tuple(K.add(a,K.scale(b,eps))for a,b in zip(base,h,strict=True));q=K.pair(y,I,y);future=K.pair(y,I,h)
  if K.sign(q)<=0 or K.sign(future)>=0:continue
  checks=[V.dot(r,list(map(V.golden,y)))for r in roots]
  if min(z[0]for z in checks)<=0:continue
  out={'status':'EXACT_STRICT_FUTURE_SPACELIKE_WITNESS','point':list(map(K.enc,y)),'norm':K.enc(q),'future_center_pairing':K.enc(future),'base_denominator':10**12,'strict_center_perturbation':str(eps),'proposal_source':F.F.PREFIX+'_adaptive_next_scout.json','proposal_source_sha256':V.sha(F.F.PREFIX+'_adaptive_next_scout.json'),'scope':'The numerical proposal is only a source of rational coordinates; the norm, orientation and all118 inequalities are independently rechecked exactly by the final verifier.'};(V.ROOT/(F.PREFIX+'_escape.json')).write_text(json.dumps(out,indent=2)+'\n');print('escape',out['norm'],'epsilon',eps,flush=True);return
 raise RuntimeError('rational reconstruction failed')
if __name__=='__main__':run()
