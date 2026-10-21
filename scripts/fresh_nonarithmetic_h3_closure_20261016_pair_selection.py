#!/usr/bin/env python3
"""Bounded compatibility search joining the three remaining H3 regions."""
import argparse,itertools,json,math
from collections import Counter
from hashlib import sha256
import mpmath as mp
import numpy as np
import fresh_nonarithmetic_h3_closure_20261016_scout as H3
P=H3.P;ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_pair_selection'
A_INPUT='results/fresh_nonarithmetic_h3_closure_20261016_nonpolar_scout.json';C_INPUT='results/fresh_nonarithmetic_h3_closure_20261016_foreign_scout.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 Q,N,H,E=P.geometry();cap=P.P.B.read(H3.CAP);sc=P.P.B.read(f'results/{H3.TAG}.json');N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));N.extend(mp.matrix([mp.mpf(x)for x in sc['proposals'][i]['point']])for i in[1,3,6]);ad=P.P.B.read(A_INPUT);cd=P.P.B.read(C_INPUT)
 W=mp.matrix([[N[i][j]for i in H3.S]for j in range(10)]);J=(W.T*Q*W)**-1;basis=[7,8,30,33];res=[N[i]-W*J*W.T*Q*N[i]for i in basis];B=mp.matrix([[(v.T*Q*w)[0]for w in res]for v in res]);BI=B**-1;BI=np.array(BI.tolist(),dtype=float)
 aa=[(i,r)for i,r in enumerate(ad['candidates'])if r['target']==0 and all(m in[2,3,4,5,6,10]for m in r['labels'])];bb=[(i,r)for i,r in enumerate(ad['candidates'])if r['target']==2];cc=list(enumerate(cd['candidates']))
 va=np.array([[r['pairings'][i]for i in basis]for _,r in aa]);vb=np.array([[r['pairings'][i]for i in basis]for _,r in bb]);vc=np.array([[r['pairings'][i]for i in basis]for _,r in cc]);ab=va@BI@vb.T;ac=va@BI@vc.T;bc=vb@BI@vc.T
 def fit(g):
  if g<-2-1e-5:return dict(kind='dotted')
  if abs(g)<1e-5:return dict(kind='right',label=2)
  if -2+1e-5<g<-1e-5:
   m=math.pi/math.acos(-g/2)
   if abs(m-round(m))<1e-4:return dict(kind='finite',label=round(m))
  return None
 pairs=[];types=Counter();tri=[]
 for i,j in itertools.product(range(len(bb)),range(len(cc))):
  rel=fit(bc[i,j])
  if not rel:continue
  types[rel['kind']]+=1;row=dict(B=bb[i][0],C=cc[j][0],pairing=float(format(bc[i,j],'.12g')),**rel);pairs.append(row)
  ka=[k for k in range(len(aa))if fit(ab[k,i])and fit(ac[k,j])]
  if ka:
   k=min(ka,key=lambda k:sum(aa[k][1]['labels']));tri.append(dict(A=aa[k][0],B=row['B'],C=row['C'],BC_relation=rel,BC_pairing=row['pairing'],label_sum=sum(aa[k][1]['labels'])+sum(bb[i][1]['labels'])+sum(cc[j][1]['labels'])))
 tri.sort(key=lambda r:(r['BC_relation']['kind']=='dotted',r['label_sum'],r['A'],r['B'],r['C']))
 paths=[f'scripts/{TAG}.py',A_INPUT,C_INPUT,f'results/{H3.TAG}.json',H3.CAP,P.CL,P.CO,P.P.BATCH]
 return dict(status='NUMERICAL_THREE_REGION_H3_CAP_COMPATIBILITY_SELECTION',sources={p:sha(p)for p in paths},counts=dict(A=len(aa),B=len(bb),C=len(cc),BC_pairs=len(pairs),**types),BC_pairs=pairs,compatible_triples=tri,scope='Numerical pairing fits only, with a conservative1e-4 effective-label tolerance. At most one low-complexity A representative is recorded for each B/C pair. Exact rebuilding, full affine inventory and complete three-face clipping are required for every selected triple.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved H3 pair selection')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],'triples',len(d['compatible_triples']));print(d['compatible_triples'][:10])
