#!/usr/bin/env python3
"""Seven exterior polars of an actual three-dimensional face, against all53 walls."""
import argparse,itertools,json
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import fresh_nonarithmetic_ridge_pair_20261016_polar as P
ROOT=P.ROOT;TAG='fresh_nonarithmetic_h3_closure_20261016_scout'
CAP='results/fresh_nonarithmetic_escape52_20261016_cap.json';CONE='results/fresh_nonarithmetic_escape52_20261016_cone.json';S=[0,1,2,3,5,6]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 Q,N,H,E=P.geometry();cap=P.P.B.read(CAP);cone=P.P.B.read(CONE);N.append(mp.matrix([P.mid(x)for x in cap['normal_interval']]));cov=[Q*n for n in N];hcov=Q*H
 face={kind:[s for s in rows if set(S)<=set(s)]for kind,rows in cone['complete_census'].items()};assert len(face['hyperideal_rays'])==7
 rows=[];V=[]
 for i,ids in enumerate(face['hyperideal_rays']):
  v=mp.matrix(P.P.null([list(cov[j])for j in ids]));v/=-(v.T*hcov)[0];q=(v.T*Q*v)[0];assert q>0;v*=mp.sqrt(2/q);V.append(v)
  rel=[dict(wall=j,pairing=mp.nstr((v.T*c)[0],45),**P.relation((v.T*c)[0]))for j,c in enumerate(cov)];bad=[r for r in rel if r['kind']not in['finite','dotted']]
  rows.append(dict(index=i,active=ids,point=[mp.nstr(x,65)for x in v],relations=rel,bad=bad));print(i,ids,'bad',bad,flush=True)
 mutual=[]
 for i,j in itertools.combinations(range(7),2):
  g=(V[i].T*Q*V[j])[0];mutual.append(dict(pair=[i,j],pairing=mp.nstr(g,45),**P.relation(g)))
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',CAP,CONE,P.CL,P.CO,P.P.BATCH]
 return dict(status='NUMERICAL_SEVEN_POLAR_ACTUAL_H3_FACE_SCREEN',sources={p:sha(p)for p in paths},white_support=S,face_census=face,proposals=rows,mutual=mutual,scope='All seven hyperideal directions of the actual H3 face are reconstructed numerically, and every candidate is checked against all53 original normals. No automatic Coxeter property is assumed for projected old walls. Exactification and complete face clipping are required before compactness claims.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved seven-polar screen')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],sum(not r['bad']for r in d['proposals']))
