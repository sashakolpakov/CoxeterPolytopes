#!/usr/bin/env python3
"""Complete new affine-circuit check after three exact mutually dotted caps."""
import argparse,json,itertools
from pathlib import Path
from hashlib import sha256
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_affine_filtered_20261016_affine_gate'
CLUSTER='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json';OLD='results/fresh_h9_persistent_cusp_screen_20261016.json';AFF='results/fresh_h9_persistent_cusp_affine_20261016.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for p,h in d.get('sources',{}).items():assert sha(p)==h,p
 return d
def run():
 c=read(CLUSTER);o=read(OLD);aff=read(AFF);assert aff['entire_old_chamber']['affine_circuits']==[]
 L=[row+[0]*3 for row in o['pair_matrix']]+[[0]*52 for _ in range(3)]
 for cap in c['caps']:
  j=cap['index'];L[j][j]=2
  for r in cap['relations']:L[j][r['wall']]=L[r['wall']][j]=r.get('label',0)
 assert all(m['kind']=='dotted'for m in c['mutual'])
 rows=[]
 theta={2:s.Integer(0),3:s.Integer(1),4:s.sqrt(2),6:s.sqrt(3)}
 for cap in c['caps']:
  j=cap['index'];ids=[i for i in range(49)if L[j][i]in theta]
  for size in range(len(ids)+1):
   for old in itertools.combinations(ids,size):
    ss=list(old)+[j]
    if any(L[a][b]not in theta for a,b in itertools.combinations(ss,2)):continue
    seen={j}
    while True:
     ext={i for i in ss if any(L[i][a]>2 for a in seen)}
     if ext<=seen:break
     seen|=ext
    if seen!=set(ss):continue
    G=s.Matrix([[2 if a==b else -theta[L[a][b]]for b in ss]for a in ss])
    det=s.simplify(G.det());assert det!=0,('singular connected standard-label subset',ss)
    rows.append(dict(support=ss,determinant=str(det)))
 paths=[f'scripts/{TAG}.py',CLUSTER,OLD,AFF]
 return dict(status='PASS_NO_AFFINE_CIRCUIT_IN_FILTERED52_CHAMBER',sources={p:sha(p)for p in paths},new_connected_subsets=len(rows),records=rows,old_affine_circuits=0,new_affine_circuits=0,
  proof='Every irreducible affine Coxeter diagram has only labels2,3,4,6, or the parallel rank-two case. No parallel pair occurs here. The old49 circuit inventory is empty. A new affine circuit contains at most one new cap since every new-new pair is dotted, and it can use only the listed nine finite neighbors of that cap. Every connected standard-label subset containing the cap is explicitly enumerated and has a nonzero exact determinant. Thus no affine circuit exists in this52-wall chamber.',
  scope='Absence of affine circuits does not prove finite volume or existence of a finite-volume completion. The complete cone and strict escape are separate.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact affine-circuit replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_connected_subsets'])
