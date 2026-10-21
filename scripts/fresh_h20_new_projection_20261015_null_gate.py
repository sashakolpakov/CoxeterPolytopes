#!/usr/bin/env python3
"""All-profile escape gate on every literal source parallel-pair null ray."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h20_b2_h18_coupling_20261014_review as A
K=A.K;dot=A.dot;mv=A.mv
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_new_projection_20261015_null_gate'
SOURCE='results/fresh_h20_new_projection_20261015_exact.json';HELPER='scripts/fresh_h20_b2_h18_coupling_20261014_review.py'
def decode(x):
 assert all(F(a)==0 for a in x[2:]);return K(F(x[0]),F(x[1]))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[decode(x)for x in row[:19]]for row in d['form'][:19]]
 N=[[decode(x)for x in row[:19]]for row in d['normals'][:40]]
 E=[decode(x)for x in d['strict_escape'][:19]];H=[decode(x)for x in d['strict_center'][:19]]
 Etrans=[decode(x)for x in d['strict_escape'][19:]]
 assert 2*dot(Etrans,Etrans)==1
 assert dot(E,mv(Q,E)).sign()<0 and dot(E,mv(Q,H)).sign()<0
 cov=[mv(Q,n)for n in N];G=[[dot(n,c)for c in cov]for n in N];q=[G[i][i]for i in range(40)]
 assert all(x==1 or x==2 for x in q)
 root=lambda x:K(1)if x==1 else K(0,1)
 rows=[];rays=[]
 for i,j in combinations(range(40),2):
  if G[i][j]*G[i][j]!=q[i]*q[j]or G[i][j].sign()>=0:continue
  c=[a/root(q[i])+b/root(q[j])for a,b in zip(N[i],N[j])]
  assert not dot(c,mv(Q,c))and dot(c,mv(Q,H)).sign()<0
  g=[dot(c,v)for v in cov];assert all(x.sign()<=0 for x in g)
  bounds=[]
  for k in range(40):
   if g[k]:
    first=K(1)if q[k]==2 else K(0,F(1,2))
    bounds.append((k,first/-g[k]))
  assert bounds
  anchor,L=bounds[0]
  for k,b in bounds[1:]:
   if(b-L).sign()>0:anchor,L=k,b
  assert L.sign()>0 and all((L-b).sign()>=0 for _,b in bounds)
  e=dot(c,mv(Q,E));margin=L*L*e*e-2
  assert e.sign()<0 and margin.sign()>0
  rows.append(dict(source_parallel_pair=[i,j],null_ray=c,source_pairings=g,incident_source_walls=[k for k,x in enumerate(g)if not x],mandatory_scale_anchor=anchor,scale_lower_bound=L,source_escape_pairing=e,squared_escape_margin=margin))
  rays.append(c)
 assert len(rows)==31
 products=[dot(c,mv(Q,e))for c,e in combinations(rays,2)]
 assert len(products)==465 and all(x.sign()<0 for x in products)
 minimum=0
 for k,r in enumerate(rows):
  if(r['squared_escape_margin']-rows[minimum]['squared_escape_margin']).sign()<0:minimum=k
 paths=[f'scripts/{TAG}.py',SOURCE,HELPER]
 return A.enc(dict(status='PASS_ALL_TRANSVERSE_PROFILES_ON31_SOURCE_PARALLEL_NULL_RAYS_RETAIN_ESCAPE',sources={p:sha(p)for p in paths},
  source_wall_count=40,source_parallel_pairs=len(rows),distinct_future_null_rays=len(rows),negative_distinct_ray_products=465,
  source_escape=E,transverse_escape=Etrans,transverse_escape_norm=1,ray_certificates=rows,minimum_margin_index=minimum,minimum_margin=rows[minimum]['squared_escape_margin'],
  theorem='Normalize a proposed wall to r²=2 and write r=t*c+w with w in the positive transverse plane. For positive t the old source wall j has nonzero product t(c,N_j); the first possible nonright Coxeter cosine is1/2, so t>=sqrt(q_j/2)/(-c,N_j). Thus t>=L. Null c implies w²=2, and Cauchy gives (w,E_trans)<=sqrt2. The stored exact signs L²(c,E_source)²>2 and (c,E_source)<0 imply (r,E)<0 for every real transverse profile, with no finite-label cutoff. At t0, acuteness to the two support roots forces w into their negative orthant, making any strict common interior impossible. Negative t violates an old source inequality. Therefore every coherent added wall with source projection on one of these31rays retains E.',
  scope='The31rays exhaust the literal parallel pairs of the retained40-wall source, not all its ideal vertices or all null directions a changed source might admit. Any finite coherent extension of the44-wall chamber using only these projection rays retains its strict spacelike escape and is infinite volume. The three explicit cusp repairs lie in this family. No obstruction to different source directions, mixed residuals, or deleting retained walls is asserted.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['minimum_margin_index'],d['minimum_margin'])
