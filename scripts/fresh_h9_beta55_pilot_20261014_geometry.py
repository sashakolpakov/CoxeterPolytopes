#!/usr/bin/env python3
"""Complete literal 12-wall cone census; no reflector extension search."""
import argparse,json
from itertools import combinations
from hashlib import sha256
from pathlib import Path
import fresh_h9_beta11_vinberg_h300_20261014_review as A
K=A.K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_beta55_pilot_20261014_geometry'
SOURCE='results/fresh_h9_beta55_pilot_20261014.json'
HELPER='scripts/fresh_h9_beta11_vinberg_h300_20261014_review.py'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def null(M):
 a=[r[:]for r in M];piv=[];i=0
 for k in range(len(a[0])):
  ids=[j for j in range(i,len(a))if a[j][k]]
  if not ids:continue
  j=ids[0];a[i],a[j]=a[j],a[i];q=a[i][k];a[i]=[x/q for x in a[i]]
  for j in range(len(a)):
   if j!=i:
    q=a[j][k];a[j]=[x-q*y for x,y in zip(a[j],a[i])]
  piv.append(k);i+=1
  if i==len(a):break
 if i!=9:return None
 free=next(k for k in range(10)if k not in piv);v=[K(k==free)for k in range(10)]
 for i,k in enumerate(piv):v[k]=-a[i][free]
 return v

def run():
 d=json.loads((ROOT/SOURCE).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 assert d['status']=='PASS_COMPLETE_BOUNDED_BETA55_PREFIX'
 G=[[A.dec(x)for x in row]for row in d['form']];N=[list(map(A.dec,r))for r in d['roots']]
 C=list(map(A.dec,d['strict_center']));H=list(map(A.dec,d['basepoint']))
 assert len(N)==12
 assert all(not N[i][j]for i in (8,9)for j in range(1,9))
 assert N[8][0]*N[9][9]-N[8][9]*N[9][0]
 NG=[A.mv(G,n)for n in N]
 rays={};supports=0;rankfail=0;infeasible=0
 for ids in combinations(range(12),9):
  supports+=1;v=null([NG[i]for i in ids])
  if v is None:rankfail+=1;continue
  pairs=[A.dot(v,p)for p in NG]
  if all(x.sign()>=0 for x in pairs):v=[-x for x in v];pairs=[-x for x in pairs]
  if not all(x.sign()<=0 for x in pairs):infeasible+=1;continue
  assert any(x.sign()<0 for x in pairs)
  first=next(x for x in v if x);scale=first if first.sign()>0 else -first;v=[x/scale for x in v]
  key=tuple((x.a,x.b)for x in v)
  if key in rays:continue
  pairs=[A.dot(v,p)for p in NG];incident=[i for i,x in enumerate(pairs)if not x]
  norm=A.pair(v,v,G);orientation=A.pair(v,C,G)
  assert orientation.sign()<0
  kind='finite'if norm.sign()<0 else'ideal'if not norm else'hyperideal'
  if kind=='finite':
   assert len(incident)==9
   M=[[A.pair(N[i],N[j],G)for j in incident]for i in incident]
   for k in range(9):
    q=M[k][k];assert q.sign()>0
    for i in range(k+1,9):
     for j in range(k+1,9):M[i][j]-=M[i][k]*M[k][j]/q
  rays[key]=dict(incident=incident,point=v,norm=norm,future_pairing=orientation,kind=kind)
 rows=list(rays.values());finite=[r for r in rows if r['kind']=='finite']
 edges=[]
 for i,j in combinations(range(len(finite)),2):
  shared=sorted(set(finite[i]['incident'])&set(finite[j]['incident']))
  if len(shared)==8:edges.append(dict(vertices=[i,j],support=shared))
 escapes=[r for r in rows if r['kind']=='hyperideal'];assert escapes
 ray=escapes[0]['point'];eps=K(1)
 while True:
  E=[x+eps*c for x,c in zip(ray,C)]
  if A.pair(E,E,G).sign()>0:break
  eps/=2
 assert A.pair(E,C,G).sign()<0 and all(A.pair(E,n,G).sign()<0 for n in N)
 counts=dict(tested_supports=supports,dependent_supports=rankfail,infeasible_supports=infeasible,
  finite_vertices=len(finite),ideal_vertices=sum(r['kind']=='ideal'for r in rows),
  hyperideal_rays=len(escapes),finite_edges=len(edges),open_finite_germs=9*len(finite)-2*len(edges))
 assert supports==220
 return A.enc(dict(status='PASS_COMPLETE_BETA55_T20_GEOMETRY_INFINITE_VOLUME',sources={p:sha(p)for p in [f'scripts/{TAG}.py',SOURCE,HELPER]},
  counts=counts,finite_vertices=finite,finite_edges=edges,all_extreme_rays=rows,
  strict_escape=dict(point=E,norm=A.pair(E,E,G),center_pairing=A.pair(E,C,G),
   wall_pairings=[A.pair(E,n,G)for n in N],ray_plus_center_epsilon=eps),
  scope='Complete polyhedral cone support census on twelve actual walls only. The strict spacelike escape proves this prefix infinite volume, not nonreflectivity of the lattice.'))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS beta55 geometry saved replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'],d['counts'])
