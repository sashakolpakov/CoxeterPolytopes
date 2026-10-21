#!/usr/bin/env python3
"""Changed H9 core: compact H8 face plus three strictly disjoint endpoint caps."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_cap_20261014'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def pd(M):
 W=[r[:]for r in M]
 for i in range(len(W)):
  d=W[i][i]
  if d.sign()<=0:return False
  for j in range(i+1,len(W)):
   for k in range(i+1,len(W)):W[j][k]-=W[j][i]*W[i][k]/d
 return True
def run():
 src=json.loads((ROOT/SOURCE).read_text())
 def sf(x):return A.E(A.Q(F(x['a'])+F(x['b'])/2,F(x['b'])/2))
 Q8=[[sf(x)for x in row]for row in src['lattice_gram']]
 U=[[sf(x)for x in row]for row in src['simple_roots']]
 H8=[sf(x)for x in src['geometry']['interior_point']]
 G8=[[A.pair(u,v,Q8)for v in U]for u in U]
 assert G8==[[sf(x)for x in row]for row in src['wall_gram']]
 assert len(U)==11 and len(Q8)==9 and A.inertia(Q8)==[8,1]
 assert A.pair(H8,H8,Q8).sign()<0 and all(A.pair(H8,u,Q8).sign()<0 for u in U)
 # Re-enumerate all165 eight-wall supports, with a bounded projective slice.
 dual=[sum((u[i]for u in U),A.Z)for i in range(9)]
 assert A.pair(dual,dual,Q8).sign()<0
 rows=[A.mv(Q8,u)for u in U];corners={};ranks=0
 for ids in combinations(range(11),8):
  v=A.kernel_one([rows[i]for i in ids])
  if v is None:ranks+=1;continue
  if A.pair(v,dual,Q8).sign()>0:v=[-x for x in v]
  if any(A.dot(v,row).sign()>0 for row in rows):continue
  assert A.pair(v,v,Q8).sign()<0 and A.pair(v,H8,Q8).sign()<0
  active=[i for i,row in enumerate(rows)if not A.dot(v,row)];assert len(active)==8
  assert pd([[G8[i][j]for j in active]for i in active])
  corners[tuple(active)]=v
 assert len(corners)==41 and set(corners)=={tuple(x)for x in src['geometry']['vertices']}
 V=[corners[tuple(ids)]for ids in src['geometry']['vertices']]
 a=A.E(A.Q(F(-1,4),F(1,4)));d=A.E(A.Q(0,F(1,4)));phi=A.PHI
 exceptional=[]
 for vi in[0,32,33]:
  v=V[vi];hh=[-A.pair(v,u,Q8)for u in U]
  pos=[h for h in hh if h.sign()>0];m=pos[0]
  for h in pos:
   if(h-m).sign()<0:m=h
  v=[x/m for x in v];hh=[h/m for h in hh]
  assert A.pair(v,v,Q8)==-a and any(h==1 for h in hh)
  assert all(not h or(h-1).sign()>=0 for h in hh)
  exceptional.append(dict(source_vertex=vi,point=v,contacts=hh))
 Q=[row+[A.Z]for row in Q8]+[[A.Z]*9+[d]]
 R=[A.Z]*9+[A.ONE];N=[u+[A.Z]for u in U]+[R]+[row['point']+[A.E(-1)]for row in exceptional]
 assert len(N)==15 and A.inertia(Q)==[9,1]
 norms=[A.pair(n,n,Q)for n in N];assert norms[:11]==[A.E(2)]*11 and norms[11]==d and norms[12:]==[A.E(F(1,4))]*3
 G=[[A.pair(n,m,Q)for m in N]for n in N];pairs=[];counts=Counter()
 for i,j in combinations(range(15),2):
  typ,m=A.relation(G[i][j],norms[i],norms[j]);assert typ in['right','finite','dotted']
  key='pi'+str(m)if typ=='finite'else typ;counts[key]+=1
  pairs.append(dict(pair=[i,j],type=typ,label=m,product=G[i][j]))
 assert dict(counts)==dict(right=79,pi3=8,pi5=2,dotted=16)
 delta=F(1)
 while True:
  H=H8+[A.E(-delta)]
  if A.pair(H,H,Q).sign()<0 and all(A.pair(H,n,Q).sign()<0 for n in N):break
  delta/=2
  assert delta>F(1,1<<30)
 facets=[]
 for i,n in enumerate(N):
  h=A.pair(H,n,Q)/norms[i];p=[x-h*y for x,y in zip(H,n)]
  assert A.pair(p,p,Q).sign()<0 and A.pair(p,H,Q).sign()<0
  vals=[A.pair(p,m,Q)for m in N]
  assert not vals[i]and all(v.sign()<0 for j,v in enumerate(vals)if j!=i)
  facets.append(dict(facet=i,point=p,norm=A.pair(p,p,Q),pairings=vals))
 vertices=[]
 for vi,v in enumerate(V):
  p=v+[A.Z];ids=src['geometry']['vertices'][vi]+[11]
  vertices.append(dict(source_vertex=vi,point=p,incident=ids))
 for k,row in enumerate(exceptional):
  v=row['point'];p=v+[-a/d];ids=src['geometry']['vertices'][row['source_vertex']]+[12+k]
  assert A.pair(p,p,Q)==-a/(4*d)
  vertices.append(dict(source_vertex=row['source_vertex'],point=p,incident=ids))
 for row in vertices:
  p,ids=row['point'],row['incident'];vals=[A.pair(p,n,Q)for n in N]
  assert [i for i,v in enumerate(vals)if not v]==sorted(ids)
  assert all(v.sign()<=0 for v in vals)
  assert A.pair(p,p,Q).sign()<0 and A.pair(p,H,Q).sign()<0
  assert pd([[G[i][j]for j in ids]for i in ids])
  row.update(norm=A.pair(p,p,Q),pairings=vals)
 assert len({A.normkey(r['point'])for r in vertices})==44
 # Complete finite-neighbor lists for each new cap give exactly its single
 # finite vertex. R has the complete41source vertices; source-only rank9
 # supports cannot be positive definite in the source Lorentz9-space.
 neighbors={}
 for i in range(12,15):
  near=[j for j in range(15)if j!=i and A.relation(G[i][j],norms[i],norms[j])[0]in['right','finite']]
  assert near==src['geometry']['vertices'][exceptional[i-12]['source_vertex']]
  neighbors[str(i)]=near
 finite_edges=[]
 for i,j in combinations(range(44),2):
  shared=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(shared)==8:
   assert pd([[G[k][l]for l in shared]for k in shared])
   finite_edges.append(dict(vertices=[i,j],support=shared))
 assert len(finite_edges)==167
 assert [r['vertices']for r in finite_edges if max(r['vertices'])>=41]==[[0,41],[32,42],[33,43]]
 # The added caps' small balls leave the middle of the compact source open.
 vsum=[sum((r['point'][i]for r in exceptional),A.Z)for i in range(9)]
 eps=F(1)
 while True:
  E=[x+A.E(eps)*h for x,h in zip(vsum,H8)]+[A.E(-10)]
  eq=A.pair(E,E,Q);ev=[A.pair(E,n,Q)for n in N]
  if eq.sign()>0 and all(v.sign()<0 for v in ev)and A.pair(E,H,Q).sign()<0:break
  eps/=2;assert eps>F(1,1<<30)
 # Arithmetic triage: this d is not square-equivalent to either earlier
 # single-coordinate lift coefficient. A field square has rational-square
 # norm, whereas these determinant ratios have norms5/64 and5/16.
 ratios=[d/(2*phi),d/(phi-1)]
 rn=[x.a.a*x.a.a-5*x.a.b*x.a.b for x in ratios]
 assert rn==[F(5,64),F(5,16)]
 cycle=4*G[11][12]*G[11][12]/(norms[11]*norms[12])
 assert cycle==A.E(A.Q(0,4))
 conj_cycle=A.E(A.Q(cycle.a.a,-cycle.a.b));assert conj_cycle.sign()<0
 sigmaQ8=[[A.E(A.Q(x.a.a,-x.a.b))for x in row]for row in Q8]
 assert A.inertia(sigmaQ8)==[9,0]
 assert A.Q(d.a.a,-d.a.b).sign()<0
 paths=[SOURCE,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',f'scripts/{TAG}.py']
 return A.enc(dict(status='PASS_CHANGED_H8_CORE_THREE_DISJOINT_CAPS_INFINITE_VOLUME',sources={p:sha(p)for p in paths},
  field='Q(sqrt5), serialized in the inherited quartic format with zero extension coefficients',
  source_support_census=dict(trials=165,rank_deficient=ranks,finite_vertices=41),
  new_coefficient=d,form=Q,normals=N,norms=norms,gram=G,pairs=pairs,pair_counts=dict(counts),
  exceptional_caps=exceptional,center=H,center_delta=str(delta),facet_points=facets,
  complete_finite_vertices=vertices,complete_finite_vertex_count=44,
  new_cap_finite_neighbors=neighbors,complete_finite_edges=finite_edges,complete_finite_edge_count=167,
  open_finite_corner_germ_count=44*9-2*167,
  strict_escape=dict(point=E,norm=eq,pairings=ev,epsilon=str(eps),center_pairing=A.pair(E,H,Q)),
  host_triage=dict(old_coefficients=[2*phi,phi-1],determinant_ratios=ratios,rational_norms=[str(x)for x in rn],
   new_host_conjugate_signature=[9,1],two_cycle=cycle,conjugate_two_cycle=conj_cycle),
  scope='Changed core; none of the mixed12 walls are retained. Exact15actualfacets, complete44finitevertices and167finiteedges, including3new vertical finiteedges;62remaining finite-corner germs are not claimed to form bounded ideal vertices. Strict escape proves infinite volume. Negative conjugate two-cycle obstructs quasi-arithmeticity of any finite-covolume completion retaining R and one new cap; current chamber is not a lattice. No compact H9construction or whole infinite-vertex census is asserted.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['complete_finite_vertex_count'],'finite vertices',flush=True)
