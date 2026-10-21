#!/usr/bin/env python3
"""Two exact affine E7 completions of the85 chamber's H13 face cusps."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_e8_interior_caps_20261014_exact_review as K
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_intrinsic_caps_20261014_cusps'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014.json';Z=K.Z
SQ2=(F(0),F(1))+(F(0),)*6
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 return x
def run():
 d=json.loads((ROOT/BASE).read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q=[[F(x)for x in row]for row in d['form']];old=list(map(K.vector,d['normals']))
 H=K.vector(d['strict_center']);E=K.vector(d['strict_escape']['point'])
 sparse=[[(j,a)for j,a in enumerate(row)if a]for row in Q];cache={}
 def dot(x,y):
  if y not in cache:cache[y]=tuple(tuple(sum(a*y[j][k]for j,a in row)for k in range(8))for row in sparse)
  out=Z
  for a,b in zip(x,cache[y]):
   if any(a)and any(b):out=K.add(out,K.multiply(a,b))
  return out
 def lin(vectors,coeff):
  out=(Z,)*21
  for v,a in zip(vectors,coeff):out=K.vadd(out,K.vscale(v,a))
  return out
 S=d['E7_support'];G=[[K.rational(dot(old[i],old[j]))for j in S]for i in S]
 assert K.pd(G);I=K.inverse(G);positions=[i for i in range(7)if I[i][i]==2]
 assert positions==[3]and S[3]==15
 weights=I[3];assert min(weights)>0
 omega=lin([old[i]for i in S],list(map(K.scalar,weights)))
 assert dot(omega,omega)==K.scalar(2)
 assert [dot(omega,old[i])for i in S]==[K.scalar(int(j==3))for j in range(7)]
 rays=[K.vector(row['ambient_point'])for row in d['face_ideal_vertices']];assert len(rays)==2
 slopes=[]
 for nu,row in zip(rays,d['face_ideal_vertices']):
  assert dot(nu,nu)==Z and dot(nu,omega)==Z and K.sign(dot(nu,H))<0
  g=[dot(nu,n)for n in old];assert all(K.sign(x)<=0 for x in g)
  assert [i for i,x in enumerate(g)if x==Z]==row['ambient_incident']
  slopes.append(g)
 assert K.sign(dot(*rays))<0
 norms0=list(map(F,d['norms']));trials=[];chosen=None
 for k in range(21):
  t=F(2**k);caps=[K.vadd(K.vscale(omega,K.scalar(-1)),K.vscale(nu,K.scalar(t)))for nu in rays]
  failures=[]
  for j,(r,sg)in enumerate(zip(caps,slopes)):
   assert dot(r,r)==K.scalar(2)
   if K.sign(dot(H,r))>=0:failures.append(['center',j])
   for i,n in enumerate(old):
    p=dot(r,n)
    if sg[i]==Z:
     expected=K.scalar(-int(i==S[3]));assert p==expected
    elif not(K.sign(p)<0 and K.sign(K.sub(K.multiply(p,p),K.scalar(2*norms0[i])))>0):failures.append(['pair',j,i])
  p=dot(*caps)
  if not(K.sign(p)<0 and K.sign(K.sub(K.multiply(p,p),K.scalar(4)))>0):failures.append(['mutual'])
  trials.append(dict(shift=t,failed_gates=failures))
  if not failures:chosen=t;new=caps;break
 assert chosen is not None
 N=old+new;norms=norms0+[F(2)]*2;pairrows=[];counts=Counter()
 for i in [85,86]:
  for j in range(i):
   g=dot(N[i],N[j]);assert K.sign(g)<=0
   if g==Z:kind='right'
   elif K.multiply(g,g)==K.scalar(norms[i]*norms[j]/4):kind='pi3'
   else:assert K.sign(K.sub(K.multiply(g,g),K.scalar(norms[i]*norms[j])))>0;kind='dotted'
   counts[kind]+=1;pairrows.append(dict(walls=[i,j],pairing=g,relation=kind))
 assert len(pairrows)==171
 assert K.sign(dot(H,H))<0 and all(K.sign(dot(H,n))<0 for n in N)
 facetpoints=[]
 for i,r in enumerate(N):
  p=K.vadd(H,K.vscale(r,K.scale(dot(H,r),-1/norms[i])))
  assert K.sign(dot(p,p))<0 and K.sign(dot(p,H))<0
  assert all(dot(p,n)==Z if j==i else K.sign(dot(p,n))<0 for j,n in enumerate(N))
  facetpoints.append(p)
 links=[]
 for k,(nu,row)in enumerate(zip(rays,d['face_ideal_vertices'])):
  incident=[i for i,n in enumerate(N)if dot(nu,n)==Z]
  assert incident==row['ambient_incident']+[85+k]
  assert all(K.sign(dot(nu,n))<0 for i,n in enumerate(N)if i not in incident)
  blocks=[]
  for block in row['intrinsic_affine_blocks']:
   ids=block['incident'];walls=[d['face_facet_walls'][i]for i in ids]
   scales=[SQ2 if i in [11,12]else K.scalar(1)for i in ids]
   kernel=[K.scale(s,F(a))for s,a in zip(scales,block['positive_kernel'])]
   assert all(K.sign(x)>0 for x in kernel)
   BG=[[K.rational(K.multiply(K.multiply(scales[i],scales[j]),dot(N[a],N[b])))for j,b in enumerate(walls)]for i,a in enumerate(walls)]
   assert BG==[[F(x)for x in rr]for rr in block['gram']]
   assert K.pd([rr[:-1]for rr in BG[:-1]])
   total=lin([N[i]for i in walls],kernel);mu=F(block['null_multiple'])
   assert mu>0 and total==K.vscale(nu,K.scalar(mu))
   assert all(dot(total,N[i])==Z for i in walls)
   blocks.append(dict(walls=walls,positive_kernel=kernel,null_multiple=mu,quotient_rank=len(walls)-1))
  walls=S+[85+k];kernel=list(map(K.scalar,weights))+[K.scalar(1)]
  BG=[[K.rational(dot(N[i],N[j]))for j in walls]for i in walls]
  assert K.pd([rr[:-1]for rr in BG[:-1]])
  total=lin([N[i]for i in walls],kernel);assert total==K.vscale(nu,K.scalar(chosen))
  assert all(dot(total,N[i])==Z for i in walls)
  blocks.append(dict(walls=walls,positive_kernel=kernel,null_multiple=chosen,quotient_rank=7))
  flat=[i for b in blocks for i in b['walls']]
  assert len(flat)==len(set(flat))and sorted(flat)==incident
  assert sum(b['quotient_rank']for b in blocks)==19
  assert all(dot(N[i],N[j])==Z for a in range(len(blocks))for b in range(a)for i in blocks[a]['walls']for j in blocks[b]['walls'])
  links.append(dict(point=nu,incident=incident,positive_affine_blocks=blocks,quotient_rank=19,bounded_link=True))
 preserved=[]
 for row in d['face_finite_vertices']:
  p=K.vector(row['ambient_point']);values=[dot(p,r)for r in new]
  assert all(K.sign(x)<0 for x in values)
  preserved.append(dict(incident=row['ambient_incident'],new_pairings=values))
 assert len(preserved)==44
 ep=[dot(E,n)for n in N]
 assert all(K.sign(x)<0 for x in ep)and K.sign(dot(E,E))>0 and K.sign(dot(E,H))<0
 # In the E7 face each added normal projects to t*nu, a future null
 # redundant inequality; the complete H13 face itself is unchanged.
 assert all(K.vadd(r,omega)==K.vscale(nu,K.scalar(chosen))for r,nu in zip(new,rays))
 paths=[BASE,f'scripts/{TAG}.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_87_WALL_BRANCH_TWO_E7_CUSP_LINKS_COMPLETED_ESCAPE_REMAINS',sources={p:sha(p)for p in paths},shift_search=trials,selected_shift=chosen,E7_support=S,weight_node=S[3],weight=omega,form=Q,normals=N,norms=norms,new_pair_counts=dict(counts),new_pairs=pairrows,strict_center=H,facet_points=facetpoints,completed_cusp_links=links,preserved_face_vertices=preserved,strict_escape=dict(point=E,norm=dot(E,E),all87_pairings=ep),scope='Retains all85 and adds precisely2 norm2 caps -omega15+t*nu. Both actual ambient ideal links have positive affine components covering their full incident sets and quotientrank19;44 face finitevertices survive and the15facet finitevolumeH13 face is unchanged. All87facets actual;171newpairs exact. Same strictspacelike escape remains, so ambient volume infinite. No complete ambient cusp/vertex census and no general completion claim.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'shift',d['selected_shift'],d['new_pair_counts'])
