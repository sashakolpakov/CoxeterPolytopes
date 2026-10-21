#!/usr/bin/env python3
"""Compressed new-corner frontier and exact new parallel-cusp diagnosis."""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from itertools import combinations
import fresh_h20_e7_residual_gate_20261014_exact_data as X
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_finite_edgewalk_20261014_frontier';PRIMARY='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(tuple,list)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
 return x

def rref(A):
 M=[list(row)for row in A];piv=[];k=0
 for j in range(len(M[0])):
  found=next((i for i in range(k,len(M))if M[i][j]),None)
  if found is None:continue
  M[k],M[found]=M[found],M[k];q=M[k][j];M[k]=[x/q for x in M[k]]
  for i in range(len(M)):
   if i!=k:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[k])]
  piv.append(j);k+=1
  if k==len(M):break
 return M,piv

def kernel(A):
 M,piv=rref(A);n=len(A[0]);out=[]
 for j in range(n):
  if j in piv:continue
  v=[F(0)]*n;v[j]=F(1)
  for i,k in enumerate(piv):v[k]=-M[i][j]
  out.append(v)
 return out

def run():
 d=json.loads((ROOT/PRIMARY).read_text());data=X.load();K=X.K;dot=data['dot']
 for p,h in d['sources'].items():assert sha(p)==h,p
 def simple(v):
  assert all(not any(map(F,a[1]))and not any(map(F,a[0][8:]))for a in v)
  return tuple(tuple(map(F,a[0][:8]))for a in v)
 phys={i:simple(v)for i,v in enumerate(d['normals'])if i!=87}
 delta=tuple(K.add(a,b)for a,b in zip(phys[9],phys[88]));assert dot(delta,delta)==K.Z and K.sign(dot(delta,data['H']))<0
 pairs={i:dot(delta,n)for i,n in phys.items()};assert all(K.sign(x)<=0 for x in pairs.values())
 D=tuple(map(F,d['field']['discriminant']));Z=B.Z;bp=[tuple(tuple(map(F,a))for a in c)for c in d['normals'][87]]
 bg=[Z,Z]
 for i,row in enumerate(data['Q']):
  for j,q in enumerate(row):
   if q:
    a=tuple(delta[i])+(F(0),)*8
    for k in range(2):bg[k]=B.add(bg[k],B.scale(B.mul(a,bp[j][k]),q))
 def sg(x):
  a,b=x;sa,sb=B.sign(a),B.sign(b)
  if not sb:return sa
  if not sa or sa==sb:return sb
  return sa*B.sign(B.sub(B.mul(a,a),B.mul(D,B.mul(b,b))))
 assert sg(bg)<0
 active=sorted(i for i,g in pairs.items()if g==K.Z);assert len(active)==20
 sq2=(F(0),F(1))+(F(0),)*6
 def scaled(i):return tuple(K.multiply(a,sq2)for a in phys[i])if i in[17,79]else phys[i]
 rows=[scaled(i)for i in active];G=[[K.rational(dot(a,b))for b in rows]for a in rows]
 rank=len(rref(G)[1]);assert rank==18
 rem=set(range(20));components=[]
 while rem:
  c={min(rem)};todo=list(c)
  for i in todo:
   for j in sorted(rem-c):
    if G[i][j]:c.add(j);todo.append(j)
  ids=sorted(c);rem-=c;M=[[G[i][j]for j in ids]for i in ids];ker=kernel(M)
  components.append(dict(walls=[active[i]for i in ids],gram=M,rank=len(rref(M)[1]),kernel=ker))
 # Exhibit positive-norm tangent lineality, not merely a rank count.
 basis=[phys[i]for i in data['S']]+[scaled(i)for i in data['parent']['face_facet_walls'][:13]]+[K.vector(data['parent']['face_basis'][-1])]
 assert len(basis)==21
 QB=[[K.rational(dot(a,b))for b in basis]for a in basis]
 M=[[K.rational(dot(a,b))for b in basis]for a in rows];ker=kernel(M);assert len(ker)==2
 def norm(v):return sum(v[i]*QB[i][j]*v[j]for i in range(21)for j in range(21))
 co=next(v for v in ker if norm(v)>0)
 def total(xs):
  a=K.Z
  for x in xs:a=K.add(a,x)
  return a
 z=tuple(total(K.scale(v[k],c)for c,v in zip(co,basis))for k in range(21))
 assert K.rational(dot(z,z))==norm(co)>0 and dot(z,delta)==K.Z and all(dot(z,n)==K.Z for n in rows)
 # Exact active-set bookkeeping: every finite endpoint of a germ that
 # retains a new cap lies in the complete40 new-cap vertex census.
 vertices=d['preserved_old_vertices']+d['new_cap_finite_vertices'];sets=[set(v['incident'])for v in vertices]
 selected={tuple(sorted(v['incident'])):v for v in d['selected_edge_endpoints']}
 for v in d['new_cap_finite_vertices']:
  if len(set(v['incident'])&{88,89})==1:assert tuple(sorted(v['incident']))in selected
 directions=[]
 for i in range(44,84):
  for j in sorted(sets[i]):
   support=sets[i]-{j};ends=[k for k,s in enumerate(sets)if k!=i and support<=s];assert len(ends)<=1
   if not support&{88,89}:assert ends and ends[0]<44
   directions.append(dict(vertex=i,released=j,support=sorted(support),other_finite_vertex=ends[0]if ends else None))
 finite=sum(v['other_finite_vertex']is not None for v in directions);assert len(directions)==800 and finite==259
 paths=[PRIMARY,X.BASE,X.PARENT,f'scripts/{TAG}.py','scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py']
 return enc(dict(status='PASS_NEW40_CORNER_FRONTIER_AND_EXACT_DEFICIENT_PARALLEL_CUSP',sources={p:sha(p)for p in paths},new_corner_directions=directions,counts=dict(directions=800,second_finite_endpoint=259,no_second_finite_endpoint=541),new_parallel_cusp=dict(walls=[9,88],point=delta,all89_basefield_pairings=pairs,bridge_pairing=bg,incident=active,scaled_incident_gram=G,quotient_normal_rank=rank,components=components,lineality_vector=z,lineality_norm=dot(z,z),lineality_basis_coefficients=co),scope='Complete finite-endpoint status only for800 germs at the40 new vertices.541 have no second finite vertex; ideal endpoint/link classification of those541 is not asserted. The explicit new parallel cusp has all90 signs and a positive-norm tangent lineality vector annihilating every incident wall, so its full horospherical link is unbounded.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'])
