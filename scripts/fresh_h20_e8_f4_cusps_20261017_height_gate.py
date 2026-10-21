#!/usr/bin/env python3
"""All-real height obstruction for the fixed affineE8/endpoint-affineA7 completion.

Both orientations of the one-dimensional missing Euclidean direction are
included. Every coherent acute/Coxeter wall in these profiles misses T.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_e8_f4_cusps_20261017_height_gate';INPUT=P.CUSP
WITNESSES={(0,-1,0):10,(0,-1,-1):3,(0,1,0):3,(0,1,-1):10,(1,-1,0):2,(1,-1,-1):9,(1,1,0):9,(1,1,-1):2}
def expr(x):return str(s.simplify(x))
def run():
 D=json.loads((ROOT/INPUT).read_text());G=s.Matrix([[s.sympify(x)for x in row]for row in D['core_gram']]);ids=D['core_ids']
 T=list(range(1,11))+[12,13,14,15,61,87,151];ix=[ids.index(i)for i in T];J=G.extract(ix,ix).inv();rows=[];e8=[]
 for ci,c in enumerate(D['cusps']):
  nu=s.Matrix(list(map(s.sympify,c['null_coefficients'])));z=s.Matrix(list(map(s.sympify,c['normalized_free_coefficients'])));np=G*nu;nt=np.extract(ix,[0])
  def polynomial(b):
   bp=G*b;bt=bp.extract(ix,[0]);return[s.simplify(2-(bt.T*J*bt)[0]),s.simplify(-2*(nt.T*J*bt)[0]),s.simplify(-(nt.T*J*nt)[0])]
  b=s.Matrix(list(map(s.sympify,c['base_new_root_coefficients'][0])));coef=polynomial(b);wall=10 if ci==0 else 9
  assert(G*b)[ids.index(wall)]==0 and np[ids.index(wall)]==-1
  shifted=[s.simplify(sum(coef)),s.simplify(coef[1]+2*coef[2]),coef[2]]
  bounds=[P.alg(x)for x in shifted];assert all(x[1]<0 for x in bounds)
  e8.append(dict(cusp=ci,forced_nonzero_product_wall=wall,product='-h',coherence_forces='h>0 from theta+r=h*nu and a common strict center',Coxeter_forces='h>=1 since two norm2 nonorthogonal Coxeter roots have negative product of magnitude at least1',projected_norm_coefficients=list(map(expr,coef)),shifted_at_h1=list(map(expr,shifted)),shifted_bounds=list(map(P.save,bounds))))
  AB=s.eye(21)[:,[ids.index(i)for i in c['A6_ids']]];AJ=(AB.T*G*AB).inv()
  for orientation in(-1,1):
   for end,sign in[(0,orientation),(-1,-orientation)]:
    b=-AB*AJ[:,end]+sign*z;coef=polynomial(b);wall=WITNESSES[ci,orientation,end];j=ids.index(wall)
    slope=s.simplify(np[j]);constant=s.simplify((G*b)[j]);assert P.alg(slope)[1]<0
    threshold=s.simplify(-constant/slope);assert P.alg(threshold)[0]>0
    shifted=[s.simplify(coef[0]+coef[1]*threshold+coef[2]*threshold**2),s.simplify(coef[1]+2*coef[2]*threshold),coef[2]]
    bounds=[P.alg(x)for x in shifted];assert all(x[1]<0 for x in bounds)
    rows.append(dict(cusp=ci,orientation=orientation,endpoint=c['A6_ids'][end],acute_witness_wall=wall,old_pair_constant=expr(constant),old_pair_slope=expr(slope),necessary_height_lower=expr(threshold),projected_norm_coefficients=list(map(expr,coef)),shifted_at_necessary_lower=list(map(expr,shifted)),strict_negative_shifted_bounds=list(map(P.save,bounds))))
 assert len(rows)==8 and len(e8)==2
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',INPUT,'results/fresh_h20_e8_exact_cluster_20261017.json']
 return dict(status='PASS_ALL_REAL_HEIGHTS_FIXED_E8_A7_COMPLETION_MISS_T',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},T=T,E8_height_cases=e8,A7_height_cases=rows,proof='For r=b+h*nu the projection toT-perp has square q(h)=2-(T.r)^T J_T(T.r). Each saved shifted polynomial has all3coefficients strictly negative. Acuteness or coherence plus the minimum nonright Coxeter product puts h at or above that shift. Hence q(h)<0 throughout the entire admissible real halfline, so the hyperplane misses the hyperbolicT-face.',completeness='The eight A7 records cover both signs of the one-dimensional positive direction and both endpoint contacts at both cusps. Any change of that direction by a multiple ofnu is absorbed into the arbitrary real height. This covers the specified norm2 affineA7 endpoint completion; it does not enumerate other affine link types or altered old walls.',scope='Exact all-real parameter obstruction within the fixed affineE8 plus endpoint-affineA7 cusp-completion profiles. Lowering heights, making heights unequal, or reversing the transverse direction cannot add a facet to the targetH3. Other cusp completion diagrams and other endpoint profiles remain open.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
