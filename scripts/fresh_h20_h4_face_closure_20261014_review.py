#!/usr/bin/env python3
"""Independent SymPy checks for the uniform H4-face obstruction.

Reads literal integer source and the corner-walk input. No producer or
repository arithmetic implementation is imported. Geometric/arithmetic-group
implications are reviewed separately in the matching note.
"""
from pathlib import Path
from itertools import combinations_with_replacement
from collections import Counter
import argparse,hashlib,json
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_h4_face_closure_20261014'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
OUT=ROOT/f'results/{TAG}_review.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def sign(x):
 x=s.simplify(s.expand(x));a=s.simplify(x.expand().coeff(s.sqrt(5),0));b=s.simplify(x.expand().coeff(s.sqrt(5),1))
 assert s.simplify(x-a-b*s.sqrt(5))==0 and a.is_Rational and b.is_Rational
 if b==0:return s.sign(a)
 if a==0:return s.sign(b)
 if s.sign(a)==s.sign(b):return s.sign(a)
 return s.sign(a)*s.sign(a*a-5*b*b)
def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=-s.Matrix(raw['basis_gram']).applyfunc(s.Rational)
 W=s.Matrix(raw['wall_coordinates_as_columns_transposed']).applyfunc(s.Rational);C=W*Q;G=W*Q*W.T
 assert W.shape==(210,22) and all(x.q==1 for x in W)and Q.det()!=0
 factors=s.factor_list(Q.charpoly().as_poly())[1]
 negative=sum(m*p.count_roots(-s.oo,0)for p,m in factors)
 positive=sum(m*p.count_roots(0,s.oo)for p,m in factors)
 assert (positive,negative)==(21,1)
 norms=[G[i,i]for i in range(210)];assert Counter(norms)=={s.Integer(2):42,s.Integer(4):168}
 assert all(G[i,j]<=0 for i in range(210)for j in range(i))
 assert all(all((2*C[i,j]/norms[i]).q==1 for j in range(22))for i in range(210))
 assert all(G[i,j]%2==0 for i in range(210)if norms[i]==4 for j in range(210))
 q=(1+s.sqrt(5))/4;gamma=s.simplify(q-s.Rational(4,5));assert sign(gamma)>0
 A=s.Matrix([[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,2]])
 assert A.inv()[0,0]==s.Rational(4,5)
 H=A-s.diag(1/q,0,0,0);assert all(sign(H[:i,:i].det())>0 for i in range(1,5))
 assert s.simplify(1/(2*H[0,0])-(3+s.sqrt(5))/8)==0
 # Complete finite-label pool whose squared cosine lies in Q(sqrt5).
 # Degree phi(m)/2<=2 bounds m<=32 using phi(m)^2>=m/2.
 small=[m for m in range(2,33)if s.totient(m)<=4]
 assert small==[2,3,4,5,6,8,10,12]
 pool={2:s.Integer(0),3:s.Rational(1,4),4:s.Rational(1,2),
       5:(3+s.sqrt(5))/8,6:s.Rational(3,4),10:(5+s.sqrt(5))/8}
 assert s.simplify(s.cos(s.pi/8)**2)==s.sqrt(2)/4+s.Rational(1,2)
 assert s.simplify(s.cos(s.pi/12)**2)==s.sqrt(3)/4+s.Rational(1,2)
 # Different quadratic fields exclude labels8,12. For each source root,
 # positive projected norm leaves precisely the following d,k possibilities.
 kinds={'rational_short':(2,0),'rational_long':(4,0),'modified_short':(2,-1),'modified_long':(4,-2)}
 for d in[2,4]:
  possibilities=[k for k in range(-5,1)if (d!=4 or k%2==0)and sign(d-q*k*k)>0]
  assert possibilities==([-1,0]if d==2 else[-2,0])
  assert sign(d-q*25)<0
 table=[];forbidden=set();counts=Counter()
 for k1,k2 in combinations_with_replacement(kinds,2):
  d1,a=kinds[k1];d2,b=kinds[k2]
  if a==b==0:continue
  n1=d1-q*a*a;n2=d2-q*b*b;prod=s.expand(n1*n2)
  # If |g|>4, then |g-qab|>sqrt(d1*d2)>=sqrt(n1*n2),
  # so no additional finite case can be missed by this finite integer table.
  for g in range(-4,1):
   if 4 in[d1,d2]and g%2:continue
   p=g-q*a*b;ratio=s.radsimp(s.expand(p*p/prod));delta=s.radsimp(p*p-prod)
   if p==0:status='RIGHT_RATIONAL_NEIGHBOR';assert (a==0)!=(b==0)
   elif sign(delta)>=0:status='DISJOINT_OR_PARALLEL'
   else:
    assert sign(p)<0 and sign(ratio)>0
    assert all(s.simplify(ratio-v)!=0 for v in pool.values())
    status='FINITE_FORBIDDEN_FOR_EVERY_INTEGER_LABEL';forbidden.add(s.simplify(ratio))
   counts[status]+=1;table.append({'types':[k1,k2],'old_pairing':g,'new_pairing':str(s.expand(p)),
     'norms':[str(s.expand(n1)),str(s.expand(n2))],'cosine_squared':str(s.radsimp(ratio)),'status':status})
 expected=[q*q/(2-q)**2,1/(2*(2-q)),1/(2-q)]
 assert len(forbidden)==3 and all(any(s.simplify(x-y)==0 for y in expected)for x in forbidden)
 # A nonorthogonal H3-face boundary must extend H3 to H4. Its far-end
 # label3 requires d-q*k^2=2*g^2. Only an unchanged short root can do so.
 h3=s.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]])-s.diag(1/q,0,0)
 assert all(sign(h3[:i,:i].det())>0 for i in range(1,4))
 attachments=[]
 for d,k in[(2,0),(2,-1),(4,0),(4,-2)]:
  for g in[-4,-3,-2,-1]:
   if d==4 and g%2:continue
   if s.simplify(d-q*k*k-2*g*g)==0:attachments.append([d,k,g])
 assert attachments==[[2,0,-1]]
 walk=json.loads((ROOT/WALK).read_text());assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'
 seen=set();components=0;components3=0
 for corner in walk['corner_orbits']:
  for edge in corner['edge_orbits']:
   E=tuple(edge['edge_walls'])
   if E in seen:continue
   seen.add(E);left=set(E)
   while left:
    cc=[min(left)];left.remove(cc[0])
    for i in cc:
     adj=[j for j in sorted(left)if G[i,j]];cc+=adj;left.difference_update(adj)
    if len(cc)in[3,4] and all(G[i,i]==2 for i in cc)and all(G[i,j]in[0,-1]for i in cc for j in cc if i!=j):
     degrees=sorted(sum(G[i,j]!=0 for j in cc if i!=j)for i in cc)
     if degrees==[1,1,2,2]:components+=1
     if degrees==[1,1,2]:components3+=1
 assert len(seen)==6492 and components==1660 and components3==2115
 primary=json.loads((ROOT/f'results/{TAG}.json').read_text())
 def fld(x):return s.Rational(x[0])+s.Rational(x[1])*s.sqrt(5)
 for row in primary['modified_facet_adjacency_table']:
  (d1,k1),(d2,k2)=[kinds[k]for k in row['types']]
  n1=d1-q*k1*k1;n2=d2-q*k2*k2;p=row['old_pairing']-q*k1*k2
  assert s.simplify(fld(row['projected_norms'][0])-n1)==0
  assert s.simplify(fld(row['projected_norms'][1])-n2)==0
  assert s.simplify(fld(row['projected_pairing'])-p)==0
  assert s.simplify(fld(row['rank2_determinant'])-(n1*n2-p*p))==0
  if 'squared_cosine'in row:assert s.simplify(fld(row['squared_cosine'])-p*p/(n1*n2))==0
 sources=[RAW,WALK,f'scripts/{TAG}_review.py',f'scripts/{TAG}.py',f'results/{TAG}.json']
 assert all((ROOT/p).exists()for p in sources)
 return {'status':'PASS_INDEPENDENT_SYMPY_H4_FACE_CLOSURE_REVIEW','sources':{p:sha(p)for p in sources},
    'all210_source_reflections_integral':True,'every_long_root_pairing_even':True,
    'source_signature_by_characteristic_polynomial_Sturm_counts':[int(positive),int(negative)],
    'norm_counts':{str(k):v for k,v in Counter(norms).items()},'positive_projected_root_types':{k:list(v)for k,v in kinds.items()},
    'finite_label_pool_in_Qsqrt5':sorted(pool),'pair_table':table,'pair_status_counts':dict(counts),
    'forbidden_cosine_squares':sorted(map(str,forbidden)),
    'literal_source_edge_supports':len(seen),'A4_components':components,'A3_components':components3,
    'continuous_A4_endpoint_families':2*components,'continuous_A3_endpoint_families':2*components3,
    'H3_to_H4_only_old_root_attachment':attachments,
    'scope':'Separate exact arithmetic and source/theorem review; no replay of152368 numerical proposals. Full proof additionally uses compact H4-face geometry, rational invariant sublattices, and Meyer isotropy/nonuniformity.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();z=run()
 if args.verify:assert z==json.loads(OUT.read_text())
 else:OUT.write_text(json.dumps(z,indent=2)+'\n')
 print(z['status']);print(z['pair_status_counts'])
