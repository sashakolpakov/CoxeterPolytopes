#!/usr/bin/env python3
"""Exact all-real one-coordinate Coxeter lifts of the fixed36 Vinberg directions.

The 24 short heights vanish by a graph/angle argument; each remaining height
has one of seven values.  All147 long-pair tests and the complete CSP use
standalone Fraction arithmetic in Q(sqrt2,sqrt5).  No numerical gate is used.
"""
from fractions import Fraction as F
from pathlib import Path
from itertools import combinations,product
from collections import Counter,deque
import argparse,copy,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/vinberg_x3_h19.json'
OUT='data/enumerations/fresh_dimension_bridge_20261013_lift_classification.json'
Z=(F(0),)*4;O=(F(1),F(0),F(0),F(0))
def num(n):return (F(n),F(0),F(0),F(0))
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def neg(x):return tuple(-a for a in x)
def sub(x,y):return add(x,neg(y))
def scale(x,a):return tuple(v*a for v in x)
def mul(x,y):
 z=[F(0)]*4
 for i in range(4):
  for j in range(4):z[i^j]+=x[i]*y[j]*(2 if i&j&1 else 1)*(5 if i&j&2 else 1)
 return tuple(z)
def qmul(x,y):return(x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def qsign(x):
 a,b=x
 if not b:return(a>0)-(a<0)
 if not a or(a>0)==(b>0):return(b>0)-(b<0)
 d=a*a-2*b*b;return 0 if not d else ((a>0)-(a<0))if d>0 else((b>0)-(b<0))
def sign(x):
 a,b=x[:2],x[2:];sa,sb=qsign(a),qsign(b)
 if not sb:return sa
 if not sa or sa==sb:return sb
 aa,bb=qmul(a,a),qmul(b,b);d=qsign(tuple(aa[i]-5*bb[i]for i in range(2)))
 return 0 if not d else sa if d>0 else sb

def inv(x):
 # Invert first over Q(sqrt2), then multiply by the sqrt5 conjugate.
 a,b=x[:2],x[2:];aa,bb=qmul(a,a),qmul(b,b);n=tuple(aa[i]-5*bb[i]for i in range(2));den=n[0]*n[0]-2*n[1]*n[1];assert den
 ni=(n[0]/den,-n[1]/den);ra=qmul(a,ni);rb=qmul(b,ni)
 return ra+tuple(-u for u in rb)
def enc(x):return list(map(str,x))
def poly_eval(p,x):
 z=Z
 for c in reversed(p):z=add(mul(z,x),num(c))
 return z
# Polynomials for cos²(pi/m), constant term first, and isolating interval.
ANGLES={2:([0,1],F(-1,1000),F(1,1000)),3:([-1,4],F(249,1000),F(251,1000)),
4:([-1,2],F(499,1000),F(501,1000)),5:([1,-12,16],F(654,1000),F(655,1000)),
6:([-3,4],F(749,1000),F(751,1000)),8:([1,-8,8],F(853,1000),F(854,1000)),
10:([5,-20,16],F(904,1000),F(905,1000)),12:([1,-16,16],F(933,1000),F(934,1000)),
15:([1,-96,416,-576,256],F(956,1000),F(957,1000)),16:([1,-32,160,-256,128],F(961,1000),F(962,1000)),
20:([1,-48,304,-512,256],F(975,1000),F(976,1000)),24:([1,-64,320,-512,256],F(982,1000),F(983,1000)),
30:([1,-32,224,-448,256],F(989,1000),F(990,1000))}
def trim(p):
 p=list(map(F,p))
 while len(p)>1 and not p[-1]:p.pop()
 return p

def padd(p,q):return trim([(p[i]if i<len(p)else 0)+(q[i]if i<len(q)else 0)for i in range(max(len(p),len(q)))])
def pmul(p,q):
 z=[F(0)]*(len(p)+len(q)-1)
 for i,a in enumerate(p):
  for j,b in enumerate(q):z[i+j]+=a*b
 return trim(z)
def pdiv(p,q):
 p=trim(p);q=trim(q);assert q!=[0];out=[F(0)]*max(1,len(p)-len(q)+1)
 while p!=[0]and len(p)>=len(q):
  i=len(p)-len(q);a=p[-1]/q[-1];out[i]=a
  for j,b in enumerate(q):p[i+j]-=a*b
  p=trim(p)
 return trim(out),p

def pgcd(p,q):
 while trim(q)!=[0]:p,q=q,pdiv(p,q)[1]
 return[x/p[-1]for x in trim(p)]
def pval(p,x):
 z=F(0)
 for a in reversed(p):z=z*x+a
 return z

def variations(chain,x):
 sig=[(v>0)-(v<0)for p in chain if(v:=pval(p,x))]
 return sum(a!=b for a,b in zip(sig,sig[1:]))
def sturm(p):
 out=[trim(p),trim([i*p[i]for i in range(1,len(p))])]
 while out[-1]!=[0]:
  r=pdiv(out[-2],out[-1])[1]
  if r==[0]:break
  out.append([-x for x in r])
 return out

def cosine_poly(m):
 a=[F(1)];b=[F(-1),F(2)]
 for i in range(2,m+1):a,b=b,padd(pmul([-2,4],b),[-x for x in a])
 b[0]-=1;return trim(b)
def angle_check():
 for m,(p,l,u)in ANGLES.items():
  ch=sturm(p);assert variations(ch,l)-variations(ch,u)==1
  assert variations(ch,u)-variations(ch,F(1))==0
  assert pdiv(cosine_poly(m),p)[1]==[0]
  for d in range(1,m):
   if m%d==0:assert len(pgcd(p,cosine_poly(d)))==1
  assert len(pgcd(p,[i*p[i]for i in range(1,len(p))]))==1
  if m>2:assert len(p)-1==sum(math.gcd(i,m)==1 for i in range(1,m))//2
 # Complete degrees for phi(m)/2 dividing4. Any odd prime p|m has p-1|8,
 # so p belongs to{3,5}; exponent bounds are 2^a,a<=4,3^b,b<=1,5^c,c<=1.
 got=[]
 for a,b,c in product(range(5),range(2),range(2)):
  m=2**a*3**b*5**c
  if m<3:continue
  degree=sum(math.gcd(i,m)==1 for i in range(1,m))//2
  if degree in(1,2,4):got.append(m)
 assert sorted(got)==[m for m in ANGLES if m>2]

def pair_label(p,qa,qb):
 if sign(p)>0:return'positive',None
 if p==Z:return'2',Z
 x=mul(mul(p,p),inv(mul(qa,qb)));s=sign(sub(x,O))
 if s>=0:return('parallel'if s==0 else'dotted'),x
 for m,(poly,l,u)in ANGLES.items():
  if poly_eval(poly,x)==Z and sign(sub(x,num(l)))>0 and sign(sub(num(u),x))>0:return str(m),x
 return'forbidden',x

def run():
 angle_check();raw=json.loads((ROOT/SOURCE).read_text());G=[[-x for x in row]for row in raw['wall_gram']]
 adj={i:{j for j in range(24)if j!=i and G[i][j]<0}for i in range(24)}
 seen={0};queue=[0]
 while queue:
  for j in adj[queue.pop()]:
   if j not in seen:seen.add(j);queue.append(j)
 assert len(seen)==24
 assert all(G[i][j]==-1 for i in range(24)for j in adj[i])
 assert not any(k in adj[i]for i in range(24)for j in adj[i]for k in adj[j]if k!=i)
 # At most2 same-positive and2 same-negative heights in this triangle-free
 # graph. Hence a zero exists, and its pi/3 edge propagates zeros everywhere.
 phi=(F(1,2),F(0),F(1,2),F(0));gold=(F(0),F(3,2),F(0),F(-1,2))
 hs=[Z,num(2),num(-2),num(1),num(-1),gold,neg(gold)]
 assert mul(gold,gold)==sub(num(10),scale(phi,6))
 norms=[add(num(6),scale(mul(h,h),3))for h in hs]
 assert all(sign(q)>0 for q in norms)
 for a,m in [(1,3),(3,4),(5,5)]:assert pair_label(num(-3),num(2),norms[a])[0]==str(m)
 tables={};records=[]
 for g in(-6,-12,-15):
  table=[set()for _ in hs]
  for a,b in product(range(7),repeat=2):
   p=add(num(g),scale(mul(hs[a],hs[b]),3));lab,x=pair_label(p,norms[a],norms[b])
   if lab not in('positive','forbidden'):table[a].add(b)
   records.append({'old_pairing':g,'height_codes':[a,b],'pairing':enc(p),'label':lab,'cosine_squared':None if x is None else enc(x)})
  tables[g]=table
 eps={24:1};queue=[24]
 while queue:
  i=queue.pop()
  for j in range(24,36):
   if G[i][j]!=-6:continue
   if j not in eps:eps[j]=-eps[i];queue.append(j)
   assert eps[j]==-eps[i]
 assert len(eps)==12
 solutions=[]
 def extend(assignment,domains):
  if len(assignment)==12:solutions.append(tuple(assignment[i]for i in range(12)));return
  i=min((j for j in range(12)if j not in assignment),key=lambda j:len(domains[j]))
  for v in sorted(domains[i]):
   nextdom=[set(x)for x in domains];ok=True
   for j in range(12):
    if j==i or j in assignment:continue
    nextdom[j]&=tables[G[i+24][j+24]][v]
    if not nextdom[j]:ok=False;break
   if ok:extend(assignment|{i:v},nextdom)
 extend({},[set(range(7))for _ in range(12)])
 rational=set()
 for phase in(1,-1):
  for bits in product((0,1),repeat=12):
   signs=[phase*eps[i+24]for i in range(12)];vals=[2*signs[i]if bits[i]==0 else-signs[i]for i in range(12)]
   rational.add(tuple(hs.index(num(x))for x in vals))
 golden=[tuple(5 if phase*eps[i+24]>0 else 6 for i in range(12))for phase in(1,-1)]
 assert set(solutions)==rational|set(golden)|{(0,)*12}
 assert len(rational)==8192 and len(solutions)==8195
 # Reconstruct a primitive positive affine D16 vector, avoiding root solvers:
 # use its known graph arms and solve linear equations by Fraction elimination.
 from fresh_dimension_bridge_20261013_verify import inverse,rank,E,scalar,mv,dot,vec
 geometry=json.loads((ROOT/'data/enumerations/fresh_dimension_bridge_20261013.json').read_text())
 ci=next(i for i,c in enumerate(raw['vertices_and_cusps']['cusps'])if'D16~'in c['component_types'])
 cusp=raw['vertices_and_cusps']['cusps'][ci];component=next(c for c in cusp['component_facets_one_based']if len(c)==17)
 J=[j-1 for j in component];A=[[F(G[i][j])for j in J[:-1]]for i in J[:-1]];iv=inverse(A)
 assert all(j<24 for j in J) # This forced affine relation is unchanged by long-wall heights.
 coef=[-sum(iv[i][j]*G[J[j]][J[-1]]for j in range(16))for i in range(16)]+[F(1)]
 assert all(x>0 for x in coef)
 assert all(sum(G[i][j]*coef[z]for z,j in enumerate(J))==0 for i in J)
 long=next(i-1 for i in cusp['facets_one_based']if i>24)
 short=next(i-1 for i in cusp['facets_one_based']if i<=24 and G[i-1][long]==-3)
 delta_pairings=[sum(G[i][j]*coef[z]for z,j in enumerate(J))for i in range(36)]
 assert all(x<=0 for x in delta_pairings)and delta_pairings[short]==delta_pairings[long]==0
 basis_inverse=inverse([row[:20]for row in G[:20]])
 old=[tuple(scalar(sum((basis_inverse[i][j]*G[j][n]for j in range(20)),F(0)))for i in range(20))+(E(),)for n in range(36)]
 Q=tuple(tuple(scalar(G[i][j])if i<20 and j<20 else scalar(3)if i==j==20 else E()for j in range(21))for i in range(21))
 h=tuple(sum((r[j]for r in old),E())for j in range(21))
 assert old==list(map(vec,geometry['old_normals']))and h==vec(geometry['strict_center'])
 delta=tuple(sum((scalar(coef[z])*old[j][k]for z,j in enumerate(J)),E())for k in range(21))
 assert dot(delta,mv(Q,delta))==E()and dot(delta,mv(Q,h)).sign()<0
 q5=norms[5];cos2=mul(num(9),inv(scale(q5,2)))
 assert cos2==scale(add(num(3),(F(0),F(0),F(1),F(0))),F(1,8))
 trace=add(num(15),scale(cos2,4));assert trace==add(num(16),phi)
 return {'status':'COMPLETE_EXACT_FIXED_PROJECTION_ONE_COORDINATE_LIFT_CLASSIFICATION',
         'source':SOURCE,'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
         'extra_coordinate_square':3,'height_field':'Q(sqrt2,sqrt5), basis 1,sqrt2,sqrt5,sqrt10',
         'height_codes':list(map(enc,hs)),'short_zero_proof':{'connected_vertices':24,'negative_pair_graph_triangle_free':True,'maximum_same_sign_support':2,'forced_zero_by_pi_over_3_propagation':True},
         'safe_finite_angle_labels':list(ANGLES),'angle_polynomials_and_isolations':{str(m):{'polynomial':p,'interval':[str(l),str(u)]}for m,(p,l,u)in ANGLES.items()},
         'long_pair_checks':records,'height_signs':[[i,eps[i]]for i in sorted(eps)],
         'labelled_assignments_including_zero':8195,'rational_pi3_pi4_assignments':8192,'golden_pi5_assignments':[list(x)for x in golden],
         'solutions_digest':hashlib.sha256(json.dumps(sorted(solutions)).encode()).hexdigest(),
         'rational_family_formula':'h_i in {2 phase eps_i, -phase eps_i}, independently for12long walls, phase=+/-1; simultaneous paired48walls for each phase',
         'golden_forced_cusp':{'source_cusp_index':ci,'affine_D16_walls':J,'positive_coefficients':list(map(str,coef)),
             'ray':[[str(x.a),str(x.b)]for x in delta],'source_wall_pairings':list(map(str,delta_pairings)),
             'incident_pi5_pair':[short,long],'lifted_long_norm':enc(q5),'cosine_squared':enc(cos2),'horosphere_reflection_product_trace':enc(trace),'trace_formula':'16+phi'},
         'scope':'Complete lift classification retains all36 original projected wall directions as Coxeter facets. The two golden lifts have no finite-volume acute-compatible completion because a forced cusp retains I2(5). The8192 rational lifts and the48/70-wall extensions are not excluded from arbitrary further completion.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();r=run()
 if a.verify:assert json.loads((ROOT/OUT).read_text())==r
 else:(ROOT/OUT).write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps({k:r[k]for k in ['status','labelled_assignments_including_zero','rational_pi3_pi4_assignments','solutions_digest']}),flush=True)
