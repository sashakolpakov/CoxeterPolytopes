#!/usr/bin/env python3
"""Exact graph/classification gate for two H17 sources and rank-three attachments."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_h17_rank_three_coupling_20261015'
P19='data/known_examples/k3_high_rank_nef_chambers.json'
P48='data/enumerations/coble_type_mii_full_ns_h17.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dot(x,y):return -x[0]*y[1]-x[1]*y[0]+sum(a*b for a,b in zip(x[2:],y[2:]))
def gram(N):return[[dot(x,y)for y in N]for x in N]
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,list):return[enc(x)for x in v]
 if isinstance(v,dict):return{k:enc(x)for k,x in v.items()}
 return v
def component_data(G):
 adj={i:{j for j in range(len(G))if i!=j and G[i][j]<0 and 4*G[i][j]**2==G[i][i]*G[j][j]}for i in range(len(G))}
 seen=set();answer=[]
 for i in range(len(G)):
  if i in seen:continue
  todo=[i];seen.add(i);C=[]
  while todo:
   a=todo.pop();C.append(a)
   for b in sorted(adj[a]-seen):seen.add(b);todo.append(b)
  C=sorted(C);clique=next((z for z in combinations(C,4)if all(G[a][b]==0 for a,b in combinations(z,2))),None)
  answer.append(dict(vertices=C,pi3_edges=[[a,b]for a,b in combinations(C,2)if b in adj[a]],orthogonal_four=None if clique is None else list(clique)))
 return answer

def run():
 a=json.loads((ROOT/P19).read_text());case=next(c for c in a['cases']if c['rho']==18)
 assert case['wall_count']==19 and case['hyperbolic_dimension']==17
 G19=[[2*int(i==j)for j in range(19)]for i in range(19)]
 for value,pairs in case['intersections'].items():
  for i,j in pairs:assert G19[i-1][j-1]==0;G19[i-1][j-1]=G19[j-1][i-1]=-int(value)
 C19=component_data(G19);assert len(C19)==1 and C19[0]['orthogonal_four']==[0,2,5,7]
 b=json.loads((ROOT/P48).read_text());labels=b['facets']['labels']
 # Reconstruct the physical48 roots without importing a source producer.
 N=[]
 for row in range(4):
  v=[0]*18;v[0]=1
  for col in range(4):v[2+4*row+col]=-1
  N.append(v)
 for col in range(4):
  v=[0]*18;v[1]=1
  for row in range(4):v[2+4*row+col]=-1
  N.append(v)
 for i in range(16):v=[0]*18;v[2+i]=1;N.append(v)
 sigmas=list(permutations(range(4)))
 for sigma in sigmas:
  v=[1,1]+[0]*16
  for row,col in enumerate(sigma):v[2+4*row+col]=-1
  N.append(v)
 assert N==b['facets']['physical_NS_wall_vectors']
 G=gram(N);assert G==[[-x for x in row]for row in b['facets']['physical_wall_gram']]
 C=component_data(G);assert[len(c['vertices'])for c in C]==[24,12,12]
 assert C[0]['orthogonal_four']==[0,1,2,3]
 assert s.Matrix(N).rank()==18
 cosets=[]
 for comp in C[1:]:
  todo=set(comp['vertices']);parts=[]
  while todo:
   i=min(todo);part={i}|{j for j in todo if G[i][j]==-2};assert len(part)==4
   todo-=part;parts.append(sorted(part))
  assert len(parts)==3
  for i,j in combinations(comp['vertices'],2):
   same=any(i in z and j in z for z in parts);assert G[i][j]==(-2 if same else-1)
  for i in comp['vertices']:
   assert G[i][i]==2
   E=[j for j in range(8,24)if G[i][j]!=0]
   assert len(E)==4 and all(G[i][j]==-1 and G[j][j]==1 for j in E)
  cosets.append(parts)
 crossright=[i,j] = next([i,j]for i in C[1]['vertices']for j in C[2]['vertices']if G[i][j]==0)
 assert all(G[i][j]in(0,-2)for i in C[1]['vertices']for j in C[2]['vertices'])
 # Exact cost-one dominant-weight classifications. A nonzero finite contact
 # has magnitude1 (pi3) or at least sqrt2; parallel/dotted magnitudes>=2.
 rt2,rt5=s.sqrt(2),s.sqrt(5);phi=(1+rt5)/2
 support_grams={
  'A3':s.Matrix([[2,-1,0],[-1,2,-1],[0,-1,2]]),
  'B3':s.Matrix([[2,-1,0],[-1,2,-rt2],[0,-rt2,2]]),
  'H3':s.Matrix([[2,-phi,0],[-phi,2,-1],[0,-1,2]]),
  'A2+A1':s.Matrix([[2,-1,0],[-1,2,0],[0,0,2]]),
  'I2(4)+A1':s.Matrix([[2,-rt2,0],[-rt2,2,0],[0,0,2]]),
  'A1^3':2*s.eye(3)}
 profiles={'A3':[[0,1,0]],'B3':[[1,0,0]],'H3':[],
           'A2+A1':[[0,0,rt2]],'I2(4)+A1':[[1,0,0],[0,1,0],[0,0,rt2]],
           'A1^3':[[rt2,0,0],[0,rt2,0],[0,0,rt2],[1,1,0],[1,0,1],[0,1,1]]}
 tables={}
 for name,B in support_grams.items():
  J=s.simplify(B.inv());assert all(x>=0 for x in J)
  for i,j in combinations(range(3),2):
   # Except A1^3 every multi-contact cost exceeds1, using just diagonals.
   if name!='A1^3':assert s.simplify(J[i,i]+J[j,j])>1
  qs=[s.simplify((s.Matrix(c).T*J*s.Matrix(c))[0])for c in profiles[name]];assert all(q==1 for q in qs)
  dots=[[s.simplify((s.Matrix(c).T*J*s.Matrix(d))[0])for d in profiles[name]]for c in profiles[name]]
  halfedges=[[i,j]for i,j in combinations(range(len(dots)),2)if dots[i][j]==s.Rational(1,2)]
  assert halfedges==([[3,4],[3,5],[4,5]]if name=='A1^3'else[])
  tables[name]=dict(gram=[[str(x)for x in row]for row in B.tolist()],inverse_diagonal=[str(J[i,i])for i in range(3)],
                   cost_one_profiles=[[str(x)for x in row]for row in profiles[name]],half_product_edges=halfedges,
                   pair_products=[[str(x)for x in row]for row in dots])
 assert s.simplify(support_grams['A3'].inv()[0,0])==s.Rational(3,4)
 assert list(support_grams['A3'].inv().diagonal())==[s.Rational(3,4),1,s.Rational(3,4)]
 assert list(support_grams['B3'].inv().diagonal())==[1,2,s.Rational(3,2)]
 assert list(support_grams['A2+A1'].inv().diagonal())==[s.Rational(2,3),s.Rational(2,3),s.Rational(1,2)]
 assert list(support_grams['I2(4)+A1'].inv().diagonal())==[1,1,s.Rational(1,2)]
 assert 1<s.Rational(4,3)<2 and 1<s.Rational(3,2)<2
 assert all(x>1 for x in support_grams['H3'].inv().diagonal())
 # I2(m)+A1, m>=5: inverse diagonal 1/(2sin²(pi/m))>=1+1/sqrt5>1;
 # therefore only the A1 pi4 profile remains, and it has no half-product pair.
 assert s.simplify(1/(2*(5-rt5)/8))==1+1/rt5
 # Construct the unique nonordinary family up to exact coordinate isometries.
 base=[v+[0,0,0]for v in N];pairs=[(0,1),(0,2),(1,2)]
 for part,pair in zip(cosets[0],pairs):
  for i in part:
   for j in pair:base[i][18+j]=-1
 for i in range(3):v=[0]*21;v[18+i]=1;base.append(v)
 H=list(map(F,[5,5]+[-1]*19))
 X0=list(map(F,[5,3,0,0,-3,0,0,0,-2,-1,-1,-2,0,0,0,-3,0,0,0,0,-3]))
 assert dot(H,H)==-31 and dot(X0,X0)==7
 E=[x+h/F(1000)for x,h in zip(X0,H)]
 assert dot(E,E)>0 and dot(E,H)<0
 def classify(G,i,j):
  g=G[i][j];q=G[i][i]*G[j][j]
  assert g<=0
  if g==0:return'right'
  value=F(4*g*g,q)
  if value==1:return'pi3'
  if value==2:return'pi4'
  if value==3:return'pi6'
  assert value>=4
  return'parallel'if value==4 else'dotted'
 GB=gram(base);hist=Counter(classify(GB,i,j)for i,j in combinations(range(51),2))
 assert all(dot(H,v)<0 for v in base)and all(dot(E,v)<0 for v in base)
 assert all(dot(X0,v)<=0 for v in base)
 assert s.Matrix(base).rank()==21
 scaled=[[F(x,2)for x in v]if i<8 or i in C[1]['vertices']else list(map(F,v))for i,v in enumerate(base)]
 GS=gram(scaled);cartan=[[2*GS[i][j]/GS[j][j]for j in range(51)]for i in range(51)]
 assert all(x.denominator==1 for row in cartan for x in row)
 # Explicit facet actuality by perpendicular projection of the common center.
 facet_points=[]
 for i,r in enumerate(base):
  q=GB[i][i];p=[h-dot(H,r)/q*x for h,x in zip(H,r)]
  assert dot(p,p)<0 and dot(p,r)==0
  assert all(dot(p,v)<0 for j,v in enumerate(base)if j!=i)
  facet_points.append(p)
 variants=[];sets=set()
 for column_swap in(False,True):
  for tp in permutations(range(3)):
   perm=list(range(21))
   if column_swap:
    for row in range(4):perm[2+4*row],perm[3+4*row]=3+4*row,2+4*row
   for j in range(3):perm[18+j]=18+tp[j]
   def transform(v):
    out=[F(0)]*21
    for j,x in enumerate(v):out[perm[j]]=x
    return out
   NN=list(map(transform,base));EE=transform(E)
   assert gram(NN)==GB and transform(H)==H
   assert all(dot(EE,n)<0 for n in NN)and dot(EE,EE)==dot(E,E)
   key=tuple(sorted(map(tuple,NN)));assert key not in sets;sets.add(key)
   # Literal projection data identify exactly one changed parity component.
   changed=[];assignment={}
   for n in NN:
    if any(n[18:]):
     if any(n[:18]):
      i=N.index([int(x)for x in n[:18]]);changed.append(i);assignment[i]=tuple(j for j in range(3)if n[18+j]==-1)
   ci=next(j for j,cc in enumerate(C[1:])if set(changed)==set(cc['vertices']))
   assert all(len({assignment[i]for i in part})==1 for part in cosets[ci])
   assert len({assignment[part[0]]for part in cosets[ci]})==3
   variants.append(dict(column_swap=column_swap,transverse_permutation=list(tp),coordinate_permutation=perm,
                        changed_component=ci,escape=EE))
 assert len(variants)==12
 paths=[P19,P48,f'scripts/{TAG}.py']
 return enc(dict(status='PASS_TWO_H17_SOURCES_ALL_RANK_THREE_DIRECT_COUPLINGS_INFINITE',sources={p:sha(p)for p in paths},
  rank_k_criterion='A pi3-connected source component with k+1 pairwiseorthogonal roots is entirely ordinary for every acute Coxeter inducing attachment to a spherical rankk support.',
  H17_19=dict(gram=G19,components=C19,one_based_orthogonal_clique=[1,3,6,8],conclusion='All19 inducing directions ordinary; transverse product escape.'),
  H17_48=dict(labels=labels,components=C,parallel_cosets=cosets,cross_parity_right_pair=crossright,
    forced_ordinary=list(range(24)),special_projection_scale='1/sqrt2',special_transverse_norm=1),
  rank_three_cost_one=tables,unbounded_I2_tail='m>=5 diagonal>=1+1/sqrt5>1; only A1 pi4 remains.',
  complete_residual_alternatives=['allordinary','A1^3; precisely one parity component special; its threeparallel4cliques assigned bijectively the threedoublepi3 weights'],
  nonordinary_representative=dict(normals=base,gram=GB,norms=[GB[i][i]for i in range(51)],pair_histogram=dict(hist),
    strict_center=H,center_norm=dot(H,H),facet_points=facet_points,boundary_spacelike=X0,boundary_norm=7,
    strict_escape=E,escape_norm=dot(E,E),escape_center_pairing=dot(E,H),escape_pairings=[dot(E,v)for v in base],
    scaled_integral_generators=scaled,cartan=cartan,ambient_rank=21),
  exact_isometric_variants=variants,
  scope='Only the spherical support walls and all walls inducing the fixed source facets. Arbitrary finite labels and arbitrary spherical rankthree support, including I2(m)+A1 tails. Extra noninducing walls, changed sourcefacet directions or omitted sourcefacets are outside. No claim of global H20 nonexistence.',
  arithmetic_triage='The12 nonordinary examples preserve the fullrank rational lattice generated by scaledroots, because the51x51 Cartan matrix is integral. No reflectivity claim; strict escape already proves these finite chambers infinitevolume.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['exact_isometric_variants']),d['nonordinary_representative']['pair_histogram'],d['nonordinary_representative']['escape_norm'])
