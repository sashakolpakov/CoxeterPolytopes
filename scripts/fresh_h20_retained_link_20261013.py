#!/usr/bin/env python3
"""Exact data and inverse-metric gates for mixed affine-factor cusp sections.

Only the standard library is used. The all-real reduction is proved in the
matching report; this rechecks every literal affine component and every
connected one-facet deletion used by that reduction.
"""
from fractions import Fraction as F
from itertools import combinations
from collections import Counter
from pathlib import Path
import argparse, hashlib, json, math

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_retained_link_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK='results/fresh_h21_source_corner_walk_20261013.json'
ABS='results/fresh_h20_cusp_sections_20261013_absorption.json'
OUT=ROOT/f'results/{TAG}.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b):return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x):return [dot(r,x) for r in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def encode(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [encode(y) for y in x]
 if isinstance(x,dict):return {str(k):encode(v) for k,v in x.items()}
 return x
def inverse(A):
 n=len(A);M=[[F(x) for x in r]+[F(i==j) for j in range(n)] for i,r in enumerate(A)]
 for j in range(n):
  k=next(i for i in range(j,n) if M[i][j]);M[j],M[k]=M[k],M[j]
  q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:
    q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:] for r in M]
def pivots(A):
 A=[[F(x) for x in r] for r in A];out=[]
 while A:
  q=A[0][0];assert q>0;out.append(q)
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))] for i in range(1,len(A))]
 return out
def components(A):
 left=set(range(len(A)));out=[]
 while left:
  todo=[min(left)];left.remove(todo[0])
  for i in todo:
   neighbors=[j for j in sorted(left) if A[i][j]];todo+=neighbors;left.difference_update(neighbors)
  out.append(sorted(todo))
 return out

def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-int(x) for x in r] for r in raw['basis_gram']]
 W=[[int(x) for x in r] for r in raw['wall_coordinates_as_columns_transposed']]
 assert len(W)==210 and all(len(w)==22 for w in W)
 C=[mv(Q,w) for w in W];G=[[dot(w,c) for c in C] for w in W]
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 assert set(G[i][i] for i in range(210))=={2,4}
 total=[sum(w[d] for w in W) for d in range(22)];gg=math.gcd(*total);h=[x//gg for x in total]
 assert pair(h,h,Q)==-14 and all(dot(h,c)<0 for c in C)
 walk=json.loads((ROOT/WALK).read_text())
 assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS' and not walk['pending_orbit_ids']
 assert walk['sources'][RAW]==sha(RAW)
 cusps=[c for c in walk['corner_orbits'] if c['type']=='ideal'];assert len(cusps)==18
 rows=[];counts=Counter();rank_one=[];cusp_rows=[];a1_bounds=[]
 for cusp in cusps:
  u=cusp['primitive_ray'];ids=cusp['active_walls'];up=[dot(u,c) for c in C]
  assert pair(u,u,Q)==0 and pair(u,h,Q)<0 and max(up)==0
  assert [i for i,z in enumerate(up) if not z]==ids
  factors=cusp['affine_components']
  assert sorted(sum([c['walls'] for c in factors],[]))==ids
  assert sum(c['rank'] for c in factors)==20
  assert [[ids[i] for i in cc] for cc in components([[G[i][j] for j in ids] for i in ids])]==[c['walls'] for c in factors]
  factor_rows=[]
  for comp in factors:
   ws=comp['walls'];n=comp['rank'];ks=list(map(F,comp['positive_kernel']));lam=F(comp['ray_multiple'])
   assert n==len(ws)-1 and min(ks)>0 and lam>0
   A=[[G[i][j] for j in ws] for i in ws]
   assert len(components(A))==1 and all(dot(r,ks)==0 for r in A)
   assert [sum(k*W[i][d] for k,i in zip(ks,ws,strict=True)) for d in range(22)]==[lam*x for x in u]
   finite_angles=Counter()
   for i,j in combinations(range(len(ws)),2):
    z=F(A[i][j]**2,A[i][i]*A[j][j]);assert z in[0,F(1,4),F(1,2),1]
    if z==1:assert n==1
    finite_angles[str(z)]+=1
   fr={'walls':ws,'rank':n,'positive_kernel':ks,'null_multiple':lam,'old_cosine_squared_counts':dict(finite_angles)}
   factor_rows.append(fr);counts['affine_factors']+=1
   if n==1:
    rank_one.append({'cusp':cusp['id'],'walls':ws});continue
   comp_bounds=[]
   for omitted in ws:
    keep=[i for i in ws if i!=omitted];B0=[[G[i][j] for j in keep] for i in keep]
    ps=pivots(B0);B=inverse(B0)
    assert all(x>=0 for r in B for x in r)
    counts['positive_definite_deletions']+=1
    if len(components(B0))!=1:
     counts['disconnected_deletions']+=1;continue
    assert all(x>0 for r in B for x in r)
    counts['connected_deletions']+=1
    attachments=[]
    for pos,i in enumerate(keep):
     neighbors=[j for j in keep if j!=i and G[i][j]<0];assert neighbors
     # Prefer the stronger pi/3 gate if one is present.
     j=min(neighbors,key=lambda j:(F(G[i][j]**2,G[i][i]*G[j][j]),j))
     old_cos=F(G[i][j]**2,G[i][i]*G[j][j]);assert old_cos in[F(1,4),F(1,2)]
     D=B[pos][pos]*G[i][i]
     bound=1/(1-old_cos);assert D>=bound
     item={'wall':i,'neighbor':j,'old_cosine_squared':old_cos,'inverse_diagonal':B[pos][pos],
           'normal_square':G[i][i],'normalized_inverse_diagonal':D,'two_node_lower_bound':bound}
     attachments.append(item);comp_bounds.append(D)
     counts['connected_attachment_profiles']+=1;counts['old_pi3_profiles' if old_cos==F(1,4) else 'old_pi4_profiles']+=1
    rows.append({'cusp':cusp['id'],'factor':ws,'omitted_wall':omitted,'retained':keep,'LDL':ps,'attachments':attachments})
   if cusp['id']==0:
    assert all(z in['0','1/4'] for z in finite_angles)
    expected=F(2) if n==5 else F(7,4)
    assert n in[5,7] and min(comp_bounds)==expected
    a1_bounds.append({'factor':ws,'rank':n,'minimum_normalized_inverse_diagonal':min(comp_bounds)})
  cusp_rows.append({'id':cusp['id'],'ray':u,'factor_ranks':[c['rank'] for c in factors],'factors':factor_rows})
 assert rank_one==[{'cusp':0,'walls':[162,164]}]
 assert sorted(x['rank'] for x in a1_bounds)==[5,7,7]
 # For delta=sqrt(5)-2: 2delta>1/3, using rational squares alone.
 assert F(5)>F(13,6)**2 and math.isqrt(5)**2!=5
 # Old pi3 -> m>=4, or old pi4 -> m>=6: inverse-metric norm cost >=2/3.
 assert (1-F(1,4)/F(1,2))/(1-F(1,4))==F(2,3)
 assert (1-F(1,2)/F(3,4))/(1-F(1,2))==F(2,3)
 # A7 + retained A1: at D>=7/4 and s>=1/2, the cross cosine is <=1/7.
 assert (1-F(7,4)*F(1,2))/(F(7,4)*(1-F(1,2)))==F(1,7)<F(1,4)
 # The only omitted-A1 exception is the separately frozen all-real theorem.
 absorption=json.loads((ROOT/ABS).read_text())
 assert absorption['status']=='PASS_CUSP0_A1_ABSORPTION_OBSTRUCTION' and absorption['cusp_id']==0
 for p,z in absorption['sources'].items():assert sha(p)==z
 assert absorption['width_squared_ratio']=='49/8' and absorption['pure_A1_plane']['actual_ridge']==[45,63]
 return encode({'status':'PASS_ALL_18_MIXED_FACTOR_CUSP_SECTION_OBSTRUCTION',
  'sources':{p:sha(p) for p in[RAW,WALK,ABS,f'scripts/{TAG}.py']},
  'cusp_count':18,'counts':dict(counts),'cusp_factors':cusp_rows,
  'unique_A1_factor':rank_one,'connected_deletion_records':rows,'cusp0_A1_attachment_bounds':a1_bounds,
  'universal_cost_gates':{'old_pi3_to_any_larger_label':'cost >= 2/3',
    'old_pi4_to_any_larger_label':'cost >= 2*(sqrt(5)-2) > 1/3',
    'old_pi4_to_label_at_least_6':'cost >= 2/3',
    'only_remaining_two_factor_case':'two pi4->pi5 attachments, total cost a nonzero rational multiple of sqrt(5)-2, not 1',
    'A1_plus_A7_cross_cosine_squared_upper':'1/7','A1_plus_D5_cost_lower':'1'},
  'scope':'Combining the report\'s all-real simplex/triangle reduction, the exact rational inverse gates, and the frozen A1-absorption theorem excludes every section with tangent normal active in at least two irreducible affine factors at any of the 18 cusp representatives. Arbitrary affine shift is included. Single-factor directions and sections avoiding all these cusps are not excluded.'})

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 ans=run()
 if args.verify:assert json.loads(OUT.read_text())==ans
 else:OUT.write_text(json.dumps(ans,indent=2)+'\n')
 print(ans['status'],ans['counts'])
