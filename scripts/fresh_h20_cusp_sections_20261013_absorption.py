#!/usr/bin/env python3
"""Exact local A1-absorption gates and an all-real actual obtuse ridge.

Standard-library replay; no search or previous verifier is imported.
The report supplies the elementary inequalities for the enumerated star cases.
"""
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import argparse, hashlib, json, math

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h20_cusp_sections_20261013_absorption'
RAW = 'data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
WALK = 'results/fresh_h21_source_corner_walk_20261013.json'
OUT = ROOT / f'results/{TAG}.json'
def sha(p): return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def dot(a,b): return sum(x*y for x,y in zip(a,b,strict=True))
def mv(A,x): return [dot(r,x) for r in A]
def pair(a,b,Q): return dot(a,mv(Q,b))
def enc(x):
 if isinstance(x,F): return str(x)
 if isinstance(x,(list,tuple)): return [enc(y) for y in x]
 if isinstance(x,dict): return {str(k):enc(v) for k,v in x.items()}
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
def face(h,ids,W,G,Q):
 ids=sorted(ids);A=[[G[i][j] for j in ids] for i in ids];ps=pivots(A);B=inverse(A)
 assert all(x>=0 for r in B for x in r)
 c=mv(B,[pair(h,W[i],Q) for i in ids])
 x=[h[d]-sum(z*W[i][d] for z,i in zip(c,ids,strict=True)) for d in range(22)]
 ys=[pair(x,w,Q) for w in W]
 assert pair(x,x,Q)<0 and pair(x,h,Q)<0
 assert [i for i,y in enumerate(ys) if not y]==ids
 assert all(ys[i]<0 for i in range(210) if i not in ids)
 return x,{'walls':ids,'point':x,'norm':pair(x,x,Q),'future_pairing':pair(x,h,Q),
           'maximum_other_pairing':max(ys[i] for i in range(210) if i not in ids),'LDL':ps}
def sign_supports(ids,G,kind):
 """All possible opposite-sign supports of an acute rank-one Gram update."""
 rows=[];counts=Counter()
 for signs in product((-1,0,1),repeat=len(ids)):
  P=[i for i,s in zip(ids,signs) if s>0];N=[i for i,s in zip(ids,signs) if s<0]
  if not P or not N: continue
  if any(G[i][j]==0 for i in P for j in N): continue
  assert min(len(P),len(N))==1
  center=P[0] if len(P)==1 else N[0];leaves=N if len(P)==1 else P
  assert all(G[center][i]==-1 for i in leaves)
  support=set(P+N)
  boundary=[i for i in support if any(G[i][j] and j not in support for j in ids)]
  if kind=='A7':
   assert len(leaves)<=2
   if len(leaves)==1:
    assert set(boundary)==support;case='adjacent_opposite_pair'
   else:
    assert set(leaves)<=set(boundary);case='two_neighbor_star'
  else:
   centers={29,48};outer=set(ids)-centers
   if len(leaves)==1:
    if set([center]+leaves)==centers:
     assert set(boundary)==centers;case='center_center_pair'
    else:
     assert len(support&centers)==1 and len(support&outer)==1
     assert support&centers<=set(boundary);case='leaf_center_pair'
   elif set(leaves)<=outer:
    assert len(leaves)==2 and center in centers and center in boundary
    case='two_outer_leaf_star'
   elif len(leaves)==2:
    assert center in centers and len(set(leaves)&centers)==1
    assert centers<=set(boundary);case='one_leaf_one_center_star'
   else:
    assert len(leaves)==3 and center in centers
    assert centers-{center}<=set(boundary);case='three_neighbor_star'
  rows.append({'positive':P,'negative':N,'boundary_changed_nodes':sorted(boundary),'case':case})
  counts[case]+=1
 return rows,dict(counts)

def run():
 raw=json.loads((ROOT/RAW).read_text())['coxeter_input']
 Q=[[-int(x) for x in r] for r in raw['basis_gram']]
 W=[[int(x) for x in r] for r in raw['wall_coordinates_as_columns_transposed']]
 assert len(W)==210 and all(len(w)==22 for w in W)
 C=[mv(Q,w) for w in W];G=[[dot(w,c) for c in C] for w in W]
 assert all(G[i][j]<=0 for i,j in combinations(range(210),2))
 total=[sum(w[k] for w in W) for k in range(22)];gg=math.gcd(*total);h=[x//gg for x in total]
 assert pair(h,h,Q)==-14 and all(pair(h,w,Q)<0 for w in W)
 walk=json.loads((ROOT/WALK).read_text());assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'
 assert walk['sources'][RAW]==sha(RAW) and not walk['pending_orbit_ids']
 ideals=[c for c in walk['corner_orbits'] if c['type']=='ideal'];assert len(ideals)==18
 cu=next(c for c in ideals if c['id']==0);u=cu['primitive_ray'];ids=cu['active_walls']
 assert pair(u,u,Q)==0 and pair(u,h,Q)==-8
 up=[pair(u,w,Q) for w in W]
 assert max(up)==0 and [i for i,z in enumerate(up) if not z]==ids
 components=cu['affine_components'];assert sorted(c['rank'] for c in components)==[1,5,7,7]
 assert sorted(sum([c['walls'] for c in components],[]))==ids
 for c in components:
  ws=c['walls'];ks=list(map(F,c['positive_kernel']));lam=F(c['ray_multiple'])
  A=[[G[i][j] for j in ws] for i in ws]
  assert min(ks)>0 and lam>0 and all(dot(r,ks)==0 for r in A)
  assert [sum(k*W[i][d] for k,i in zip(ks,ws,strict=True)) for d in range(22)]==[lam*x for x in u]
  pivots([r[:-1] for r in A[:-1]])
 assert [c['walls'] for c in components if c['rank']==1]==[[162,164]]
 assert W[162]==[2*u[d]-W[164][d] for d in range(22)]
 assert G[162][162]==G[164][164]==4 and G[162][164]==-4
 only_A1=[(c['id'],z['walls']) for c in ideals for z in c['affine_components'] if z['rank']==1]
 assert only_A1==[(0,[162,164])]
 # The A7 edge update is forced by the all-label elementary proof in the report.
 updates=[];supports=[]
 for c in components:
  if c['rank']!=7:continue
  ws=c['walls'];assert all(G[i][i]==2 for i in ws)
  assert all(sum(G[i][j]<0 for j in ws if i!=j)==2 for i in ws)
  assert all(F(x)==1 for x in c['positive_kernel']) and F(c['ray_multiple'])==1
  rows,counts=sign_supports(ws,G,'A7');supports.append({'walls':ws,'counts':counts,'rows':rows})
  keep=ws[:-1];B=inverse([[G[i][j] for j in keep] for i in keep])
  for a,b in combinations(ws,2):
   if G[a][b]!=-1:continue
   values={i:F(i==a)-F(i==b) for i in ws}
   coeff=mv(B,[values[i] for i in keep])+[F(0)]
   e=[sum(k*W[i][d] for k,i in zip(coeff,ws,strict=True)) for d in range(22)]
   assert [pair(e,W[i],Q) for i in ws]==[values[i] for i in ws]
   assert pair(e,e,Q)==F(7,8) and pair(e,u,Q)==0
   assert max(coeff)-min(coeff)==F(7,8)
   H=[[G[i][j]-values[i]*values[j] for j in ws] for i in ws]
   assert all(sum(r)==0 for r in H)
   pivots([r[:-1] for r in H[:-1]])
   angles=Counter()
   for i,j in combinations(range(8),2):
    assert H[i][j]<=0
    z=F(H[i][j]**2,H[i][i]*H[j][j]);assert z in[0,F(1,4),F(1,2)]
    angles[str(z)]+=1
   assert angles[str(F(1,2))]==2 and angles[str(F(1,4))]==5
   updates.append({'walls':ws,'opposite_pair':[a,b],'coefficients':coeff,'vector':e,'square':F(7,8),
                   'functional_width':F(7,8),'projected_cosine_squared_counts':dict(angles)})
 assert len(updates)==16
 # e_B=(sqrt(2)/8)*(r162-delta). Its q is1/8 and its functional width squared1/8.
 tangent=[W[162][d]-u[d] for d in range(22)]
 assert pair(tangent,tangent,Q)==4 and pair(tangent,u,Q)==0
 assert F(2,64)*4==F(1,8)
 width_A_squared=F(49,64);width_B_squared=F(1,8)
 assert width_A_squared>width_B_squared and width_A_squared/width_B_squared==F(49,8)
 # D5 topology and its only boundary equality case.
 D=[4,9,14,26,29,48]
 assert next(c['walls'] for c in components if c['rank']==5)==D
 assert [sum(G[i][j]<0 for j in D if i!=j) for i in D]==[1,1,1,1,3,3]
 marks=[1,1,1,1,2,2];assert all(sum(G[i][j]*m for j,m in zip(D,marks))==0 for i in D)
 ds,dc=sign_supports(D,G,'D5')
 keep=D[:-1];b=[F(i==29) for i in keep];coeff=mv(inverse([[G[i][j] for j in keep] for i in keep]),b)
 eq=[sum(k*W[i][d] for k,i in zip(coeff,keep,strict=True)) for d in range(22)]
 assert [pair(eq,W[i],Q) for i in D]==[0,0,0,0,1,-1] and pair(eq,eq,Q)==1
 # Remaining pure-A1 plane: all |t|<1/2 meet the strict interior.
 e0=[F(x,2) for x in tangent]
 assert pair(e0,e0,Q)==1 and pair(e0,u,Q)==pair(e0,h,Q)==0
 assert [e0[d]+F(u[d],2) for d in range(22)]==[F(x,2) for x in W[162]]
 assert [e0[d]-F(u[d],2) for d in range(22)]==[-F(x,2) for x in W[164]]
 interior=[]
 for k in[164,162]:
  _,cert=face(h,[k],W,G,Q);interior.append(cert)
 i,j=45,63
 assert G[i][i]==G[j][j]==2 and G[i][j]==0
 assert [pair(e0,W[k],Q) for k in[i,j]]==[F(-1,2),F(1,2)]
 assert [up[k] for k in[i,j]]==[-1,-1]
 zL,L=face(h,[164,i,j],W,G,Q);zR,R=face(h,[162,i,j],W,G,Q)
 A=-pair(u,zL,Q);B=-pair(u,zR,Q)
 assert A>0 and B>0 and pair(zL,zR,Q)<0
 assert set(L['walls'])&set(R['walls'])=={i,j}
 assert pair(e0,zL,Q)==F(1,2)*pair(u,zL,Q)
 assert pair(e0,zR,Q)==F(-1,2)*pair(u,zR,Q)
 # Projected pairing1/4-t²>0 and determinant3-4t²>0 on |t|<1/2.
 # No finite cosine approximation or label cutoff is used.
 return enc({'status':'PASS_CUSP0_A1_ABSORPTION_OBSTRUCTION',
  'sources':{p:sha(p) for p in[RAW,WALK,f'scripts/{TAG}.py']},
  'cusp_id':0,'cusp':u,'affine_factors':['D5','A7','A7','A1'],'active_walls':ids,
  'only_rank_one_component_in_18_representatives':only_A1,
  'A7_opposite_supports':supports,'A7_unique_update_recipes':updates,
  'A7_width_squared':width_A_squared,'A1_available_width_squared':width_B_squared,
  'width_squared_ratio':F(49,8),'D5_star_counts':dc,'D5_star_supports':ds,
  'D5_center_pair_equality':{'walls':keep,'coefficients':coeff,'vector':eq,'square':F(1)},
  'pure_A1_plane':{'base_normal':e0,'base_square':F(1),'direction':u,
   'interior_parameter_interval':'-1/2<t<1/2','interior_endpoint_faces':interior,
   'actual_ridge':[i,j],'left_face':L,'right_face':R,
   'interpolation':'(1/2-t)*B*zL+(t+1/2)*A*zR','A':A,'B':B,
   'projected_pairing_coefficients':[F(1,4),0,-1],
   'projected_Gram_determinant_coefficients':[3,0,-4]},
  'scope':'Every interior section through this literal cusp in which both A1 inequalities are redundant in its Euclidean link. The combinatorial star/width proof is in the report. Not all cusp sections, not all sections retaining one or both A1 link facets, and not a general H20 nonexistence theorem.'})

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
 result=run()
 if args.verify:assert json.loads(OUT.read_text())==result
 else:OUT.write_text(json.dumps(result,indent=2)+'\n')
 print(result['status'],len(result['A7_unique_update_recipes']),result['D5_star_counts'])
