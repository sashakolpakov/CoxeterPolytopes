#!/usr/bin/env python3
"""Exact conditional two-cusp completion for the proposed E8 triangle branch.

Only the literal21-root Coxeter core is used for exact algebra. Compatibility
of the ambient162-wall recipe remains a hypothesis, not an output of this
certificate. No numerical root coordinates or producer arithmetic imported.
"""
import argparse,hashlib,json
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_cusp_completion_20261017'
RECIPE='results/fresh_h20_e8_triangle_cluster_20261016_recipe.json'
T=list(range(1,11))+[12,13,14,15,61,87,151];IDS=T+[158,159,160,161]
E8=[158,12,13,14,15,61,87,151]
PATHS=[[2,1,159,161,160],[1,159,161,160,10]]
A6S=[list(range(4,10)),list(range(3,9))]
def mat(A):return [[str(s.simplify(x))for x in row]for row in A.tolist()]
def vec(v):return [str(s.simplify(x))for x in v]
def simp(v):return v.applyfunc(s.simplify)
def run():
 G=s.eye(21)*2
 edges=[(i,i+1,3)for i in range(1,10)]+[(12,13,3),(13,14,3),(14,15,3),(61,87,3),(158,12,3),(158,61,3),(158,151,3),(159,1,3),(160,10,3),(159,161,4),(160,161,3)]
 for i,j,m in edges:G[IDS.index(i),IDS.index(j)]=G[IDS.index(j),IDS.index(i)]=-2*s.cos(s.pi/m)
 determinant=s.factor(G.det());assert determinant.is_negative
 good=[i for i in range(21)if IDS[i]!=161];_,dd=G.extract(good,good).LDLdecomposition(hermitian=False);assert all(x>0 for x in dd.diagonal())
 def emb(ids):return s.eye(21)[:,[IDS.index(i)for i in ids]]
 def dot(u,v):return s.simplify((u.T*G*v)[0])
 EB=emb(E8);EG=EB.T*G*EB;assert EG.det()==1
 theta=EG.inv()[:,E8.index(15)];assert all(x.is_Integer and x>0 for x in theta)and (theta.T*EG*theta)[0]==2
 # A convenient auxiliary timelike vector is the sum of the distinct
 # positive affine null rays. It need not itself be a strict chamber center.
 initial_nus=[]
 for path in PATHS:
  FB=emb(path);kernel=(FB.T*G*FB).nullspace()[0]
  if kernel[0]<0:kernel=-kernel
  assert all(x>0 for x in kernel)
  initial_nus.append(simp(FB*kernel))
 H=initial_nus[0]+initial_nus[1]
 assert dot(initial_nus[0],initial_nus[1])<0 and dot(H,H)<0
 cusps=[];nus=[];bases=[]
 for index,(path,a6)in enumerate(zip(PATHS,A6S)):
  FB=emb(path);FG=FB.T*G*FB;kernel=FG.nullspace()[0]
  if kernel[0]<0:kernel=-kernel
  assert all(x>0 for x in kernel)and FG.rank()==4
  nu=simp(FB*kernel);nus.append(nu);assert dot(nu,nu)==0 and dot(H,nu)<0
  pairs=simp(G*nu);incident=[IDS[i]for i,x in enumerate(pairs)if x==0]
  assert sorted(incident)==sorted(path+E8+a6)
  assert all(x<=0 for x in pairs)
  AB=emb(a6);AG=AB.T*G*AB;AJ=AG.inv();first=AB*AJ[:,0];last=AB*AJ[:,-1]
  assert dot(first,first)==s.Rational(6,7)and dot(last,last)==s.Rational(6,7)and dot(first,last)==s.Rational(1,7)
  assert first+last==AB*s.ones(6,1)
  # Ground the affine block at its last node. The resulting18 roots are PD.
  positive=path[:-1]+E8+a6;PB=emb(positive);PG=PB.T*G*PB;PJ=PG.inv()
  _,dd=PG.LDLdecomposition(hermitian=False);assert all(s.simplify(x)>0 for x in dd.diagonal())
  Hp=simp(H-PB*PJ*PB.T*G*H);hn=dot(Hp,Hp);an=dot(Hp,nu);assert hn<0 and an<0
  chosen=None
  for seed in IDS:
   V=emb([seed]);Vp=simp(V-PB*PJ*PB.T*G*V)
   # Orthogonal projection off the nondegenerate span(nu,Hp).
   a=dot(Vp,nu);b=dot(Vp,Hp)
   z=simp(Vp-((-hn*a/an**2+b/an)*nu+(a/an)*Hp))
   zz=dot(z,z)
   if zz==0:continue
   assert zz.is_positive
   assert simp(PB.T*G*z)==s.zeros(18,1)and dot(z,nu)==0 and dot(z,Hp)==0
   z=simp(s.sqrt(s.Rational(8,7)/zz)*z);assert dot(z,z)==s.Rational(8,7)
   chosen=(seed,zz,z);break
  assert chosen is not None
  seed,zz,z=chosen
  bE=-EB*theta;b1=simp(-first+z);b2=simp(-last-z);B=[bE,b1,b2];bases.append(B)
  BG=s.Matrix([[dot(u,v)for v in B]for u in B]);assert BG==s.Matrix([[2,0,0],[0,2,-1],[0,-1,2]])
  assert all(dot(b,nu)==0 for b in B)
  # Check every incident contact and the literal affine component ranks.
  allvec=FB.row_join(EB).row_join(AB).row_join(s.Matrix.hstack(*B))
  big=simp(allvec.T*G*allvec)
  assert all(x<=0 for i in range(big.rows)for j,x in enumerate(big.row(i))if i!=j)
  eaff=EG.row_join(EB.T*G*bE).col_join((bE.T*G*EB).row_join(s.Matrix([[2]])))
  aaff=AG.row_join(AB.T*G*s.Matrix.hstack(b1,b2)).col_join((s.Matrix.hstack(b1,b2).T*G*AB).row_join(s.Matrix([[2,-1],[-1,2]])))
  assert eaff.rank()==8 and aaff.rank()==7
  assert simp(EB*theta+bE)==s.zeros(21,1)
  assert simp(AB*s.ones(6,1)+b1+b2)==s.zeros(21,1)
  assert sorted([dot(bE,emb([i]))for i in E8])==[-1]+[0]*7
  assert simp(AB.T*G*b1)==-s.eye(6)[:,0]and simp(AB.T*G*b2)==-s.eye(6)[:,-1]
  # At height1 each of the two new affine positive relations is a positive
  # multiple of the same existing null ray; all mutual products stay fixed.
  raised=[b+nu for b in B]
  assert simp(EB*theta+raised[0]-nu)==s.zeros(21,1)
  assert simp(AB*s.ones(6,1)+raised[1]+raised[2]-2*nu)==s.zeros(21,1)
  cusps.append(dict(index=index,affine_F4_path=path,affine_F4_gram=mat(FG),positive_kernel=vec(kernel),null_coefficients=vec(nu),core_pairings=vec(pairs),incident_core_ids=incident,incident_types=['affineF4','E8','A6'],incident_Euclidean_rank=18,E8_ids=E8,A6_ids=a6,E8_highest_root=vec(theta),E8_highest_root_square='2',A6_endpoint_weight_squares=['6/7','6/7'],A6_endpoint_weight_product='1/7',grounded_positive_ids=positive,free_direction_seed=seed,unnormalized_free_square=str(zz),normalized_free_coefficients=vec(z),normalized_free_square='8/7',base_new_root_coefficients=[vec(b)for b in B],base_new_Gram=mat(BG),completed_component_ranks=[4,8,7],completed_quotient_rank=19,completed_incident_roots=22,positive_affine_relations=['original F4 kernel gives nu','theta_E8 + newE8 = hE*nu','sum_A6 + newA + newB = (hA+hB)*nu']))
 cross=dot(nus[0],nus[1]);assert cross<0
 crossrows=[]
 for i,b in enumerate(bases[0]):
  for j,c in enumerate(bases[1]):crossrows.append(dict(pair=[i,j],quadratic=str(cross),linear=str(s.simplify(dot(b,nus[1])+dot(nus[0],c))),constant=str(dot(b,c))))
 # Exact blocker establishing that the easy large-height walls miss T.
 blockers=[]
 for ci,nu in enumerate(nus):
  ids=[j for j in T if dot(nu,emb([j]))<0]
  assert ids==([3,10]if ci==0 else[2,9])
  blockers.append(dict(cusp=ci,nonincident_T_roots=ids,reason='Sufficient large-height inequalities make every new root dotted to both these retained T roots; hence it cannot meet the actual T-face.'))
 src=f'scripts/{TAG}.py'
 return dict(status='PASS_EXACT_CONDITIONAL_TWO_CUSP_COMPLETION_FORMULAS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in [src,RECIPE]},core_ids=IDS,core_gram=mat(G),core_determinant=str(determinant),core_inertia=[20,1,0],cusps=cusps,distinct_null_pairing=str(cross),cross_cusp_pair_polynomials=crossrows,T_face_blockers=blockers,sufficient_height_inequalities=dict(old_nonincident='h > max_j (b.N_j + sqrt(2*N_j^2) + 1)/(-nu.N_j)',common_future_center='h > (b.H + 1)/(-nu.H)',cross_cusp_common_height='Writing pair=-a*h^2+b*h+c with a>0, h>=1+(abs(b)+abs(c+3))/a implies pair<-3.',positive_relations='All three heights at each cusp are strictly positive.'),hypotheses=['The ambient162-wall recipe is exactly promoted to a coherent acute Coxeter chamber with these literal core products.','Each of159,160,161 is strictly dotted to every retained old wall outside T, as numerically checked but not exactly certified by the inherited recipe.','A common future timelike strict interior point H exists.'],conclusion='Conditionally, six further norm2 walls simultaneously complete the two cusp links to affineF4+affineE8+affineA7. Sufficiently large integer heights give strict dotted old nonincident and cross-cusp contacts and preserve coherence. This existence theorem does not certify an actual168-wall numerical/exact realization.',scope='Exact core algebra and conditional finite-height construction only. The inherited162 branch remains numerical. No finite-volume claim; the simple sufficiently-large-height completion misses T and leaves its unfinished H3 region unchanged.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['distinct_null_pairing'],flush=True)
