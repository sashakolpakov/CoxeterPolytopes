#!/usr/bin/env python3
"""Exact illustrative nonreflective D19-genus cusp via Niemeier E6^4.

This is not a reconstruction of Esselmann's full high-rank classification.
No project arithmetic or lattice tables are imported.
"""
import argparse,hashlib,itertools,json,math
from pathlib import Path
import sympy as s
from sympy.matrices.normalforms import hermite_normal_form,smith_normal_form
from sympy.polys.domains import ZZ
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_esselmann_d19_countercusp_20261017'
def mat(A):return [[str(x)for x in row]for row in A.tolist()]
def vec(A):return [str(x)for x in A]
def integer(A):assert all(x.q==1 for x in A);return A.applyfunc(int)
def lattice_basis(columns,denominator):
 H=hermite_normal_form(integer(columns*denominator))
 return H/denominator
def dgram(n):
 G=2*s.eye(n)
 for i,j in [(0,2),(1,2)]+[(i,i+1)for i in range(2,n-1)]:G[i,j]=G[j,i]=-1
 return G
def order(v):return s.ilcm(*[x.q for x in v])
def disc_generator(G,wanted):
 I=G.inv()
 for i in range(G.rows):
  v=I[:,i]
  if order(v)==4 and I[i,i]%2==wanted:return dict(column=i,coordinates=vec(v),order=4,square=str(I[i,i]),square_mod_2=str(wanted))
 raise AssertionError('No requested discriminant generator')
def run():
 G=2*s.eye(6)
 for i,j in [(0,1),(0,2),(2,3),(0,4),(4,5)]:G[i,j]=G[j,i]=-1
 I=G.inv();om=I[:,3];w=3*om;assert G.det()==3 and (om.T*G*om)[0]==s.Rational(4,3)
 assert all(x>0 for x in I) and min(I.diagonal())==s.Rational(4,3)
 assert all(x.q==1 for x in w) and math.gcd(*map(int,w))==1 and (w.T*G*w)[0]==12
 ids=[0,1,2,4,5];D5=G.extract(ids,ids);assert D5.det()==4 and all((G*w)[i]==0 for i in ids)
 # Full finite Weyl root orbit, checked literally.
 roots={tuple(s.eye(6)[:,i])for i in range(6)};todo=list(roots)
 while todo:
  r=todo.pop();v=s.Matrix(r)
  for i in range(6):
   z=list(r);z[i]-=(G*v)[i];z=tuple(z)
   if z not in roots:roots.add(z);todo.append(z)
 assert len(roots)==72 and all((s.Matrix(r).T*G*s.Matrix(r))[0]==2 for r in roots)
 assert s.Matrix(list(roots)).rank()==6
 # A norm2 vector moved to the dominant chamber has at most one weight
 # coordinate: two nonzero coordinates already cost at least8/3.
 dominant2=[i for i in range(6)if I[i,i]==2 and all(x.q==1 for x in I[:,i])]
 assert dominant2==[1] and tuple(I[:,1])in roots
 # Tetracode: every nonzero word has exactly3 nonzero coordinates.
 gens=[(1,0,1,1),(0,1,1,-1)]
 code=sorted({tuple((a*x+b*y)%3 for x,y in zip(*gens))for a,b in itertools.product(range(3),repeat=2)})
 assert len(code)==9 and all(sum(x!=0 for x in c)==3 for c in code if any(c))
 assert [c for c in code if c[1:]==(0,0,0)]==[(0,0,0,0)]
 Q=s.diag(G,G,G,G);g=[]
 for c in gens:g.append(s.Matrix.vstack(*(x*om for x in c)))
 assert (g[0].T*Q*g[0])[0]==4 and (g[1].T*Q*g[1])[0]==4 and (g[0].T*Q*g[1])[0]==0
 NB=lattice_basis(s.eye(24).row_join(g[0]).row_join(g[1]),3);NG=integer(NB.T*Q*NB)
 assert NG.det()==1 and all(x%2==0 for x in NG.diagonal())
 assert abs(NB.det())==s.Rational(1,9)
 embedding=s.eye(24)[:,ids];coords=integer(NB.inv()*embedding)
 sn=smith_normal_form(coords,domain=ZZ);assert all(abs(sn[i,i])==1 for i in range(5))
 # Complement: <w> + three E6 blocks, with the SAME two glue vectors.
 Q0=s.diag(s.Matrix([[12]]),G,G,G)
 a=s.zeros(19,1);a[0]=s.Rational(1,3);a[7:13,0]=om;a[13:19,0]=om
 b=s.zeros(19,1);b[1:7,0]=om;b[7:13,0]=om;b[13:19,0]=-om
 MB=lattice_basis(s.eye(19).row_join(a).row_join(b),3);MG=integer(MB.T*Q0*MB)
 assert Q0.det()==324 and abs(MB.det())==s.Rational(1,9) and MG.det()==4
 assert all(x%2==0 for x in MG.diagonal())
 _,LD=MG.LDLdecomposition(hermitian=False);assert all(x>0 for x in LD.diagonal())
 # Literal embedding M -> N and orthogonality to the primitive D5.
 E=s.zeros(24,19);E[:6,0]=w;E[6:24,1:19]=s.eye(18)
 assert E.T*Q*E==Q0 and embedding.T*Q*E==s.zeros(5,19)
 integral_M_in_N=integer(NB.inv()*E*MB)
 D19=dgram(19);assert D19.det()==4
 dm=disc_generator(MG,s.Rational(3,4));d5=disc_generator(D5,s.Rational(5,4));d19=disc_generator(D19,s.Rational(3,4))
 assert (-s.Rational(d5['square']))%2==s.Rational(dm['square_mod_2'])
 assert all((4*x).q==1 for x in MG.inv())
 possible=[q for q in range(2,9,2)if 8%q==0];assert possible==[2,4,8]
 bounds=[]
 for q in [4,8]:
  lower=s.Rational(q*q,4)*s.Rational(4,3);assert lower>q
  bounds.append(dict(root_square=q,nonzero_E6_projection_lower_bound=str(lower),strictly_greater_than_root_square=True))
 # The first-only intersection is E6, and its D5 perpendicular is Zw.
 assert math.gcd(*map(int,w))==1
 root_embeddings=[]
 for block in range(3):
  for r in sorted(roots):
   v=s.zeros(19,1);v[1+6*block:7+6*block,0]=s.Matrix(r)
   root_embeddings.append(integer(MB.inv()*v))
 RR=s.Matrix.hstack(*root_embeddings);assert RR.rank()==18 and RR.cols==216
 assert all((r.T*MG*r)[0]==2 for r in root_embeddings)
 # Illustrative contrast ONLY: D4 in the same first E6 leaves two norm4
 # reflection directions, rather than the nonreflecting norm12 D5 axis.
 D4ids=[0,1,2,4];D4=G.extract(D4ids,D4ids);assert D4.det()==4
 x=2*I[:,3]-I[:,5];y=-I[:,3]+2*I[:,5];XY=s.Matrix.hstack(x,y)
 assert all(z.q==1 for z in XY) and XY.T*G*XY==s.Matrix([[4,-2],[-2,4]])
 assert G.extract(D4ids,list(range(6)))*XY==s.zeros(4,2)
 R20=s.diag(XY.T*G*XY,G,G,G)
 a20=s.zeros(20,1);a20[:2,0]=s.Matrix([s.Rational(2,3),s.Rational(1,3)])
 a20[8:14,0]=om;a20[14:20,0]=om
 b20=s.zeros(20,1);b20[2:8,0]=om;b20[8:14,0]=om;b20[14:20,0]=-om
 B20=lattice_basis(s.eye(20).row_join(a20).row_join(b20),3)
 G20=integer(B20.T*R20*B20);assert G20.det()==4 and all(z%2==0 for z in G20.diagonal())
 E20=s.zeros(24,20);E20[:6,:2]=XY;E20[6:24,2:20]=s.eye(18)
 assert E20.T*Q*E20==R20
 integer(NB.inv()*E20*B20)
 D4coords=integer(NB.inv()*s.eye(24)[:,D4ids]);D4smith=smith_normal_form(D4coords,domain=ZZ)
 assert all(abs(D4smith[i,i])==1 for i in range(4))
 reflections20=[]
 for i in range(20):
  r=integer(B20.inv()*s.eye(20)[:,i]);q=(r.T*G20*r)[0]
  reflection=integer(s.eye(20)-s.Rational(2,q)*r*(r.T*G20))
  assert reflection.T*G20*reflection==G20
  reflections20.append(dict(base_index=i,square=str(q),integral_reflection=True))
 S20=smith_normal_form(G20,domain=ZZ)
 assert sorted(abs(S20[i,i])for i in range(20))==[1]*18+[2,2]
 cosets20={tuple(z%1 for z in G20.inv()[:,i])for i in range(20)}
 cosets20|={tuple(s.Rational(0)for _ in range(20))}
 while True:
  more={tuple((a+b)%1 for a,b in zip(u,v))for u in cosets20 for v in cosets20}
  if more<=cosets20:break
  cosets20|=more
 assert len(cosets20)==4
 ds20=sorted((s.Matrix(v).T*G20*s.Matrix(v))[0]%2 for v in cosets20 if any(v))
 assert ds20==[1,1,1]
 source=f'scripts/{TAG}.py'
 return dict(status='PASS_EXACT_ILLUSTRATIVE_D19_GENUS_NONREFLECTIVE_CUSP',sources={source:hashlib.sha256((ROOT/source).read_bytes()).hexdigest()},references=dict(esselmann='https://doi.org/10.1006/jnth.1996.0141',vinberg_genus_lemma='https://www.numdam.org/item/AST_1992__209__279_0.pdf',niemeier_glue_table='https://www.jstage.jst.go.jp/article/jmath/72/2/72_599/_pdf',nikulin_genus='https://www.mathnet.ru/eng/im1677'),E6=dict(gram=mat(G),inverse=mat(I),omega=vec(om),omega_square='4/3',primitive_axis=vec(w),axis_square=12,dual_minimum='4/3',root_count=72,dominant_norm2_weight_index=1,roots=[vec(r)for r in sorted(roots)]),D5=dict(indices=ids,gram=mat(D5),discriminant=d5,primitive_embedding_smith_diagonal=[1]*5),Niemeier=dict(tetracode_generators=[list(x)for x in gens],tetracode_words=[list(x)for x in code],base_gram=mat(Q),glue_vectors=[vec(x)for x in g],basis=mat(NB),gram=mat(NG),index=9,determinant=1,nonzero_glue_coset_minimum=4,first_component_intersection='E6 exactly'),complement=dict(base_gram=mat(Q0),glue_vectors=[vec(a),vec(b)],basis=mat(MB),gram=mat(MG),embedding_in_N_basis=mat(integral_M_in_N),rank=19,index=9,determinant=4,discriminant=dm,D19_discriminant=d19,LDL_pivots=vec(LD.diagonal()),root_count=216,root_span_rank=18,primitive_reflector_possible_squares=possible,long_root_exclusions=bounds,axis_minimum=12),D4_contrast=dict(indices=D4ids,axis_basis=mat(XY),axis_gram=mat(XY.T*G*XY),base_gram=mat(R20),glue_vectors=[vec(a20),vec(b20)],basis=mat(B20),gram=mat(G20),rank=20,determinant=4,discriminant_group='(Z/2)^2',nonzero_discriminant_squares=[str(z)for z in ds20],spanning_reflections=reflections20,scope='This one D4 complement is reflective; no assertion that every lattice in genusD20 is reflective.'),proof_inputs=['Finite Weyl group moves any dual vector to the dominant weight cone; positive inverse entries and diagonal minimum4/3 prove min(E6*)=4/3.','Primitive complement discriminant is opposite; alternatively the explicit19Gram directly gives Z/4 with q=3/4.','Even-lattice genus is determined by signature and discriminant form (Nikulin Corollary1.9.4).','Scharlau–Walhorn Lemma1.3: reflectivity of U+M requires every lattice in genus(M) reflective.'],scope='Exact explicit countercusp for U+D19 through one nonreflective positive lattice in genusD19. NOT the complete Esselmann exclusion of all rank21 Lorentzian lattices; no novelty claim.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['complement']['rank'],d['complement']['root_span_rank'],flush=True)
