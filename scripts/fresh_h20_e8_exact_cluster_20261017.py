#!/usr/bin/env python3
"""Exact algebraic four-wall E8 recipe with outward rational sign certificates.

The exact vectors are defined by projection, inversion and positive square
roots, not by interval endpoints. Every prescribed equality follows from that
construction; interval arithmetic proves its nondegeneracy and strict gates.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
import fresh_h20_s61_multicap_20261014_exact_review as R
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e8_exact_cluster_20261017'
GEOM=f'results/{TAG}_geometry.json';RECIPE='results/fresh_h20_e8_triangle_cluster_20261016_recipe.json'
CUSP='results/fresh_h20_e8_cusp_completion_20261017.json'
BASE='results/fresh_h20_compact_quad_cap_20261016.json'
# Rational open windows select the branches of the frozen numerical recipe.
# They add no equality constraint and are verified with outward bounds.
BRANCH_WINDOWS={159:(145,F(-11989,1000),F(-11988,1000)),160:(145,F(-10659,1000),F(-10658,1000)),161:(145,F(-15067,1000),F(-15066,1000))}
pt,add,sub,neg,mul,div,sqrt=R.pt,R.plus,R.minus,R.neg,R.times,R.div,R.root
def save(x):
 d=1<<160;return[str(F(R.floor(x[0]*d),d)),str(F(R.ceil(x[1]*d),d))]
def vec(v):return[save(x)for x in v]
def total(xs):
 z=pt(0)
 for x in xs:z=add(z,x)
 return z
def inner(a,b):return total(mul(x,y)for x,y in zip(a,b))
def mv(M,v):return[inner(row,v)for row in M]
def inv(M):
 n=len(M);A=[row[:]+[pt(i==j)for j in range(n)]for i,row in enumerate(M)];piv=[]
 for j in range(n):
  k=next(i for i in range(j,n)if A[i][j][0]>0 or A[i][j][1]<0)
  A[j],A[k]=A[k],A[j];d=A[j][j];piv.append(save(d));A[j]=[div(x,d)for x in A[j]]
  for i in range(n):
   if i==j:continue
   d=A[i][j];A[i]=[sub(x,mul(d,y))for x,y in zip(A[i],A[j])]
 assert all(R.overlap(A[i][j],pt(i==j))for i in range(n)for j in range(n))
 return[row[n:]for row in A],piv
def alg(x):
 x=sp.sympify(x)
 if x.is_Rational:return pt(F(int(x.p),int(x.q)))
 if x.is_Add:return total(alg(t)for t in x.args)
 if x.is_Mul:
  a=pt(1)
  for t in x.args:a=mul(a,alg(t))
  return a
 if x.is_Pow:
  b,e=x.args
  if e in(sp.Rational(1,2),sp.Rational(-1,2)):
   v=sqrt(alg(b));return v if e>0 else div(pt(1),v)
  if e.is_Integer:
   v=pt(1)
   for _ in range(abs(int(e))):v=mul(v,alg(b))
   return div(pt(1),v)if e<0 else v
 raise AssertionError(x)
def cycle_deletions():
 nodes=list(range(1,11))+[159,160,161]
 edges=[(i,i+1,F(1))for i in range(1,10)]+[(1,159,F(1)),(10,160,F(1)),(159,161,F(2)),(160,161,F(1))]
 out=[]
 for omit in nodes:
  d={i:F(2)for i in nodes if i!=omit};e={(i,j):q for i,j,q in edges if omit not in(i,j)};piv=[];ok=True
  while d:
   i=min(i for i in d if sum(i in p for p in e)<=1);v=d.pop(i);piv.append([i,str(v)])
   if v<=0:ok=False;break
   p=next((p for p in e if i in p),None)
   if p is not None:
    j=next(j for j in p if j!=i);d[j]-=e.pop(p)/v
  out.append(dict(omitted=omit,PD=ok,pivots=piv))
 assert[r['omitted']for r in out if r['PD']]==[1,159,160,161]
 return out
def run():
 D=json.loads((ROOT/GEOM).read_text());recipe=json.loads((ROOT/RECIPE).read_text());cusp=json.loads((ROOT/CUSP).read_text());base=json.loads((ROOT/BASE).read_text())
 assert base['known_finite_vertices_at_least']==184
 assert D['status']=='PASS_REPLAYED_INDEPENDENT_158_PHYSICAL_GEOMETRY_EXPORT'
 for doc in(D,recipe,cusp):
  for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 Q=[list(map(F,row))for row in D['form']];qr=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 N=[[tuple(map(F,x))for x in v]for v in D['normals']];qs=list(map(F,D['norms']));Hold=[tuple(map(F,x))for x in D['center']]
 assert len(N)==158
 def pair(v,w):return total(mul(pt(q),mul(v[i],w[j]))for i,row in enumerate(qr)for j,q in row)
 def lin(vs,cs):return[total(mul(c,v[i])for v,c in zip(vs,cs))for i in range(21)]
 def plus(v,w):return[add(a,b)for a,b in zip(v,w)]
 def scale(v,c):return[mul(c,x)for x in v]
 T=list(range(1,11))+[12,13,14,15,61,87,151];TG=sp.eye(17)*2
 for path in[list(range(1,11)),[12,13,14,15],[61,87]]:
  for i,j in zip(path,path[1:]):TG[T.index(i),T.index(j)]=TG[T.index(j),T.index(i)]=-1
 TJ=TG.inv();J=[[pt(F(x))for x in row]for row in TJ.tolist()]
 for i,a in enumerate(T):
  for j,b in enumerate(T):assert R.overlap(pair(N[a],N[b]),pt(F(TG[i,j])))
 def projectT(v):return plus(v,scale(lin([N[i]for i in T],mv(J,[pair(N[i],v)for i in T])),pt(-1)))
 HT=projectT(Hold);oldres=[projectT(v)for v in N]
 # Work in a well-conditioned four-root residual basis. Iterating naive
 # physical-coordinate intervals loses correlations through large old
 # coordinates; this is a basis change, not a change of exact construction.
 basis_ids=[95,142,145,31];B=[oldres[j]for j in basis_ids]
 M=[[pair(v,w)for w in B]for v in B];MI,_=inv(M)
 def coordinates(v):return mv(MI,[pair(b,v)for b in B])
 OC=[coordinates(v)for v in oldres]
 for k,j in enumerate(basis_ids):OC[j]=[pt(i==k)for i in range(4)]
 HC=coordinates(HT)
 def rp(v,w):return inner(v,mv(M,w))
 def rl(vs,cs):return[total(mul(c,v[i])for v,c in zip(vs,cs))for i in range(4)]
 oldT=[[pair(N[j],N[i])for j in T]for i in range(158)]
 newres=[];newprojections=[];newT=[];steps=[]
 for index,rec in enumerate(recipe['steps']):
  targetT=[pt(0)for _ in T]
  for j,m in rec['support_contacts']:assert m==3;targetT[T.index(j)]=pt(-1)
  coeff=mv(J,targetT);fixed=lin([N[i]for i in T],coeff);fixedq=inner(targetT,coeff)
  frame=[];f=[];definitions=[]
  for c in rec['free_and_previous_constraints']:
   j=c['index'];target=alg(c['full_pair'])
   if j<158:frame.append(OC[j])
   else:frame.append(newres[j-158])
   f.append(sub(target,inner(coeff,oldT[j]if j<158 else newT[j-158])));definitions.append(c)
  assert len(frame)==3
  FG=[[rp(v,w)for w in frame]for v in frame]
  # Reuse exact identities already imposed by earlier definitions instead
  # of interval-evaluating the same large cancelling scalar products again.
  frame_ids=[c['index']for c in definitions]
  for a,i in enumerate(frame_ids):
   for b,j in enumerate(frame_ids):
    if max(i,j)<158:continue
    if i==j:full=pt(2)
    elif min(i,j)>=158:full=alg(recipe['target_gram'][i-158][j-158])
    else:
     ni,oi=max(i,j),min(i,j)
     previous=next(c for c in recipe['steps'][ni-158]['free_and_previous_constraints']if c['index']==oi)
     full=alg(previous['full_pair'])
    ti=oldT[i]if i<158 else newT[i-158];tj=oldT[j]if j<158 else newT[j-158]
    exact=sub(full,inner(ti,mv(J,tj)));assert R.overlap(FG[a][b],exact);FG[a][b]=exact
  JI,piv=inv(FG);p=rl(frame,mv(JI,f));pq=inner(f,mv(JI,f))
  chosen=None
  for label,seed in[(f'projected_old_{j}',OC[j])for j in(31,150,145,142,95)]+[('projected_old_center',HC)]:
   hs=[rp(v,seed)for v in frame]
   k=plus(seed,scale(rl(frame,mv(JI,hs)),pt(-1)));kq=sub(rp(seed,seed),inner(hs,mv(JI,hs)))
   if kq[0]>0 or kq[1]<0:chosen=(label,k,kq);break
  assert chosen,(index,list(map(float,kq)),piv)
  label,k,kq=chosen;rad=div(sub(sub(pt(2),fixedq),pq),kq);assert rad[0]>0
  print('FRAME',index,label,'kq',list(map(float,kq)),'rad',list(map(float,rad)),flush=True)
  branches=[];accepted=[]
  for sign in(-1,1):
   y=plus(p,scale(k,mul(pt(sign),sqrt(rad))));r=plus(fixed,lin(B,y));pairs=[pair(r,n)for n in N]
   bad=[j for j in range(158)if j not in T and not(pairs[j][1]<0 and pairs[j][1]**2>2*qs[j])]
   window=BRANCH_WINDOWS.get(158+index)
   match=window is None or(window[1]<pairs[window[0]][0]and pairs[window[0]][1]<window[2])
   branches.append(dict(sign=sign,failed_old_dotted=bad,old_center_pair=save(pair(r,Hold)),branch_window_match=match,window_pair=None if window is None else save(pairs[window[0]])))
   if not bad and match:accepted.append((sign,y,r,pairs))
  assert len(accepted)==1,(index,branches)
  sign,y,r,pairs=accepted[0]
  assert R.overlap(pair(r,r),pt(2))
  for j,t in zip(T,targetT):assert R.overlap(pairs[j],t)
  for c in definitions:assert R.overlap(pairs[c['index']],alg(c['full_pair']))
  contacts=[]
  prescribed={j:'pi/3'for j,m in rec['support_contacts']}
  for j,g in enumerate(pairs):
   if j in T:rel=prescribed.get(j,'right')
   elif j<158:rel='dotted'
   else:
    target=next(c['full_pair']for c in definitions if c['index']==j)
    rel={'0':'right','-1':'pi/3','-sqrt(2)':'pi/4'}[target]
   contacts.append(dict(wall=j,relation=rel,pair=save(g)))
  N.append(r);qs.append(F(2));newres.append(y);newprojections.append(fixed);newT.append(targetT)
  steps.append(dict(root=158+index,T_target_products=vec(targetT),T_projection_square=save(fixedq),linear_constraints=definitions,residual_frame_gram=[vec(row)for row in FG],inverse_pivots=piv,kernel_seed=label,kernel_square=save(kq),particular_square=save(pq),positive_radical_square=save(rad),selected_sign=sign,branch_window=None if 158+index not in BRANCH_WINDOWS else list(map(str,BRANCH_WINDOWS[158+index])),branches=branches,physical_normal=vec(r),all_prior_contacts=contacts,contact_counts=dict(Counter(c['relation']for c in contacts))))
  print('PASS E8 exact wall',158+index,'old center',float(pair(r,Hold)[1]),flush=True)
 # Usually the inherited strict center works; this is certified rather than
 # inferred from pairwise acuteness.
 H=Hold;hp=[pair(H,n)for n in N];assert pair(H,H)[1]<0 and all(g[1]<0 for g in hp)
 facets=[sub(pair(H,H),div(mul(g,g),pt(q)))for g,q in zip(hp,qs)];assert all(x[1]<0 for x in facets)
 # Every full finite support containing a new root lies in the literal core,
 # because all new/old-outside-T products are dotted. E8 deletion cannot
 # repair the indefinite13-cycle; its thirteen exact path deletions were
 # independently tabulated by the frozen rational recipe.
 ids=cusp['core_ids'];assert ids==T+[158,159,160,161]
 deletions=cycle_deletions();CG=[[alg(x)for x in row]for row in cusp['core_gram']]
 for i,a in enumerate(ids):
  for j,b in enumerate(ids):assert R.overlap(pair(N[a],N[b]),CG[i][j])
 corners=[]
 for active in recipe['new_finite_corner_supports']:
  inds=[ids.index(j)for j in active];SG=[[CG[i][j]for j in inds]for i in inds];SI,piv=inv(SG)
  # Exact finite classification supplies PD; certify positive leading
  # Schur pivots in this literal ordering too.
  W=[row[:]for row in SG];ldl=[]
  for k in range(20):
   v=W[k][k];assert v[0]>0;ldl.append(save(v))
   for i in range(k+1,20):
    for j in range(i,20):W[i][j]=W[j][i]=sub(W[i][j],div(mul(W[i][k],W[k][j]),v))
  v=plus(H,scale(lin([N[j]for j in active],mv(SI,[pair(H,N[j])for j in active])),pt(-1)))
  vp=[pair(v,n)for n in N];vq=pair(v,v);assert vq[1]<0 and pair(v,H)[1]<0
  for j,g in enumerate(vp):
   if j in active:assert R.overlap(g,pt(0))
   else:assert g[1]<0,(active,j,save(g))
  corners.append(dict(active=active,positive_pivots=ldl,point=vec(v),norm=save(vq),all162_pairs=vec(vp)))
 # Exact defective cusp gives a canonical strict spacelike escape. The
 # frozen conditional formula expresses both nu and its transverse z in
 # the literal core; it becomes unconditional for this certified branch.
 cusps=[];escape=None
 for ci,c in enumerate(cusp['cusps']):
  nu=lin([N[j]for j in ids],[alg(x)for x in c['null_coefficients']]);z=lin([N[j]for j in ids],[alg(x)for x in c['normalized_free_coefficients']])
  vals=[pair(nu,n)for n in N];incident=c['incident_core_ids']
  assert R.overlap(pair(nu,nu),pt(0))and pair(nu,H)[1]<0
  for j,g in enumerate(vals):
   if j in incident:assert R.overlap(g,pt(0))
   else:assert g[1]<0,(ci,j,save(g))
  cusps.append(dict(index=ci,incident=incident,point=vec(nu),all162_pairs=vec(vals),Euclidean_rank=18))
  if ci==0:
   assert R.overlap(pair(z,z),pt(F(8,7)))and R.overlap(pair(nu,z),pt(0))
   for exponent in range(2,80):
    t=F(1,2**exponent);E=plus(plus(nu,scale(z,pt(t))),scale(H,pt(t**3)))
    eq=pair(E,E);eh=pair(E,H);es=[pair(E,n)for n in N]
    if eq[0]>0 and eh[1]<0 and all(g[1]<0 for g in es):
     escape=dict(definition='nu0+t*z0+t^3*H',t=str(t),point=vec(E),norm=save(eq),center_pair=save(eh),all162_pairs=vec(es));break
   assert escape
 oldcusps=[]
 for raw in D['old_cusps']:
  v=[tuple(map(F,x))for x in raw];p=[pair(v,N[j])for j in range(158,162)];assert all(g[1]<0 for g in p);oldcusps.append(vec(p))
 paths=[f'scripts/{TAG}.py',GEOM,RECIPE,CUSP,BASE,'scripts/fresh_h20_s61_multicap_20261014_exact_review.py']
 return dict(status='PASS_EXACT_162_WALL_E8_TRIANGLE_CLUSTER',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},exact_definition='Sequential orthogonal T projection, exact3x3 Gram inversion in the residual basis[95,142,145,31], one-dimensional kernel projection and its positive square root; branch selected by strict old dotted gates and explicit rational windows. Intervals enclose these algebraic operations, not substitute rounded vectors.',steps=steps,actual_facets=162,strict_center=vec(H),center_norm=save(pair(H,H)),all_center_pairs=vec(hp),facet_point_norms=vec(facets),new_finite_vertices=corners,complete_new_finite_count=4,cycle_deletions=deletions,finite_support_classification='All new/old-outside-T pairs are dotted. In the21-root core every E8 deletion retains the indefinite13-cycle; the13 exact cycle deletions have exactly4PD cases, precisely the four saved supports.',old_finite_vertices_preserved=True,old_finite_persistence_proof='At an old finite vertex v, expand a new wall in its20 spherical simple normals plus v. The spherical inverse is entrywise nonnegative and all root products nonpositive. A nonpositive v coefficient would pair nonnegatively with the common strict center, impossible. Hence every new wall pairs strictly negatively with every old finite vertex.',known_finite_vertices_at_least=base['known_finite_vertices_at_least']+4,five_old_cusp_new_pairs=oldcusps,defective_cusps=cusps,strict_escape=escape,compact_new_triangle_support=T+[158],compact_new_triangle_sides=[159,160,161],scope='Exact coherent162-wall Coxeter chamber, all actual facets, all4new finite vertices, compact new triangular two-face, five retained proper cusps and two defective cusps, and an explicit strict spacelike escape. No finite-volume or complete ambient vertex census claim. Conditional cusp completion is now applicable but its six added walls are not constructed here.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
