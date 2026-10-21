#!/usr/bin/env python3
"""Two exact nonright-S6 escape cutters; native wall norms are retained."""
import argparse,itertools,json
from collections import Counter
from functools import lru_cache
from hashlib import sha256
import sympy as sy
import fresh_nonarithmetic_escape52_20261016_cap as C53
import fresh_h9_persistent_cusp_affine_20261016 as AF
U=C53.U;B,D,A,R,C,K,Z,O=U.B,U.D,U.A,U.R,U.C,U.K,U.Z,U.O
I=U.I;ROOT=C53.ROOT;TAG='fresh_nonarithmetic_nonright_s6_20261016_caps'
OLD='results/fresh_nonarithmetic_escape52_20261016_cap.json'
CONE='results/fresh_nonarithmetic_escape52_20261016_cone.json'
CL=f'results/{U.TAG}.json'
S=[0,1,2,3,6,7,8,10]
RECIPES=[dict(scout_index=844,anchor=34,labels={6:5,8:3,34:5}),dict(scout_index=877,anchor=12,labels={3:3,10:10,12:4})]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 old=B.read(U.CAP);cl=B.read(CL);c53=B.read(OLD);co=B.read(CONE)
 Q=[[A.B.db(x)for x in row]for row in old['form']]
 N0=[tuple(U.dec(x)for x in row)for row in old['normals']]
 H0=tuple(U.dec(x)for x in old['center']);E0=tuple(U.dec(x)for x in co['strict_escape']['point'])
 QI=[[U.iv(D.lift(x))for x in row]for row in Q]
 NI=[[U.iv(x)for x in n]for n in N0]+[[I(*x)for x in r['normal_interval']]for r in cl['caps']]+[[I(*x)for x in c53['normal_interval']]]
 norms=[U.iv(D.pair(n,n,Q))for n in N0]+[I(2)]*4
 H=[U.iv(x)for x in H0];E=[U.iv(x)for x in E0];hc=U.imv(QI,H);ec=U.imv(QI,E)
 cov=[U.imv(QI,n)for n in NI];hh=U.idot(H,hc)
 inv=B.read('results/fresh_h9_persistent_cusp_screen_20261016.json')
 L=[row+[0]*6 for row in inv['pair_matrix']]+[[0]*55 for _ in range(6)]
 for cap in cl['caps']+[dict(index=52,relations=c53['relations'])]:
  j=cap['index'];L[j][j]=2
  for rel in cap['relations']:L[j][rel['wall']]=L[rel['wall']][j]=rel.get('label',0)
 caps=[]
 for j,recipe in enumerate(RECIPES,53):
  frame=S+[recipe['anchor']];labels=recipe['labels']
  G=[[U.iv(D.pair(N0[i],N0[k],Q))for k in frame]for i in frame]
  J,pivots=C53.inverse(G)
  # A norm-two new root has native product -sqrt(q_old/2)*2cos(pi/m).
  contact=[-U.sqrt(norms[i]/2)*U.theta(labels.get(i,2))for i in frame]
  pc=U.imv(J,contact);p=C53.combine(pc,[NI[i]for i in frame])
  hp=[U.idot(H,cov[i])for i in frame];hcoeff=U.imv(J,hp)
  k=[U.rnd(x-y)for x,y in zip(H,C53.combine(hcoeff,[NI[i]for i in frame]))]
  pp=U.idot(contact,pc);qk=U.rnd(hh-U.idot(hp,hcoeff));delta=U.rnd((2-pp)/qk);lam=U.sqrt(delta)
  branches=[];chosen=None
  for sign in[-1,1]:
   n=[U.rnd(x+sign*lam*y)for x,y in zip(p,k)];he=U.idot(n,hc);ee=U.idot(n,ec)
   acute=all(i in frame or U.idot(n,cov[i]).h<0 for i in range(53))
   branches.append(dict(sign=sign,center_pairing=he.enc(),escape_pairing=ee.enc(),strict_off_frame_acute=acute))
   if acute and he.h<0<ee.l:assert chosen is None;chosen=(sign,n,he,ee)
  assert chosen is not None
  sign,n,he,ee=chosen;relations=[]
  for i in range(j):
   g=U.idot(n,cov[i])
   if i in frame:
    m=labels.get(i,2);relations.append(dict(wall=i,label=m,kind='right'if m==2 else'finite',native_product=contact[frame.index(i)].enc(),pairing_interval=g.enc(),exact_by_recipe=True));L[j][i]=L[i][j]=m
   else:
    gap=U.rnd(g*g-2*norms[i]);assert g.h<0<gap.l,('not strict dotted',j,i,g.enc(),gap.enc())
    relations.append(dict(wall=i,kind='dotted',pairing_interval=g.enc(),squared_product_minus_norm_product=gap.enc()))
  caps.append(dict(index=j,scout_index=recipe['scout_index'],frame=frame,finite_contacts=[[i,labels.get(i,2)]for i in frame],native_contact_vector=[x.enc()for x in contact],frame_gram=[[x.enc()for x in row]for row in G],inverse_pivots=pivots,particular=[x.enc()for x in p],particular_norm=pp.enc(),kernel=[x.enc()for x in k],kernel_norm=qk.enc(),radical_square=delta.enc(),selected_sign=sign,branches=branches,normal_interval=[x.enc()for x in n],norm_exact=2,center_pairing=he.enc(),escape_pairing=ee.enc(),relations=relations))
  NI.append(n);cov.append(U.imv(QI,n));norms.append(I(2));L[j][j]=2
 # The only possible affine circuits containing a new root use standard
 # crystallographic bonds. The two new roots are dotted to each other.
 affcounts=[]
 for cap in caps:
  new=cap['index'];ids=[r['wall']for r in cap['relations']if r.get('label')in AF.VALUES]
  queue=[(new,)];seen=set(queue);counts=Counter()
  while queue:
   ss=queue.pop();M=tuple(tuple(AF.scale(AF.ONE,2)if i==j else AF.VALUES[L[i][j]]for j in ss)for i in ss)
   kind,rank,nullity=AF.classify(M);counts[kind]+=1;assert kind!='positive_semidefinite',('new affine',ss)
   if kind!='positive_definite':continue
   for k in ids:
    if k in ss:continue
    ms=[L[k][i]for i in ss]
    if any(m not in AF.VALUES for m in ms)or all(m==2 for m in ms):continue
    nxt=tuple(sorted(ss+(k,)))
    if nxt not in seen:seen.add(nxt);queue.append(nxt)
  affcounts.append(dict(wall=new,connected_standard_label_subset_counts=dict(counts)))
 # Complete finite-vertex census on the two new facets: their nine finite
 # neighbors are all old walls, and their mutual pair is dotted.
 theta={2:sy.Integer(0),3:sy.Integer(1),4:sy.sqrt(2),5:(1+sy.sqrt(5))/2,6:sy.sqrt(3),10:sy.sqrt((5+sy.sqrt(5))/2)}
 supports=[];vertices=[]
 for cap in caps:
  new=cap['index']
  for ids in itertools.combinations(cap['frame'],8):
   ss=list(ids)+[new];labels=[[L[i][j]if i!=j else 2 for j in ss]for i in ss]
   if any(labels[i][j]==0 for i in range(9)for j in range(i)):
    supports.append(dict(support=ss,kind='contains_dotted_pair'));continue
   M=sy.Matrix([[2 if i==j else-theta[labels[i][j]]for j in range(9)]for i in range(9)])
   leading=[sy.simplify(M[:i,:i].det())for i in range(1,10)];pd=all(z.is_positive is True for z in leading)
   assert pd or any(z.is_nonpositive is True for z in leading),('undecided finite support',ss,leading)
   supports.append(dict(support=ss,leading_minors=[str(z)for z in leading],positive_definite=pd))
   if not pd:continue
   VN=[NI[i]for i in ss]
   GI=[[norms[a]if a==b else-U.sqrt(norms[a]*norms[b])/2*U.theta(L[a][b])for b in ss]for a in ss]
   JI,vp=C53.inverse(GI);hs=[U.idot(v,hc)for v in VN];cs=U.imv(JI,hs)
   v=[U.rnd(x-y)for x,y in zip(H,C53.combine(cs,VN))];qv=U.rnd(hh-U.idot(hs,cs));vh=U.idot(v,hc);assert qv.h<0 and vh.h<0
   nonincident=[]
   for j,c in enumerate(cov):
    if j not in ss:
     g=U.idot(v,c);assert g.h<0,('finite witness off wall',ss,j,g.enc());nonincident.append(dict(wall=j,pairing=g.enc()))
   vertices.append(dict(support=ss,point_interval=[x.enc()for x in v],norm_interval=qv.enc(),future_chart_interval=vh.enc(),inverse_pivots=vp,nonincident_pairings=nonincident))
 paths=[f'scripts/{TAG}.py',f'scripts/{C53.TAG}.py',f'scripts/{U.TAG}.py',f'scripts/{AF.TAG}.py',U.CAP,CL,OLD,CONE,'results/fresh_h9_persistent_cusp_screen_20261016.json','results/fresh_h9_persistent_cusp_affine_20261016.json','results/fresh_nonarithmetic_affine_filtered_20261016_affine_gate.json','results/fresh_nonarithmetic_nonright_s6_20261016_scout.json','results/fresh_nonarithmetic_nonright_s6_20261016_selection.json']
 return dict(status='PASS_EXACT_TWO_NONRIGHT_S6_ESCAPE_CUTTERS',sources={p:sha(p)for p in paths},retained_walls=53,resulting_walls=55,form=old['form'],center=old['center'],cut_escape=co['strict_escape']['point'],caps=caps,affine_checks=affcounts,finite_support_census=supports,new_finite_vertices=vertices,actual_facets=55,
 exact_recipe='For each saved native nine-root frame F, use its exact Gram G and c_i=-sqrt(q_i/2)*2cos(pi/m_i). Define p=NF G^-1c,k=H-NF G^-1(<H,n_i>),r=p+sign*sqrt((2-p²)/k²)k. Certified inverse pivots and positive radical squares establish the exact roots; projection identities establish native contacts and norm2. Every other pairing is interval-certified strictly dotted, including the mutual pair.',
 facet_proof='All55 spacelike acute normals share strict timelikeH. qH-<H,n>n is a future timelike relative-interior point on each facet, strict against every other normal. All old finite vertices persist by the PD Z-Gram argument. Each new facet has exactly nine finite neighbors, so its finite-vertex census is complete by the eighteen listed eight-neighbor supports.',
 no_affine_proof='Retained53 has no affine circuits. All connected standard-label subdiagrams containing each new root are explicitly grown from that root until non-PD. No affine one occurs; labels5/10 cannot occur in irreducible affine diagrams. The new mutual pair is dotted, so no affine diagram can contain both.',
 scope='Exact coherent55-wall local geometry and complete new-facet finite-vertex census. Both roots cut the renewed53 escape. Complete ambient cone/finite volume are separate; no completion or volume claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved two-cap replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'new finite vertices',len(d['new_finite_vertices']))
