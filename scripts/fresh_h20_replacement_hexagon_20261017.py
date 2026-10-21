#!/usr/bin/env python3
"""Exact171-wall alternative: delete31,158; add168..172; compact hexagon.

Exact vectors are defined by algebraic projection and quadratic completion.
Outward rational intervals certify the strict gates and branch choices.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
import sympy as s
import fresh_h20_e8_exact_cluster_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_replacement_hexagon_20261017'
pt,add,sub,neg,mul,div,sqrt=P.pt,P.add,P.sub,P.neg,P.mul,P.div,P.sqrt
save,vec,total,inner,mv,alg,inv=P.save,P.vec,P.total,P.inner,P.mv,P.alg,P.inv
REMOVED={31,158}
WINDOWS={169:(145,F(-5484,1000),F(-5481,1000)),170:(95,F(-2025,1000),F(-2022,1000)),171:(142,F(-5180,1000),F(-5177,1000)),172:(95,F(-3160,1000),F(-3157,1000))}
CENTER=['-578503/1000000','598849/1000000','-185489/1000000','623055/1000000']
ESCAPE=['259264/1000000','-434859/1000000','263921/1000000','635828/1000000']
RECIPES=[dict(root=169,T_contacts={61:3},anchors=[(95,2),(142,5),(152,5)]),dict(root=170,T_contacts={},anchors=[(142,2),(152,2),(169,2)]),dict(root=171,T_contacts={},anchors=[(95,4),(169,3),(161,3)]),dict(root=172,T_contacts={},anchors=[(152,2),(169,2),(171,2)])]

def run():
 G=json.loads((ROOT/P.GEOM).read_text());D=json.loads((ROOT/f'results/{P.TAG}.json').read_text());U=json.loads((ROOT/'results/fresh_h20_e8_f4_cusps_20261017.json').read_text());C=json.loads((ROOT/P.CUSP).read_text());RE=json.loads((ROOT/'results/fresh_h20_prism_replacement_20261017.json').read_text())
 for doc in(G,D,U,C,RE):
  for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 N=[[tuple(map(F,x))for x in v]for v in G['normals']]+[[tuple(map(F,x))for x in z['physical_normal']]for z in D['steps']]+[[tuple(map(F,x))for x in z['physical_normal']]for z in U['new_walls']]
 qs=list(map(F,G['norms']))+[F(2)]*10;Hold=[tuple(map(F,x))for x in D['strict_center']];assert len(N)==168
 Q=[list(map(F,row))for row in G['form']];nz=[(i,j,q)for i,row in enumerate(Q)for j,q in enumerate(row)if q]
 def pair(v,w):return total(mul(pt(q),mul(v[i],w[j]))for i,j,q in nz)
 def plus(v,w):return[add(x,y)for x,y in zip(v,w)]
 def scale(v,c):return[mul(c,x)for x in v]
 def lin(vs,cs):return[total(mul(c,v[i])for v,c in zip(vs,cs))for i in range(len(vs[0]))]
 T=C['core_ids'][:17];TG=s.eye(17)*2
 for path in[list(range(1,11)),[12,13,14,15],[61,87]]:
  for i,j in zip(path,path[1:]):TG[T.index(i),T.index(j)]=TG[T.index(j),T.index(i)]=-1
 TJ=TG.inv();J=[[pt(F(x))for x in row]for row in TJ.tolist()]
 def projectT(v):return plus(v,scale(lin([N[j]for j in T],mv(J,[pair(N[j],v)for j in T])),pt(-1)))
 B=[projectT(N[j])for j in(95,142,145,31)];M=[[pair(v,w)for w in B]for v in B];MI,_=inv(M)
 def coords(v):return mv(MI,[pair(b,v)for b in B])
 OC=[coords(projectT(v))for v in N];oldT=[[pair(N[j],v)for j in T]for v in N]
 for k,j in enumerate((95,142,145,31)):OC[j]=[pt(i==k)for i in range(4)]
 def rp(v,w):return inner(v,mv(M,w))
 def norm_target(label,q):return pt(0)if label==2 else neg(mul(sqrt(pt(2*q)),alg(s.cos(s.pi/label))))
 known={};steps=[];counts=Counter()
 def add_step(index,r,y,tt,prescribed,definition):
  assert index==len(N)
  assert P.R.overlap(pair(r,r),pt(2)),('norm',index,save(pair(r,r)))
  rows=[]
  for j,n in enumerate(N):
   if j in REMOVED:continue
   g=pair(r,n)
   if j in prescribed:
    label=prescribed[j];target=norm_target(label,qs[j]);assert P.R.overlap(g,target),(index,j,save(g),save(target));rel='right'if label==2 else f'pi/{label}';known[index,j]=known[j,index]=target
   else:
    assert g[1]<0 and g[1]**2>2*qs[j],('not dotted',index,j,save(g),str(qs[j]));rel='dotted'
   counts[rel]+=1;rows.append(dict(wall=j,relation=rel,pair=save(g)))
  N.append(r);qs.append(F(2));OC.append(y);oldT.append(tt);known[index,index]=pt(2)
  for j,t in zip(T,tt):known[index,j]=known[j,index]=t
  steps.append(dict(root=index,definition=definition,T_target_products=vec(tt),physical_normal=vec(r),residual_coordinates=vec(y),all_retained_prior_contacts=rows,contact_counts=dict(Counter(z['relation']for z in rows))))
  print('EXACT HEXAGON WALL',index,steps[-1]['contact_counts'],flush=True)
 # Exact branch1 from the already certified core, now used on a different
 # retained set. Its failed wall158 is explicitly absent.
 rc=RE['exact_recipe'];cs=[add(alg(p),mul(alg(rc['amplitude']),alg(u)))for p,u in zip(rc['particular'],rc['axis'])];r=lin([N[j]for j in rc['core_ids']],cs);y=lin([OC[j]for j in rc['core_ids']],cs);tt=[pt(0)]*17
 add_step(168,r,y,tt,{**{j:2 for j in T},159:2,160:3,161:2},dict(type='exact core p+sqrt((2-p^2)/u^2)*u',source='results/fresh_h20_prism_replacement_20261017.json',core_coefficients=vec(cs)))
 for rec in RECIPES:
  i=rec['root'];tt=[pt(-1 if j in rec['T_contacts']else 0)for j in T];tc=mv(J,tt);fixed=lin([N[j]for j in T],tc);fixedq=inner(tt,tc);ids=[j for j,m in rec['anchors']];frame=[OC[j]for j in ids];targets=[norm_target(lab,qs[j])for j,lab in rec['anchors']];f=[sub(t,inner(tc,oldT[j]))for j,t in zip(ids,targets)]
  FG=[[rp(v,w)for w in frame]for v in frame]
  for a,j in enumerate(ids):
   for b,k in enumerate(ids):
    if j==k:full=pt(qs[j])
    elif(j,k)in known:full=known[j,k]
    else:continue
    exact=sub(full,inner(oldT[j],mv(J,oldT[k])));assert P.R.overlap(FG[a][b],exact);FG[a][b]=exact
  FI,piv=inv(FG);p=lin(frame,mv(FI,f));pq=inner(f,mv(FI,f));chosen=None
  for j in(31,145,95,142,168):
   seed=OC[j];hs=[rp(v,seed)for v in frame];k=plus(seed,scale(lin(frame,mv(FI,hs)),pt(-1)));kq=sub(rp(seed,seed),inner(hs,mv(FI,hs)))
   if kq[0]>0 or kq[1]<0:chosen=(j,k,kq);break
  assert chosen;i_seed,k,kq=chosen;rad=div(sub(sub(pt(2),fixedq),pq),kq);assert rad[0]>0,(i,save(rad))
  window=WINDOWS[i];matches=[];br=[]
  for sign in(-1,1):
   yy=plus(p,scale(k,mul(pt(sign),sqrt(rad))));rr=plus(fixed,lin(B,yy));g=pair(rr,N[window[0]]);hit=window[1]<g[0]and g[1]<window[2];br.append(dict(sign=sign,window_pair=save(g),matches=hit))
   if hit:matches.append((sign,yy,rr))
  assert len(matches)==1,(i,br);sign,yy,rr=matches[0];pres={j:rec['T_contacts'].get(j,2)for j in T};pres.update(dict(rec['anchors']))
  definition=dict(type='T projection plus3-anchor Gram inverse and1D quadratic completion',anchors=[dict(wall=j,label=m,target=save(t))for(j,m),t in zip(rec['anchors'],targets)],frame_gram=[vec(v)for v in FG],inverse_pivots=piv,kernel_seed=i_seed,kernel_square=save(kq),particular_square=save(pq),radical_square=save(rad),selected_sign=sign,window=[window[0],str(window[1]),str(window[2])],branches=br)
  add_step(i,rr,yy,tt,pres,definition)
 ret=[j for j in range(173)if j not in REMOVED];assert len(ret)==171 and sum(counts.values())==840
 dual=lin([N[j]for j in T],[total(row)for row in J]);Hface=lin(B,[pt(F(x))for x in CENTER]);assert pair(Hface,Hface)[1]<0 and pair(Hface,Hold)[1]<0
 H=None
 for exponent in range(8,100):
  eps=F(1,2**exponent);v=plus(Hface,scale(dual,pt(-eps)));hp=[pair(v,N[j])for j in ret]
  if pair(v,v)[1]<0 and all(g[1]<0 for g in hp):H=v;break
 assert H is not None and pair(H,Hold)[1]<0;hq=pair(H,H);facets=[sub(hq,div(mul(g,g),pt(qs[j])))for j,g in zip(ret,hp)];assert all(g[1]<0 for g in facets)
 # S is exactly A10+A4+A3+A1. Its inverse is rational despite the new
 # nonright contacts of the six side roots to the seed169.
 S=T+[169];SG=s.eye(18)*2;SG[:17,:17]=TG;SG[17,T.index(61)]=SG[T.index(61),17]=-1;SJ=SG.inv();SI=[[pt(F(x))for x in row]for row in SJ.tolist()]
 st={j:[pair(N[k],N[j])for k in S]for j in ret};projected={j:sub(pt(qs[j]),inner(st[j],mv(SI,st[j])))for j in ret};sides=[];negative=[]
 for j in ret:
  if j in S:assert P.R.overlap(projected[j],pt(0));continue
  if projected[j][0]>0:sides.append(j)
  else:assert projected[j][1]<0,(j,save(projected[j]));negative.append(j)
 assert sides==[95,142,152,170,171,172],sides
 cycle=[95,142,170,152,172,171];corners=[];pair_tests=[]
 for a in range(6):
  for b in range(a):
   j,k=sides[a],sides[b];off=sub(pair(N[j],N[k]),inner(st[j],mv(SI,st[k])));det=sub(mul(projected[j],projected[k]),mul(off,off));adj=abs(cycle.index(j)-cycle.index(k))in(1,5)
   if adj:assert det[0]>0
   else:assert off[1]<0 and det[1]<0,(j,k,save(off),save(det))
   pair_tests.append(dict(pair=[k,j],adjacent=adj,projected_product=save(off),determinant=save(det)))
 for a,j in enumerate(cycle):
  k=cycle[(a+1)%6];active=S+[j,k];W=[[pair(N[x],N[y])for y in active]for x in active]
  # Use the defining spherical identities on the18-root support.
  for x in range(18):
   for y in range(18):assert P.R.overlap(W[x][y],pt(F(SG[x,y])));W[x][y]=pt(F(SG[x,y]))
  WI,_=inv(W);Z=[row[:]for row in W];piv=[]
  for z in range(20):
   v=Z[z][z];assert v[0]>0;piv.append(save(v))
   for x in range(z+1,20):
    for y in range(x,20):Z[x][y]=Z[y][x]=sub(Z[x][y],div(mul(Z[x][z],Z[z][y]),v))
  hv=[pair(H,N[x])for x in active];coeff=mv(WI,hv);v=plus(H,scale(lin([N[x]for x in active],coeff),pt(-1)));vq=sub(hq,inner(hv,coeff));assert vq[1]<0
  # Evaluate products through the defining projection identity. Direct
  # coordinate norms lose correlations at these very shallow corners.
  vp=[sub(pair(H,N[x]),inner(coeff,[pair(N[y],N[x])for y in active]))for x in ret]
  assert P.R.overlap(pair(v,v),vq)and P.R.overlap(pair(v,H),vq)
  for x,g in zip(ret,vp):
   if x in active:assert P.R.overlap(g,pt(0))
   else:assert g[1]<0,(active,x,save(g))
  corners.append(dict(active=active,positive_pivots=piv,point=vec(v),norm=save(vq),all171_pairs=vec(vp),product_identity='H.n - (G_corner^-1 h).(corner.n); v.H=v^2=H^2-h^T G_corner^-1 h'))
 Eface=lin(B,[pt(F(x))for x in ESCAPE]);assert pair(Eface,Eface)[0]>0;escape=None
 for exponent in range(8,120):
  epsE=F(1,2**exponent);e=plus(Eface,scale(dual,pt(-epsE)));eq=pair(e,e);eh=pair(e,H);ep=[pair(e,N[j])for j in ret]
  if eq[0]>0 and eh[1]<0 and all(g[1]<0 for g in ep):escape=dict(basis=[95,142,145,31],rational_coefficients=ESCAPE,support_perturbation=str(epsE),point=vec(e),norm=save(eq),center_pair=save(eh),all171_pairs=vec(ep));break
 assert escape
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.GEOM,f'results/{P.TAG}.json',P.CUSP,'results/fresh_h20_e8_f4_cusps_20261017.json','results/fresh_h20_prism_replacement_20261017.json','results/fresh_h20_h4_hexagon_closure_20261017.json']
 return dict(status='PASS_EXACT_171_WALL_REPLACEMENT_WITH_COMPACT_HEXAGON',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},removed_old_wall_ids=sorted(REMOVED),retained_old_wall_count=166,new_wall_labels=list(range(168,173)),actual_wall_labels=ret,actual_facets=171,old_base='Frozen168 retained except31 and158; this is an alternative branch, not a retaining extension.',exact_definitions=steps,contact_counts=dict(counts),new_prior_pair_count=840,center=dict(residual_basis=[95,142,145,31],rational_coefficients=CENTER,support_perturbation=str(eps),point=vec(H),norm=save(hq),all171_pairs=vec(hp)),facet_point_norms=vec(facets),compact_hexagon=dict(support=S,support_gram=[[str(x)for x in row]for row in SG.tolist()],cyclic_sides=cycle,all_other_projected_norms={str(j):save(projected[j])for j in ret if j not in S},pair_tests=pair_tests,finite_vertices=corners,complete_side_count=6,complete_vertex_count=6),strict_escape=escape,scope='Exact coherent171-wall Coxeter chamber on the explicit alternative retaining166 old roots and five new roots. All840 new/prior gates, common strict future center, every actual facet, complete compact six-sided codim18 face and six actual finite vertices, plus strict spacelike escape. No whole ambient finite-vertex/cusp census and no finite-volume claim; the strict escape proves infinite volume.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
