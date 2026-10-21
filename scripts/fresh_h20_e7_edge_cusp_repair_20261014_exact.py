#!/usr/bin/env python3
"""Five contact-producing affine cusp repairs on the frozen90 chamber."""
import argparse,json
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from hashlib import sha256
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
import fresh_h20_e7_finite_edgewalk_20261014_frontier as G
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_e7_edge_cusp_repair_20261014_exact'
BASE='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json';FRONT='results/fresh_h20_e7_finite_edgewalk_20261014_frontier.json';OLD='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json';Z=B.Z
add,sub,mul,scale,sg=B.add,B.sub,B.mul,B.scale,B.sign

def sc(q):return(F(q),)+(F(0),)*15

def inv(a):
 if len(a)==1:assert a[0];return(1/a[0],)
 h=len(a)//2;x,y=a[:h],a[h:];den=sub(mul(x,x),scale(mul(y,y),B.P[h.bit_length()-1]));v=inv(den)
 return mul(x,v)+scale(mul(y,v),-1)

def total(xs):
 a=Z
 for x in xs:a=add(a,x)
 return a

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return[enc(y)for y in x]
 if isinstance(x,dict):return{str(k):enc(y)for k,y in x.items()}
 return x

def run():
 d,f,old=[json.loads((ROOT/p).read_text())for p in[BASE,FRONT,OLD]]
 for doc in[d,f,old]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 def bv(v):return tuple(tuple(map(F,x))+(F(0),)*(16-len(x))for x in v)
 def ep(x):return tuple(tuple(map(F,a))for a in x)
 def ev(v):return tuple((x,Z)for x in v)
 EN=[tuple(map(ep,v))for v in d['normals']];N={i:tuple(a for a,b in v)for i,v in enumerate(EN)if all(b==Z for a,b in v)};assert len(N)==89
 Q=[[F(x)for x in row]for row in d['form']];sp=[[(j,q)for j,q in enumerate(row)if q]for row in Q];cache={}
 def dot(v,w):
  if w not in cache:cache[w]=tuple(tuple(sum(q*w[j][k]for j,q in row)for k in range(16))for row in sp)
  return total(mul(a,b)for a,b in zip(v,cache[w]))
 def lin(vs,cs):return tuple(total(mul(c,v[k])for c,v in zip(cs,vs))for k in range(21))
 H=tuple(a for a,b in d['strict_center']);H=tuple(tuple(map(F,x))for x in H)
 D=tuple(map(F,d['field']['discriminant']));EZ=(Z,Z);ecache={}
 def ea(a,b):return(add(a[0],b[0]),add(a[1],b[1]))
 def es(a,q):return(scale(a[0],q),scale(a[1],q))
 def em(a,b):return(add(mul(a[0],b[0]),mul(D,mul(a[1],b[1]))),add(mul(a[0],b[1]),mul(a[1],b[0])))
 def sign(a):
  x,y=a;sx,sy=sg(x),sg(y)
  if not sy:return sx
  if not sx or sx==sy:return sy
  return sx*sg(sub(mul(x,x),mul(D,mul(y,y))))
 def et(xs):
  a=EZ
  for x in xs:a=ea(a,x)
  return a
 def edot(v,w):
  if all(b==Z for a,b in v)and all(b==Z for a,b in w):return(dot(tuple(a for a,b in v),tuple(a for a,b in w)),Z)
  if w not in ecache:ecache[w]=tuple((tuple(sum(q*w[j][0][k]for j,q in row)for k in range(16)),tuple(sum(q*w[j][1][k]for j,q in row)for k in range(16)))for row in sp)
  return et(em(a,b)for a,b in zip(v,ecache[w]))
 sq2=tuple(F(i==1)for i in range(16));sq3=tuple(F(i==8)for i in range(16));delta=bv(f['new_parallel_cusp']['point']);L=bv(f['new_parallel_cusp']['lineality_vector']);assert dot(delta,delta)==Z and dot(L,L)==sc(F(98,27))and dot(delta,L)==Z
 A5=lin([N[i]for i in[12,13,18,31,65]],[sc(-1)]*5);A1=lin([N[15]],[sc(-1)]);C3=lin([N[17],N[79],N[83]],[scale(sq2,-1),scale(sq2,-1),sc(-1)]);line=lin([L],[scale(sq3,F(3,7))])
 bases=[A5,A1,C3,line,lin([line],[sc(-1)])];shifts=[sc(1),sc(2),sc(1),add(sc(1),scale(sq3,F(4,7))),scale(sq3,F(10,7))];names=['A5','A1','C3_highest_short','line_plus','line_minus']
 new=[];recipes=[]
 for name,b,k in zip(names,bases,shifts):
  assert sg(k)>0 and dot(b,b)==sc(2)and dot(b,delta)==Z
  r=lin([b,delta],[sc(1),k]);assert dot(r,r)==sc(2)and dot(r,delta)==Z
  new.append(r);recipes.append(dict(name=name,base=b,shift=k,normal=r))
 for i,r in enumerate(new):N[90+i]=r;EN.append(ev(r))
 norms=list(map(F,d['norms']))+[F(2)]*5;counts=Counter();pairs=[]
 for i in range(90,95):
  for j in range(i):
   g=edot(EN[i],EN[j]);prod=norms[i]*norms[j];assert sign(g)<=0
   if sign(g)==0:kind='right'
   elif em(g,g)==(sc(prod/4),Z):kind='pi3'
   elif em(g,g)==(sc(prod/2),Z):kind='pi4'
   elif em(g,g)==(sc(3*prod/4),Z):kind='pi6'
   elif em(g,g)==(sc(prod),Z):kind='parallel'
   else:assert sign(ea(em(g,g),(sc(-prod),Z)))>0;kind='dotted'
   counts[kind]+=1;pairs.append(dict(walls=[i,j],pairing=g,relation=kind))
 assert len(pairs)==460 and all(sign(edot(ev(H),n))<0 for n in EN)
 facet_points=[]
 for i,r in enumerate(new):
  p=lin([H,r],[sc(1),scale(dot(H,r),F(-1,2))]);assert dot(p,r)==Z and sg(dot(p,p))<0 and sg(dot(p,H))<0
  facet_points.append(dict(wall=90+i,point=p,norm=dot(p,p)))
 # Every original facet remains actual by the same strict-center
 # projection lemma, since all new pairings with it are nonpositive.
 vertices=d['preserved_old_vertices']+d['new_cap_finite_vertices'];pres=[]
 for i,row in enumerate(vertices):
  v=bv(row['point']);g=[dot(v,r)for r in new];assert all(sg(x)<0 for x in g);pres.append(dict(vertex=i,new_pairings=g))
 old_cusps=[]
 for row in old['completed_cusp_links']:
  nu=bv(row['point']);g=[dot(nu,r)for r in new];assert all(sg(x)<0 for x in g);old_cusps.append(dict(point=nu,new_pairings=g,quotient_rank=19))
 cpairs=[edot(ev(delta),n)for n in EN];assert all(sign(x)<=0 for x in cpairs);active=[i for i,x in enumerate(cpairs)if sign(x)==0];assert active==f['new_parallel_cusp']['incident']+list(range(90,95))
 def scaled(i):return lin([N[i]],[sq2])if i in[17,79]else N[i]
 rows=[scaled(i)for i in active]
 def rational(a):assert not any(a[1:]);return a[0]
 Gram=[[rational(dot(a,b))for b in rows]for a in rows];assert len(G.rref(Gram)[1])==19
 rem=set(range(len(rows)));components=[];coord=next(i for i,x in enumerate(delta)if sg(x))
 while rem:
  c={min(rem)};todo=list(c)
  for i in todo:
   for j in sorted(rem-c):
    if Gram[i][j]:c.add(j);todo.append(j)
  ids=sorted(c);rem-=c;M=[[Gram[i][j]for j in ids]for i in ids];kernel=G.kernel(M);assert len(kernel)==1
  w=kernel[0]
  if all(x<0 for x in w):w=[-x for x in w]
  assert all(x>0 for x in w)and len(G.rref(M)[1])==len(ids)-1
  from fresh_h20_e8_interior_caps_20261014_exact_review import pd
  assert pd([row[:-1]for row in M[:-1]])
  s=lin([rows[i]for i in ids],list(map(sc,w)));mu=mul(s[coord],inv(delta[coord]));assert sg(mu)>0 and all(a==mul(mu,b)for a,b in zip(s,delta))
  components.append(dict(walls=[active[i]for i in ids],gram=M,positive_kernel=w,quotient_rank=len(ids)-1,positive_null_multiple=mu))
 assert len(components)==6 and sum(c['quotient_rank']for c in components)==19
 # Explicit actual finite contacts outside the old cusp link.
 ridges=[]
 for i,j in [(90,0),(90,14),(91,14),(92,0),(92,10),(93,8),(94,0),(94,10)]:
  g=dot(N[i],N[j]);prod=norms[i]*norms[j];det=sub(sc(prod),mul(g,g));assert sg(det)>0
  hi,hj=dot(H,N[i]),dot(H,N[j]);ci=mul(sub(scale(hi,norms[j]),mul(g,hj)),inv(det));cj=mul(sub(scale(hj,norms[i]),mul(g,hi)),inv(det));p=lin([H,N[i],N[j]],[sc(1),scale(ci,-1),scale(cj,-1)])
  assert sg(dot(p,p))<0 and sg(dot(p,H))<0
  assert all(sign(edot(ev(p),n))==0 if k in[i,j]else sign(edot(ev(p),n))<0 for k,n in enumerate(EN))
  ridges.append(dict(walls=[i,j],point=p,norm=dot(p,p)))
 # Bounded continuation test on the541 previously unresolved germs at
 # the40 new corners. First new intersection must be timelike and strict
 # on every other wall; no complete new-cap vertex census is claimed.
 byvertex={i:[]for i in range(44,84)}
 for row in f['new_corner_directions']:
  if row['other_finite_vertex']is None:byvertex[row['vertex']].append(row)
 edge_endpoints=[];finite_trials=0
 for vi,directions in byvertex.items():
  if not directions:continue
  row=vertices[vi];v=bv(row['point']);ids=row['incident'];rr=[scaled(i)for i in ids];GM=[[rational(dot(a,b))for b in rr]for a in rr];GI=B.inv(GM);v2=dot(v,v);rv=[dot(v,r)for r in new]
  for edge in directions:
   j=edge['released'];pos=ids.index(j);s=sq2 if j in[17,79]else sc(1);dual=lin(rr,[scale(s,q)for q in GI[pos]]);d2=dot(dual,dual)
   assert dot(v,dual)==Z and dot(dual,N[j])==sc(1)
   choices=[]
   for k,r in enumerate(new):
    den=dot(r,dual)
    if sg(den)<0:
     t=mul(rv[k],inv(den));assert sg(t)>0;choices.append((t,k))
   if not choices:continue
   t,k=choices[0]
   for tt,kk in choices[1:]:
    if sg(sub(tt,t))<0:t,k=tt,kk
   norm=add(v2,mul(mul(t,t),d2))
   if sg(norm)>=0:continue
   finite_trials+=1;p=lin([v,dual],[sc(1),scale(t,-1)]);assert dot(p,p)==norm
   inc=edge['support']+[90+k];signs=[sign(edot(ev(p),n))for n in EN]
   if any(x!=0 if i in inc else x>=0 for i,x in enumerate(signs)):continue
   assert sg(dot(p,H))<0
   edge_endpoints.append(dict(old_vertex=vi,released=j,new_wall=90+k,parameter=t,point=p,norm=norm,incident=inc))
 assert finite_trials==12 and len(edge_endpoints)==12
 assert len({tuple(sorted(r['incident']))for r in edge_endpoints})==12
 assert Counter(r['new_wall']for r in edge_endpoints)=={92:6,93:6}
 proper_ideal_edges=[]
 for vi,row in enumerate(edge_endpoints):
  common=sorted(set(row['incident'])&set(active))
  if len(common)!=19:continue
  rr=[scaled(i)for i in common];MM=[[rational(dot(a,b))for b in rr]for a in rr]
  assert pd(MM)and sg(dot(row['point'],delta))<0
  # Both endpoints are feasible; every wall outside the common support
  # is strict at at least one end, hence strict on the open segment.
  for i,n in enumerate(EN):
   gp=edot(ev(row['point']),n);gn=cpairs[i]
   assert sign(gp)==sign(gn)==0 if i in common else sign(gp)<0 or sign(gn)<0
  proper_ideal_edges.append(dict(finite_endpoint=vi,support=common,finite_null_pairing=dot(row['point'],delta)))
 assert len(proper_ideal_edges)==6 and all(edge_endpoints[r['finite_endpoint']]['new_wall']==93 for r in proper_ideal_edges)
 escape=tuple(map(ep,d['strict_escape']['point']));epairs=[edot(escape,n)for n in EN];assert all(sign(g)<0 for g in epairs)and sign(edot(escape,escape))>0 and sign(edot(escape,ev(H)))<0
 paths=[BASE,FRONT,OLD,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','scripts/fresh_h20_e7_finite_edgewalk_20261014_frontier.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_95_WALL_CONTACT_CUSP_REPAIR_STILL_INFINITE',sources={p:sha(p)for p in paths},field=d['field'],form=Q,normals=EN,norms=norms,recipes=recipes,new_pair_counts=dict(counts),new_pairs=pairs,strict_center=ev(H),new_facet_points=facet_points,preserved_vertices=pres,preserved_old_bounded_cusps=old_cusps,completed_cusp=dict(point=delta,all95_pairings=cpairs,incident=active,components=components,quotient_rank=19),actual_nonincident_ridges=ridges,edge_continuation_scope=dict(tested_germs=541,timelike_first_new_intersections=finite_trials,actual_new_endpoints=len(edge_endpoints),complete_new_vertex_census=False),new_edge_endpoints=edge_endpoints,finite_to_proper_ideal_edges=proper_ideal_edges,strict_escape=dict(point=escape,norm=edot(escape,escape),all95_pairings=epairs),scope='Retains all90 and all84 specified finitevertices. Five fixed quotient roots at contact-producing shifts complete the full rank19 cusp link. All460 new pairs and95facet actuality follow exact gates/commoncenter. Oldtwo bounded cusps preserved. Bounded541germ continuation test supplies actual endpoints only; no complete newvertex or globaledge census. Existing strictescape survives.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['new_pair_counts'],d['edge_continuation_scope'])
