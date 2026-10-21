#!/usr/bin/env python3
"""Complete 1728 edge germs of one repaired cusp; no new cap search."""
import argparse, json, itertools, time
from fractions import Fraction as F
from pathlib import Path
from hashlib import sha256
from collections import Counter
import fresh_h20_e7_causal_residuals_20261014_bridge_review as B
from fresh_h20_e7_edge_cusp_repair_20261014_exact import inv as binv
ROOT=Path(__file__).resolve().parents[1]; TAG='fresh_h20_e7_repaired_cusp_walk_20261014'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PREV='results/fresh_h20_e7_finite_edgewalk_20261014_exact.json'
Z=B.Z; add,sub,mul,scale,sg=B.add,B.sub,B.mul,B.scale,B.sign
EZ=(Z,Z)
def sc(q):return (F(q),)+(F(0),)*15
def total(xs):
 a=Z
 for x in xs:a=add(a,x)
 return a
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,(list,tuple)):return [enc(y) for y in x]
 if isinstance(x,dict):return {str(k):enc(v) for k,v in x.items()}
 return x
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 start=time.monotonic(); d=json.loads((ROOT/BASE).read_text()); old=json.loads((ROOT/PREV).read_text())
 for doc in [d,old]:
  for p,h in doc['sources'].items():assert sha(p)==h,p
 D=tuple(map(F,d['field']['discriminant']));assert sg(D)>0
 def ea(a,b):return (add(a[0],b[0]),add(a[1],b[1]))
 def es(a,q):return (scale(a[0],q),scale(a[1],q))
 def em(a,b):
  if a[1]==Z:return (mul(a[0],b[0]),mul(a[0],b[1]))
  if b[1]==Z:return (mul(a[0],b[0]),mul(a[1],b[0]))
  return(add(mul(a[0],b[0]),mul(D,mul(a[1],b[1]))),add(mul(a[0],b[1]),mul(a[1],b[0])))
 def sign(a):
  x,y=a;sx,sy=sg(x),sg(y)
  if not sy:return sx
  if not sx or sx==sy:return sy
  return sx*sg(sub(mul(x,x),mul(D,mul(y,y))))
 def et(xs):
  a=EZ
  for x in xs:a=ea(a,x)
  return a
 def eparse(v):return tuple(tuple(tuple(map(F,a)) for a in x) for x in v)
 def bv(v):return tuple(tuple(map(F,a)) for a in v)
 N=list(map(eparse,d['normals']));H=eparse(d['strict_center']);nu=tuple((a,Z) for a in bv(d['completed_cusp']['point']))
 Q=[[F(x) for x in row] for row in d['form']];sparse=[[(j,q) for j,q in enumerate(row) if q] for row in Q];cache={}
 def dot(v,w):
  if w not in cache:cache[w]=tuple(et(es(w[j],q) for j,q in row) for row in sparse)
  return et(em(a,b) for a,b in zip(v,cache[w]))
 def lin(vs,cs):return tuple(et(em(c,v[k]) for c,v in zip(cs,vs)) for k in range(21))
 def rational(x):assert x[1]==Z and not any(x[0][1:]);return x[0][0]
 Hn=[dot(H,n) for n in N];nun=[dot(nu,n) for n in N];H2=dot(H,H);Hnu=dot(H,nu)
 assert all(sign(x)<0 for x in Hn) and sign(H2)<0 and sign(Hnu)<0 and dot(nu,nu)==EZ
 assert H2[1]==Hnu[1]==Z
 incident=[i for i,g in enumerate(nun) if sign(g)==0];assert incident==d['completed_cusp']['incident']
 noninc=[i for i in range(95) if i not in incident];assert all(sign(nun[i])<0 for i in noninc)
 den=[es(x,-1) for x in nun]
 sq2=(tuple(F(i==1) for i in range(16)),Z)
 scaled={i:lin([N[i]],[sq2]) if i in [17,79] else N[i] for i in incident}
 roots=[scaled[i] for i in incident];G=[[rational(dot(a,b)) for b in roots] for a in roots]
 components=d['completed_cusp']['components'];assert sorted(sum((c['walls'] for c in components),[]))==sorted(incident)
 for c in components:
  ids=c['walls'];g=[[G[incident.index(i)][incident.index(j)] for j in ids] for i in ids]
  assert g==[[F(x) for x in row] for row in c['gram']]
  weights=list(map(F,c['positive_kernel']));assert all(x>0 for x in weights)
  assert all(sum(g[i][j]*weights[j] for j in range(len(ids)))==0 for i in range(len(ids)))
  assert all(G[incident.index(i)][incident.index(j)]==0 for i in ids for j in incident if j not in ids)
 # Only25 positive-definite component projections are needed.
 rootpairs={i:[dot(scaled[i],n) for n in N] for i in incident};hroot={i:dot(H,scaled[i]) for i in incident}
 from fresh_h20_e8_interior_caps_20261014_exact_review import pd
 choices=[]
 for ci,c in enumerate(components):
  rows=[]
  for omitted in c['walls']:
   ids=[i for i in c['walls'] if i!=omitted];g=[[rational(dot(scaled[i],scaled[j])) for j in ids] for i in ids];assert pd(g)
   gi=B.inv(g);coef=[et(es(hroot[j],v) for j,v in zip(ids,row)) for row in gi]
   P=lin([scaled[i] for i in ids],coef);p2=et(em(a,hroot[i]) for i,a in zip(ids,coef));pn=[et(em(a,rootpairs[i][j]) for i,a in zip(ids,coef)) for j in range(95)]
   assert all(a[1]==Z for a in P) and p2[1]==Z and dot(P,nu)==EZ and dot(P,P)==p2
   assert all(dot(P,scaled[i])==hroot[i] for i in ids)
   rows.append(dict(omitted=omitted,retained=ids,projection=P,norm=p2,pairings=pn))
  choices.append(rows)
 print('cached',sum(map(len,choices)),'component projections',round(time.monotonic()-start,2),flush=True)
 # Numerical ordering proposes a maximizing wall only. Every threshold is
 # then compared exactly, and every final endpoint sign is checked exactly.
 rad=[float(x)**.5 for x in d['field']['base_square_root_basis']]
 def val(a):return sum(float(q)*r for q,r in zip(a,rad))
 dv=val(D)**.5
 def ev(a):return val(a[0])+val(a[1])*dv
 denv={i:ev(den[i]) for i in noninc};assert all(x>0 for x in denv.values())
 records=[];vertices={};ideal={};counts=Counter();winners=Counter();ties=Counter();nullpoints={}
 known=old['preserved_old_vertices']+old['new_cap_finite_vertices']+d['new_edge_endpoints']
 known_sets={tuple(sorted(r['incident'])):i for i,r in enumerate(known)}
 for idx,selected in enumerate(itertools.product(*choices)):
  omitted=[r['omitted'] for r in selected];support=sorted(sum((r['retained'] for r in selected),[]));assert len(support)==19
  w2=ea(H2,es(et(r['norm'] for r in selected),-1));assert sign(w2)<0
  wn=[ea(Hn[j],es(et(r['pairings'][j] for r in selected),-1)) for j in range(95)]
  assert all(sign(wn[i])==0 for i in support) and all(sign(wn[i])<0 for i in omitted)
  win=max(noninc,key=lambda i:ev(wn[i])/denv[i])
  for j in noninc:
   if sign(ea(em(wn[j],den[win]),es(em(wn[win],den[j]),-1)))>0:win=j
  dd,tt=den[win],wn[win]
  xp=[ea(em(dd,x),em(tt,y)) for x,y in zip(wn,nun)]
  signs=[sign(x) for x in xp];assert all(s<=0 for s in signs)
  active=[i for i,s in enumerate(signs) if s==0];assert set(support)<=set(active) and all(signs[i]<0 for i in omitted)
  nn=ea(em(em(dd,dd),w2),es(em(em(dd,tt),Hnu),2));nh=ea(em(dd,w2),em(tt,Hnu));nnu=em(dd,Hnu)
  sn=sign(nn);kind='finite' if sn<0 else 'ideal_contact' if sn==0 else 'open_opposite_end'
  counts[kind]+=1;winners[win]+=1;ties[len(active)-19]+=1
  rec=dict(omitted=omitted,support=support,first_contact=win,contact_walls=[i for i in active if i not in support],type=kind,endpoint_norm=nn)
  if sn<=0:
   assert sign(nh)<0 and sign(nnu)<0
   key=tuple(active);pool=vertices if sn<0 else ideal
   if key not in pool:
    W=lin([H]+[r['projection'] for r in selected],[(sc(1),Z)]+[(sc(-1),Z)]*6);assert dot(W,W)==w2 and dot(W,nu)==Hnu
    X=lin([W,nu],[dd,tt]);assert dot(X,X)==nn and dot(X,H)==nh and dot(X,nu)==nnu
    assert all(dot(X,n)==g for n,g in zip(N,xp))
    if sn<0:
     assert len(active)==20 # actual finite Coxeter vertices are simple
     pool[key]=dict(incident=active,point=X,norm=nn,future_pairing=nh,cusp_pairing=nnu,known_vertex=known_sets.get(key),germs=[])
    else:pool[key]=dict(incident=active,point=X,norm=nn,future_pairing=nh,cusp_pairing=nnu,germs=[])
   pool[key]['germs'].append(idx);rec['endpoint_key']=active
  else:
   # At the opposite light-cone point Z=(-2Hnu)W+W²nu, every
   # nonincident old wall must still be strict: exact unbounded-end proof.
   zp=[ea(em(es(Hnu,-2),x),em(w2,y)) for x,y in zip(wn,nun)]
   assert all(sign(zp[i])==0 if i in support else sign(zp[i])<0 for i in range(95))
   rec['opposite_null_incident']=support
  records.append(rec)
  if (idx+1)%128==0:print(idx+1,dict(counts),round(time.monotonic()-start,2),flush=True)
 assert len(records)==1728
 # Distinct finite vertices are separated by their exact full active sets;
 # each has precisely20 independent roots because the19 PD support plus
 # first wall has a nonzero edge derivative, and X is timelike.
 finite=list(vertices.values());ideals=list(ideal.values())
 newcount=sum(r['known_vertex'] is None for r in finite)
 paths=[BASE,PREV,f'scripts/{TAG}.py','scripts/fresh_h20_e7_causal_residuals_20261014_bridge_review.py','scripts/fresh_h20_e7_edge_cusp_repair_20261014_exact.py','scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
 return enc(dict(status='PASS_COMPLETE_REPAIRED_CUSP_EDGE_CENSUS',sources={p:sha(p) for p in paths},field=d['field'],cusp_point=nu,component_sizes=list(map(len,choices)),component_projection_count=25,counts=dict(counts),first_contact_counts=dict(winners),contact_tie_counts=dict(ties),distinct_finite_vertices=len(finite),finite_vertices_new_to_known96=newcount,distinct_ideal_contact_points=len(ideals),component_projections=choices,edge_germs=records,finite_vertices=finite,ideal_contact_points=ideals,scope='Complete1728 edges from this one bounded cusp, not a whole95-wall vertex/edge census. Exact positive scaling avoids field division. Every first hit is certified against all95 walls; open opposite ends have exactly19 PD active roots and no bounded horospherical link. The inherited global strict escape still proves infinite volume.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],d['distinct_finite_vertices'],d['finite_vertices_new_to_known96'])
