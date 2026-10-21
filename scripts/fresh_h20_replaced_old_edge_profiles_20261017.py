#!/usr/bin/env python3
"""Old finite-edge spherical extensions on the one-wall replacement base.

Exact finite-diagram deletion census from numerically recognized old labels;
native geometry and endpoint events remain numerical discovery.
"""
import argparse,hashlib,itertools,json,math
from collections import Counter
import networkx as nx
import numpy as np
import sympy as s
from scipy.linalg import null_space
import fresh_h20_opposite_prism_cap_20261017 as P
import fresh_h20_prism_replacement_20261017 as R
import fresh_h20_e8_partner_endpoints_20261017_bridges as B
ROOT=P.ROOT;TAG='fresh_h20_replaced_old_edge_profiles_20261017'

def census(W,gram):
 old=nx.Graph();old.add_nodes_from(W);labels=[]
 for i in range(19):
  for j in range(i):
   co=-gram[i,j]/math.sqrt(gram[i,i]*gram[j,j]);lab=R.classify(co)
   assert isinstance(lab,int),(W[i],W[j],co,lab)
   if lab>2:old.add_edge(W[i],W[j],label=lab);labels.append([W[i],W[j],lab])
 comps=sorted((sorted(c)for c in nx.connected_components(old)),key=lambda c:(-len(c),c));GG=s.eye(19)*2
 for i,j,lab in labels:GG[W.index(i),W.index(j)]=GG[W.index(j),W.index(i)]=-2*s.cos(s.pi/lab)
 JJ=GG.inv();rows={};matches=0
 for name,g in P.P.models():
  for deleted in g:
   h=g.copy();h.remove_node(deleted)
   for bits in itertools.product([0,1],repeat=len(comps)):
    ids=sum((c for c,b in zip(comps,bits)if b),[])
    if len(ids)!=len(h):continue
    gm=nx.algorithms.isomorphism.GraphMatcher(h,old.subgraph(ids),edge_match=lambda a,b:a['label']==b['label'])
    for mp in gm.isomorphisms_iter():
     cs=tuple(sorted((mp[j],g[deleted][j]['label'])for j in g.neighbors(deleted)));matches+=1
     if len(cs)==1 and old.degree(cs[0][0])==0:continue
     z=rows.setdefault(cs,dict(contacts=[list(x)for x in cs],models=[]));z['models'].append(dict(type=name,deleted=deleted,mapping={str(k):v for k,v in mp.items()}))
 out=[]
 for cs,z in sorted(rows.items(),key=lambda x:(-len(x[0]),x[0])):
  d=s.zeros(19,1)
  for i,lab in cs:d[W.index(i)]=2*s.cos(s.pi/lab)
  cost=s.simplify((d.T*JJ*d)[0]);res=s.simplify(2-cost);assert res.is_positive
  z.update(id=len(out),normalized_cost=str(cost),residual_norm=str(res));out.append(z)
 return dict(support=W,old_edges=labels,components=comps,profiles=out,parametric_singletons=[j for j in W if old.degree(j)==0],raw_matches=matches)

def conic(cen,row,G,Y,CT,J,q):
 T=G['T'];m=G['m'];edge=cen['edge'];U=Y[edge];inv=np.linalg.inv(U@m@U.T);dt=np.zeros(17);de=np.zeros(2)
 for j,lab in row['contacts']:
  val=-math.sqrt(2*q[j])*math.cos(math.pi/lab)
  if j in T:dt[T.index(j)]=val
  else:de[edge.index(j)]=val
 offset=CT@J@dt;f=de-offset[edge];p=(inv@f)@U;h=2-dt@J@dt-p@m@p;assert h>1e-9
 ht=np.array([1.,0,0,0]);ht-=((ht@m@U.T)@inv)@U;ht/=math.sqrt(-(ht@m@ht))
 es=null_space(np.r_[U@m,[ht@m]])[:,0];es/=math.sqrt(es@m@es)
 return [(sgn,p,math.sqrt(h)/2*(ht+sgn*es),math.sqrt(h)/2*(-ht+sgn*es),dt,offset)for sgn in(-1,1)]

def run():
 G,Y,CT,J,q,V,C=P.geometry();m=G['m'];T=G['T'];rep=P.P.load('results/fresh_h20_prism_replacement_20261017.json');Y=np.r_[Y,[np.array(rep['residual4'][1])]];CT=np.r_[CT,np.zeros((1,17))];q=np.r_[q,2];ret=[j for j in range(169)if j!=158]
 full=Y@m@Y.T+CT@J@CT.T;pos=q-np.einsum('ij,jk,ik->i',CT,J,CT);facets=[j for j in ret if j not in T and pos[j]>1e-7]
 bv=rep['individual_replacements'][1]['face']['vertices'];V=np.array([[1.]+z['point']for z in bv if z['radius_squared']<1-1e-6]);V/=np.sqrt(-np.einsum('ij,jk,ik->i',V,m,V))[:,None]
 ledger=[];censuses=[]
 for edge in itertools.combinations(facets,2):
  if edge[0]>=159 or edge[1]==168:continue # New partner-edge and base sides are already covered.
  W=T+list(edge)
  if np.linalg.eigvalsh(full[np.ix_(W,W)]).min()<1e-8:continue
  ends=[k for k,v in enumerate(V)if max(abs(Y[list(edge)]@m@v))<1e-7]
  ledger.append(dict(edge=list(edge),finite_endpoint_indices=ends,open=len(ends)<2))
  if len(ends)>=2:continue
  cen=census(W,full[np.ix_(W,W)]);cen.update(edge=list(edge),finite_endpoint_indices=ends);censuses.append(cen)
  print('OLD_EDGE',edge,'profiles',len(cen['profiles']),'singletons',cen['parametric_singletons'],flush=True)
 stats=Counter();out=[];near=[];trials=[]
 for ci,cen in enumerate(censuses):
  for row in cen['profiles']:
   for sg,p,a,b,dt,offset in conic(cen,row,G,Y,CT,J,q):
    aa=a@m@Y.T;bb=b@m@Y.T;cc=p@m@Y.T+offset;hh=2-dt@J@dt
    for j in facets:
     if j in cen['edge']:continue
     isolated=max(abs(dt))<1e-9 and max(abs(CT[j]))<1e-8;labs=range(2,65)if isolated else range(2,6)
     for lab in labs:
      stats['equations']+=1;target=-math.sqrt(2*q[j])*(0 if lab==2 else math.cos(math.pi/lab))
      for k in B.roots(aa[j],cc[j]-target,bb[j]):
       stats['real_branches']+=1;r=p+k*a+b/k
       if max(r@m@V.T)>=-1e-8:stats['old_vertex_not_strict']+=1;continue
       ps=r@m@Y.T+offset;co=-ps/np.sqrt(2*q);native=[R.classify(x)for x in co]
       if hh*pos[j]-(ps[j]-offset[j])**2<=1e-8:stats['extra_not_finite']+=1;continue
       bad=[u for u in ret if native[u]in('obtuse','angle_gap')];stats['intrinsic_events']+=1;trials.append([ci,row['id'],sg,j,lab,k,bad])
       rec=dict(census=ci,edge=cen['edge'],profile=row['id'],contacts=row['contacts'],branch=sg,anchor=j,label=lab,k=k,residual4=r.tolist(),T_pairings=dt.tolist(),bad_wall_ids=bad,bad_gate_details=[dict(wall=u,kind=native[u],cosine=float(co[u]))for u in bad])
       if len(bad)<=3:near.append(rec)
       if not bad:
        if any(max(abs(np.array(z['residual4'])-r))<1e-6 and max(abs(np.array(z['T_pairings'])-dt))<1e-6 for z in out):continue
        rec['all_retained_labels']=[native[u]for u in ret];rec['face']=R.poly(np.r_[Y,[r]],ret+[169]);out.append(rec);stats['accepted']+=1
  print('SOLVED_OLD_EDGE',ci,dict(stats),'accepted',len(out),'near',len(near),flush=True)
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',f'scripts/{P.P.TAG}.py',f'scripts/{R.TAG}.py',f'scripts/{B.TAG}.py','results/fresh_h20_prism_replacement_20261017.json']
 return dict(status='FINITE_OLD_EDGE_DELETION_CENSUS_WITH_BOUNDED_NATIVE_EVENT_SCOUT',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},retained_old_ids=ret,actual_T_facets=facets,old_edge_ledger=ledger,endpoint_censuses=censuses,statistics=dict(stats),candidates=out,at_most_three_bad_gates=near,all_intrinsic_event_trials=trials,scope='Alternative base retains all168 old walls except158 and adds the T-right replacement. Fixed endpoint profiles are an exact finite-diagram deletion census conditional on numerically recognized old native labels. Both conic branches and all actual-old-facet finite events are solved numerically; nonisolated event labels2..5 and isolated2..64. Explicit singleton endpoint families, unbounded isolated event labels, parabolic endpoints, and intermediate facets are not exhausted. No exact ambient promotion or all-family exclusion.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['statistics'],flush=True)
