#!/usr/bin/env python3
"""Exact local scout for joining two affine cusp factors in an H20 section.

Each quotient direction pairs nontrivially with two facet normals in one
affine factor. Unlike the earlier single-component bridge cut, the chosen
wall pair need not be adjacent, and two distinct factors are mixed.
"""
from __future__ import annotations
from itertools import combinations,product
from collections import Counter
import json
import numpy as np
import sympy as sp
import fresh_h20_new_route_20261013 as old
ROOT=old.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_joint_cusp_20261013.json'
ALLOWED=(sp.Rational(0),sp.Rational(1,4),sp.Rational(1,2),sp.Rational(3,4),sp.Rational(1))

def two_node_directions(Q,W,G,h,u,comp):
 M=G[np.ix_(comp,comp)];marks=old.positive_kernel(M);out=[]
 for a,b in combinations(range(len(comp)),2):
  values=sp.zeros(len(comp),1);values[a]=marks[b];values[b]=-marks[a]
  sol,params=sp.Matrix(M.tolist()).gauss_jordan_solve(values);sol=sol.subs({p:0 for p in params})
  e=sum((c*W[i] for c,i in zip(sol,comp)),np.zeros(22,dtype=object));e-=sp.Rational(h@Q@e,h@Q@u)*u
  out.append({'e':e,'q':e@Q@e,'v':W@Q@e,'bridge':(comp[a],comp[b]),'comp':comp})
 return out

def local_ratios(A,B,G,active):
 # qA+x*qB is the unscaled normal norm. Each edge from an endpoint to
 # a zero-pairing wall imposes one of three rational squared cosines.
 qa,qb=A['q'],B['q'];a,b=A['bridge'];c,d=B['bridge']
 for i,other,branch in ((a,b,0),(c,d,1)):
  n=next((j for j in active if j!=i and G[i,j] and (A['v'][j]==B['v'][j]==0)),None)
  if n is None:continue
  roots=set()
  for z in ALLOWED[1:-1]:
   delta=sp.Rational(G[i,i])-sp.Rational(G[i,n]**2,G[n,n])/z
   if branch==0:
    if delta:roots.add(sp.cancel((A['v'][i]**2/delta-qa)/qb))
   else:
    den=B['v'][i]**2-delta*qb
    if den:roots.add(sp.cancel(delta*qa/den))
  return sorted(x for x in roots if x>0),'endpoint_zero_edge'
 # This branch is retained for completeness of the declared bridge pattern.
 x=sp.Symbol('x');q=qa+x*qb
 N=x*A['v'][a]**2*B['v'][c]**2
 den=(G[a,a]*q-A['v'][a]**2)*(G[c,c]*q-x*B['v'][c]**2)
 roots=set()
 for z in ALLOWED[1:-1]:
  for r in sp.solve(N-z*den,x):
   if r.is_Rational and r>0:roots.add(r)
 return sorted(roots),'cross_edge_rational_roots'

def gram_join(A,B,G,active,x):
 q=A['q']+x*B['q'];bc=set(B['comp']);v=A['v']+B['v']
 return np.array([[sp.Rational(G[i,j])/(x if i in bc and j in bc else 1)-v[i]*v[j]/q
                    for j in active] for i in active],dtype=object)

def run():
 Q,W,G,h=old.base.load();raw=json.loads(old.CUSP_SOURCE.read_text());counts=Counter();survivors=[];rows=[]
 for ci,row in enumerate(raw['cusps']):
  u=np.array(row['vector'],dtype=object)
  if max(W@Q@u)>0:continue
  ids=[i for i,p in enumerate(W@Q@u) if p==0];comps=old.components(G,ids)
  dirs=[two_node_directions(Q,W,G,h,u,c) for c in comps]
  for ca,cb in combinations(range(len(comps)),2):
   for A0,B0 in product(dirs[ca],dirs[cb]):
    for sa,sb in product((1,-1),repeat=2):
     A={**A0,'v':sa*A0['v'],'e':sa*A0['e'],'bridge':A0['bridge'] if sa==1 else A0['bridge'][::-1]}
     B={**B0,'v':sb*B0['v'],'e':sb*B0['e'],'bridge':B0['bridge'] if sb==1 else B0['bridge'][::-1]}
     removed=[A['bridge'][1],B['bridge'][1]];active=[i for i in ids if i not in removed]
     xs,mode=local_ratios(A,B,G,active);counts['oriented_component_two_wall_pairs']+=1;counts[mode]+=1
     record={'cusp_index':ci,'oriented_bridges_one_based':[[i+1 for i in A['bridge']],[i+1 for i in B['bridge']]],'ratios':[str(x) for x in xs],'outcomes':[]}
     for x in xs:
      H=gram_join(A,B,G,active,x);bad=None
      for i in range(len(active)):
       if H[i,i]<=0:bad={'kind':'nonpositive_norm','wall':active[i]+1};break
      if bad is None:
       for a,b in combinations(range(len(active)),2):
        z=H[a,b]**2/(H[a,a]*H[b,b])
        if H[a,b]>0 or z not in ALLOWED:
         bad={'kind':'pair_conflict','walls':[active[a]+1,active[b]+1],'inner_product':str(H[a,b]),'cosine_squared':str(z)};break
      if bad is None:
       cs=[]
       for js in old.components(H):
        marks=old.positive_kernel(H[np.ix_(js,js)])
        if marks is None:bad={'kind':'no_positive_affine_kernel','walls':[active[i]+1 for i in js]};break
        cs.append({'walls_one_based':[active[i]+1 for i in js],'marks':marks,'rank':len(js)-1})
      if bad is None:
       assert sum(c['rank'] for c in cs)==19
       saved={'cusp_index':ci,'u':list(map(int,u)),'eA':list(map(str,A['e'])),'eB':list(map(str,B['e'])),
              'qA':str(A['q']),'qB':str(B['q']),'s_squared':str(x),'active':list(i+1 for i in active),
              'removed':[i+1 for i in removed],'components':cs,'oriented_bridges_one_based':record['oriented_bridges_one_based'],
              'scaled_gram':[[str(v) for v in row] for row in H]}
       survivors.append(saved);print('LOCAL SURVIVOR',ci,record['oriented_bridges_one_based'],x,[c['rank'] for c in cs],flush=True)
      counts['ratio_candidates']+=1;counts['local_affine_survivors' if bad is None else bad['kind']]+=1
      record['outcomes'].append({'s_squared':str(x),'rejection':bad})
     rows.append(record)
 print(dict(counts),flush=True)
 OUTPUT.write_text(json.dumps({'schema':'joint-affine-cusp-local-v2','counts':dict(counts),'cases':rows,'local_survivors':survivors,'scope':'Local scout mixing two arbitrary two-node directions in distinct recorded cusp factors, retaining all incident facets except two prescribed endpoints. Endpoint-zero edges force rational ratios; exceptional cross-only cases use rational roots only. No general section or facet-subset exclusion.'},indent=2)+'\n')
if __name__=='__main__':run()
