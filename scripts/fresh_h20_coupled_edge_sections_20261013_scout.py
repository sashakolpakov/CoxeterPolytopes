#!/usr/bin/env python3
"""Bounded numerical proposals for genuinely coupled section pairings."""
from pathlib import Path
from collections import Counter,defaultdict
from itertools import combinations
import json, math, hashlib
import numpy as np
import fresh_h20_h4_edge_sections_20261013_scout as S
ROOT=S.ROOT;TAG='fresh_h20_coupled_edge_sections_20261013'
def chain(cc,G):
 if len(cc)<2 or any(G[i,i]!=2 for i in cc):return None
 if any(G[i,j]not in[0,-1] for i,j in combinations(cc,2)):return None
 adj={i:[j for j in cc if i!=j and G[i,j]] for i in cc}
 if sorted(map(len,adj.values()))!=[1,1]+[2]*(len(cc)-2):return None
 out=[min(i for i in cc if len(adj[i])==1)]
 while len(out)<len(cc):out.append(next(j for j in adj[out[-1]] if j not in out))
 return out
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 Q,W,G,h=S.load();C=W@Q;walk=json.loads((ROOT/S.WALK).read_text());corners=walk['corner_orbits']
 buckets=defaultdict(list);seen=set()
 for c in corners:
  for ei,edge in enumerate(c['edge_orbits']):
   E=edge['edge_walls'];key=tuple(E)
   if key in seen:continue
   seen.add(key);components=S.comps(E,G);chains=[z for cc in components if (z:=chain(cc,G))]
   singles=[cc[0] for cc in components if len(cc)==1]
   base={'corner':c['id'],'edge_index':ei,'edge_walls':E,'endpoint_orbits':[c['id'],edge['endpoint_orbit']],
         'endpoint_types':[c['type'],corners[edge['endpoint_orbit']]['type']]}
   for cc in chains:
    # Each expression is coefficient * sqrt(rational radicand).
    n=len(cc);buckets[f'A{n}_to_A1_B{n-1}'].append(dict(base,template=f'A{n}_to_A1_B{n-1}',
       component=cc,pairings=[[cc[0],'1','1'],[cc[1],'-1','1']]))
    if n in[2,3] and singles:
     j=singles[0];target='B3' if n==2 else 'F4'
     buckets[f'A1_A{n}_to_{target}'].append(dict(base,template=f'A1_A{n}_to_{target}',
        component=[j]+cc,pairings=[[j,'1',str(int(G[j,j]))+'/5'],[cc[0],'1','1']]))
    if n==3:
     buckets['A3_to_A2_A1_three_contacts'].append(dict(base,template='A3_to_A2_A1_three_contacts',
       component=cc,pairings=[[cc[0],'1/3','6'],[cc[1],'-1/2','6'],[cc[2],'1/3','6']]))
 selected=[]
 for typ,items in sorted(buckets.items()):
  seen_types=set()
  for r in items:
   key=tuple(sorted(r['endpoint_types']))
   if key in seen_types:continue
   seen_types.add(key);selected.append(r)
   if len(seen_types)==3:break
 from fractions import Fraction as F
 rows=[];allowed=np.array([math.cos(math.pi/m)**2 for m in range(2,129)])
 for si,seed in enumerate(selected):
  c=corners[seed['corner']];edge=c['edge_orbits'][seed['edge_index']];E=seed['edge_walls'];A=G[np.ix_(E,E)].astype(float);I=np.linalg.inv(A)
  a=np.zeros(20)
  for j,coef,rad in seed['pairings']:a[E.index(j)]=float(F(coef))*math.sqrt(float(F(rad)))
  v=I@a;w=v@W[E];spatial=a@v;gamma=1-spatial
  assert gamma>1e-8
  H=A-np.outer(a,a);assert np.linalg.eigvalsh(H)[0]>1e-7
  for i,j in combinations(range(20),2):
   assert H[i,j]<1e-8
   ratio=H[i,j]**2/(H[i,i]*H[j,j]);assert np.min(abs(allowed-ratio))<1e-8
  left=np.array(c['primitive_ray'],float);left/=-(left@Q@h)
  right=np.array(edge['endpoint_primitive_ray'],float);right/=-(right@Q@h)
  for ts in ['1/3','1/2','2/3']:
   t=float(F(ts));P=(1-t)*left+t*right;qp=P@Q@P
   z0=right-(right@Q@P/qp)*P;qz=z0@Q@z0;assert qp<0 and qz>0
   pp=C@P;assert np.max(abs(pp[E]))<1e-7 and max(pp[i]for i in range(210)if i not in E)<-1e-8
   for sign in [-1,1]:
    z=sign*math.sqrt(gamma/qz)*z0;normal=w+z;ep=C@normal
    assert abs(normal@Q@normal-1)<1e-7 and abs(normal@Q@P)<1e-7
    duals=I@W[E]-np.outer(v/gamma,z)
    bad=[];neighbors=[]
    for k,release in enumerate(E):
     dp=C@duals[k];bounds=[(pp[j]/dp[j],j)for j in range(210)if j not in E and dp[j]<-1e-9]
     if not bounds:continue
     step,j=min(bounds);x=P-step*duals[k];xp=C@x;xx=x@Q@x
     active=np.flatnonzero(abs(xp)<1e-7).tolist()
     if xx>=-1e-8 or len(active)!=20 or max(xp)>1e-7:continue
     B=G[np.ix_(active,active)]-np.outer(ep[active],ep[active])
     if np.linalg.eigvalsh(B)[0]<1e-7:continue
     neighbors.append({'release':release,'new_wall':j,'active':active,'norm':float(xx)})
     defects=[]
     for ia,ib in combinations(range(20),2):
      p=B[ia,ib];ratio=p*p/(B[ia,ia]*B[ib,ib])
      if p>1e-7:defects.append((p,'obtuse',active[ia],active[ib],ratio))
      elif p<-1e-7 and np.min(abs(allowed-ratio))>1e-5:
       if ratio<.25:defects.append((min(ratio,.25-ratio),'between_0_and_1/4',active[ia],active[ib],ratio))
       elif ratio<.5:defects.append((min(ratio-.25,.5-ratio),'between_1/4_and_1/2',active[ia],active[ib],ratio))
     if defects:
      p,kind,ia,ib,ratio=max(defects);bad.append({'release':release,'new_wall':j,'active':active,
         'bad_pair':[ia,ib],'kind':kind,'cosine_squared':float(ratio),'score':float(p)})
    row={'seed':si,'template':seed['template'],'parameter':ts,'radical_sign':sign,
      'spatial_norm':float(spatial),'normal':normal.tolist(),'checked_finite_neighbors':len(neighbors)}
    if bad:row.update(status='NUMERICAL_ACTUAL_VERTEX_DEFECT',witness=max(bad,key=lambda z:z['score']))
    else:row.update(status='UNRESOLVED',neighbors=neighbors)
    rows.append(row)
  print(si,seed['template'],seed['endpoint_types'],'cases',len(rows),'unresolved',sum(r['status']=='UNRESOLVED' for r in rows),flush=True)
 out={'status':'BOUNDED_NUMERICAL_COUPLED_TWO_AND_THREE_PAIRING_SECTIONS',
  'sources':{p:sha(p)for p in[S.RAW,S.WALK,'scripts/'+Path(__file__).name,'scripts/fresh_h20_h4_edge_sections_20261013_scout.py']},
  'literal_edge_supports_examined':len(seen),'template_proposal_counts':{k:len(v)for k,v in buckets.items()},
  'selected_edges':selected,'cases':rows,'counts':dict(Counter(r['status']for r in rows)),
  'scope':'One literal source edge per template and available unordered endpoint-type pair; positions1/3,1/2,2/3 and both radical signs. No complete parameter or edge-profile claim.'}
 (ROOT/'results'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n');print(out['counts'])
if __name__=='__main__':run()
