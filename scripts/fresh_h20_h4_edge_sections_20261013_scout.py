#!/usr/bin/env python3
"""Bounded numerical A4->H4 interior-edge section pilot; no exact predicates."""
from pathlib import Path
import collections,hashlib,itertools,json,math
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_h4_edge_sections_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';WALK='results/fresh_h21_source_corner_walk_20261013.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def load():
 r=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=-np.array(r['basis_gram'],dtype=np.float64);W=np.array(r['wall_coordinates_as_columns_transposed'],dtype=np.float64)
 G=np.rint(W@Q@W.T).astype(int);h=W.sum(axis=0).astype(int);h=h//math.gcd(*map(int,h));return Q,W,G,h.astype(float)
def comps(E,G):
 left=set(E);out=[]
 while left:
  cc=[min(left)];left.remove(cc[0])
  for i in cc:
   nn=[j for j in sorted(left)if G[i,j]];cc+=nn;left.difference_update(nn)
  out.append(sorted(cc))
 return out
def A4(cc,G):
 return len(cc)==4 and all(G[i,i]==2 for i in cc) and sorted(sum(G[i,j]!=0 for j in cc)-1 for i in cc)==[1,1,2,2] and all(G[i,j]in[0,-1]for i in cc for j in cc if i!=j)
def run():
 Q,W,G,h=load();C=W@Q;walk=json.loads((ROOT/WALK).read_text());corners=walk['corner_orbits']
 catalog=[];seen=set();literal=set();component_count=0
 for c in corners:
  for ei,edge in enumerate(c['edge_orbits']):
   E=edge['edge_walls'];cs=comps(E,G);a4=[cc for cc in cs if A4(cc,G)]
   if tuple(E)not in literal:literal.add(tuple(E));component_count+=len(a4)
   key=tuple(sorted([c['id'],edge['endpoint_orbit']]))
   if not a4 or key in seen:continue
   seen.add(key);catalog.append({'corner':c['id'],'edge_index':ei,'endpoint_orbits':list(key),'types':[c['type'],corners[edge['endpoint_orbit']]['type']], 'component':a4[0],'component_sizes':sorted(map(len,cs))})
 # Nine invariant-distinct endpoint-orbit pairs, including all endpoint types.
 selected=[];signatures=set()
 for types in [('ideal','ideal'),('ideal','finite'),('finite','finite')]:
  picks=[r for r in catalog if tuple(r['types'])==types]
  for r in picks:
   sig=(types,tuple(r['component_sizes']))
   if sig in signatures:continue
   signatures.add(sig);selected.append(r)
   if sum(tuple(x['types'])==types for x in selected)>= (1 if types[1]=='ideal' else 4):break
 assert len(selected)==9
 qe=(1+math.sqrt(5))/4;gamma=qe-.8;rows=[]
 allowed=np.array([math.cos(math.pi/m)**2 for m in range(2,101)])
 for si,seed in enumerate(selected):
  c=corners[seed['corner']];edge=c['edge_orbits'][seed['edge_index']];E=edge['edge_walls'];cc=seed['component'];A=G[np.ix_(E,E)].astype(float);I=np.linalg.inv(A)
  v=np.array(c['primitive_ray'],float);v/=-(v@Q@h);b=np.array(edge['endpoint_primitive_ray'],float);b/=-(b@Q@h)
  ends=[i for i in cc if sum(G[i,j]!=0 for j in cc)-1==1];assert len(ends)==2
  for anchor in ends:
   ai=E.index(anchor);w=I[ai]@W[E];assert abs(w@Q@w-.8)<1e-8
   B=A.copy();B[ai,ai]-=1/qe;BI=np.linalg.inv(B);assert np.linalg.eigvalsh(B)[0]>1e-5
   for t in [1/3,1/2,2/3]:
    P=(1-t)*v+t*b;qp=P@Q@P;z=b-(b@Q@P/qp)*P;qz=z@Q@z;assert qp<0 and qz>0
    pp=C@P;assert np.max(np.abs(pp[E]))<1e-8 and max(pp[i]for i in range(210)if i not in E)<-1e-7
    for sign in [-1,1]:
     normal=w+sign*math.sqrt(gamma/qz)*z;ep=C@normal;assert abs(normal@Q@normal-qe)<1e-7 and abs(normal@Q@P)<1e-8
     projected=W[E]-np.outer(np.eye(20)[ai]/qe,normal);D=BI@projected
     row={'seed':si,'anchor':anchor,'parameter':str({1/3:'1/3',1/2:'1/2',2/3:'2/3'}[t]),'radical_sign':sign,'normal':normal.tolist(),'initial_point':P.tolist(),'initial_norm':float(qp)}
     bads=[];neighbors=[]
     for k,release in enumerate(E):
      dp=C@D[k];bounds=[(pp[j]/dp[j],j)for j in range(210)if j not in E and dp[j]<-1e-9]
      if not bounds:continue
      tt,j=min(bounds);x=P-tt*D[k];xp=C@x;xx=x@Q@x
      act=np.flatnonzero(np.abs(xp)<1e-7).tolist()
      if xx>=-1e-8 or len(act)!=20 or max(xp)>1e-7:continue
      HA=G[np.ix_(act,act)]-np.outer(ep[act],ep[act])/qe
      if np.linalg.eigvalsh(HA)[0]<1e-7:continue
      pairs=[]
      for ia,ib in itertools.combinations(range(20),2):
       p=HA[ia,ib];cs=p*p/(HA[ia,ia]*HA[ib,ib])
       if p>1e-7:pairs.append((p,'obtuse',act[ia],act[ib],cs))
       elif p<-1e-7 and np.min(abs(allowed-cs))>1e-5:
        if cs<.25:pairs.append((min(cs,.25-cs),'between_0_and_1/4',act[ia],act[ib],cs))
        elif cs<.5:pairs.append((min(cs-.25,.5-cs),'between_1/4_and_1/2',act[ia],act[ib],cs))
        else:pairs.append((float(np.min(abs(allowed-cs))),'other_numeric_bad',act[ia],act[ib],cs))
      neighbors.append({'release':release,'new_wall':j,'active':act,'norm':float(xx),'maximum_other_pairing':float(max(xp[i]for i in range(210)if i not in act))})
      if pairs:
       best=max(pairs);bads.append({'release':release,'new_wall':j,'active':act,'norm':float(xx),'bad_pair':list(best[2:4]),'kind':best[1],'cosine_squared':float(best[4]),'score':float(best[0])})
     row['checked_finite_neighbors']=len(neighbors)
     if bads:row['status']='NUMERICAL_ACTUAL_FINITE_VERTEX_DEFECT';row['witness']=max(bads,key=lambda x:x['score'])
     else:row['status']='UNRESOLVED_AFTER_FIRST_FINITE_NEIGHBORS';row['neighbors']=neighbors
     rows.append(row)
  print('seed',si,'cases',len(rows),'unresolved',sum(r['status'].startswith('UNRESOLVED')for r in rows),flush=True)
 out={'status':'BOUNDED_NUMERICAL_A4_TO_H4_EDGE_SECTION_PILOT','sources':{p:sha(p)for p in[RAW,WALK,'scripts/'+Path(__file__).name]},'catalog_summary':{'unique_literal_edge_supports':len(literal),'A4_components':component_count,'distinct_endpoint_orbit_pairs_with_A4':len(catalog)},'selected_edges':selected,'cases':rows,'counts':dict(collections.Counter(r['status']for r in rows)),'scope':'Nine distinct unordered source endpoint-orbit pairs; one A4 component each; both endpoints; positions1/3,1/2,2/3; both radical signs. Numerical actual-vertex proposals only. No all-parameter or complete edge-orbit theorem; no finite-volume Coxeter candidate claimed.'}
 (ROOT/'results'/f'{TAG}_scout.json').write_text(json.dumps(out,indent=2)+'\n');print(out['counts'])
if __name__=='__main__':run()
