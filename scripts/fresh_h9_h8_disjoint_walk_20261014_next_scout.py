#!/usr/bin/env python3
"""Bounded71-germ scout using100-digit physical pairings before float events."""
import argparse,gzip,json,math
from collections import Counter
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import numpy as np
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_walk_20261014_next_scout'
FILES=['results/fresh_h9_h8_disjoint_walk_20261014_exact.json','results/fresh_h9_h8_disjoint_walk_20261014_next_frontier.json']
ANCHORS=[2,3,4,5,6,7,8,9,10,11,12,16,20,24,32,64];TAILS=[6,7,8,10,12,16,24,32,64]
mp.mp.dps=100;S5=mp.sqrt(5);RAD=6+16*S5/5

def scalar(x):return mp.mpf(x[0][0])+mp.mpf(x[0][1])*S5+(mp.mpf(x[1][0])+mp.mpf(x[1][1])*S5)*mp.sqrt(RAD)
def qscalar(x):return mp.mpf(x[0])+mp.mpf(x[1])*S5

def intervals(coefs):
 cuts=[1.]
 for a,b,c in coefs:
  assert a<0
  disc=b*b-4*a*c
  if disc>=0:
   for z in[(-b-math.sqrt(disc))/(2*a),(-b+math.sqrt(disc))/(2*a)]:
    if z>1 and math.isfinite(z):cuts.append(z)
 cuts=sorted(set(cuts));out=[]
 for lo,hi in zip(cuts,cuts[1:]+[math.inf]):
  if hi-lo<1e-8:continue
  z=(lo+hi)/2 if math.isfinite(hi)else max(2*lo,lo+1)
  if all(a*z*z+b*z+c< -1e-8 for a,b,c in coefs):out.append((lo,hi,z))
 return out

def run():
 d,f=[json.loads((ROOT/p).read_text())for p in FILES]
 Q=mp.matrix([[scalar(x)for x in row]for row in d['form']]);N=[mp.matrix([scalar(x)for x in row])for row in d['normals']]
 E=mp.matrix([scalar(x)for x in d['strict_escape']['point']]);norms=[scalar(x)for x in d['norms']]
 scales=[mp.sqrt(2/x)for x in norms];G=np.array([[float(scalar(x)*scales[i]*scales[j])for j,x in enumerate(row)]for i,row in enumerate(d['gram'])])
 NE=np.array([float((n.T*Q*E)[0]*scales[i])for i,n in enumerate(N)])
 fvals=-2*np.cos(np.pi/np.arange(2,257));counts=Counter();found={};stats=[]
 for row in f['queue']:
  eid=row['id'];S=row['support'];outside=[i for i in range(19)if i not in S];pd=row['profile_data']
  J=np.array([[float(qscalar(x))for x in rr]for rr in pd['normalized_inverse']])
  P=mp.matrix([scalar(x)for x in row['start_point']]);D=mp.matrix([scalar(x)for x in row['direction']])
  pq=scalar(row['norm_polynomial'][0]);dq=scalar(row['norm_polynomial'][2]);assert pq<0 and dq>0
  PN=[(P.T*Q*n)[0]*scales[i]for i,n in enumerate(N)];DN=[(D.T*Q*n)[0]*scales[i]for i,n in enumerate(N)]
  up=[PN[i]/mp.sqrt(-pq)+DN[i]/mp.sqrt(dq)for i in range(19)]
  vm=[DN[i]/mp.sqrt(dq)-PN[i]/mp.sqrt(-pq)for i in range(19)]
  assert all(up[i]<0 for i in outside)
  up=np.array([float(x)for x in up]);vm=np.array([float(x)for x in vm]);pn=np.array([float(x)for x in PN]);dn=np.array([float(x)for x in DN])
  PE=(P.T*Q*E)[0];DE=(D.T*Q*E)[0]
  ep=float(PE/mp.sqrt(-pq)+DE/mp.sqrt(dq));em=float(DE/mp.sqrt(dq)-PE/mp.sqrt(-pq));nullt=float(mp.sqrt(-pq/dq))
  profiles=[r['contacts']for r in pd['finite_profiles']]+[[[r['wall'],m]]for r in pd['dihedral_tail_profiles']for m in TAILS]
  before=len(found);accepted=0
  for contacts in profiles:
   c=np.zeros(8)
   for i,m in contacts:c[S.index(i)]=-2*math.cos(math.pi/m)
   rho=2-c@J@c;assert rho>0
   b=c@J@G[S];a=math.sqrt(rho)/2*up;cc=math.sqrt(rho)/2*vm;pe=float(c@J@NE[S]);sr=math.sqrt(rho)/2
   def record(z,recipe):
    nonlocal accepted
    if not math.isfinite(z)or z<=1+1e-9:return
    vals=a*z+b+cc/z
    err=np.min(abs(vals[:,None]-fvals[None,:]),axis=1)
    if not np.all((vals< -2-1e-7)|(err<1e-8)):return
    t=nullt*(z*z-1)/(z*z+1)
    ev=pn+t*dn
    if max(ev[outside])>=-1e-8:return
    # q(P+tD)<0 and forward time are exact family facts for1<z<infinity.
    en=float(pq)+t*t*float(dq)
    if en>=-1e-8:counts['near_null_endpoint_skips']+=1;return
    ee=pe+sr*(ep*z+em/z);accepted+=1
    rec=dict(edge_id=eid,contacts=contacts,**recipe,exponential_parameter=float(z),
       normalized_wall_pairings=vals.tolist(),endpoint_forward_parameter=t,endpoint_norm=en,
       endpoint_other_pairings=ev[outside].tolist(),escape_pairing=float(ee),cuts_escape=bool(ee>1e-8),
       finite_neighbors=[int(i)for i in range(19)if vals[i]>-2+1e-7])
    key=tuple(np.round(vals,8))
    if key not in found:found[key]=rec
   ints=intervals(list(zip(a[outside],b[outside]+2,cc[outside])))
   counts['all_dotted_parameter_intervals']+=len(ints)
   for lo,hi,z in ints:record(z,dict(method='all_dotted_interval',interval=[lo,hi if math.isfinite(hi)else None]))
   for anchor in outside:
    aa=a[anchor];assert aa<0
    for m in ANCHORS:
     bb=b[anchor]+2*math.cos(math.pi/m);c0=cc[anchor];dis=bb*bb-4*aa*c0
     if dis<0:continue
     for branch in[-1,1]:
      z=(-bb+branch*math.sqrt(dis))/(2*aa);counts['real_anchor_roots']+=1
      record(z,dict(method='finite_anchor',anchor=anchor,anchor_label=m,quadratic_branch=branch))
  counts['sampled_profiles']+=len(profiles)
  stats.append(dict(edge_id=eid,profiles=len(profiles),accepted_events=accepted,new_distinct_caps=len(found)-before))
  print(stats[-1],'total',len(found),'cutters',sum(x['cuts_escape']for x in found.values()),flush=True)
 bank=sorted(found.values(),key=lambda x:(not x['cuts_escape'],-x['escape_pairing'],x['edge_id']))
 return dict(status='NUMERICAL_NEXT_NINE_VERTEX_ENDPOINT_PROPOSALS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in FILES+[f'scripts/{TAG}.py']},
   counts=dict(counts),finite_anchor_labels=ANCHORS,sampled_dihedral_tail_labels=TAILS,edges=stats,candidate_count=len(bank),escape_cutter_count=sum(r['cuts_escape']for r in bank),candidates=bank,
   scope='Bounded numerical events on71exact open germs from9new vertices; previous56open germs not repeated. Physical scalar pairings computed at100digits before float quadratic events. Global finite-angle gate labels2..256 with1e-8tolerance; dotted gap1e-7. Symbolic tails sampled only at listed labels. Near-null endpoints with float norm margin below1e-8are skipped and counted. No all-label/all-parameter exclusion or exact cap claim.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();bank=d.pop('candidates')
 bp=ROOT/f'data/enumerations/{TAG}_bank.json.gz';blob=gzip.compress(json.dumps(bank,separators=(',',':')).encode(),mtime=0)
 d['candidate_bank']=dict(path=str(bp.relative_to(ROOT)),sha256=sha256(blob).hexdigest(),count=len(bank));d['best_escape_pairing_samples']=bank[:24]
 p=ROOT/f'results/{TAG}.json'
 if args.verify:assert blob==bp.read_bytes()and d==json.loads(p.read_text())
 else:bp.write_bytes(blob);p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['candidate_count'],d['escape_cutter_count'],flush=True)
