#!/usr/bin/env python3
"""Independent SymPy reconstruction of the two vertex-free-face caps."""
import argparse,json
from hashlib import sha256
from pathlib import Path
from itertools import combinations
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_vertex_free_caps_20261015_review'
BASE='results/fresh_h9_quartet_polar_second_20261014.json'
CAPS='results/fresh_h9_vertex_free_caps_20261015.json'
FACES='results/fresh_h9_vertex_free_edges_20261015.json'
K,Z,O,C=R.K,R.Z,R.O,R.C
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def scalar(x):return list(map(str,x.to_list()))
def vec(x):return[scalar(v)for v in x]
def pd(M):
 a=[r[:]for r in M]
 for k in range(len(a)):
  p=a[k][k]
  if R.sg(p)<=0:return False
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/p
 return True

def run():
 paths=[BASE,CAPS,FACES,'results/fresh_h9_quartet_polar_extension_20261014_review.json']
 b,d,f,oldreview=[json.loads((ROOT/p).read_text())for p in paths]
 for doc in(b,d,f,oldreview):
  for p,h in doc.get('sources',{}).items():assert sha(p)==h,p
 assert len(d['caps'])==2
 Q,N,H,E=R.dm(b['form']),R.dm(b['normals']),R.dv(b['center']),R.dv(b['strict_escape']['point'])
 pair=lambda x,y:C.pair(x,y,Q)
 NC=[C.mv(Q,n)for n in N];norms=[pair(n,n)for n in N]
 phi=(O+R.BASIS[1])/2;targets={2:Z,3:-O,5:-phi}
 caps=[];records=[];paircounts={}
 for row in d['caps']:
  contacts=dict(row['finite_contacts']);basis=[]
  # Reverse the supplied contact order, using independent row elimination.
  for i in reversed(list(contacts)):
   assert norms[i]==2*O
   trial=basis+[i]
   if len(C.eliminate([NC[j]for j in trial])[1])==len(trial):basis=trial
   if len(basis)==10:break
  assert len(basis)==10
  n=C.mv(C.inverse([NC[i]for i in basis]),[targets[contacts[i]]for i in basis])
  assert n==R.dv(row['normal'])and pair(n,n)==2*O
  vals=[C.dot(n,c)for c in NC]
  assert vals==R.dv(row['pairings'])and all(vals[i]==targets[m]for i,m in contacts.items())
  for i,g in enumerate(vals):
   t,m,u=R.relation(g,2*O,norms[i]);assert t!='forbidden'
   saved=row['relations'][i];assert(t,m)==(saved['kind'],saved['label'])
   assert u==(None if saved['invariant']is None else R.dec(saved['invariant']))
   key=t+(str(m)if m else'');paircounts[key]=paircounts.get(key,0)+1
  assert pair(n,H)==R.dec(row['center_pairing'])and R.sg(pair(n,H))<0
  assert pair(n,E)==R.dec(row['old_escape_pairing'])and R.sg(pair(n,E))>0
  S=f['faces'][row['source_face']]['support'];saved=row['new_finite_endpoint']
  assert saved['support']==S and len(S)==8
  # Direct nullspace reconstruction uses the opposite row order to the producer.
  V=C.null_one([C.mv(Q,n)]+[NC[i]for i in reversed(S)])
  if R.sg(pair(V,H))>0:V=[-x for x in V]
  sv=R.dv(saved['point']);assert R.canon(V)==R.canon(sv)
  assert R.sg(pair(sv,sv))<0 and R.sg(pair(sv,H))<0 and not pair(sv,n)
  vp=[C.dot(sv,c)for c in NC];assert vp==R.dv(saved['old_wall_pairings'])
  assert all(R.sg(x)==(0 if i in S else-1)for i,x in enumerate(vp))
  support=[N[i]for i in S]+[n];sg=[[pair(x,y)for y in support]for x in support]
  assert pd(sg)
  inverse=C.inverse(sg);assert all(R.sg(x)>=0 for rr in inverse for x in rr)
  caps.append(n);records.append(dict(source_face=row['source_face'],support=S,
   independent_reconstruction_basis=basis,normal=vec(n),endpoint=vec(sv),
   endpoint_norm=scalar(pair(sv,sv)),old_escape_pairing=scalar(pair(n,E)),
   physical_endpoint_zero_set=S+[25+len(caps)-1]))
 assert not pair(caps[0],caps[1])
 assert R.dec(d['mutual_pair']['pairing'])==Z and d['mutual_pair']['kind']=='right'and d['mutual_pair']['label']==2
 paircounts['right']=paircounts.get('right',0)+1
 allN=N+caps;allcov=NC+[C.mv(Q,n)for n in caps]
 for k,row in enumerate(d['caps']):
  V=R.dv(row['new_finite_endpoint']['point']);values=[C.dot(V,c)for c in allcov]
  ids=row['new_finite_endpoint']['support']+[25+k]
  assert[i for i,x in enumerate(values)if not x]==sorted(ids)
  assert all(R.sg(x)==(0 if i in ids else-1)for i,x in enumerate(values))
  records[k]['mate_pairing']=scalar(values[26-k])
 assert R.canon(R.dv(d['caps'][0]['new_finite_endpoint']['point']))!=R.canon(R.dv(d['caps'][1]['new_finite_endpoint']['point']))
 # Explicit actual-facet checks for the full27-wall extension.
 assert R.sg(pair(H,H))<0 and all(R.sg(C.dot(H,c))<0 for c in allcov)
 for i,n in enumerate(allN):
  q=pair(n,n);h=C.dot(H,allcov[i]);P=[a-h/q*b for a,b in zip(H,n)]
  assert R.sg(pair(P,P))<0 and R.sg(pair(P,H))<0
  assert all(R.sg(C.dot(P,c))==(0 if j==i else-1)for j,c in enumerate(allcov))
 paths += [f'scripts/{TAG}.py','scripts/fresh_h9_vertex_free_caps_20261015.py',
  'scripts/fresh_h9_quartet_polar_extension_20261014_review.py',
  'scripts/fresh_h9_h8_all_high_20261014_review.py','scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
 return dict(status='PASS_INDEPENDENT_TWO_VERTEX_FREE_ESCAPE_CUTTERS',sources={p:sha(p)for p in paths},
  arithmetic='Separate SymPy ANP field/192bit outward signs; reversed contact frames and physical nullspaces, no producer arithmetic imports.',
  caps=records,new_pair_counts=paircounts,new_pairs=51,actual_facets=27,facet_sign_checks=729,
  actual_distinct_new_endpoints=2,endpoint_sign_checks=54,
  scope='Both norm2 caps retain the old25wall chamber and common center, are mutuallyright, cut the oldescape, and give distinct actual finite endpoints on the specified formerly vertex-free lines. Full27wall cone geometry and any renewedescape remain separate tests; finitevolume is not inferred.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved independent two-cap replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
