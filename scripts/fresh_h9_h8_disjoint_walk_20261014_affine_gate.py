#!/usr/bin/env python3
"""All-real/all-label escape exclusion for14proper ideal endpoint profiles."""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
import fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate as B
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_disjoint_walk_20261014_affine_gate'
FILES=['results/fresh_h9_h8_disjoint_cap_20261014.json','results/fresh_h9_h8_disjoint_walk_20261014.json','results/fresh_h9_h8_disjoint_walk_20261014_affine.json']
I,ei,qi=B.I,B.ei,B.qi
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d,f,p=[json.loads((ROOT/x).read_text())for x in FILES]
 for doc in[d,f,p]:
  for path,h in doc['sources'].items():assert sha(path)==h,path
 Q,N,H,E=A.dm(d['form']),A.dm(d['normals']),A.dv(d['center']),A.dv(d['strict_escape']['point'])
 norms=[A.pair(n,n,Q)for n in N];scales=[(2/ei(x)).sqrt()for x in norms]
 G=[[ei(A.pair(n,m,Q))*scales[i]*scales[j]for j,m in enumerate(N)]for i,n in enumerate(N)]
 GE=[ei(A.pair(E,n,Q))*scales[i]for i,n in enumerate(N)]
 assert A.pair(E,E,Q).sign()>0 and A.pair(E,H,Q).sign()<0 and all(x.hi<0 for x in GE)
 records=[];total=Counter()
 for pr in p['supports']:
  if not pr['affine_profiles']:continue
  row=f['queue'][pr['germ']];S=pr['support'];assert S==row['support']
  J=[[A.Q(*x)for x in rr]for rr in row['profile_data']['normalized_inverse']]
  M=[[A.Q(*x)for x in rr]for rr in row['profile_data']['normalized_Gram']]
  for i,ii in enumerate(S):
   for j,jj in enumerate(S):
    assert sum((M[i][k]*J[k][j]for k in range(8)),A.Q())==A.Q(i==j)
    raw=A.pair(N[ii],N[jj],Q);want=A.E(M[i][j]);assert raw.sign()==want.sign()and 4*raw*raw==want*want*norms[ii]*norms[jj]
  outside=[j for j in range(15)if j not in S]
  X,Y=A.dv(row['start_point']),A.dv(row['direction']);xx=-A.pair(X,X,Q);yy=A.pair(Y,Y,Q)
  assert xx.sign()>0 and yy.sign()>0 and not A.pair(X,Y,Q)
  assert A.pair(X,H,Q).sign()<0 and all(not A.pair(Z,N[i],Q)for Z in[X,Y]for i in S)
  dx,dy=ei(xx).sqrt(),ei(yy).sqrt()
  aa=[ei(A.pair(X,N[j],Q))*scales[j]/dx for j in outside]+[ei(A.pair(X,E,Q))/dx]
  cc=[ei(A.pair(Y,N[j],Q))*scales[j]/dy for j in outside]+[ei(A.pair(Y,E,Q))/dy]
  assert aa[-1].hi<0<cc[-1].lo and(aa[-1]+cc[-1]).hi<0
  assert all(a.hi<=0 and(a+c).hi<0 for a,c in zip(aa,cc))
  cols=[[G[i][j]for i in S]for j in outside]+[[GE[i]for i in S]]
  edge=dict(edge_id=row['id'],outside=outside,a=aa,c=cc)
  for prof in pr['affine_profiles']:
   contact={S.index(i):m for i,m in prof['contacts']};w0=[];w1=[]
   for i in range(8):
    w0.append(sum((J[i][j]for j,m in contact.items()if m==3),A.Q()))
    w1.append(sum((J[i][j]for j,m in contact.items()if m==4),A.Q()))
    assert w0[-1].sign()>=0 and w1[-1].sign()>=0 and(w0[-1]or w1[-1])
    assert [w0[-1].enc(),w1[-1].enc()]==prof['positive_kernel_first8_as_a_plus_b_sqrt2'][i]
   cost0=sum((w0[i]if m==3 else 2*w1[i]for i,m in contact.items()),A.Q())
   cost1=sum((w1[i]if m==3 else w0[i]for i,m in contact.items()),A.Q())
   assert cost0==2 and not cost1
   weights=[qi(a)+B.SQ2*qi(b)for a,b in zip(w0,w1)]
   bb=[-B.isdot(weights,col)for col in cols]
   # The inherited one-parameter engine handles rho=0without alteration.
   result=B.certify(edge,dict(kind='finite',contacts=prof['contacts'],b=bb,rho=I(0)))
   result['kind']='affine';total.update(result['counts']);total['profiles']+=1
   result.update(germ=row['id'],outside=outside,a=[x.enc()for x in aa],b=[x.enc()for x in bb],c=[x.enc()for x in cc])
   records.append(result);print('germ',row['id'],'contacts',prof['contacts'],'boxes',result['counts']['boxes'],flush=True)
 assert len(records)==14
 paths=FILES+[f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py']
 return dict(status='PASS_ALL_REAL_AFFINE_ENDPOINT_PROFILES_MISS_ESCAPE',sources={p:sha(p)for p in paths},
   precision_bits=B.BITS,counts=dict(total),records=records,maximum_depth=max(r['maximum_depth']for r in records),
   parameterization='r=p+h(U+V),h>=0,p²=2,p=minus support dual of normalized nonright contacts;U²=-1,V²=1,UV=0.',
   allowed_pairing_relaxation='Labels2..15isolated; all integer labels>=16,parallel,dotted retained in(-infinity,-2cos(pi/16)].',
   cosine_enclosures=dict(pi=B.PI.enc(),isolated=[x.enc()for x in B.ALLOWED_POINTS],tail_threshold=B.HIGH_THRESHOLD.enc()),
   scope='Every one-wall bounded ideal endpoint completion of one of the original15-wall seed62open finite-corner germs strictly misses its literal escape. Complete14local affine profiles are inherited, with weights/cost/positivekernel rechecked. All real nonnegative null coefficients and all integer old-wall labels are covered. Does not exclude later71supports, intermediate walls, ideal ends with unbounded links, or other cap directions. Interval engine is reused, not independently implemented.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert d==json.loads(p.read_text())
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['counts'],flush=True)
