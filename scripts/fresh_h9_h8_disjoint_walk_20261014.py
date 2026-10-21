#!/usr/bin/env python3
"""Exact 62 open finite-corner germs of the changed H8-based H9 seed."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
import fresh_h9_mixed_d4_endpoint_profiles_20261014 as P
ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_h9_h8_disjoint_cap_20261014';TAG='fresh_h9_h8_disjoint_walk_20261014'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
 d=json.loads((ROOT/f'results/{BASE}.json').read_text())
 for p,h in d['sources'].items():assert sha(p)==h,p
 Q,N,H=A.dm(d['form']),A.dm(d['normals']),A.dv(d['center'])
 G=A.dm(d['gram']);V=d['complete_finite_vertices'];edges=d['complete_finite_edges']
 closed={(v,tuple(row['support']))for row in edges for v in row['vertices']}
 labels={tuple(row['pair']):row['label'] for row in d['pairs'] if row['type']!='dotted'}
 queue=[];counts=Counter();profile_total=tail_total=0;types=Counter()
 for vi,row in enumerate(V):
  ids=sorted(row['incident']);v=A.dv(row['point']);openj=[j for j in ids if (vi,tuple(k for k in ids if k!=j))not in closed]
  if not openj:continue
  GI=A.inverse([[G[i][j]for j in ids]for i in ids])
  for j in openj:
   S=[i for i in ids if i!=j];col=ids.index(j)
   D=[-sum((GI[a][col]*N[k][b]for a,k in enumerate(ids)),A.Z)for b in range(10)]
   assert not A.pair(v,D,Q)
   dq=A.pair(D,D,Q);vq=A.pair(v,v,Q);assert dq.sign()>0 and vq.sign()<0
   assert [A.pair(D,N[k],Q)for k in ids]==[A.E(-1 if k==j else 0)for k in ids]
   vp=[A.pair(v,n,Q)for n in N];dp=[A.pair(D,n,Q)for n in N]
   upper=[]
   for k in range(15):
    if k in S:assert not vp[k]and not dp[k];continue
    if k==j:assert not vp[k]and dp[k]==-1;continue
    assert vp[k].sign()<0
    if dp[k].sign()>0:
     t=-vp[k]/dp[k];assert (t*t*dq+vq).sign()>0
     upper.append((t,k))
   if upper:
    tmax,k=upper[0]
    for t,k0 in upper[1:]:
     if (t-tmax).sign()<0:tmax,k=t,k0
    for p in range(1,80):
     t=tmax*A.E(F(2**p-1,2**p))
     if (vq+t*t*dq).sign()>0:break
    else:raise AssertionError('no spacelike interior')
   else:
    tmax=None;k=None;t=A.ONE
    while (vq+t*t*dq).sign()<=0:t*=2
   w=[a+t*b for a,b in zip(v,D)];vals=[A.pair(w,n,Q)for n in N]
   assert all(x.sign()==(0 if i in S else -1)for i,x in enumerate(vals))
   assert A.pair(w,H,Q).sign()<0 and A.pair(w,w,Q).sign()>0
   ll={(a,b):labels[tuple(sorted((S[a],S[b])))]for a in range(8)for b in range(a+1,8)}
   profiles=P.profiles(list(range(8)),ll)
   for rowp in profiles['finite_profiles']:
    rowp['contacts']=[[S[i],m]for i,m in rowp['contacts']]
   for rowp in profiles['dihedral_tail_profiles']:rowp['wall']=S[rowp['wall']]
   profiles['support']=S
   profiles['support_components']=[[S[i]for i in c]for c in profiles['support_components']]
   desc='+'.join(sorted(profiles['support_component_types']));types[desc]+=1
   profile_total+=len(profiles['finite_profiles']);tail_total+=len(profiles['dihedral_tail_profiles'])
   queue.append(dict(id=len(queue),start_vertex=vi,released_wall=j,support=S,start_point=v,
      direction=D,norm_polynomial=[vq,A.Z,dq],first_old_wall_parameter=tmax,first_old_wall=k,
      spacelike_parameter=t,spacelike_point=w,spacelike_norm=A.pair(w,w,Q),
      spacelike_pairings=vals,profile_data=profiles))
   counts['bottom' if vi<41 else 'top']+=1
   print('germ',len(queue),'vertex',vi,'release',j,'type',desc,flush=True)
 assert len(queue)==62 and dict(counts)=={'bottom':38,'top':24}
 assert len({tuple(row['support'])for row in queue})==62
 paths=[f'results/{BASE}.json',f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py','scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014.py']
 return A.enc(dict(status='PASS_EXACT_62_OPEN_FINITE_CORNER_GERMS',sources={p:sha(p)for p in paths},
   complete_finite_vertex_dependency=44,complete_finite_edge_dependency=167,
   open_count=62,origin_counts=dict(counts),support_types=dict(types),
   finite_profile_count=profile_total,symbolic_dihedral_tail_count=tail_total,queue=queue,
   scope='Complete open edge germs incident to the44finite vertices. Every opposite light-cone limit has exactly its8PD support walls active, hence an unbounded link. No census of edges without finite endpoints. Complete local spherical profiles plus symbolic I2 tails; global cap feasibility and affine endpoint extensions remain separate.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/f'results/{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['finite_profile_count'],out['symbolic_dihedral_tail_count'],flush=True)
