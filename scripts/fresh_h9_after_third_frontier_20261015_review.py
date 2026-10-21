#!/usr/bin/env python3
"""Independent SymPy review of the28-wall new germs and vertex-free lines."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations,product
from pathlib import Path
import fresh_h9_quartet_polar_extension_20261014_review as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_after_third_frontier_20261015_review'
BASE='results/fresh_h9_replacement_escape_extension_20261015.json'
LINES='results/fresh_h9_after_third_vertex_free_edges_20261015.json'
FRONT='results/fresh_h9_after_third_endpoint_frontier_20261015.json'
PRIOR='results/fresh_h9_replacement_escape_extension_20261015_review.json'
K,Z,O,C=R.K,R.Z,R.O,R.C
SQ5=R.BASIS[1]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def q5(v):return K.convert(F(v[0]))+K.convert(F(v[1]))*SQ5
def small(x):return K.convert(F(x))
def scalar(v):return list(map(str,v.to_list()))
def pd(M):
 a=[r[:]for r in M]
 for k in range(len(a)):
  p=a[k][k]
  if R.sg(p)<=0:return False
  for i in range(k+1,len(a)):
   for j in range(k+1,len(a)):a[i][j]-=a[i][k]*a[k][j]/p
 return True

def components(nodes,es):
 adj={i:set()for i in nodes}
 for (i,j),m in es.items():
  if m>2:adj[i].add(j);adj[j].add(i)
 out=[];left=set(nodes)
 while left:
  todo=[min(left)];seen=set(todo)
  for i in todo:
   for j in sorted(adj[i]-seen):seen.add(j);todo.append(j)
  out.append(sorted(seen));left-=seen
 return out,adj

def typ(nodes,es):
 n=len(nodes)
 if n==1:return'A1'
 ed={p:m for p,m in es.items()if m>2 and set(p)<=set(nodes)}
 cs,adj=components(nodes,ed);assert len(cs)==1 and len(ed)==n-1
 if n==2:return{3:'A2',4:'B2'}.get(next(iter(ed.values())),'I2(%s)'%next(iter(ed.values())))
 branch=[i for i in nodes if len(adj[i])>2]
 if branch:
  assert len(branch)==1 and len(adj[branch[0]])==3 and set(ed.values())=={3}
  lengths=[]
  for j in adj[branch[0]]:
   prev,cur,k=branch[0],j,1
   while adj[cur]-{prev}:
    nxt=next(iter(adj[cur]-{prev}));prev,cur,k=cur,nxt,k+1
   lengths.append(k)
  lengths.sort()
  if lengths[:2]==[1,1]:return'D%d'%n
  assert(n,tuple(lengths))in[(6,(1,2,2)),(7,(1,2,3)),(8,(1,2,4))]
  return'E%d'%n
 cur=min(i for i in nodes if len(adj[i])==1);prev=None;ls=[]
 while adj[cur]-({prev}if prev is not None else set()):
  nxt=next(iter(adj[cur]-({prev}if prev is not None else set())));ls.append(ed[tuple(sorted((cur,nxt)))]);prev,cur=cur,nxt
 if set(ls)=={3}:return'A%d'%n
 heavy=[(i,m)for i,m in enumerate(ls)if m!=3];assert len(heavy)==1
 i,m=heavy[0]
 if m==4 and i in(0,n-2):return'B%d'%n
 if m==4 and n==4 and i==1:return'F4'
 assert m==5 and n in(3,4)and i in(0,n-2)
 return'H%d'%n

def run():
 paths=[BASE,LINES,FRONT,PRIOR];docs=[]
 for path in paths:
  d=json.loads((ROOT/path).read_text())
  for p,h in d.get('sources',{}).items():assert sha(p)==h,p
  docs.append(d)
 b,line,front,prior=docs
 assert prior['status']=='PASS_INDEPENDENT_FULL_28_WALL_CONE'
 assert prior['sources'][BASE]==sha(BASE)
 assert b['field']['new_positive_radical_square']==['6','16/5']
 Q,N,H=R.dm(b['form']),R.dm(b['normals']),R.dv(b['center'])
 assert len(N)==28
 NC=[C.mv(Q,n)for n in N]
 G=[[C.dot(n,c)for c in NC]for n in N];assert G==R.dm(b['gram'])
 pair=lambda x,y:C.pair(x,y,Q)
 labels={}
 for i,j in combinations(range(28),2):
  t,m,_=R.relation(G[i][j],G[i][i],G[j][j]);assert t!='forbidden'
  if t in('right','finite'):labels[i,j]=m
 @lru_cache(None)
 def support_pd(S):return pd([[G[i][j]for j in S]for i in S])
 rays=[dict(active=a)for a in b['complete_census']['hyperideal_rays']]
 assert len(rays)==727 and b['counts']['cusps']==0 and b['counts']['opposite_sheet_rays']==0
 masks=[sum(1<<i for i in r['active'])for r in rays]
 @lru_cache(None)
 def exterior_ray(i):
  ids=rays[i]['active'];v=C.null_one([NC[j]for j in ids])
  vals=[C.dot(v,c)for c in NC]
  if any(R.sg(z)>0 for z in vals):v=[-z for z in v];vals=[-z for z in vals]
  assert all(R.sg(z)<=0 for z in vals)and[k for k,z in enumerate(vals)if not z]==ids
  assert R.sg(pair(v,v))>0
  return v
 expected=[];counts=Counter()
 for i,j in combinations(range(len(rays)),2):
  common=masks[i]&masks[j]
  if common.bit_count()<8:continue
  counts['shared_rank8_candidates']+=1
  ids=[k for k in range(28)if common>>k&1]
  candidates=[S for S in combinations(ids,8)if all(p in labels for p in combinations(S,2))]
  if not candidates:counts['not_spherical_clique']+=1;continue
  positive=[S for S in candidates if support_pd(S)]
  if not positive:counts['not_positive_definite']+=1;continue
  assert len(ids)==8 and positive==[tuple(ids)]
  # Exactly two complete-cone rays generate this 2D face, independently
  # checked against both exterior and finite ray active sets.
  assert[k for k,m in enumerate(masks)if m&common==common]==[i,j]
  assert not any(set(ids)<=set(v['incident'])for v in b['finite_vertices'])
  u,v=exterior_ray(i),exterior_ray(j)
  qu,qv,uv=pair(u,u),pair(v,v),pair(u,v)
  assert R.sg(qu)>0 and R.sg(qv)>0 and R.sg(uv)<0 and R.sg(qu*qv-uv*uv)<0
  p=[a-uv/qv*c for a,c in zip(u,v)];pv=pair(p,v);pn=pair(p,p)
  assert not pv and R.sg(pn)<0 and R.sg(pair(p,H))<0
  vals=[C.dot(p,c)for c in NC]
  assert[i for i,x in enumerate(vals)if not x]==ids and all(R.sg(x)<=0 for x in vals)
  # The negative null parameter lies strictly beyond the first exterior ray:
  # (uv/qv)^2-(-pn/qv)=qu/qv>0, with uv/qv<0.
  assert (uv/qv)**2+pn/qv==qu/qv and R.sg(qu/qv)>0
  saved=line['faces'][len(expected)]
  assert saved['support']==ids and saved['hyperideal_endpoints']==[i,j]
  assert R.dv(saved['first_ray'])==u and R.dv(saved['second_ray'])==v
  assert R.dv(saved['start_point'])==p and R.dv(saved['direction'])==v
  assert R.dv(saved['norm_polynomial'])==[pn,Z,qv]and R.dv(saved['start_pairings'])==vals
  expected.append(dict(support=ids,exterior_rays=[i,j],negative_norm=scalar(pn),lorentz_determinant=scalar(qu*qv-uv*uv)))
  counts['complete_geodesic_faces']+=1
 assert dict(counts)==line['counts']==dict(shared_rank8_candidates=3261,not_spherical_clique=1666,
    not_positive_definite=1524,complete_geodesic_faces=71)and len(expected)==71
 assert len({tuple(r['support'])for r in expected})==71
 print('PASS independent71vertex-free lines',flush=True)
 # Derive all open germs at the three new cap vertices from the complete87-corner list.
 vertices=b['finite_vertices'];closed=set();assert len(vertices)==87
 edges=[]
 for i,j in combinations(range(len(vertices)),2):
  common=sorted(set(vertices[i]['incident'])&set(vertices[j]['incident']))
  if len(common)==8:
   edges.append(dict(vertices=[i,j],support=common));closed.add((i,tuple(common)));closed.add((j,tuple(common)))
 assert edges==b['finite_edges']and len(edges)==249
 new=[i for i,v in enumerate(vertices)if any(k>=27 for k in v['incident'])]
 assert len(new)==3
 endpoint_types=Counter(sum(i in new for i in row['vertices'])for row in edges)
 assert endpoint_types=={0:247,1:1,2:1} and 9*84-2*247-1==261
 wanted=[(i,j,tuple(k for k in vertices[i]['incident']if k!=j))for i in new for j in vertices[i]['incident']
         if(i,tuple(k for k in vertices[i]['incident']if k!=j))not in closed]
 assert wanted==[(r['start_vertex'],r['released_wall'],tuple(r['support']))for r in front['queue']]
 assert len(wanted)==len(set(S for _,_,S in wanted))==24
 iv={};profiles_total=tails_total=inert_components=0;types=Counter();germs=[];profile_summaries=[]
 squares={3:O,4:2*O,5:(3*O+SQ5)/2}
 phi=(O+SQ5)/2
 for row in front['queue']:
  vi,j,S=row['start_vertex'],row['released_wall'],tuple(row['support']);ids=vertices[vi]['incident']
  assert support_pd(S)
  if vi not in iv:iv[vi]=C.inverse([[G[a][b]for b in ids]for a in ids])
  col=ids.index(j);J=iv[vi];v=R.dv(vertices[vi]['point'])
  D=[-sum((J[a][col]*N[k][b]for a,k in enumerate(ids)),Z)for b in range(10)]
  assert D==R.dv(row['direction'])and v==R.dv(row['start_point'])
  vq,dq=pair(v,v),pair(D,D)
  assert R.sg(vq)<0 and R.sg(dq)>0 and not pair(v,D)
  vp=[C.dot(v,c)for c in NC];dp=[C.dot(D,c)for c in NC]
  assert[k for k,x in enumerate(vp)if not x]==ids and all(R.sg(x)<=0 for x in vp)
  assert [dp[k]for k in ids]==[-O if k==j else Z for k in ids]
  assert R.dv(row['norm_polynomial'])==[vq,Z,dq]
  upper=[]
  for k in range(28):
   if k in S:assert not vp[k]and not dp[k];continue
   if k==j:continue
   assert R.sg(vp[k])<0
   if R.sg(dp[k])>0:
    t=-vp[k]/dp[k];assert R.sg(t)>0 and R.sg(vq+t*t*dq)>0;upper.append((t,k))
  if upper:
   t,k=upper[0]
   for x,i in upper[1:]:
    if R.sg(x-t)<0:t,k=x,i
   assert R.dec(row['first_old_wall_parameter'])==t and row['first_old_wall']==k
  else:assert row['first_old_wall_parameter']is None and row['first_old_wall']is None
  t=R.dec(row['spacelike_parameter']);assert R.sg(t)>0
  w=[a+t*b for a,b in zip(v,D)];vals=[C.dot(w,c)for c in NC]
  assert w==R.dv(row['spacelike_point'])and vals==R.dv(row['spacelike_pairings'])
  assert all(R.sg(x)==(0 if i in S else-1)for i,x in enumerate(vals))
  assert R.sg(pair(w,H))<0 and pair(w,w)==R.dec(row['spacelike_norm'])and R.sg(pair(w,w))>0
  # Independent normalized blocks, including the exact true inverse diagonal
  # of the inert I2(10); its off-diagonal is strictly positive, not zero.
  es={p:labels[p]for p in combinations(S,2)};cs,_=components(S,es);ts=[typ(c,es)for c in cs]
  pf=row['profile_data'];assert cs==pf['support_components']and ts==pf['support_component_types']
  types['+'.join(sorted(ts))]+=1
  Jphys=C.inverse([[G[a][b]for b in S]for a in S]);diag={i:Jphys[S.index(i)][S.index(i)]*G[i][i]/2 for i in S}
  effective=[[Z for _ in S]for _ in S];blocked=[]
  for cc,tt in zip(cs,ts):
   if tt=='I2(10)':
    assert len(cc)==2;blocked+=cc;inert_components+=1
    i,j=cc;a=3*O+SQ5
    assert diag[i]==diag[j]==a and R.sg(a-2*O)>0
    assert 4*G[i][j]*G[i][j]/(G[i][i]*G[j][j])==(5*O+SQ5)/2
    offsq=(25*O+11*SQ5)/2
    assert a*a-offsq==a/2 and R.sg(offsq)>0
    for i in cc:effective[S.index(i)][S.index(i)]=a
   else:
    M=[[2*O if i==j else {2:Z,3:-O,5:-phi}[es[tuple(sorted((i,j)))]]for j in cc]for i in cc]
    assert pd(M);inv=C.inverse(M)
    for a,i in enumerate(cc):
     assert inv[a][a]==diag[i]
     for bb,jj in enumerate(cc):effective[S.index(i)][S.index(jj)]=inv[a][bb]
  assert effective==[[q5(x)for x in rr]for rr in pf['effective_projection_inverse']]
  assert sorted(blocked)==pf['forbidden_attachment_nodes']
  assert [cc for cc,tt in zip(cs,ts)if tt=='I2(10)']==pf['inert_I2_10_components']
  assert all(R.sg(x)>=0 for rr in effective for x in rr)
  assert all(R.sg(diag[i]-small(F(1,2)))>=0 for i in S)
  where={i:k for k,c in enumerate(cs)for i in c};allow=[i for i in reversed(S)if i not in blocked]
  accepted={};reject=Counter();trials=0
  for n in range(4):
   for touched in combinations(allow,n):
    if len({where[i]for i in touched})!=n:continue
    for ms in(product((5,4,3),repeat=n)if n<3 else[(3,3,3)]):
     trials+=1;cost=sum((squares[m]*diag[i]for i,m in zip(touched,ms)),Z);rho=2*O-cost
     if R.sg(rho)<=0:reject['affine'if not rho else'indefinite']+=1;continue
     key=tuple(sorted(zip(touched,ms)));ee=dict(es);ee.update({tuple(sorted((i,28))):m for i,m in key})
     parts,_=components(list(S)+[28],ee);tt=[typ(c,ee)for c in parts]
     accepted[key]=(cost,rho,sorted(tt),tt[next(i for i,c in enumerate(parts)if 28 in c)])
  actual={tuple(map(tuple,r['contacts'])):r for r in pf['finite_profiles']}
  assert set(actual)==set(accepted)and len(actual)==len(pf['finite_profiles'])
  for key,(cost,rho,tt,nt)in accepted.items():
   rr=actual[key];assert cost==q5(rr['projection_square'])and rho==q5(rr['residual_norm'])
   assert tt==sorted(rr['endpoint_components'])and nt==rr['new_node_component']
  tails=[cc[0]for cc in cs if len(cc)==1]
  assert tails==[r['wall']for r in pf['dihedral_tail_profiles']]
  assert all(diag[i]==small(F(1,2))for i in tails)
  assert trials==pf['finite_trial_count']and dict(reject)==pf['rejected_trial_counts']
  profiles_total+=len(actual);tails_total+=len(tails)
  germs.append(dict(id=row['id'],support=list(S),start_vertex=vi,released_wall=j,spacelike_norm=scalar(pair(w,w))))
  profile_summaries.append(dict(id=row['id'],finite=len(actual),tails=len(tails),blocked=sorted(blocked),trials=trials))
  if len(germs)%20==0:print('checked germs',len(germs),flush=True)
 assert profiles_total==front['finite_profile_count']==247
 assert tails_total==front['symbolic_dihedral_tail_count']==18
 assert dict(types)==front['support_types']
 paths += [f'scripts/{TAG}.py','scripts/fresh_h9_after_third_vertex_free_edges_20261015.py',
  'scripts/fresh_h9_after_third_endpoint_frontier_20261015.py','scripts/fresh_h9_polar_endpoint_profiles_20261015.py',
  'scripts/fresh_h9_quartet_polar_extension_20261014_review.py','scripts/fresh_h9_h8_all_high_20261014_review.py',
  'scripts/fresh_h9_h8_facet_cover_20261014_line_review.py']
 return dict(status='PASS_INDEPENDENT_AFTER_THIRD_71_LINES_AND_24_GERMS',sources={p:sha(p)for p in paths},
  arithmetic='Separate SymPy ANP field and192bit rational outward signs, no producer quartic or profile imports.',
  line_counts=dict(counts),complete_geodesic_faces=expected,new_cap_vertices=len(new),open_germs=germs,
  finite_edge_types_by_new_endpoints={str(k):v for k,v in endpoint_types.items()},older_open_germs=261,
  finite_profiles=profiles_total,dihedral_tails=tails_total,profile_summaries=profile_summaries,
  inert_I2_10_component_occurrences=inert_components,profile_types=dict(types),
  conceptual_checks=['Two exterior generators and PD8 support put both future null limits strictly inside the feasible2D cone face.',
   'Finite-start and strict-spacelike same-support witnesses show each opposite null limit has precisely8PD incident walls and unbounded link.',
   'The inherited profile engine permits inert I2(10) components; none occurs in these24 supports.',
   'Cycles, four-contact stars, heavy three-contact stars and label>=6 paths exclude every omitted finite profile; isolated I2 tails remain symbolic.'],
  scope='Independent complete3261candidate-pair bookkeeping and71line geometry; all24new cap-vertex germs and247finite/18symbolic local profiles. The prior independent28wall cone census is inherited, not repeated. Global cap feasibility, affine endpoint extensions and finite-volume completion remain open.')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved independent after-third line/frontier review')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
