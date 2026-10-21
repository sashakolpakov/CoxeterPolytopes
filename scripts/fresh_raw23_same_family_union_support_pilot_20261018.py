#!/usr/bin/env python3
"""Exact three retained raw23-pair unions; no occurrence inequality yet."""
import argparse,json
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import fresh_raw23_same_family_relative_pilot_20261018 as S
C=S.C;K=C.K;P=C.P;V=C.V;R=S.R;A=S.A;L=S.L;T=S.T;H=S.H;ROOT=S.ROOT
TAG='fresh_raw23_same_family_union_support_pilot_20261018'

def run():
 prior=S.read(S.TAG)
 for p,h in prior['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 assert prior['counts']['relative_cases']==175 and prior['counts']['retained_relative_cases']==3
 g=[[tuple(e)for e in row]for row in prior['geometry']['native_gram']]
 tau=K('7/2','1/2');cases=[]
 for index,row in enumerate(prior['relative_cases']):
  if not row['result']['allowed']:continue
  kinds=row['kinds'];assert kinds in([1,1],[1,6],[6,6])
  frame={int(i):{int(j):m for j,m in r.items()}for i,r in row['finite_frame'].items()}
  stars=[row['anchor'],row['other']]
  active=sorted({v for st in stars for _,v in st['active_embedding']});assert len(active)==9
  name={v:i for i,v in enumerate(active)}
  covs=[tuple((i,tuple(e))for i,e in cv)for st in stars for cv in st['native']['contacts']]
  common=set(covs[:2])&set(covs[2:]);assert len(common)==1
  common=next(iter(common));assert len(common)==1 and common[0][1][0]==0
  qroot=common[0][0];assert common[0][1][1:]==('7/2','1/2')
  other=[next(cv for cv in covs[2*i:2*i+2]if cv!=common)for i in range(2)]
  cv=[common]+other;assert len(set(cv))==3
  # Nodes0..8 are the contacted finite frame;9 is the common exterior
  # wall;10,11 are the two distinct exterior walls.
  ex={c:9+i for i,c in enumerate(cv)};scale=[K(1)]*12
  if kinds==[1,1]:scale[9]=tau
  elif kinds==[6,6]:scale[name[qroot]]=tau
  else:
   scale[9]=scale[11]=tau
   mp=dict(stars[1]['active_embedding'])
   for i in mp.values():
    if i!=qroot:scale[name[i]]=tau
  # scale means the squared positive scaling of the physical norm2 root.
  G=[[K()for j in range(12)]for i in range(12)]
  for i in range(12):G[i][i]=2*scale[i]
  for a,b in combinations(active,2):
   m=frame[a].get(b,2)
   if m==2:continue
   assert scale[name[a]]==scale[name[b]]and m in(3,5)
   x=K(-1)if m==3 else K('-1/2','-1/2')
   G[name[a]][name[b]]=G[name[b]][name[a]]=x*scale[name[a]]
  for i,c in enumerate(cv):
   for base,e in c:
    j=name[base];k=9+i
    if e[0]==0:
     assert i==0 and base==qroot and scale[j]*scale[k]==tau;z=-tau
    else:
     assert e==(3,'0','0')and scale[j]==scale[k];z=-scale[j]
    G[j][k]=G[k][j]=z
  for i,kind in enumerate(kinds):
   z=K()if kind==1 else -scale[9]
   assert kind==1 or scale[9]==scale[10+i]
   G[9][10+i]=G[10+i][9]=z
  G[10][11]=G[11][10]=K(-1,-1)if kinds[0]==kinds[1]else K(-3,-2)
  # Check each guessed scaled cross product against the full native
  # projections and unique future heights, in its own immutable field.
  J,_,_,tens=T.projection_data(frame);assert not tens
  cross=[]
  for i,j in combinations(range(3),2):
   z=G[9+i][9+j];squared=z*z/(scale[9+i]*scale[9+j])
   expected=(2,'0','0')if z==K()else(0,str(squared.a),str(squared.b))
   f=A.pair_field((cv[i],cv[j]),(expected,));pr=R.projection(cv[i],cv[j],J,f)
   hs=[]
   for c in(cv[i],cv[j]):
    q=R.projection(c,c,J,f);assert not any(q[1:]);hs.append(q[0]-2);assert hs[-1].sign()>0
   native=f.negative_contact(expected)
   assert f.signs_difference(tuple(x-y for x,y in zip(pr,native)),hs[0]*hs[1])==0
   cross.append(dict(nodes=[9+i,9+j],scaled_product=z.out(),physical_squared_product=squared.out()))
  inertia=P.congruence(G)['inertia'];conjugate=P.congruence([[K(x.a,-x.b)for x in r]for r in G])['inertia']
  assert inertia==[9,1,2]
  vertices=[]
  for ids in combinations(range(12),9):
   it=P.congruence(P.sub(G,ids))['inertia']
   if it==[9,0,0]:vertices.append(list(ids))
  known=[];prism_maps=[]
  for k,st in enumerate(stars):
   active_map=dict(st['active_embedding']);mp={i:name[j]for i,j in st['active_embedding']};om=row['omissions'][k]
   for i in om:
    cov=tuple(sorted((active_map[j],tuple(g[i][j]))for j in active_map if g[i][j]!=A.RIGHT))
    mp[i]=ex[cov]
   assert set(mp)==set(range(7))and len(set(mp.values()))==7
   passive=sorted(set(range(9))-set(mp.values()))
   vv=[sorted(passive+[mp[i]for i in range(7)if i not in s['omitted']])for s in prior['geometry']['finite_vertex_types']]
   assert all(v in vertices for v in vv)and len(vv)==10
   known+=vv;prism_maps.append(dict(raw_map=sorted(mp.items()),passive=passive,vertices=vv))
  known={tuple(v)for v in known};extra=[v for v in vertices if tuple(v)not in known]
  cases.append(dict(relative_index=index,kinds=kinds,finite_base_native_nodes=active,
   full_frame=row['finite_frame'],squared_scales=[s.out()for s in scale],scaled_gram=P.out(G),
   cross_product_checks=cross,inertia=inertia,conjugate_inertia=conjugate,
   all_positive_definite_nine_supports=vertices,prism_maps=prism_maps,
   inherited_union_vertex_count=len(known),extra_positive_definite_supports=extra,
   passive_rank=4,actual_vertex_argument_inherited=True))
 paths={f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',f'results/{S.TAG}.json',
  'results/fresh_sqrt5_prism_simplex_base_incidence_20261017.md',
  'results/fresh_raw357_extra_vertex_incidence_pilot_20261018.md'}
 paths|={str(Path(m.__file__).relative_to(ROOT))for m in(C,P,V,R,A,L,T,H)}
 return S.V.normal(dict(status='PASS_RAW23_THREE_RETAINED_UNION_SUPPORT_PILOT',
  sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
  cases=cases,complete_support_decisions=3*220,optimizer_invocations=0,exported_columns=[],
  scope='Three retained relative union cores only. Exact conjugate/physical Gram and complete nine-support decisions; positive supports are actual vertices conditional on occurrence, via the frozen acute supporting-wall projection theorem. No new incidence quota, recovery multiplicity, capacity, matrix column or bound is claimed.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/f'results/{TAG}.json'
 assert args.verify or not p.exists();out=run()
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],[(r['kinds'],len(r['all_positive_definite_nine_supports']),r['inherited_union_vertex_count'],len(r['extra_positive_definite_supports']),r['conjugate_inertia'])for r in out['cases']],flush=True)
