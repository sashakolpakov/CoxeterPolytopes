#!/usr/bin/env python3
"""Seven complete sole-A1 cap profiles at one compact simplex core.

Exact native pair exclusions and typed coface coefficients only. No broad
245-profile census, vertex graph edit, optimizer, or new global bound.
"""
import argparse,json
from collections import Counter,defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_core1_mixed_coface_capacity_pilot_20261018'
RAW=(23,24,25,33,36,38,39)
STEMS=dict(incidence='fresh_sqrt5_prism_simplex_base_incidence_20261017',
 raw='fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017',
 forced='fresh_sqrt5_forced_base_simplex_domain_20261017',
 actual='fresh_sqrt5_actual_provenance_dual_20261017',
 continuation='fresh_all_marker_families_dual_continuation_20261018',
 current='fresh_raw357_coface_full_model_20261018',
 witness='fresh_raw23_surviving_coface_witness_20261018')
Z=(F(0),F(0));ONE=(F(1),F(0))
def q(a=0,b=0):return(F(a),F(b))
def add(a,b):return(a[0]+b[0],a[1]+b[1])
def neg(a):return(-a[0],-a[1])
def sub(a,b):return add(a,neg(b))
def mul(a,b):return(a[0]*b[0]+5*a[1]*b[1],a[0]*b[1]+a[1]*b[0])
def inv(a):
 d=a[0]**2-5*a[1]**2;assert d
 return(a[0]/d,-a[1]/d)
def sg(a):
 x,y=a
 if not y:return(x>0)-(x<0)
 if not x:return(y>0)-(y<0)
 if(x>0)==(y>0):return(x>0)-(x<0)
 return((x>0)-(x<0))*((x*x>5*y*y)-(x*x<5*y*y))
def inverse(g):
 n=len(g);a=[list(r)+[ONE if i==j else Z for j in range(n)]for i,r in enumerate(g)]
 for i in range(n):
  k=next(k for k in range(i,n)if a[k][i]!=Z);a[i],a[k]=a[k],a[i]
  s=inv(a[i][i]);a[i]=[mul(v,s)for v in a[i]]
  for k in range(n):
   if k!=i:
    s=a[k][i];a[k]=[sub(u,mul(s,v))for u,v in zip(a[k],a[i])]
 return[r[n:]for r in a]
def bil(x,J,y):return sumq(mul(mul(a,J[i][j]),b)for i,a in enumerate(x)for j,b in enumerate(y))
def sumq(xs):
 out=Z
 for x in xs:out=add(out,x)
 return out
def enc(v):
 if isinstance(v,F):return str(v)
 if isinstance(v,dict):return{k:enc(x)for k,x in v.items()}
 if isinstance(v,(list,tuple)):return[enc(x)for x in v]
 return v
def read(s):return json.loads((ROOT/f'results/{s}.json').read_text())

def run():
 ds={k:read(s)for k,s in STEMS.items()};inc=ds['incidence'];raw=ds['raw']['records']
 core=next(c for c in inc['cores']if c['index']==1)
 assert core['nodes']==5 and core['simplex_dimension']==4 and core['spherical_deletions']==list(range(5))
 G=[[q(2)if i==j else Z for j in range(5)]for i in range(5)]
 for i,j,e in core['edges']:
  assert e[0]in(3,5);G[i][j]=G[j][i]=q(-1)if e[0]==3 else q('-1/2','-1/2')
 autos=[p for p in permutations(range(5))if all(G[i][j]==G[p[i]][p[j]]for i in range(5)for j in range(5))]
 assert autos==[tuple(range(5))];J=inverse(G)
 profiles=[];selected={}
 for t in inc['prism_templates']:
  for b in t['bases']:
   if b['base_core']!=1 or b['finite_extra']!=['A1']:continue
   r=raw[t['source_core']];mp=dict(b['core_isomorphism']);om=b['deleted_normal'];cq=b['opposite_base_facet']
   if set(range(r['nodes']))!={om,cq}|set(mp):continue
   assert sorted(mp.values())==list(range(5)) and om not in mp and cq not in mp
   es=dict(zip(combinations(range(r['nodes']),2),r['row']['labels']))
   label=lambda i,j:es[tuple(sorted((i,j)))]
   assert all(label(cq,i)==2 for i in mp)and label(om,cq)==0
   for i,j in combinations(mp,2):
    m=label(i,j);assert m in(2,3,5)
    assert G[mp[i]][mp[j]]==(Z if m==2 else q(-1)if m==3 else q('-1/2','-1/2'))
   x=[Z]*5
   for i,v in mp.items():
    m=label(om,i);assert m in(2,3);x[v]=Z if m==2 else q(-1)
   tau=tuple(map(F,r['parameters'][','.join(map(str,sorted((om,cq))))]))
   norm=bil(x,J,x);assert add(norm,mul(q('1/2'),tau))==q(2)
   z=dict(raw=t['source_core'],canonical_core=t['canonical_core'],prism_dimension=t['prism_dimension'],
    deleted=om,coface=cq,base_map=b['core_isomorphism'],contacts=x,tau=tau,projection_norm=norm)
   profiles.append(z);selected[t['canonical_core']]=z
 profiles.sort(key=lambda z:z['raw']);assert tuple(z['raw']for z in profiles)==RAW
 pairs=[]
 for a,b in combinations(profiles,2):
  dot=bil(a['contacts'],J,b['contacts'])
  square_gap=sub(mul(a['tau'],b['tau']),mul(q(4),mul(dot,dot)))
  assert sg(dot)<0 and sg(square_gap)>0
  # dot+sqrt(tau_a*tau_b)/2 >0: forbidden for two outward normals.
  pairs.append(dict(raw=[a['raw'],b['raw']],projection_product=dot,
    fourfold_positive_product_gap=square_gap,same_actual_A1_coface_excluded=True))
 # Distinct A1 normals are orthogonal. Retain every unexcluded pair.
 different=[]
 for i,a in enumerate(profiles):
  for b in profiles[i:]:
   dot=bil(a['contacts'],J,b['contacts']);gap=sg(dot)<0 and sg(add(dot,ONE))>0
   if gap:assert a['raw']==b['raw']==23
   else:assert sg(add(dot,ONE))<=0
   different.append(dict(raw=[a['raw'],b['raw']],product=dot,excluded_by_gap=gap))
 # Each source context is one maximal prism/quad aggregate, never an
 # orientation. Exactly one deleted2 base extends to this core1 simplex.
 provenance=ds['actual']['provenance_masses']
 explicit={provenance[pi]['source_context']:f for pi,f in ds['continuation']['prepared']['provenance_to_face']}
 source_by_face=defaultdict(list)
 for ctx,f,ci,k,bs in inc['contexts']:source_by_face[explicit.get(ctx,f)].append(ctx)
 grouped=defaultdict(list);chosen=[]
 for ctx,f,ci,k,bs in inc['contexts']:
  if ci not in selected:continue
  hits=[(t,m)for t,m in bs if inc['base_types'][t]['core']==1]
  assert len(hits)==1 and hits[0][1]==1
  t=hits[0][0];f=explicit.get(ctx,f)
  assert source_by_face[f]==[ctx]
  row=dict(context=ctx,face=f,raw=selected[ci]['raw'],canonical_core=ci,base_type=t,multiplicity=1)
  grouped[t].append(row);chosen.append(row)
 assert len(chosen)==1645 and len(grouped)==235 and len({r['face']for r in chosen})==1645
 forced=dict(ds['forced']['base_to_new_face'])
 alpha={f:F(x)for f,x in ds['current']['candidate']['face_multipliers']}
 witness={f:F(x)for f,x in ds['witness']['candidate']['face_multipliers']}
 columns=[];absent=[]
 for t,rows in sorted(grouped.items()):
  typ=inc['base_types'][t];assert sorted(r['raw']for r in rows)==list(RAW)
  cap=typ['spectator_forest'].count('A1');assert cap>=1
  sf=typ['old_face']if typ['represented_in_lp']else forced.get(t)
  z=dict(base_type=t,simplex_face=sf,forest=typ['spectator_forest'],capacity=cap,sources=rows,
    positive_face_coefficients=sorted([[r['face'],1]for r in rows]),vertex_cost=0)
  if sf is None:
   assert not any(witness.get(r['face'],F())for r in rows)
   z['current_source_mass']=sum((alpha.get(r['face'],F())for r in rows),F())
   absent.append(z);continue
  load=sum((alpha.get(r['face'],F())for r in rows),F());cost=cap*alpha.get(sf,F())
  newload=sum((witness.get(r['face'],F())for r in rows),F());newcost=cap*witness.get(sf,F())
  assert newload<=newcost
  z.update(integer_face_coefficients=sorted(z['positive_face_coefficients']+[[sf,-cap]]),
    current_load=load,current_capacity=cost,current_excess=load-cost,
    raw23_witness_load=newload,raw23_witness_capacity=newcost)
  columns.append(z)
 assert len(columns)==62 and len(absent)==173
 violations=[c for c in columns if c['current_excess']>0];assert len(violations)==3
 # Qualified singleton screen: each of the50 surviving raw23 source
 # masses is bounded by EACH actual simplex-base mass even at capacity1.
 singleton=[]
 bycontext={r[0]:r for r in inc['contexts']};assert len(bycontext)==len(inc['contexts'])
 for r in ds['witness']['kept_Q_sources']:
  ctx,f,ci,k,bs=bycontext[r['source_context']]
  assert ctx==r['source_context']and ci==19 and explicit.get(ctx,f)==r['face']and len(bs)==2
  for t,m in bs:
   typ=inc['base_types'][t];sf=typ['old_face']if typ['represented_in_lp']else forced.get(t)
   assert m==1 and sf is not None and witness[r['face']]<=witness.get(sf,F())
   singleton.append(dict(face=r['face'],base_type=t,simplex_face=sf,mass=witness[r['face']],simplex_mass=witness[sf]))
 assert len(singleton)==100
 paths={f'scripts/{TAG}.py'}|{f'results/{s}.json'for s in STEMS.values()}|{
  'results/fresh_sqrt5_prism_simplex_base_incidence_20261017.md',
  'results/fresh_sqrt5_cross_coface_cap_capacity_20261017.md'}
 return enc(dict(status='PASS_CORE1_SEVEN_PROFILE_MIXED_COFACE_CAPACITY_PILOT',
  sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in sorted(paths)},
  core=core,gram=G,inverse=J,automorphisms=autos,profiles=profiles,same_coface_pairs=pairs,
  different_coface_pairs=different,typed_sources=chosen,existing_target_columns=columns,
  absent_target_rows=absent,violating_existing_columns=violations,raw23_capacity1_base_checks=singleton,
  counts=dict(profiles=7,same_coface_exclusions=21,different_coface_tests=28,
   source_classes=1645,simplex_types=235,existing_targets=62,absent_targets=173,
   positive_absent_targets=sum(c['current_source_mass']>0 for c in absent),current_violations=3),
  full_matrix_exports_installed=False,optimizer_invocations=0,
  scope='Complete seven sole-A1 profiles at simplex core1 only. New mixed coface capacities are valid necessary face-incidence constraints; three current-dual violations are not an optimized improvement. The separate raw23-only full-matrix witness satisfies every62 represented column and both raw23 base-capacity1 tests. No general cross-simplex impossibility or complete245-profile census is claimed.'))

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();p=ROOT/f'results/{TAG}.json'
 assert args.verify or not p.exists();out=run()
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],out['counts'],flush=True)
