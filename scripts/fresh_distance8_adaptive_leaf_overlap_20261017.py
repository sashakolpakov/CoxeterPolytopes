#!/usr/bin/env python3
"""Nine literal high-core overlap contexts: all-label adaptive finite payments.

This is a local certificate, not a complete distance-eight overlap census.
Elementary rational algebra is inherited; the profile/domain checks below are new.
"""
import argparse,json,hashlib
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,product
from pathlib import Path
import fresh_selective_distance9_dimension108_20261017_review as R
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_adaptive_leaf_overlap_20261017'
INPUTS=['results/fresh_high_rank_quasilanner_census_20261017.json',
 'scripts/fresh_selective_distance9_dimension108_20261017_review.py',
 'scripts/fresh_selective_distance10_dimension116_20261017_review.py',
 'results/fresh_distance9_supporting_extensions_review_20261017.md',
 'results/fresh_finite8_a17_d11_coupon_budget_20261017.md']
CASES=[
 (22,(1,2),25,(1,2),((5,9),(6,8),(7,7),(8,6))),
 (25,(0,1),25,(1,2),((3,9),)),
 (25,(0,1),25,(1,2),((5,9),(6,8),(7,7),(8,6),(9,5))),
 (25,(0,1),25,(1,4),((9,3),)),
 (25,(1,2),22,(1,2),((6,8),(7,7),(8,6),(9,5))),
 (25,(1,2),25,(0,1),((9,3),)),
 (25,(1,2),25,(0,1),((5,9),(6,8),(7,7),(8,6),(9,5))),
 (25,(1,2),25,(1,2),((4,9),(5,8),(6,7),(7,6),(8,5),(9,4))),
 (25,(1,4),25,(0,1),((3,9),)),
]
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
 if isinstance(x,(tuple,list)):return[enc(v)for v in x]
 return x

def profiles(g,iv,corners,du,dv):
 n=len(g);ci=[]
 for c in corners:
  ns=sorted(i for b in c['blocks']for i in b['nodes'])
  h=R.extract(g,ns);R.psd(h,True);ci.append((ns,R.inverse(h)))
 forbidden={i for ns,j in ci for k,i in enumerate(ns)if j[k][k]>2}
 allowed=[i for i in range(n)if i not in forbidden]
 out=[];reasons=Counter();cert=[]
 for size in range(1,len(allowed)+1):
  for ss in combinations(allowed,size):
   short=next(((a,b)for a,b in combinations(ss,2)if min(du[a]+dv[b],du[b]+dv[a])+2<=8),None)
   if short:reasons['one_root_shortcut']+=1;continue
   vals=[sum(j[a][b]for a,i in enumerate(ns)for b,k in enumerate(ns)if i in ss and k in ss)for ns,j in ci]
   if max(vals)>2:reasons['corner_monotone_projection']+=1;continue
   sub=R.extract(iv,ss)
   if any(x<0 for row in sub for x in row):
    assert len(ss)==2
    choices=[(k,1-k,sub[1-k][1-k]-sub[k][1-k]**2/sub[k][k])for k in(0,1)if sub[k][k]>0]
    pivot,other,lower=max(choices,key=lambda x:x[2]);assert lower>2
    cert.append(dict(support=ss,pivot=pivot,remaining_coordinate=other,schur_lower=lower))
    reasons['all_real_quadratic_lower_bound']+=1;continue
   q=sum(sum(row)for row in sub)
   if q>2:reasons['nonnegative_core_projection']+=1;continue
   if size>1:
    assert q==2 and all(x>0 for row in sub for x in row)
    out.append(dict(support=ss,scale_square=1,projection=q,residual=2-q))
   else:
    assert q>0 and 2/q<=F(5,2) # label5 squared=(3+sqrt5)/2 > 5/2
    for scale in(1,2):
     if scale*q<=2 and all(scale*v<=2 for v in vals):
      out.append(dict(support=ss,scale_square=scale,projection=scale*q,residual=2-scale*q))
 return out,dict(allowed=allowed,forbidden=sorted(forbidden),support_rejections=dict(reasons),negative_projection_certificates=cert,
                  label_bound='All surviving contacts have label3 or4; multi-contact profiles are unit contacts.')

def mutual(iv,a,b):return sum(iv[i][j]for i in a['support']for j in b['support'])
def pair_short(iv,a,b,du,dv):
 s=mutual(iv,a,b);rho=a['residual']*b['residual'];sq=a['scale_square']*b['scale_square']*s*s
 dist=min(min(du[i]+dv[j],du[j]+dv[i])+3 for i in a['support']for j in b['support'])
 assert s<0 and sq>rho and dist<=8
 return dict(left=a,right=b,projection_base=s,projection_square=sq,residual_product=rho,shortcut=dist)

def leaf_record(g,c,leaf,untouched):
 ns=c['marked_component'];assert ns and leaf in c['preferred_corner']
 neighbors={i:[j for j in ns if j!=i and g[i][j]]for i in ns}
 assert len(neighbors[leaf])==1
 neighbor=neighbors[leaf][0]
 assert leaf in untouched and neighbor in untouched and len(neighbors[neighbor])==2
 assert g[leaf][neighbor]==-1 and len(ns)>=10
 p,q=c['preferred_corner'];dist=R.graph_dist(R.extract(g,ns),[ns.index(p)])[ns.index(q)]
 assert dist==9
 return dict(corner=c['preferred_corner'],leaf=leaf,neighbor=neighbor,base_order=len(ns),distance=dist,
             conclusion='Whole corner is finite; the marked simple leaf remains a leaf, hence finite leaf9 pays.')

def run():
 ds={r['index']:r for r in json.loads((ROOT/INPUTS[0]).read_text())['diagrams']};rows=[]
 for ix,(di,u,dj,v0,match)in enumerate(CASES):
  native,g,sc,mp=R.combined(ds,di,dj,match);assert set(sc)=={1}
  v=[mp[i]for i in v0];u=list(u);iv=R.inverse(g);sig=R.inertia(g);assert sig==[len(g)-1,1,0]or tuple(sig)==(len(g)-1,1,0)
  du,dv=R.graph_dist(native,u),R.graph_dist(native,v);assert min(du[i]for i in v)==9
  corners=[R.inspect_corner(g,sc,u+v,p,q)for p,q in product(u,v)]
  ps,domain=profiles(g,iv,corners,du,dv);contacted={i for p in ps for i in p['support']}
  untouched=set(range(len(g)))-contacted
  row=dict(index=ix,left=[di,u],right=[dj,list(v0)],match=match,right_copy_map=mp,gram=g,inverse=iv,signature=sig,
           all_corners=corners,profiles=ps,domain=domain)
  if ix in(1,3,5,8):
   options=[]
   for c in corners:
    ns=c['marked_component']
    if not ns:continue
    for leaf in c['preferred_corner']:
     try:options.append(leaf_record(g,c,leaf,untouched))
     except AssertionError:pass
   assert options;row.update(method='sealed_marked_leaf',payment=options[0])
  elif ix in(0,4,7):
   # Choose the actual core corner in which both marked roots are simple leaves.
   options=[]
   for c in corners:
    ns=c['marked_component']
    if not ns:continue
    a,b=c['preferred_corner'];ng={i:[j for j in ns if j!=i and g[i][j]]for i in ns}
    if not(len(ng[a])==len(ng[b])==1):continue
    if not(all(ng[x][0]in untouched and len(ng[ng[x][0]])==2 for x in(a,b))):continue
    left=[p for p in ps if a in p['support']];right=[p for p in ps if b in p['support']]
    assert all(a not in p['support']or b not in p['support']for p in ps)
    tests=[pair_short(iv,x,y,du,dv)for x,y in product(left,right)]
    leaves=[leaf_record(g,c,x,untouched|{x})for x in(a,b)]
    options.append(dict(corner=c['preferred_corner'],leaf_alternatives=leaves,pair_exclusions=tests))
   assert options;row.update(method='adaptive_untouched_marked_leaf',payment=options[0])
  else:
   assert ix in(2,6)
   # Exchange the two copies for the reversed orientation; these roles are literal.
   far=[1,2,3]if ix==2 else[11,12,13]
   good=11 if ix==2 else 1;bad=13 if ix==2 else 3
   mark_other=0 if ix==2 else 10
   c=next(c for c in corners if set(c['preferred_corner'])=={mark_other,good})
   gs=[p for p in ps if good in p['support']];fs=[p for p in ps if set(p['support'])&set(far)]
   assert len(gs)==1 and gs[0]['support']==(good,)and gs[0]['scale_square']==1
   tests=[pair_short(iv,x,gs[0],du,dv)for x in fs]
   leaf=leaf_record(g,c,good,untouched|{good})
   # If no far-side root occurs, the opposite base endpoint and neighbor stay fixed.
   ns=c['marked_component'];ng={i:[j for j in ns if j!=i and g[i][j]]for i in ns}
   sealed=set(range(len(g)))-{i for p in ps if not(set(p['support'])&set(far))for i in p['support']}
   fixed=[x for x in ns if len(ng[x])==1 and x in sealed and ng[x][0]in sealed and len(ng[ng[x][0]])==2]
   assert fixed and len(ns)==12
   active=[p for p in ps if not(set(p['support'])&set(far))]
   relevant=[p for p in active if set(p['support'])&set(ns)]
   assert relevant==gs
   a=iv[good][good];b=iv[bad][bad];ab=iv[good][bad]
   assert(a,b,ab)==(F(5,4),F(1),F(1,2))
   res=[[2-a,-ab],[-ab,2-b]];ji=R.inverse(res)
   assert ji==[[F(2),F(1)],[F(1),F(3,2)]]
   # Every remaining profile either is this unit outside root or is incompatible with G.
   other_checks=[]
   for p in active:
    if p==gs[0]:continue
    if p['support']==(bad,)and p['scale_square']==1:continue
    s=mutual(iv,p,gs[0]);sq=p['scale_square']*s*s;rp=p['residual']*gs[0]['residual']
    assert s>0 and sq>rp
    other_checks.append(dict(profile=p,projection_square=sq,residual_product=rp))
   assert (ab+1)**2>(2-a)*(2-b)
   row.update(method='adaptive_leaf_or_bounded_relevant_tail',payment=dict(corner=c['preferred_corner'],
      far_support_nodes=far,if_far_present=leaf,pair_exclusions=tests,if_far_absent_fixed_simple_leaves=fixed,
      good_node=good,outside_node=bad,residual_gram=res,residual_inverse=ji,other_profile_exclusions=other_checks,
      alone_tail_rank_bound=3,joint_tail_rank_bound=1,whole_component_bound=16,
      conclusion='Finite: either marked leaf9, or whole component order at most16.'))
  rows.append(row)
 sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in INPUTS}
 sources['scripts/'+TAG+'.py']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
 return enc(dict(status='PASS',scope='Nine literal contexts only; no full overlap census or dimension claim.',
   hypotheses=['actual quadrilateral face star','no marked path of length at most8','all unspecified cross-core contacts already excluded'],
   method_counts=dict(Counter(r['method']for r in rows)),cases=rows,sources=sources))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run();path=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print('PASS',out['method_counts'],len(out['cases']),flush=True)
