#!/usr/bin/env python3
"""Separate exact arithmetic replay of the82-wall polar extension, no producer."""
import argparse,copy,hashlib,json
from collections import Counter
from pathlib import Path
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv,label,check as base_check
ROOT=Path(__file__).resolve().parents[1]
SOURCE='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
OUT='data/enumerations/fresh_dimension_bridge_20261013_polar_verify.json'
def rank(a):
 a=[list(r)for r in a];n=len(a);w=len(a[0]);r=0
 for j in range(w):
  k=next((i for i in range(r,n)if a[i][j]!=Z),None)
  if k is None:continue
  a[r],a[k]=a[k],a[r];p=a[r][j];a[r]=[x/p for x in a[r]]
  for i in range(r+1,n):
   p=a[i][j]
   if p!=Z:a[i]=[x-p*y for x,y in zip(a[i],a[r])]
  r+=1
  if r==n:break
 return r

def check(d,base):
 assert d['status']=='EXACT_POLAR_EXTENSIONS_WITH_STRICT_ESCAPE'
 assert hashlib.sha256((ROOT/d['source']).read_bytes()).hexdigest()==d['source_sha256']
 Q=tuple(map(vec,base['form']));roots=list(map(vec,base['normals']));cov=[mv(Q,r)for r in roots]
 ns=[dot(r,c)for r,c in zip(roots,cov)];h=vec(base['strict_center']);hc=mv(Q,h);old=vec(base['strict_escape']['vector'])
 alltypes=Counter(row[2]for row in base['pair_checks']);addedpairs=0;actualnorms=[]
 for cut in d['cuts']:
  idx=cut['index'];assert idx==len(roots)
  r=vec(cut['vector']);rc=mv(Q,r);q=dot(r,rc);ch=dot(r,hc)
  assert q.sign()>0 and ch.sign()<0
  assert q==parse(cut['norm'])and ch==parse(cut['center_pairing'])
  ids=cut['active_indices'];assert len(ids)==len(set(ids))==20 and all(0<=i<idx for i in ids)
  assert rank([cov[i]for i in ids])==20
  ps=[dot(r,c)for c in cov];assert all(ps[i]==Z for i in ids)
  assert ps==list(map(parse,cut['pairings']))and all(p.sign()<=0 for p in ps)
  labs=[label(p,q,n)for p,n in zip(ps,ns)];assert labs==cut['pair_kinds']
  if idx==70:assert dot(r,mv(Q,old))==parse(cut['first_old_escape_pairing'])and parse(cut['first_old_escape_pairing']).sign()>0
  else:assert cut['first_old_escape_pairing']is None
  # Each polar truncates its own feasible spacelike ray, since r²>0.
  # Existing bounded cusp links cannot acquire a new acute incident facet;
  # explicitly confirm strict separation at all107 inherited rays.
  assert all(dot(r,mv(Q,vec(row['ray']))).sign()<0 for row in base['cusp_rows'])
  alltypes.update(labs);addedpairs+=idx;actualnorms.append(q.enc());roots.append(r);cov.append(rc);ns.append(q)
 assert len(roots)==d['wall_count']==82 and len(d['cuts'])==12
 assert sum(alltypes.values())==82*81//2
 esc=d['strict_escape'];y=vec(esc['vector']);yc=mv(Q,y);q=dot(y,yc);ch=dot(y,hc)
 assert q.sign()>0 and ch.sign()<0 and q==parse(esc['norm'])and ch==parse(esc['center_pairing'])
 ps=[dot(r,yc)for r in roots];assert all(p.sign()<0 for p in ps)and ps==list(map(parse,esc['all_wall_pairings']))
 return {'status':'PASS_INDEPENDENT_82_WALL_H20_POLAR_EXTENSION_WITH_STRICT_ESCAPE',
         'source':SOURCE,'source_sha256':hashlib.sha256((ROOT/SOURCE).read_bytes()).hexdigest(),
         'base_source':d['source'],'base_source_sha256':d['source_sha256'],'walls':82,'new_polars':12,
         'new_previous_pair_checks':addedpairs,'total_pair_checks':sum(alltypes.values()),
         'total_pair_type_counts':dict(sorted(alltypes.items())),'new_norms':actualnorms,
         'all_inherited_cusps_preserved_strictly':107,'strict_escape_norm':q.enc(),'strict_escape_center_pairing':ch.enc(),
         'scope':'Independent quadratic arithmetic validates every added normal, norm, all relations to every preceding wall, strict center,20independent orthogonal active walls, positive cutting of its own feasible spacelike ray, preservation of all107 inherited cusps, and final strict spacelike escape. All82 walls are true facets by acute projection. No completeness claim for the remaining frontier or arbitrary further additions.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
 d=json.loads((ROOT/SOURCE).read_text());base=json.loads((ROOT/d['source']).read_text());base_check(base);out=check(d,base)
 bad=copy.deepcopy(d);bad['cuts'][0]['norm']=['0','0']
 try:check(bad,base)
 except AssertionError:pass
 else:raise AssertionError('forged polar norm accepted')
 bad=copy.deepcopy(d);bad['strict_escape']['vector']=base['strict_center']
 try:check(bad,base)
 except AssertionError:pass
 else:raise AssertionError('timelike escape accepted')
 out['forged_norm_and_timelike_escape_rejected']=True
 if a.verify:assert out==json.loads((ROOT/OUT).read_text())
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out),flush=True)
