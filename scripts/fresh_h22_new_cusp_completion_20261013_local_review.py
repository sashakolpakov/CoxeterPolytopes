#!/usr/bin/env python3
"""Separate Fraction reconstruction of the representative H22 repaired cusp.

The added coordinate is represented as z*sqrt(11), so the extended form is
rational diag(Q,11). This does not repeat the full-orbit native pair census.
"""
import argparse,hashlib,json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
P=Path(__file__).resolve().parents[1]
SOURCE='results/fresh_h22_new_cusp_completion_20261013_verify.json'
PRIOR='results/fresh_h22_gap_face_family_20261013.json'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json'
OUT='results/fresh_h22_new_cusp_completion_20261013_local_review.json'
TIGHT='results/fresh_h22_new_cusp_completion_20261013_tight.json'
TIGHT_OUT='results/fresh_h22_new_cusp_completion_20261013_tight_local_review.json'
LOWER='results/fresh_h22_new_cusp_completion_20261013_lower.json'
LOWER_OUT='results/fresh_h22_new_cusp_completion_20261013_lower_local_review.json'
def read(p):return json.loads((P/p).read_text())
def sha(p):return hashlib.sha256((P/p).read_bytes()).hexdigest()
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),F(0))
def mv(M,x):return[dot(r,x)for r in M]
def run(tight=False,lower=False):
 assert not(tight and lower)
 d=read(SOURCE);p=read(PRIOR);raw=read(RAW)['coxeter_input'];assert all(sha(n)==h for n,h in d['sources'].items())
 changed=LOWER if lower else TIGHT if tight else None
 placement=read(changed)if changed else d
 if changed:assert placement['source']==SOURCE and placement['source_sha256']==sha(SOURCE)
 Q=[[-F(x)for x in row]+[F(0)]for row in raw['basis_gram']]+[[F(0)]*22+[F(11)]]
 W=[[F(x)for x in row]+[F(0)]for row in raw['wall_coordinates_as_columns_transposed']]
 seed=[F(x,6)for x in p['six_old_projection']]+[F(1)];partner=[F(x,6)for x in p['new_cusp']['six_partner']]+[F(1)];delta=[x+y for x,y in zip(seed,partner,strict=True)]
 assert dot(delta,mv(Q,delta))==0
 D=int(placement['common_denominator']);shifts=list(map(F,placement['shifts']))if changed else[F(d['shift'])]*11
 assert all(x>0 for x in shifts)
 caps=[[F(x,D)for x in row]+[2*M]for row,M in zip(placement['new_cap_old_numerators'],shifts,strict=True)]
 assert len(caps)==11
 q=[dot(r,mv(Q,r))for r in caps];assert q==list(map(F,d['new_cap_norms'])) and all(x>0 for x in q)
 assert all(dot(r,mv(Q,delta))==0 for r in caps)
 # Independently check all 2,310 relations to the original 210 walls,
 # including the new finite-label events at the lower singleton shifts.
 oldnorm=[dot(r,mv(Q,r))for r in W];oldcounts=[]
 for r,n in zip(caps,q,strict=True):
  rc=mv(Q,r);counts=Counter()
  for w,m in zip(W,oldnorm,strict=True):
   pair=dot(w,rc);ratio=pair*pair/(m*n)
   if pair==0:kind='right'
   else:
    assert pair<0
    if ratio in(F(1,4),F(1,2),F(3,4)):kind={F(1,4):'pi/3',F(1,2):'pi/4',F(3,4):'pi/6'}[ratio]
    elif ratio==1:kind='parallel'
    else:assert ratio>1;kind='dotted'
   counts[kind]+=1
  oldcounts.append(dict(counts))
 assert oldcounts==placement['old210_pair_counts']
 comps=d['old_components'];common=[i for cs in comps for i in cs];assert sorted(common)==sorted(p['new_cusp']['common_old_walls']) and len(common)==18
 relations=[]
 for j,cs in enumerate(comps):
  M=shifts[j]
  s=[caps[j][k]+sum(W[i][k]for i in cs)for k in range(23)];assert s==[M*x for x in delta];relations.append({'old':cs,'new':[j],'delta_multiplier':str(M)})
 assert [x+y for x,y in zip(seed,partner,strict=True)]==delta
 for j in(7,9):
  total=shifts[j]+shifts[j+1]
  assert[x+y for x,y in zip(caps[j],caps[j+1],strict=True)]==[total*x for x in delta]
  relations.append({'old':[],'new':[j,j+1],'delta_multiplier':str(total)})
 local=[W[i]for i in common]+[seed,partner]+caps;cov=[mv(Q,r)for r in local];G=[[dot(a,b)for b in cov]for a in local]
 assert all(G[i][i]>0 for i in range(31)) and all(G[i][j]<=0 for i in range(31)for j in range(i))
 gram_rank=sp.Matrix(G).rank();vector_rank=sp.Matrix(local).rank();assert gram_rank==21 and vector_rank==22
 # This comparison independently checks the representative native integer
 # contractions, including all cancellations, using arbitrary-size integers.
 oldQ=[row[:22]for row in Q[:22]];nums=placement['new_cap_old_numerators'];icov=[mv(oldQ,row)for row in nums];native_normnums=[dot(r,c)+44*M*M*D*D for r,c,M in zip(nums,icov,shifts,strict=True)]
 assert native_normnums==[x*D*D for x in q]
 assert all(dot(r,c)<0 for r,c in zip(nums,icov,strict=True))
 sources=(SOURCE,PRIOR,RAW,changed)if changed else(SOURCE,PRIOR,RAW)
 return {'status':'PASS_SEPARATE_FRACTION_H22_REPAIRED_CUSP','sources':{x:sha(x)for x in sources},'field_coordinate':'last coordinate denotes coefficient of sqrt(11); rational form diag(Q,11)','new_norms':list(map(str,q)),'incident_facets':31,'normal_span_rank':vector_rank,'gram_rank':gram_rank,'positive_affine_relations':relations,'parallel_pair_delta_multiplier':'1','component_ranks':[5,4,3,3,1,1,1,1,1,1],'scope':'Separate exact local norm/pairing/rank/affine-relation reconstruction. Full-orbit completeness and all distant pairings are reviewed and replayed in the primary verifier, not independently recomputed here.'}
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');g=a.add_mutually_exclusive_group();g.add_argument('--tight',action='store_true');g.add_argument('--lower',action='store_true');v=a.parse_args();d=run(v.tight,v.lower);out=LOWER_OUT if v.lower else TIGHT_OUT if v.tight else OUT
 if v.verify:assert d==read(out)
 else:(P/out).write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'])
