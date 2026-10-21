#!/usr/bin/env python3
"""Exact polar candidates and a mutually compatible twelve-cap H9 batch."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_nonarithmetic_escape35_20261016_cone as D
A,R,C,K,Z,O=D.A,D.R,D.C,D.K,D.Z,D.O
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_joint_regions_20261016_batch'
PRIMARY=D.PRIMARY;CONE=f'results/{D.TAG}.json';SCOUT='results/fresh_nonarithmetic_joint_regions_20261016_polar.json'
RAY_IDS=[834,452,382,973,12,161,232,14,2,233,313,1,594,930,13]
BATCH=list(range(12))
I=R.L.I;ZERO=D.ZERO
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):
 d=json.loads((ROOT/p).read_text())
 for f,h in d['sources'].items():assert sha(f)==h,f
 return d
def rnd(z,bits=160):
 den=1<<bits
 return I(F(z.l.numerator*den//z.l.denominator,den),F(-((-z.h.numerator*den)//z.h.denominator),den))
@lru_cache(maxsize=100000)
def iv(a):
 z=I()
 for x,b in zip(a,D.BIV):
  if x:z=z+D.biv(x)*b
 return rnd(z)
def dot_sign(v,w):
 z=I()
 for a,b in zip(v,w):z=z+iv(a)*iv(b)
 if z.l>0:return 1
 if z.h<0:return-1
 return D.sign(D.dot(v,w))
def vecadd(v,w):return tuple(D.add(a,b)for a,b in zip(v,w))
def vecscale(v,s):return tuple(D.mul(a,s)for a in v)
def vecsub(v,w):return tuple(D.sub(a,b)for a,b in zip(v,w))
def geometry():
 local=read(PRIMARY);cone=read(CONE);old=read(D.START)
 Q=[[A.B.db(x)for x in row]for row in local['form']]
 H=tuple(D.lift(A.B.db(x))for x in local['center']);N=[tuple(D.flatten_raw(x)for x in row)for row in local['normals']]
 cov=[tuple(D.mv(Q,n))for n in N];baseN=R.m4(old['normals']);basecov=[C.mv(Q,n)for n in baseN]
 initial={tuple(s)for rows in old['complete_census'].values()for s in rows}
 crossings={s['wall']:{tuple(x['new_support']):x for x in s['crossings']}for s in cone['clipping_stages']}
 @lru_cache(maxsize=None)
 def ray(stage,ids):
  if stage==31:
   assert ids in initial
   v=C.null_one([basecov[i]for i in ids]);vals=[C.dot(v,c)for c in basecov];g=next(x for x in vals if x)
   if R.sg(g)>0:v=[-x for x in v];g=-g
   return tuple(D.lift(x/-g)for x in v)
  if stage not in ids:return ray(stage-1,ids)
  rec=crossings[stage][ids];v=ray(stage-1,tuple(rec['negative_support']));w=ray(stage-1,tuple(rec['positive_support']))
  gv,gw=D.dot(v,cov[stage]),D.dot(w,cov[stage]);assert D.sign(gv)<0<D.sign(gw)
  out=vecsub(vecscale(v,gw),vecscale(w,gv));assert D.sign(D.dot(out,cov[stage]))==0
  return out
 return local,cone,Q,H,N,cov,ray
def relation(g,qn):
 sg=D.sign(g)
 if not sg:return dict(kind='right',label=2)
 if sg>0:return dict(kind='obtuse')
 square=D.mul(g,g)
 if D.sign(D.sub(square,qn))>0:return dict(kind='dotted')
 for m,c in[(3,O),(4,2*O),(5,(3*O+R.BASIS[1])/2),(6,3*O),(10,(5*O+R.BASIS[1])/2)]:
  if D.sign(D.sub(D.scale(square,4*O),D.scale(qn,c)))==0:return dict(kind='finite',label=m)
 lo=D.sub(D.scale(square,4*O),qn);hi=D.sub(D.scale(qn,2*O),D.scale(square,4*O))
 if D.sign(lo)>0 and D.sign(hi)>0:return dict(kind='forbidden_finite_gap',neighbor_labels=[3,4])
 raise AssertionError('Unresolved relation; no alphabet assumption is permitted')
def encode(x):
 if isinstance(x,dict):return{k:encode(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):
  if len(x)==8 and all(hasattr(a,'to_list')for a in x):return D.encfield(x)
  return[encode(a)for a in x]
 if hasattr(x,'to_list'):return R.enc(x)
 return x
def run():
 scout=read(SCOUT);assert [r['ray_index']for r in scout['proposals']]==RAY_IDS
 local,cone,Q,H,N,cov,ray=geometry();hcov=tuple(D.mv(Q,H));hh=D.dot(H,hcov)
 caps=[];norms=[];capcov=[];records=[]
 for index,ray_id in enumerate(RAY_IDS):
  ids=cone['complete_census']['hyperideal_rays'][ray_id];v=ray(35,tuple(ids));vc=tuple(D.mv(Q,v));q=D.dot(v,vc);hp=D.dot(v,hcov)
  assert D.sign(q)>0 and D.sign(hp)<0
  vals=[D.dot(v,c)for c in cov];rels=[relation(g,D.mul(q,D.dot(n,c)))for g,n,c in zip(vals,N,cov)]
  assert all(x['kind']in['right','finite','dotted']for x in rels)
  assert [i for i,g in enumerate(vals)if D.sign(g)==0]==ids
  caps.append(v);norms.append(q);capcov.append(vc)
  records.append(dict(index=index,ray_index=ray_id,right_support=ids,normal=v,norm=q,center_pairing=hp,old_pairings=vals,old_relations=rels,counts=dict(Counter(r['kind']for r in rels))))
  print('exact polar',index,ray_id,records[-1]['counts'],flush=True)
 pairs=[];bad=[]
 for i,j in combinations(range(15),2):
  g=D.dot(caps[i],capcov[j]);rel=relation(g,D.mul(norms[i],norms[j]));pairs.append(dict(pair=[i,j],pairing=g,**rel))
  if rel['kind']not in['right','finite','dotted']:bad.append([i,j,rel['kind']])
 assert bad==[[1,12,'obtuse'],[1,13,'obtuse'],[4,14,'forbidden_finite_gap']]
 assert all(p['kind']in['right','finite','dotted']for p in pairs if max(p['pair'])<12)
 # All remaining old finite vertices are retained; the geometric preservation
 # theorem also follows from acuteness, but these signs are evaluated directly.
 oldpoints=[tuple(D.flatten_raw(x)for x in v['point'])for v in local['finite_vertices']]
 assert len(oldpoints)==100
 for i in BATCH:assert all(dot_sign(p,capcov[i])<0 for p in oldpoints)
 N48=N+[caps[i]for i in BATCH];C48=cov+[capcov[i]for i in BATCH]
 q48=[D.dot(v,c)for v,c in zip(N48,C48)];hp48=[D.dot(H,c)for c in C48]
 assert all(D.sign(q)>0 for q in q48) and all(D.sign(g)<0 for g in hp48)
 facets=[]
 for i,(n,q,hp)in enumerate(zip(N48,q48,hp48)):
  V=vecsub(vecscale(H,q),vecscale(n,hp));vv=D.sub(D.mul(D.mul(q,q),hh),D.mul(q,D.mul(hp,hp)))
  assert D.sign(vv)<0 and D.sign(D.mul(q,hh))<0
  # q*H-hp*n is future and has exactly one zero; the off-wall sign
  # follows from q>0,hp_j<0,hp_i<0 and every g_ij<=0.
  assert D.sign(D.dot(V,C48[i]))==0
  for j in range(48):
   if j!=i:assert dot_sign(n,C48[j])<=0
  facets.append(dict(wall=i,point=V,norm=vv))
 # The first uncut base-field ray gives a strict escape against the entire
 # fifteen-candidate bank, including its mutually incompatible proposals.
 escape_id=16;ids=cone['complete_census']['hyperideal_rays'][escape_id]
 assert ids==[0,1,2,4,5,6,8,17,27]
 V=ray(35,tuple(ids));qh=D.dot(V,hcov);assert not any(qh[1:]);V=tuple(D.scale(x,-O/qh[0])for x in V)
 assert D.sign(D.pair(V,V,Q))>0
 allcov=cov+capcov;vp=[D.dot(V,c)for c in allcov];assert all(D.sign(x)<=0 for x in vp)
 for power in range(13,31):
  E=vecadd(V,vecscale(H,D.lift(K.convert(F(1,1<<power)))));eq=D.pair(E,E,Q)
  if D.sign(eq)>0:break
 else:raise AssertionError('Escape perturbation failed')
 ep=[D.dot(E,c)for c in allcov];assert all(D.sign(x)<0 for x in ep) and D.sign(D.dot(E,hcov))<0
 oldE=tuple(D.lift(A.B.db(x))for x in local['strict_escape']['point']);oldcuts=[D.sign(D.dot(oldE,c))for c in capcov];assert oldcuts[4]==1
 # Exhaustive exact sign scoring on the already complete old cone.
 allrays=[ray(35,tuple(ids))for ids in cone['complete_census']['hyperideal_rays']]
 cuts=[]
 for i,c in enumerate(capcov):
  signs=[dot_sign(v,c)for v in allrays];positive=[j for j,s in enumerate(signs)if s>0]
  assert positive==scout['proposals'][i]['cut_rays']
  cuts.append(dict(candidate=i,positive_rays=positive,counts={str(k):v for k,v in sorted(Counter(signs).items())}))
  print('exact coverage',i,len(positive),flush=True)
 union=sorted(set().union(*(set(cuts[i]['positive_rays'])for i in BATCH)))
 allunion=sorted(set().union(*(set(x['positive_rays'])for x in cuts)))
 assert len(union)==541 and len(allunion)==545
 paths=[f'scripts/{TAG}.py',f'scripts/{D.TAG}.py',PRIMARY,CONE,SCOUT,D.START,D.START_REVIEW]
 return encode(dict(status='PASS_EXACT_FIFTEEN_POLARS_AND_COMPATIBLE_TWELVE_CAP_BATCH',sources={p:sha(p)for p in paths},
  field_encoding='Eight coefficients over frozen octic field; bit0=sqrt(1+sqrt5),bit1=sqrt(24+72sqrt5/5),bit2=sqrt(10+68sqrt5/15).',form=Q,center=H,candidates=records,mutual_pairs=pairs,incompatible_pairs=bad,batch_indices=BATCH,batch_normals=N48,batch_norms=q48,actual_facets=48,facet_witnesses=facets,strict_old_finite_vertex_preservation_count=100,
  old_escape_candidate_signs=oldcuts,coverage=cuts,batch_cut_union=union,all_fifteen_cut_union=allunion,
  common_strict_escape=dict(point=E,norm=eq,pairings_all_fiftyone=ep,ray=V,ray_index=escape_id,active=ids,epsilon=f'1/{1<<power}'),
  scope='Exact15 individually compatible polar walls and one48-facet chamber retaining12 mutuallycompatible caps. This bank is initially numerically selected; failed polar directions need a separate all-label gap proof. All15 share a strict escape, so every subset and the48-facet batch have infinite volume. Complete new cone is separate.'))
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved exact polar batch replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],len(d['batch_cut_union']),flush=True)
