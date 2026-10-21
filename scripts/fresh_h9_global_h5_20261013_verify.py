#!/usr/bin/env python3
"""Rigorous common escape for exact formula counterparts of saved H5 banks."""
import argparse,json,sys
from pathlib import Path
from fractions import Fraction as Q
from functools import lru_cache
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_h5_closure_next_20261013_verify as C
from scripts import fresh_h9_global_h5_20261013_schema as SC
F=C.F;V=C.V;K=C.K;N=C.N
PRE='data/enumerations/fresh_h9_global_h5_20261013'
FIXED=['seed2','compact1','compact2','compact3']
def atan_bounds(q,n):
 s=sum((Q(-1)**j)*q**(2*j+1)/Q(2*j+1)for j in range(n));tail=q**(2*n+1)/Q(2*n+1)
 return (s-tail,s)if n%2 else(s,s+tail)
@lru_cache(None)
def cos2(m):
 try:return F.cos2(m)
 except ValueError:pass
 a,b=atan_bounds(Q(1,5),100);c,d=atan_bounds(Q(1,239),32);pi=(V.rat(16*a-4*d)[0],V.rat(16*b-4*c)[1]);x=V.div(pi,V.rat(m));xx=V.mul(x,x);term=V.rat(1);out=term
 # Alternating Taylor remainder, valid throughout 0<=x<=pi/2.
 for k in range(1,61):
  term=V.div(V.neg(V.mul(term,xx)),V.rat((2*k-1)*(2*k)));out=V.add(out,term)
 rem=Q(2)**122
 for k in range(1,123):rem/=k
 eps=V.rat(rem)[1];return V.mul(V.rat(2),(out[0]-eps,out[1]+eps))
def run(verify=False):
 G,I,h,g,R,names,center,pairs,old,sources=C.old_state();latest=V.read(C.PRE+'_verify.json')
 for p,s in latest['sources'].items():assert V.sha(p)==s,p
 for row in latest['new_caps']:
  name=row['specification']['name']
  for x in row['new_pairs']:pairs[F.key(name,x['wall'])]=x['relation']
  R.append([C.O.frozen_interval(x)for x in row['normal']]);names.append(name)
 assert len(R)==125
 data=V.read(PRE+'_recipes.json')
 for p,s in data['sources'].items():assert V.sha(p)==s,p
 y=tuple(map(K.parse,data['exact_escape_proposal']));q=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(future)<0
 iy=list(map(V.golden,y));oldps=[V.dot(r,iy)for r in R];assert all(p[0]>0 for p in oldps)
 x=[V.golden(K.neg(z))for z in K.mv(I,y)];fixed=[R[names.index(n)]for n in FIXED]
 for i,a in enumerate(FIXED):
  for b in FIXED[:i]:assert pairs[F.key(a,b)]=={'kind':'finite','label':2}
 coefficients=[V.div(V.pair(x,g,e),V.rat(2))for e in fixed];assert all(t[1]<0 for t in coefficients)
 xf=[V.sub(a,b)for a,b in zip(x,N.combine(coefficients,fixed),strict=True)];fq=V.pair(xf,g,xf);face_future=V.dot(xf,list(map(V.golden,h)));assert fq[0]>0 and face_future[0]>0
 faceps=[V.pair(xf,g,r)for n,r in zip(names,R,strict=True)if n not in FIXED];assert all(p[1]<0 for p in faceps)
 extras=[]
 anchorbank=V.read('data/enumerations/fresh_h9_h5_closure_next_20261013_scout.json')['candidates']
 assert [s['source_row']for s in data['extra_anchors']]==[7315,4818,9678]
 assert [s['name']for s in data['extra_anchors']]==['alternate_endpoint0','alternate_endpoint1','alternate_endpoint2']
 for spec in data['extra_anchors']:
  literal=anchorbank[spec['source_row']];assert spec['incident']==FIXED+literal['edge']+[literal['anchor']] and spec['profile']==[2]*4+literal['profile']+[literal['label']] and spec['branch']==literal['branch']
  r,proof=C.make(spec,g,R,names,center,pairs);assert V.dot(r,list(map(V.golden,h)))[0]>0
  ff=dict(zip(spec['incident'],spec['profile'],strict=True))
  for n in names:
   if n in ff:pairs[F.key(spec['name'],n)]={'kind':'finite','label':ff[n]}
   else:pairs[F.key(spec['name'],n)]={'kind':'dotted'}
  R.append(r);names.append(spec['name']);extras.append({'specification':spec,'projection':proof})
 lookup=dict(zip(names,R));cache={};counts=Counter();maxpair={};branch_counts=Counter();selected_bound=-10**1000;matching=0
 source_rows={};seen={};bankkeys=dict(SC.BANKS)
 for row in data['recipes']:
  tag=row['source']
  if tag not in source_rows:
   dd=V.read('data/enumerations/'+tag+'.json');source_rows[tag]=dd[bankkeys[tag]];seen[tag]=set()
  assert row['row'] not in seen[tag];seen[tag].add(row['row'])
  literal=source_rows[tag][row['row']];assert (row['incident'],row['profile'])==SC.canonical(tag,row['row'],literal)
  ss=row['incident'];ll=row['profile'];sg=row['branch'];assert len(ss)==len(ll)==9 and len(set(ss))==9 and sg in(-1,1)
  key=tuple(ss)
  if key not in cache:
   E=[lookup[n]for n in ss];M=[]
   for i,a in enumerate(ss):
    rr=[]
    for j,b in enumerate(ss):
     rel=pairs.get(F.key(a,b))
     rr.append(V.rat(2)if i==j else V.neg(cos2(rel['label']))if rel and rel['kind']=='finite'else V.pair(E[i],g,E[j]))
    M.append(rr)
   MI=F.matinv(M);cp=[V.pair(center,g,a)for a in E];mc=V.mv(MI,cp);qz=V.sub(V.pair(center,g,center),V.dot(cp,mc));assert qz[0]>0 or qz[1]<0
   ex=[V.neg(V.dot(a,iy))for a in E];mx=V.mv(MI,ex);zx=V.sub(V.neg(V.dot(center,iy)),V.dot(cp,mx));proj=N.combine(mc,E);z=[V.sub(a,b)for a,b in zip(center,proj,strict=True)]
   cache[key]=(E,MI,cp,mc,qz,mx,zx,z)
  E,MI,cp,mc,qz,mx,zx,z=cache[key];b=[V.neg(cos2(m))for m in ll];co=V.mv(MI,b);qp=V.dot(b,co);rad=V.div(V.sub(V.rat(2),qp),qz);assert rad[0]>0,(tag,row['row'],'radicand',V.enc(rad));rt=V.sqrt(rad);rt=rt if sg==1 else V.neg(rt)
  pc=V.add(V.dot(b,mc),V.mul(rt,qz));assert pc[1]<0,(tag,row['row'],'orientation',V.enc(pc))
  px=V.add(V.dot(b,mx),V.mul(rt,zx));assert px[1]<0,(tag,row['row'],'target',V.enc(px))
  # The decimal proposal fixes a branch by the sign of its component
  # perpendicular to the nine-contact span. This is a rigorous rational
  # sign check, not a claim that a rounded vector itself has exact norm two.
  rounded=[V.rat(Q(str(t)))for t in literal['normal']]
  rz=V.sub(V.pair(rounded,g,center),V.dot([V.pair(rounded,g,e)for e in E],mc))
  selector=V.div(rz,qz);assert selector[0]>0 if sg==1 else selector[1]<0,(tag,row['row'],'branch selector')
  part=N.combine(co,E);r=[V.add(a,V.mul(rt,t))for a,t in zip(part,z,strict=True)]
  matching+=1
  # Most profiles prescribe all four finite contacts. The earliest roof scout
  # frees one of these walls; certify its remaining pairing is negative too.
  for n,e in zip(FIXED,fixed,strict=True):
   if n in ss:assert ll[ss.index(n)]>=2
   else:assert V.pair(r,g,e)[1]<0
  counts[tag]+=1;branch_counts[str(sg)]+=1;maxpair[tag]=max(maxpair.get(tag,-10**1000),px[1]);selected_bound=max(selected_bound,px[1])
  if sum(counts.values())%5000==0:print('verified',sum(counts.values()),'groups',len(cache),flush=True)
 assert sum(counts.values())==38209 and matching==38209
 assert set(counts)==set(bankkeys)
 assert all(seen[tag]==set(range(len(source_rows[tag])))for tag in counts)
 out={'status':'PASS_EXACT_COMMON_ESCAPE_ALL_38209_LITERAL_FORMULA_CAPS','sources':{p:V.sha(p)for p in list(data['sources'])+sources+[C.PRE+'_verify.json',PRE+'_recipes.json','scripts/'+Path(__file__).name,'scripts/'+Path(SC.__file__).name]},'retained_walls':125,'literal_caps':sum(counts.values()),'distinct_incidence_groups':len(cache),'formula_branch_definition':'Sign of the literal decimal proposal projected into the one-dimensional orthogonal complement; certified using rational interval arithmetic.','formula_branch_selectors_checked':matching,'extra_anchor_recipes':extras,'bank_counts':dict(counts),'branch_counts':dict(branch_counts),'maximum_target_pairing_upper':str(selected_bound),'bank_target_pairing_upper_bounds':{k:str(v)for k,v in maxpair.items()},'strict_escape':{'covector':data['exact_escape_proposal'],'norm':K.enc(q),'future_pairing':K.enc(future),'all125covector_pairings':list(map(V.enc,oldps))},'H5_escape':{'point':list(map(V.enc,xf)),'norm':V.enc(fq),'future_covector_pairing':V.enc(face_future),'old_nonfixed_pairings':list(map(V.enc,faceps)),'candidate_feasibility_proof':'Every candidate is negative at the strict H9 escape and nonpositive at the four fixed walls. Its pairing only decreases under the saved projection, since all projection coefficients are negative.'},'scope':'Finite literal formula bank only. Conflicts are ignored and optional deletions only enlarge the intersection. Every exact norm-two formula branch in the bank has a common strict spacelike escape. This neither excludes other angle profiles nor all completions of the retained seed.'}
 if verify:assert out==V.read(PRE+'_verify.json')
 else:(V.ROOT/(PRE+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');run(a.parse_args().verify)
