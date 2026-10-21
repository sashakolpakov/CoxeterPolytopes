#!/usr/bin/env python3
"""Four exact escapes and all-label conflicts exclude the enlarged H5 bank."""
import argparse,json,sys
from pathlib import Path
from collections import Counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_global_h5_20261013_verify as B
C=B.C;F=B.F;V=B.V;K=B.K;N=B.N;SC=B.SC;PRE=B.PRE
def run(verify=False):
 prior=V.read(PRE+'_verify.json');assert prior['status']=='PASS_EXACT_COMMON_ESCAPE_ALL_38209_LITERAL_FORMULA_CAPS'
 for p,s in prior['sources'].items():assert V.sha(p)==s,p
 data=V.read(PRE+'_enlarged_recipes.json')
 for p,s in data['sources'].items():assert V.sha(p)==s,p
 original=V.read(PRE+'_recipes.json');G,I,h,g,R,names,center,pairs,_,sources=C.old_state();last=V.read(C.PRE+'_verify.json')
 for p,s in last['sources'].items():assert V.sha(p)==s,p
 for row in last['new_caps']:
  name=row['specification']['name']
  for x in row['new_pairs']:pairs[F.key(name,x['wall'])]=x['relation']
  R.append([C.O.frozen_interval(z)for z in row['normal']]);names.append(name)
 assert len(R)==125;oldR=list(R);oldnames=list(names)
 for spec in original['extra_anchors']:
  r,_=C.make(spec,g,R,names,center,pairs);ff=dict(zip(spec['incident'],spec['profile'],strict=True))
  for n in names:pairs[F.key(spec['name'],n)]={'kind':'finite','label':ff[n]}if n in ff else{'kind':'dotted'}
  R.append(r);names.append(spec['name'])
 lookup=dict(zip(names,R));groups={}
 def group(ss):
  key=tuple(ss)
  if key not in groups:
   E=[lookup[n]for n in ss];M=[]
   for i,a in enumerate(ss):
    rr=[]
    for j,b in enumerate(ss):
     rel=pairs.get(F.key(a,b));rr.append(V.rat(2)if i==j else V.neg(B.cos2(rel['label']))if rel and rel['kind']=='finite'else V.pair(E[i],g,E[j]))
    M.append(rr)
   MI=F.matinv(M);cp=[V.pair(center,g,e)for e in E];mc=V.mv(MI,cp);qz=V.sub(V.pair(center,g,center),V.dot(cp,mc));assert qz[0]>0 or qz[1]<0
   z=[V.sub(a,b)for a,b in zip(center,N.combine(mc,E),strict=True)];groups[key]={'E':E,'MI':MI,'cp':cp,'mc':mc,'qz':qz,'z':z,'targets':{}}
  return groups[key]
 def expression(row):
  ss=row['incident'];gr=group(ss);b=[V.neg(B.cos2(m))for m in row['profile']];co=V.mv(gr['MI'],b);qp=V.dot(b,co);rad=V.div(V.sub(V.rat(2),qp),gr['qz']);assert rad[0]>0
  rt=V.sqrt(rad);sg=row['branch'];assert sg in(-1,1);rt=rt if sg==1 else V.neg(rt);pc=V.add(V.dot(b,gr['mc']),V.mul(rt,gr['qz']));assert pc[1]<0
  return{'spec':row,'group':gr,'b':b,'co':co,'rt':rt,'center':pc}
 def root(obj):
  if 'root'in obj:return obj['root']
  gr=obj['group'];obj['root']=[V.add(a,V.mul(obj['rt'],b))for a,b in zip(N.combine(obj['co'],gr['E']),gr['z'],strict=True)];return obj['root']
 entries=[{'root':r,'name':n}for n,r in zip(oldnames,oldR,strict=True)];source_ids={}
 for i,row in enumerate(original['recipes']):
  entries.append(expression(row));source_ids[(row['source'],row['row'])]=125+i
 offset=len(entries);assert offset==38334
 assert len(data['new_recipes'])==119
 for row in data['new_recipes']:
  old=original['recipes'][source_ids[(row['source'],row['row'])]-125];assert row['incident']==old['incident']and row['profile']==old['profile']and row['branch']==-old['branch'];entries.append(expression(row))
 assert len(entries)==38453
 targets=[]
 for i,w in enumerate(data['selected']):
  y=tuple(map(K.parse,w['exact_escape_proposal']));q=K.pair(y,I,y);future=K.pair(y,I,h);assert K.sign(q)>0 and K.sign(future)<0;iy=list(map(V.golden,y));x=[V.golden(K.neg(z))for z in K.mv(I,y)];fixed=[lookup[n]for n in B.FIXED];coef=[V.div(V.pair(x,g,e),V.rat(2))for e in fixed];assert all(z[1]<0 for z in coef)
  xf=[V.sub(a,b)for a,b in zip(x,N.combine(coef,fixed),strict=True)];fq=V.pair(xf,g,xf);ft=V.dot(xf,list(map(V.golden,h)));assert fq[0]>0 and ft[0]>0
  targets.append({'iy':iy,'record':{'covector':w['exact_escape_proposal'],'norm':K.enc(q),'future_pairing':K.enc(future),'H5_projected_norm':V.enc(fq),'H5_future_covector_pairing':V.enc(ft),'fixed_projection_coefficients':list(map(V.enc,coef))},'cutters':set(w['cutters']),'max_noncut':-10**1000,'min_cut':10**1000})
 assert len(targets)==4
 def target_pair(obj,ti):
  if 'name'in obj:return V.neg(V.dot(obj['root'],targets[ti]['iy']))
  gr=obj['group']
  if ti not in gr['targets']:
   ex=[V.neg(V.dot(e,targets[ti]['iy']))for e in gr['E']];mx=V.mv(gr['MI'],ex);zx=V.sub(V.neg(V.dot(center,targets[ti]['iy'])),V.dot(gr['cp'],mx));gr['targets'][ti]=(mx,zx)
  mx,zx=gr['targets'][ti];return V.add(V.dot(obj['b'],mx),V.mul(obj['rt'],zx))
 for ei,obj in enumerate(entries):
  for ti,t in enumerate(targets):
   px=target_pair(obj,ti)
   if ei in t['cutters']:assert px[0]>0,('cutter',ti,ei,V.enc(px));t['min_cut']=min(t['min_cut'],px[0])
   else:assert px[1]<0,('noncutter',ti,ei,V.enc(px));t['max_noncut']=max(t['max_noncut'],px[1])
  if ei%5000==0:print('sign rows',ei,flush=True)
 # Every finite Coxeter angle has pairing -2cos(pi/m); an interval strictly
 # between two consecutive such values excludes all integer labels at once.
 def conflict(p):
  if p[0]>0:return{'kind':'positive_pairing','pairing':V.enc(p)}
  for m in range(2,513):
   upper=V.neg(B.cos2(m));lower=V.neg(B.cos2(m+1))
   if lower[1]<p[0]and p[1]<upper[0]:return{'kind':'between_consecutive_finite_angles','labels':[m,m+1],'pairing':V.enc(p),'lower_angle_pairing':V.enc(lower),'upper_angle_pairing':V.enc(upper)}
  raise AssertionError(('uncertified conflict',V.enc(p)))
 def pair_to_old(a,n):
  if n in a['spec']['incident']:return a['b'][a['spec']['incident'].index(n)]
  cache=a.setdefault('old_pair_cache',{})
  if n not in cache:
   gr=a['group'];ex=[]
   for nm,e in zip(a['spec']['incident'],gr['E'],strict=True):
    rel=pairs.get(F.key(nm,n));ex.append(V.rat(2)if nm==n else V.neg(B.cos2(rel['label']))if rel and rel['kind']=='finite'else V.pair(e,g,lookup[n]))
   mx=V.mv(gr['MI'],ex);zx=V.sub(V.pair(center,g,lookup[n]),V.dot(gr['cp'],mx));cache[n]=V.add(V.dot(a['b'],mx),V.mul(a['rt'],zx))
  return cache[n]
 def wallpair(a,b):
  if 'name'in b:return pair_to_old(a,b['name'])
  gr=a['group'];ex=[pair_to_old(b,n)for n in a['spec']['incident']];mx=V.mv(gr['MI'],ex);zx=V.sub(b['center'],V.dot(gr['cp'],mx));return V.add(V.dot(a['b'],mx),V.mul(a['rt'],zx))
 assignment={int(k):v for k,v in data['assignment'].items()};assert set(assignment)==set(range(119));proof=[]
 for j in range(119):
  ti=assignment[j];assert ti in range(4);assert offset+j not in targets[ti]['cutters'];obj=entries[offset+j];checks=[]
  for ci in sorted(targets[ti]['cutters']):
   try:check=conflict(wallpair(obj,entries[ci]))
   except AssertionError as exc:raise AssertionError(('branch/cutter',j,ci,*exc.args))from exc
   checks.append({'cutter_entry':ci,**check})
  proof.append({'required_new_branch':j,'witness':ti,'cutter_conflicts':checks})
  if j%20==0:print('conflict classes',j,flush=True)
 out={'status':'PASS_EXACT_38209_PLUS119_BANK_EXCLUSION_WITH_OPTIONAL_OLD_WALL_DELETION','sources':{p:V.sha(p)for p in list(data['sources'])+[PRE+'_enlarged_recipes.json','scripts/'+Path(__file__).name,'scripts/'+Path(B.__file__).name]},'old_optional_walls':125,'original_formula_caps':38209,'new_opposite_formula_caps':119,'total_dictionary_entries':len(entries),'complete_sign_checks':len(entries)*4,'new_branch_coverage':119,'conflict_checks':sum(len(r['cutter_conflicts'])for r in proof),'witnesses':[{**t['record'],'cutters':sorted(t['cutters']),'maximum_noncutter_pairing_upper':str(t['max_noncut']),'minimum_cutter_pairing_lower':str(t['min_cut'])}for t in targets],'coverage':proof,'argument':'If no new opposite branch is selected, the inherited common escape survives. Otherwise choose any selected new branch and its assigned witness. All possible cutters are exactly listed and conflict with that required branch; therefore every selected halfspace is strictly negative at the witness. This applies even after arbitrary optional old-wall deletions.','scope':'The finite literal dictionary, with exact algebraic formula branches, admits no finite-volume pairwise Coxeter-compatible selection. The H5 conclusion keeps the four named fixed walls as the face and uses nonpositive pairings with them. It does not exclude other normals, other finite profiles, or other geometric seeds.'}
 if verify:assert out==V.read(PRE+'_enlarged_verify.json')
 else:(V.ROOT/(PRE+'_enlarged_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'conflicts',out['conflict_checks'],flush=True);return out
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--verify',action='store_true');run(a.parse_args().verify)
