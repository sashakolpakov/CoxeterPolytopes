#!/usr/bin/env python3
"""Explicit simultaneous completion of both F4 cusp links of exact162.

Select the first certified positive integer common height. Every old gate,
cross-cusp gate, center, affine relation and remaining escape is certified.
"""
import argparse,hashlib,json
from fractions import Fraction as F
from pathlib import Path
import sympy as sp
import fresh_h20_e8_exact_cluster_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_e8_f4_cusps_20261017'
INPUT='results/fresh_h20_e8_exact_cluster_20261017.json';GEOM=P.GEOM;CUSP=P.CUSP
pt,add,sub,neg,mul,div,sqrt=P.pt,P.add,P.sub,P.neg,P.mul,P.div,P.sqrt
save,vec,total,inner,mv,alg=P.save,P.vec,P.total,P.inner,P.mv,P.alg
SEED=['-2638991277/1250000000','5464275411/2500000000','-16560862307/10000000000','14082107/9765625']
def run():
 D=json.loads((ROOT/INPUT).read_text());G=json.loads((ROOT/GEOM).read_text());C=json.loads((ROOT/CUSP).read_text())
 assert D['status']=='PASS_EXACT_162_WALL_E8_TRIANGLE_CLUSTER'
 for doc in(D,G,C):
  for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
 N=[[tuple(map(F,x))for x in v]for v in G['normals']]+[[tuple(map(F,x))for x in r['physical_normal']]for r in D['steps']]
 qs=list(map(F,G['norms']))+[F(2)]*4;H=[tuple(map(F,x))for x in D['strict_center']]
 Q=[list(map(F,row))for row in G['form']];qr=[[(j,q)for j,q in enumerate(row)if q]for row in Q]
 def pair(v,w):return total(mul(pt(q),mul(v[i],w[j]))for i,row in enumerate(qr)for j,q in row)
 def lin(vs,cs):return[total(mul(c,v[i])for v,c in zip(vs,cs))for i in range(21)]
 def plus(v,w):return[add(a,b)for a,b in zip(v,w)]
 def scale(v,c):return[mul(c,x)for x in v]
 def negative_dotted(g,q):return g[1]<0 and g[1]**2>2*q
 ids=C['core_ids'];roots=[N[j]for j in ids];nus=[];bases=[];domains=[];upper=F(1)
 for ci,c in enumerate(C['cusps']):
  nu=lin(roots,[alg(x)for x in c['null_coefficients']]);nus.append(nu)
  B=[lin(roots,[alg(x)for x in row])for row in c['base_new_root_coefficients']];bases.append(B)
  inc=c['incident_core_ids'];np=[pair(nu,n)for n in N];nh=pair(nu,H);assert nh[1]<0
  for bi,b in enumerate(B):
   bp=[pair(b,n)for n in N];bh=pair(b,H);bounds=[]
   for j in range(162):
    if j in inc:continue
    assert np[j][1]<0
    v=div(add(add(bp[j],sqrt(pt(2*qs[j]))),pt(1)),neg(np[j]));upper=max(upper,v[1]);bounds.append(dict(wall=j,sufficient_height=save(v)))
   h=div(add(bh,pt(1)),neg(nh));upper=max(upper,h[1])
   domains.append(dict(cusp=ci,base=bi,incident=inc,old_base=bp,old_slope=np,center_base=bh,center_slope=nh,height_bounds=bounds,center_height_bound=save(h)))
 assert P.R.overlap(pair(nus[0],nus[1]),pt(-1))
 cross=[]
 for i,b in enumerate(bases[0]):
  for j,c in enumerate(bases[1]):
   linear=add(pair(b,nus[1]),pair(nus[0],c));constant=pair(b,c)
   bound=1+max(abs(x)for x in linear)+max(abs(x)for x in add(constant,pt(3)))
   upper=max(upper,bound);cross.append(dict(i=i,j=j,linear=linear,constant=constant,sufficient_height=str(bound)))
 guaranteed=P.R.ceil(upper)+1;trials=[];height=None
 for h in range(1,guaranteed+1):
  failures=[]
  for k,d in enumerate(domains):
   hc=add(d['center_base'],mul(pt(h),d['center_slope']))
   if hc[1]>=0:failures.append(['center',k,save(hc)])
   for j,(b,n)in enumerate(zip(d['old_base'],d['old_slope'])):
    if j in d['incident']:continue
    g=add(b,mul(pt(h),n))
    if not negative_dotted(g,qs[j]):failures.append(['old',k,j,save(g)])
  for c in cross:
   g=add(pt(-h*h),add(mul(pt(h),c['linear']),c['constant']))
   if g[1]>=-2:failures.append(['cross',c['i'],c['j'],save(g)])
  trials.append(dict(height=h,failure_count=len(failures),first_failure=failures[:1]))
  if not failures:height=h;break
 assert height is not None
 added=[plus(b,scale(nus[ci],pt(height)))for ci,B in enumerate(bases)for b in B]
 contacts=[];counts={}
 for k,r in enumerate(added):
  ci,bi=divmod(k,3);inc=set(C['cusps'][ci]['incident_core_ids']);one=[15,C['cusps'][ci]['A6_ids'][0],C['cusps'][ci]['A6_ids'][-1]][bi]
  assert P.R.overlap(pair(r,r),pt(2))
  rows=[]
  for j,n in enumerate(N+added[:k]):
   g=pair(r,n)
   if j<162:
    if j in inc:rel='pi/3'if j==one else'right';assert P.R.overlap(g,pt(-1 if j==one else 0))
    else:rel='dotted';assert negative_dotted(g,qs[j])
   else:
    cj,bj=divmod(j-162,3)
    if cj==ci:rel='pi/3'if{bi,bj}=={1,2}else'right';assert P.R.overlap(g,pt(-1 if rel=='pi/3'else 0))
    else:rel='dotted';assert g[1]<-2
   counts[rel]=counts.get(rel,0)+1;rows.append(dict(wall=j,relation=rel,pair=save(g)))
  contacts.append(dict(root=162+k,cusp=ci,base=bi,height=height,physical_normal=vec(r),all_prior_contacts=rows))
 full=N+added;norms=qs+[F(2)]*6;hp=[pair(H,n)for n in full];assert all(g[1]<0 for g in hp)
 facets=[sub(pair(H,H),div(mul(g,g),pt(q)))for g,q in zip(hp,norms)];assert all(g[1]<0 for g in facets)
 completed=[]
 for ci,c in enumerate(C['cusps']):
  nu=nus[ci];inc=sorted(c['incident_core_ids']+list(range(162+3*ci,165+3*ci)));ps=[pair(nu,n)for n in full]
  for j,g in enumerate(ps):
   if j in inc:assert P.R.overlap(g,pt(0))
   else:assert g[1]<0
  # The identities are inherited exact vector formulas; overlap checks
  # separately validate the physical implementation of those formulas.
  theta=lin([N[j]for j in c['E8_ids']],[alg(x)for x in c['E8_highest_root']])
  re=plus(theta,added[3*ci]);ra=plus(lin([N[j]for j in c['A6_ids']],[pt(1)]*6),plus(added[3*ci+1],added[3*ci+2]))
  assert all(P.R.overlap(x,y)for x,y in zip(re,scale(nu,pt(height))))
  assert all(P.R.overlap(x,y)for x,y in zip(ra,scale(nu,pt(2*height))))
  completed.append(dict(cusp=ci,point=vec(nu),incident=inc,types=['affineF4','affineE8','affineA7'],Euclidean_rank=19,normal_span_rank=20,all168_pairs=vec(ps),positive_relations={'F4':'nu','E8':f'{height}*nu','A7':f'{2*height}*nu'}))
 oldcusps=[]
 for raw in G['old_cusps']:
  v=[tuple(map(F,x))for x in raw];ps=[pair(v,r)for r in added];assert all(g[1]<0 for g in ps);oldcusps.append(vec(ps))
 oldE=[tuple(map(F,x))for x in D['strict_escape']['point']];cuts=[pair(oldE,r)for r in added];assert any(g[0]>0 for g in cuts)
 # A rational residual seed comes from the already published numerical
 # exterior region, but all following signs are exact interval checks.
 T=list(range(1,11))+[12,13,14,15,61,87,151];TG=sp.eye(17)*2
 for path in[list(range(1,11)),[12,13,14,15],[61,87]]:
  for i,j in zip(path,path[1:]):TG[T.index(i),T.index(j)]=TG[T.index(j),T.index(i)]=-1
 J=[[pt(F(x))for x in row]for row in TG.inv().tolist()]
 def projectT(v):return plus(v,scale(lin([N[j]for j in T],mv(J,[pair(N[j],v)for j in T])),pt(-1)))
 HT=projectT(H);assert pair(HT,HT)[1]<0;HT=scale(HT,div(pt(1),sqrt(neg(pair(HT,HT)))))
 seed=lin([projectT(N[j])for j in(95,142,145,31)],[pt(F(x))for x in SEED]);assert pair(seed,seed)[0]>0
 for j in range(162):
  g=pair(seed,N[j])
  if j in T:assert P.R.overlap(g,pt(0))
  else:assert g[1]<0
 lo,hi=F(0),F(1);Eface=None;bisect=[]
 for iteration in range(120):
  e=plus(scale(HT,pt(1-hi)),scale(seed,pt(hi)));eq=pair(e,e);ps=[pair(e,n)for n in full]
  assert eq[0]>0
  if all(g[1]<0 for j,g in enumerate(ps)if j not in T):Eface=e;break
  mid=(lo+hi)/2;m=plus(scale(HT,pt(1-mid)),scale(seed,pt(mid)));qm=pair(m,m)
  assert qm[0]>0 or qm[1]<0
  bisect.append(dict(t=str(mid),norm=save(qm)))
  if qm[0]>0:hi=mid
  else:lo=mid
 assert Eface is not None
 dualsum=lin([N[j]for j in T],[total(row)for row in J]);escape=None
 for exponent in range(4,120):
  eps=F(1,2**exponent);e=plus(Eface,scale(dualsum,pt(-eps)));eq=pair(e,e);eh=pair(e,H);ps=[pair(e,n)for n in full]
  if eq[0]>0 and eh[1]<0 and all(g[1]<0 for g in ps):
   escape=dict(seed_coefficients=SEED,basis=[95,142,145,31],face_mix_t=str(hi),support_perturbation=str(eps),norm=save(eq),center_pair=save(eh),point=vec(e),all168_pairs=vec(ps),norm_bisection=bisect);break
 assert escape
 for ci,c in enumerate(C['T_face_blockers']):
  for r in added[3*ci:3*ci+3]:assert all(negative_dotted(pair(r,N[j]),qs[j])for j in c['nonincident_T_roots'])
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',INPUT,GEOM,CUSP,'results/fresh_h20_e8_partner_endpoints_20261017.json']
 return dict(status='PASS_EXACT_168_WALL_TWO_F4_CUSP_COMPLETION',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},selected_common_height=height,guaranteed_sufficient_integer=guaranteed,integer_trials=trials,selection_scope='First positive common integer certified to preserve the strict inherited center and make every old nonincident/cross-cusp contact strictly dotted; no minimum-real-height or unequal-height claim.',height_bounds=[dict(cusp=d['cusp'],base=d['base'],old=d['height_bounds'],center=d['center_height_bound'])for d in domains],new_walls=contacts,contact_counts=counts,actual_facets=168,center=vec(H),all_center_pairs=vec(hp),facet_norms=vec(facets),completed_cusps=completed,five_old_cusp_new_pairs=oldcusps,proper_cusps_at_least=7,old_defective_escape_new_pairs=vec(cuts),new_finite_vertices=0,new_finite_exclusion='A new finite vertex could use only its own22-root affineF4+affineE8+affineA7 diagram: all other walls are dotted to the selected new wall. Any PD subset omits at least one root from each affine component, hence has at most19roots, below the20 needed. Cross-cusp new roots are dotted. Therefore no new finite vertex occurs.',known_finite_vertices_at_least=D['known_finite_vertices_at_least'],strict_escape=escape,unchanged_T_face=True,T_face_blockers=C['T_face_blockers'],scope='Six explicit compatible walls complete both deficient links and cut the old defective-cusp escape. All168facets and at least7proper cusps are certified. No new finite vertex; the targetT-face is unchanged and a new strict spacelike escape proves infinite volume. This is an actual local cusp completion, not a global finite-volume construction.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],'height',d['selected_common_height'],'escape',d['strict_escape']['norm'],flush=True)
