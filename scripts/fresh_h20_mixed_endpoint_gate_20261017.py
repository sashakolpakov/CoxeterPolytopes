#!/usr/bin/env python3
"""All-label mixed spherical/parabolic endpoint neighbor exclusion.

Inherits the complete endpoint domains and their two same-kind gates.
"""
import argparse,hashlib,json,math
from fractions import Fraction as F
import sympy as s
import fresh_h20_parabolic_endpoint_gate_20261017 as A
P,E=A.P,A.E;ROOT=A.ROOT;TAG='fresh_h20_mixed_endpoint_gate_20261017'
pt,add,sub,neg,mul,div,sqrt=A.pt,A.add,A.sub,A.neg,A.mul,A.div,A.sqrt
scale,plus,sx=A.scale,A.plus,A.sx

def run():
 B=A.build();G=B['G'];gp=B['gp'];dot=B['dot'];ids=B['ids']
 old=json.loads((ROOT/'results/fresh_h20_e8_endpoint_belt_20261017_exact_gate.json').read_text())
 low=json.loads((ROOT/'results/fresh_h20_all_endpoint_separation_20261017.json').read_text())
 pg=json.loads((ROOT/f'results/{A.TAG}.json').read_text())
 assert pg['status']=='PASS_EXACT_PARABOLIC_ENDPOINT_CONTINUOUS_HEIGHT_SEPARATION'
 reject={(r['edge'],r['profile'],r['label158'])for r in old['rejected_discretes']}
 lreject={(r['edge'],r['singleton_label'],r['label158'])for r in low['low_discrete_rejections']}
 assert P.alg(s.sqrt(2))[1]<F(29,20) and P.alg((1+s.sqrt(5))/2)[0]>F(29,20)
 # The Coxeter magnitudes are 0,1,sqrt2,phi,sqrt3,...,2 and then >2.
 # Combine this exact order with the inherited continuous necessary bounds.
 starts={(0,2):(1+s.sqrt(5))/2,(1,0):s.Integer(0),(1,1):s.Integer(1),(1,7):s.Integer(1),(2,0):s.sqrt(2)}
 para=[];domain=[]
 for f in B['families']:
  start=starts[f['ci'],f['profile']];lo=plus(f['base'],scale(f['nu'],P.alg(start)))
  dis=[]
  if f['ci']==2:dis=[dict(label=3,r=plus(f['base'],f['nu']))]
  para.append(dict(ci=f['ci'],id=f['profile'],contacts=f['contacts'],lo=lo,nu=f['nu'],dis=dis))
  domain.append(dict(edge=f['ci'],profile=f['profile'],continuous_start=str(start),separate_discrete158_labels=[r['label']for r in dis]))
 finite=[];same_edge_inverses=[];ee=s.eye(21)[:,ids.index(158)]
 xlabels=[(2,s.Integer(0)),(3,s.Integer(1)),(4,s.sqrt(2)),(5,(1+s.sqrt(5))/2),(6,s.sqrt(3))]
 for edge in B['edge_data']:
  ci=edge['ci'];cen=edge['cen'];SB=edge['SB'];WJ=edge['WJ'];sigma=edge['sigma'];V=edge['V'];unit=edge['unit']
  assert all(x.is_nonnegative for x in WJ)
  same_edge_inverses.append(dict(edge=ci,inverse_entrywise_nonnegative=True))
  groups=[('fixed',cen['profiles'],s.Rational(9,5)),('low',[dict(id=m,contacts=[[151,m]])for m in(3,4,5)],s.Rational(9,5)),('high',[dict(id=6,contacts=[[151,6]])],s.Integer(0))]
  for kind,rows,threshold in groups:
   for row in rows:
    dd=s.zeros(19,1)
    for j,m in row['contacts']:dd[cen['support'].index(j)]=2*s.cos(s.pi/m)
    pp=(-SB*WJ*dd).applyfunc(s.simplify);p=sx(pp);h=s.simplify(2-(pp.T*G*pp)[0]);assert h.is_positive
    off=s.simplify((pp.T*G*ee)[0]);t=div(P.alg(threshold+off),sqrt(P.alg(sigma*h)));assert t[0]>1
    k=add(t,sqrt(sub(mul(t,t),pt(1))))
    a=scale(plus(V,scale(unit,pt(-1))),div(sqrt(P.alg(h)),pt(2)))
    b=scale(plus(scale(V,pt(-1)),scale(unit,pt(-1))),div(sqrt(P.alg(h)),pt(2)))
    aa=scale(a,k);bb=scale(b,div(pt(1),k));nums=[plus(plus(p,aa),bb),plus(p,scale(aa,pt(2))),aa]
    dis=[]
    if kind!='high':
     for label,x in xlabels:
      if kind=='fixed'and(ci,row['id'],label)in reject:continue
      if kind=='low'and(ci,row['id'],label)in lreject:continue
      rad=s.simplify((x+off)**2/sigma-h);assert P.alg(rad)[0]>0
      r=plus(plus(p,scale(unit,neg(div(P.alg(x+off),sqrt(P.alg(sigma)))))),scale(V,sqrt(P.alg(rad))))
      dis.append(dict(label=label,r=r,gr=gp(r)))
    finite.append(dict(ci=ci,id=row['id'],kind=kind,nums=nums,gnums=[gp(v)for v in nums],dis=dis))
  print('FINITE_EDGE',ci,flush=True)
 stats={};maxima={};gaps=[];ex=[];stream=hashlib.sha256()
 def check(kind,co,record):
  stats[kind]=stats.get(kind,0)+1;upper=max(x[1]for x in co);maxima[kind]=max(maxima.get(kind,F(-10**99)),upper)
  stream.update(json.dumps([kind,record,[P.save(x)for x in co]],sort_keys=True).encode())
  if upper>=0:ex.append(dict(kind=kind,record=record,coefficients=[P.save(x)for x in co]))
 isolated=0
 for p in para:
  for f in finite:
   if p['ci']==f['ci']:continue
   if f['kind']=='high'and any(j==151 for j,m in p['contacts']):isolated+=1;continue
   target=0 if f['kind']=='high' else 2
   rec=[p['ci'],p['id'],f['ci'],f['kind'],f['id']]
   co=[add(P.inner(v,w),pt(target if a==0 and b<2 else 0))for a,v in enumerate([p['lo'],p['nu']])for b,w in enumerate(f['gnums'])]
   check(f['kind']+'_tail_tail',co,rec)
   for r in f['dis']:
    check(f['kind']+'_tail_discrete',[add(P.inner(p['lo'],r['gr']),pt(2)),P.inner(p['nu'],r['gr'])],rec+[r['label']])
   for d in p['dis']:
    check(f['kind']+'_discrete_tail',[add(P.inner(d['r'],w),pt(target if b<2 else 0))for b,w in enumerate(f['gnums'])],rec+[d['label']])
    for r in f['dis']:
     z=P.inner(d['r'],r['gr']);stats['discrete_discrete']=stats.get('discrete_discrete',0)+1
     if z[1]<-2:continue
     assert z[1]<0 and z[0]>-2
     co=mul(pt(F(-1,2)),z);lab=max(2,int(math.pi/math.acos(float(sum(co)/2))))
     lo,hi=E.cosine(lab),E.cosine(lab+1);assert lo[1]<co[0] and co[1]<hi[0]
     gaps.append(dict(record=rec+[d['label'],r['label']],labels=[lab,lab+1],pair=P.save(z),cosine=P.save(co),lower_cosine=P.save(lo),upper_cosine=P.save(hi)))
  print('PARABOLIC',p['ci'],p['id'],'exceptions',len(ex),'gaps',len(gaps),flush=True)
 assert not ex,ex[:1]
 assert len(gaps)==2
 paths=[f'scripts/{TAG}.py',f'scripts/{A.TAG}.py',f'results/{A.TAG}.json',A.CENSUS,A.END,P.CUSP,P.GEOM,
  'results/fresh_h20_e8_exact_cluster_20261017.json','results/fresh_h20_e8_f4_cusps_20261017.json',
  'results/fresh_h20_e8_endpoint_belt_20261017_exact_gate.json','results/fresh_h20_all_endpoint_separation_20261017.json',
  f'scripts/{P.TAG}.py',f'scripts/{E.TAG}.py']
 return dict(status='PASS_EXACT_ALL_LABEL_MIXED_FINITE_PARABOLIC_ENDPOINT_GATE',
  sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},parabolic_domains=domain,
  same_edge_checks=same_edge_inverses,high_nonorthogonal151_exclusions=isolated,case_counts=stats,
  strict_coefficient_maxima={k:str(v)for k,v in maxima.items()},coefficient_stream_sha256=stream.hexdigest(),
  finite_angle_gaps=gaps,exceptions=ex,
  scope='All complete finite and parabolic rank19-support endpoint profiles, retaining old168 and requiring every unprescribed old168 contact other than158 strictly dotted. Cross-edge mixed actual-neighbor pairs are excluded at all labels and heights; same-edge mixed pairs are obtuse by the positive-kernel/future-ray argument. Other contact patterns or intermediate facets are not covered.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
