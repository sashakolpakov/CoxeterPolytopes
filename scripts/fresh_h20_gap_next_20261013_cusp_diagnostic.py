#!/usr/bin/env python3
"""Exact all98 incidence and explicit recession at the new parallel pair18/96."""
import argparse,json
from fresh_h20_gap_next_20261013_exact import setup,read,sha,ROOT,OUT as SOURCE,kernel
from fresh_dimension_bridge_20261013_verify import E,Z,O,F,dot,mv,parse,vec
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,rank,fdot,fmv,same
from fresh_h20_edge_family_completion_20261013_cusp_repairs import base_root
from fresh_h20_edge_family_completion_20261013_exact import positive_solve
OUT='data/enumerations/fresh_h20_gap_next_20261013_cusp_diagnostic.json'

def pair(f,xc,g,r):
 if g is B:return f,fdot(f,xc,tuple(map(f.lift,r)))
 if g.d==f.d:return f,fdot(f,xc,r)
 H=Extension(f,f.lift(g.d));rr=tuple((f.lift(a),f.lift(b))for a,b in r)
 return H,fdot(H,tuple(map(H.lift,xc)),rr)

def run():
 Q,h,hc,roots=setup();d=read(SOURCE);assert d['wall_count']==98 and d['retained_walls']==96
 r18=base_root(*roots[18]);c1,c2=d['new_caps'];f=Extension(B,E(F(10)))
 assert c1['field']['radicand']==['10','0']and c2['field']['radicand']==['1','0']
 r1=tuple(map(f.parse,c1['normal']));assert all(parse(x[1])==Z for x in c2['normal'])
 r2=tuple(parse(x[0])for x in c2['normal']);roots=roots+[(f,r1),(B,r2)]
 Qf=tuple(tuple(map(f.lift,row))for row in Q);r1c=fmv(f,Qf,r1)
 assert same(f,fdot(f,r1,r1c),f.rational(2))and dot(r18,mv(Q,r18))==E(F(2))
 assert same(f,fdot(f,r1c,tuple(map(f.lift,r18))),f.rational(-2))
 delta=tuple(f.add(x,f.lift(y))for x,y in zip(r1,r18,strict=True));dc=fmv(f,Qf,delta)
 assert f.sign(fdot(f,delta,dc))==0
 center=fdot(f,delta,tuple(map(f.lift,hc)));assert f.sign(center)<0
 pairs=[];incident=[]
 for i,(g,r)in enumerate(roots):
  H,p=pair(f,dc,g,r);assert H.sign(p)<=0
  if H.sign(p)==0:incident.append(i)
  pairs.append({'wall':i,'pairing':H.encode(p),'sign':H.sign(p),'outer_radicand':None if H is f else g.d.enc()})
 ids=[0,2,4,5,7,8,10,11,12,13,14,16,19,20,26,36]
 assert incident==sorted(ids+[18,96])
 rs=[base_root(*roots[i])for i in ids];G=[[dot(a,mv(Q,b))for b in rs]for a in rs]
 assert rank(G)==16
 pd=positive_solve(G,[Z]*16);assert pd is not None
 assert all(dot(r18,mv(Q,r))==Z and f.sign(fdot(f,r1c,tuple(map(f.lift,r))))==0 for r in rs)
 # Construct an actual nonzero Euclidean recession direction. Requiring
 # orthogonality to p1 and a+b separately keeps this vector in Qsqrt6.
 p1=vec(c1['polar_certificate']['vector']);a=base_root(*roots[3]);b=tuple(x/E(F(3))for x in base_root(*roots[42]))
 rows=[mv(Q,r)for r in rs+[r18,p1,tuple(x+y for x,y in zip(a,b,strict=True)),h]]
 v,piv=kernel(rows);vc=mv(Q,v);vq=dot(v,vc)
 assert vq.sign()>0 and dot(v,hc)==Z and f.sign(fdot(f,delta,tuple(map(f.lift,vc))))==0
 rp=[]
 for i in incident:
  g,r=roots[i];z=fdot(g,r,vc if g is B else tuple(map(g.lift,vc)));assert g.sign(z)==0
  rp.append({'wall':i,'pairing':g.encode(z)})
 # v is not a multiple of delta because v.h=0 whereas delta.h<0.
 # Both signs of v are recession directions of the incident horospherical
 # inequalities; a strict interior comes from the full chamber's center.
 return {'status':'EXACT_NEW_IDEAL_RAY_HAS_UNBOUNDED_HOROSPHERICAL_LINK','source':SOURCE,'source_sha256':sha(SOURCE),'ray_formula':'r18+r96','ray':list(map(f.encode,delta)),'ray_norm':f.encode(fdot(f,delta,dc)),'ray_center_pairing':f.encode(center),
         'all98_ray_pairings':pairs,'incident_walls':incident,'other_incident_gram':[[x.enc()for x in row]for row in G],'other_incident_rank':16,'other_incident_positive_pivots':[x.enc()for x in pd[1]],
         'parallel_block':[[2,-2],[-2,2]],'parallel_positive_kernel':[1,1],'parallel_kernel_sum':'delta','parallel_rank':1,'horosphere_dimension':19,'total_quotient_rank':17,
         'recession':{'vector':[x.enc()for x in v],'norm':vq.enc(),'center_pairing':dot(v,hc).enc(),'ray_pairing':f.encode(fdot(f,delta,tuple(map(f.lift,vc)))),'incident_pairings':rp,'construction_rank':len(piv)},
         'scope':'All98 null-ray pairings are explicitly recomputed in exact composita; the complete incident Gram has positive blockrank16 plus one affineA1 block. An explicit positive-norm vector orthogonal to the center, nullray and every incident normal gives a nonzero horospherical recession direction, proving the link unbounded. No assertion about other ideal rays.'}
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
