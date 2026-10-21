#!/usr/bin/env python3
"""Two exact coupled endpoint caps and a compact H2 face in the 839-wall branch.

Fraction projection identities prove all equalities. Outward integer intervals
prove strict signs. No numerical proposal or optimizer is imported.
"""
import argparse,json
from fractions import Fraction as F
from pathlib import Path
import fresh_h22_incidence_prism_20261013_verify as V
P=V.P;TAG='fresh_h22_pi7_coupled_caps_20261013';OUT='results/'+TAG+'.json'
OLD='results/fresh_h22_pi7_cusp_pair_20261013.json';GEO='results/'+V.PRE+'_geometry.json'
EDGE='results/fresh_h22_pi7_edgewalk_20261013.json';UNCAPPED='results/'+V.PRE+'_uncapped.json'
def enc(x):return list(map(str,x))
def inverse(M):
 n=len(M);A=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(M)]
 for i in range(n):
  k=next(j for j in range(i,n)if A[j][i]);A[i],A[k]=A[k],A[i];q=A[i][i];A[i]=[x/q for x in A[i]]
  for j in range(n):
   if j!=i:
    c=A[j][i];A[j]=[x-c*y for x,y in zip(A[j],A[i])]
 return [r[n:]for r in A]
def pd(M):
 A=[[F(x)for x in r]for r in M];piv=[]
 for i in range(len(A)):
  z=A[i][i];assert z>0;piv.append(z)
  for j in range(i+1,len(A)):
   for k in range(i+1,len(A)):A[j][k]-=A[j][i]*A[i][k]/z
 return piv
def div(a,b):
 assert b[0]>0
 return V.mul(a,(V.S*V.S//b[1],-((-V.S*V.S)//b[0])))
def ivec(v):return[V.rat(x)for x in v]
def iadd(a,b):return[V.add(x,y)for x,y in zip(a,b,strict=True)]
def iscale(a,c):return[V.scale(x,c)for x in a]
def imul(a,c):return[V.mul(x,c)for x in a]
def ipair(M,a,b):
 cov=[V.isum(V.scale(x,y)for x,y in zip(b,row,strict=True))for row in M]
 return V.isum(V.mul(x,y)for x,y in zip(a,cov,strict=True))
def run(verify=False):
 cert=V.read(V.OUT);old=V.read(OLD);geo=V.read(GEO);frontier=V.read(EDGE);uncapped=V.read(UNCAPPED);recipes=V.read(V.REC)['recipes']
 assert cert['status']=='PASS_EXACT_835_FACET_PI7_INCIDENCE_BRANCH' and old['status']=='PASS_EXACT_837_FACETS_AND_TWO_BOUNDED_CUSP_LINKS'
 for obj in(cert,old,geo,frontier,uncapped):
  for p,z in obj['sources'].items():assert V.sha(p)==z
 M,R,h=V.setup();caps=[]
 for rec in recipes:
  u=list(map(F,rec['u']));v=list(map(F,rec['v']));q=F(rec['norm']);qu=V.pair(M,u,u);qv=V.pair(M,v,v)
  assert qv>0 and V.pair(M,u,v)==0
  rad=(q/qv,F(0),-q*q*qu/(4*qv));root=V.scale(V.sqrt(V.ki(rad)),rec['radical_branch']);p=V.scale(V.T,-q/2)
  caps.append([V.add(V.scale(p,a),V.scale(root,b))for a,b in zip(u,v,strict=True)])
 walls=[ivec(r)for r in R]+caps
 for rec in old['new_caps']:
  k=rec['reflected_cap_index'];normal=caps[k][:-1]+[V.neg(caps[k][-1])]
  assert [enc(x)for x in normal]==rec['normal_intervals'];walls.append(normal)
 assert len(walls)==837
 vr=geo['new_finite_vertices'][15];S=vr['base_walls'];k=vr['new_cap_index'];assert k==17 and vr['candidate']==246
 assert S==[4,6,9,11,14,16,18,21,22,23,24,25,27,28,40,46,48,49,50,63,196]
 assert all(cert['new_caps'][k]['finite_contacts'][str(i)]==2 for i in S) and F(recipes[k]['norm'])==2
 A=[[V.pair(M,R[i],R[j])for j in S]for i in S];piv=pd(A);AI=inverse(A);assert all(x>=0 for r in AI for x in r)
 b=[V.pair(M,R[i],h)for i in S];lam=V.mv(AI,b);h0=[F(h[j])-sum(c*R[i][j]for c,i in zip(lam,S,strict=True))for j in range(23)]
 assert all(V.pair(M,h0,R[i])==0 for i in S);q0=V.pair(M,h0,h0);assert q0<0
 cap=caps[k];hp=V.isum(V.scale(x,y)for x,y in zip(cap,V.mv(M,h),strict=True));pp=iadd(ivec(h0),imul(cap,V.scale(hp,F(-1,2))))
 qp=V.add(V.rat(q0),V.neg(V.scale(V.mul(hp,hp),F(1,2))));assert qp[1]<0
 # p=h0-(h,cap)cap/2 is exactly orthogonal to all22 incident walls.
 # Hence p.h=p²=q0-(h,cap)²/2, and Z=sqrt(2/(-p²))*p has square -2.
 scale=V.sqrt(div(V.rat(2),V.neg(qp)));Z=imul(pp,scale)
 dual=[]
 for release in(27,9):
  row=AI[S.index(release)];u=[sum(c*R[i][j]for c,i in zip(row,S,strict=True))for j in range(23)]
  assert all(V.pair(M,u,R[i])==int(i==release)for i in S);dual.append(u)
 u,v=dual;assert(V.pair(M,u,u),V.pair(M,v,v),V.pair(M,u,v))==(2,1,1)
 # u,v lie in span(S), so each is exactly perpendicular to Z.
 D=F(29281);rootD=V.sqrt(V.rat(D));rA=iscale(iadd(iscale(Z,209),ivec([-241*x for x in u])),F(1,120))
 rB=imul(iadd(iscale(Z,241),ivec([-418*x for x in v])),div(V.rat(1),rootD))
 assert F(-2*209**2+2*241**2,120**2)==2 and F(-2*241**2+418**2,D)==2 and -2*209*241+241*418==0
 I=S+[812+k];edges=[[i for i in I if i!=j]for j in(27,9)];new=[rA,rB];new_records=[]
 assert frontier['status']=='PASS_EXACT_835_WALL_FINITE_VERTEX_FRONTIER'
 for ei,edge in zip((342,332),edges,strict=True):
  assert frontier['rows'][ei]['status']=='CERTIFIED_UNCLOSED_EDGE_RAY' and frontier['rows'][ei]['edge']==edge
 # The old null limits can be scaled to Z-u and Z-sqrt2*v. Both are cut.
 assert F(-418+482,120)==F(8,15)
 second_null_cut=div(V.add(V.rat(-482),V.scale(V.sqrt(V.rat(2)),418)),rootD);assert second_null_cut[0]>0
 norms=[V.pair(M,r,r)for r in R]+[F(rec['norm'])for rec in recipes]+[F(2),F(2)]
 for nr,edge in zip(new,edges,strict=True):
  hpnew=ipair(M,nr,ivec(h));assert hpnew[1]<0;pairings=[]
  for i,w in enumerate(walls):
   val=ipair(M,nr,w)
   if i in edge:assert val[0]<=0<=val[1];kind='right_by_projection_identity'
   else:V.dotted(val,2,norms[i]);kind='dotted'
   pairings.append({'wall':i,'kind':kind,'interval':enc(val)})
  new_records.append({'norm':2,'right_neighbors':edge,'center_pairing':enc(hpnew),'old_pairings':pairings,'normal_intervals':[enc(x)for x in nr]})
 mutual=ipair(M,rA,rB);assert mutual[0]<=0<=mutual[1] # equality proved above, not by overlap.
 walls+=new;allnorms=norms+[F(2),F(2)];common=[i for i in I if i not in(27,9)];assert len(common)==20
 # Compact quadrilateral vertices in the exact orthogonal frame (Z,u,v).
 triples=[(1,0,0),(241,-209,0),(50369,-29281,-28800),(209,0,-241)]
 active=[I,edges[0]+[837],common+[837,838],edges[1]+[838]];vertices=[]
 for co,inc in zip(triples,active,strict=True):
  a,b,c=co;point=iadd(iscale(Z,a),ivec([b*x+c*y for x,y in zip(u,v,strict=True)]));q=-2*a*a+2*b*b+c*c+2*b*c;assert q<0
  future=ipair(M,point,ivec(h));assert future[1]<0;ps=[]
  for i,w in enumerate(walls):
   val=ipair(M,point,w)
   if i in inc:assert val[0]<=0<=val[1]
   else:assert val[1]<0,(co,i,val)
   ps.append(enc(val))
  # All required zeros follow from the dual identities and the following scalar equations.
  if 837 in inc:assert -418*a-241*(2*b+c)==0
  if 838 in inc:assert -482*a-418*(b+c)==0
  assert len(inc)==22
  vertices.append({'frame_coefficients':co,'norm':q,'incident':inc,'future_center_pairing':enc(future),'all_wall_pairings':ps})
 # The released-normal Gram in this H2 slice is the inverse of its dual Gram.
 assert inverse([[F(2),F(1)],[F(1),F(1)]])==[[F(1),F(-1)],[F(-1),F(2)]]
 # Every potentially affine diagram containing a new cap consists of a
 # principal PD edge block and one/two orthogonal new norm-two roots.
 assert set(edges[0])&set(edges[1])==set(common)
 for edge in edges:
  base=[i for i in edge if i<812];pd([[V.pair(M,R[i],R[j])for j in base]for i in base])
 # The previous strict escape survives the two new normals as well.
 esc=geo['strict_escape'];x=[F(a,esc['denominator'])for a in esc['numerator']];qx=V.pair(M,x,x);hx=V.pair(M,x,h);assert qx>0 and hx<0
 eps=[]
 for w in walls:
  val=ipair(M,w,ivec(x));assert val[1]<0;eps.append(enc(val))
 # Both repaired old cusps remain strict for these additions.
 cusp=[]
 for row in old['bounded_cusp_links']:
  ray=list(map(F,row['primitive_integer_null_ray']));ps=[ipair(M,r,ivec(ray))for r in new];assert all(p[1]<0 for p in ps);cusp.append({'ray':enc(ray),'new_pairings':[enc(p)for p in ps]})
 dr=list(map(F,uncapped['old_null_ray']))+[F(0)];assert V.pair(M,dr,dr)==0 and V.pair(M,dr,h)<0
 dps=[ipair(M,r,ivec(dr))for r in new];assert all(p[1]<0 for p in dps)
 sx=list(map(F,uncapped['strict_escape']['old_coordinates']))+[F(uncapped['strict_escape']['last_coordinate_in_sqrt2_units'])];assert V.pair(M,sx,sx)==1586 and V.pair(M,sx,h)<0
 small_pairs=[ipair(M,w,ivec(sx))for w in walls];assert all(p[1]<0 for p in small_pairs)
 out={'status':'PASS_EXACT_839_FACETS_COUPLED_ENDPOINTS_COMPACT_H2_AND_STRICT_ESCAPE','sources':{p:V.sha(p)for p in(V.OUT,V.REC,GEO,OLD,EDGE,UNCAPPED,'scripts/'+V.PRE+'_verify.py','scripts/'+TAG+'.py')},'interval_denominator':str(V.S),'facets':839,'old_vertex':15,'old_pi7_cap':812+k,'old_incident':I,'released_walls':[27,9],'old21_positive_pivots':enc(piv),'old21_inverse':[enc(r)for r in AI],'rational_center_projection':enc(h0),'rational_projection_norm':str(q0),'full_center_projection_norm_interval':enc(qp),'rational_duals':[enc(r)for r in dual],'dual_gram':[[2,1],[1,1]],'normalized_timelike_formula':'Z=sqrt(2/(-p²))*p, p=h0-(h,cap829)*cap829/2; Z²=-2','new_formulas':['(209*Z-241*u)/120','(241*Z-418*v)/sqrt29281'],'new_caps':new_records,'new_mutual_pairing':'0','new_pair_counts':{'right':43,'dotted':1632},'old_open_edges_closed':[342,332],'null_limit_cut_pairings':['8/15',enc(second_null_cut)],'common20_walls':common,'compact_H2_vertices':vertices,'intrinsic_quadrilateral_angles':['pi/4','pi/2','pi/2','pi/2'],'new_positive_affine_diagrams':0,'preserved_bounded_cusps':cusp,'strict_escape':{'numerator':esc['numerator'],'denominator':esc['denominator'],'norm':str(qx),'center_pairing':str(hx),'all839_pairings':eps},'untouched_old_cusp':{'ray':enc(dr),'new_pairings':[enc(p)for p in dps]},'small_integer_strict_escape':{'vector':enc(sx),'norm':1586,'center_pairing':str(V.pair(M,sx,h)),'all839_pairings':[enc(p)for p in small_pairs]},'scope':'Two genuinely new all-right endpoint caps meet old pi7 cap829 at right angles and each other at a right angle. They retain all837 prior walls, close two old open edges, and give three new finite vertices forming an actual compact H2 quadrilateral with the old vertex. No new positive affine diagram is introduced. Both prior strict escapes survive and a specified old cusp remains untouched. Infinite volume is proved; no complete global vertex/frontier census is claimed.'}
 if verify:assert json.loads(json.dumps(out))==V.read(OUT)
 else:(P/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'new_pairs',1675,'escape',str(qx),flush=True)
 return out
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
