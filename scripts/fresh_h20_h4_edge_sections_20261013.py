#!/usr/bin/env python3
"""Nine exact actual-ridge failures of a bounded A4->H4 edge-section pilot."""
from dataclasses import dataclass
from fractions import Fraction as F
from pathlib import Path
import argparse,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_h4_edge_sections_20261013'
RAW='data/enumerations/borcherds_d4_h21_og10_carrier_obstruction.json';WALK='results/fresh_h21_source_corner_walk_20261013.json';SCOUT=f'results/{TAG}_scout.json'
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def sg(x):return (x>0)-(x<0)
@dataclass(frozen=True)
class K:
 a:F=F(0);b:F=F(0)
 def __post_init__(self):object.__setattr__(self,'a',F(self.a));object.__setattr__(self,'b',F(self.b))
 def __add__(self,o):o=kk(o);return K(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return K(-self.a,-self.b)
 def __sub__(self,o):return self+-kk(o)
 def __rsub__(self,o):return kk(o)+-self
 def __mul__(self,o):o=kk(o);return K(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):o=kk(o);d=o.a*o.a-5*o.b*o.b;assert d;return self*K(o.a/d,-o.b/d)
 def __rtruediv__(self,o):return kk(o)/self
 def __bool__(self):return bool(self.a or self.b)
 def sign(self):
  a,b=sg(self.a),sg(self.b)
  if not a:return b
  if not b or a==b:return a
  return a*sg(self.a*self.a-5*self.b*self.b)
 def enc(self):return [str(self.a),str(self.b)]
def kk(x):return x if isinstance(x,K)else K(x)
RAD=K(1)
@dataclass(frozen=True)
class L:
 a:K=K();b:K=K()
 def __post_init__(self):object.__setattr__(self,'a',kk(self.a));object.__setattr__(self,'b',kk(self.b))
 def __add__(self,o):o=ll(o);return L(self.a+o.a,self.b+o.b)
 __radd__=__add__
 def __neg__(self):return L(-self.a,-self.b)
 def __sub__(self,o):return self+-ll(o)
 def __rsub__(self,o):return ll(o)+-self
 def __mul__(self,o):o=ll(o);return L(self.a*o.a+RAD*self.b*o.b,self.a*o.b+self.b*o.a)
 __rmul__=__mul__
 def __truediv__(self,o):o=ll(o);d=o.a*o.a-RAD*o.b*o.b;assert d;return self*L(o.a/d,-o.b/d)
 def __rtruediv__(self,o):return ll(o)/self
 def __bool__(self):return bool(self.a or self.b)
 def sign(self):
  a,b=self.a.sign(),self.b.sign()
  if not a:return b
  if not b or a==b:return a
  return a*(self.a*self.a-RAD*self.b*self.b).sign()
 def enc(self):return [self.a.enc(),self.b.enc()]
def ll(x):return x if isinstance(x,L)else L(x)
def dot(a,b):return sum((x*y for x,y in zip(a,b,strict=True)),0)
def mv(A,x):return [dot(r,x)for r in A]
def pair(a,b,Q):return dot(a,mv(Q,b))
def inv(A):
 n=len(A);M=[[F(x)for x in r]+[F(i==j)for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  i=next(i for i in range(j,n)if M[i][j]);M[i],M[j]=M[j],M[i];q=M[j][j];M[j]=[x/q for x in M[j]]
  for i in range(n):
   if i!=j:q=M[i][j];M[i]=[x-q*y for x,y in zip(M[i],M[j],strict=True)]
 return [r[n:]for r in M]
def pd(A):
 A=[[kk(x)for x in r]for r in A];out=[]
 while A:
  q=A[0][0];assert q.sign()>0;out.append(q.enc())
  A=[[A[i][j]-A[i][0]*A[0][j]/q for j in range(1,len(A))]for i in range(1,len(A))]
 return out
def enc(v):return [x.enc()if isinstance(x,(K,L))else str(x)for x in v]
def run():
 global RAD
 r=json.loads((ROOT/RAW).read_text())['coxeter_input'];Q=[[-int(x)for x in rr]for rr in r['basis_gram']];W=[[int(x)for x in rr]for rr in r['wall_coordinates_as_columns_transposed']]
 C=[mv(Q,w)for w in W];G=[[dot(w,c)for c in C]for w in W]
 assert len(W)==210 and all(len(w)==22 for w in W)
 assert all(G[i][j]<=0 for i in range(210)for j in range(i))
 for i in range(210):
  for j in range(i):
   cs=F(G[i][j]**2,G[i][i]*G[j][j]);assert cs>=1 or cs in[F(0),F(1,4),F(1,2),F(3,4)]
 h=[sum(w[j]for w in W)for j in range(22)];gg=math.gcd(*h);h=[x//gg for x in h]
 assert pair(h,h,Q)==-14 and all(dot(h,c)<0 for c in C)
 walk=json.loads((ROOT/WALK).read_text());assert walk['status']=='PASS_COMPLETE_CORNER_ORBIT_CENSUS'and not walk['pending_orbit_ids'];assert walk['sources'][RAW]==sha(RAW)
 scout=json.loads((ROOT/SCOUT).read_text());assert all(sha(p)==z for p,z in scout['sources'].items())
 qe=K(F(1,4),F(1,4));gamma=qe-F(4,5);assert gamma==K(F(-11,20),F(1,4))and gamma.sign()>0
 rows=[]
 for si,seed in enumerate(scout['selected_edges']):
  caseid,recipe=next((i,x)for i,x in enumerate(scout['cases'])if x['seed']==si and x['parameter']=='1/2'and x['radical_sign']==1)
  c=walk['corner_orbits'][seed['corner']];edge=c['edge_orbits'][seed['edge_index']];E=edge['edge_walls'];cc=seed['component'];anchor=recipe['anchor'];ai=E.index(anchor)
  assert len(E)==20 and len(set(E))==20 and anchor in cc
  assert len(cc)==4 and all(G[i][i]==2 for i in cc)
  assert sorted(sum(G[i][j]!=0 for j in cc)-1 for i in cc)==[1,1,2,2]
  assert all(G[i][j]in[0,-1]for i in cc for j in cc if i!=j)
  assert sum(G[anchor][j]!=0 for j in cc)==2
  assert all(G[i][j]==0 for i in cc for j in E if j not in cc)
  v=list(map(F,c['primitive_ray']));b=list(map(F,edge['endpoint_primitive_ray']))
  for vec in[v,b]:
   assert pair(vec,vec,Q)<=0 and pair(vec,h,Q)<0 and max(dot(vec,co)for co in C)==0
  v=[x/-pair(v,h,Q)for x in v];b=[x/-pair(b,h,Q)for x in b]
  P=[(x+y)/2 for x,y in zip(v,b,strict=True)];qp=pair(P,P,Q);pp=[dot(P,co)for co in C]
  assert qp<0 and [i for i,x in enumerate(pp)if not x]==E and all(x<=0 for x in pp)
  A=[[G[i][j]for j in E]for i in E];oldpiv=pd(A);I=inv(A)
  w=[sum(I[ai][j]*W[i][d]for j,i in enumerate(E))for d in range(22)]
  assert pair(w,w,Q)==F(4,5)and all(dot(w,C[i])==int(i==anchor)for i in E)
  z=[x-pair(b,P,Q)/qp*y for x,y in zip(b,P,strict=True)];qz=pair(z,z,Q)
  assert qz>0 and pair(P,z,Q)==0 and all(dot(z,C[i])==0 for i in E)
  RAD=gamma/qz
  # Negative field norm excludes squares in K and validates all field divisions.
  assert RAD.sign()>0 and RAD.a*RAD.a-5*RAD.b*RAD.b<0
  e=[L(x,y)for x,y in zip(w,z,strict=True)];ep=[dot(e,co)for co in C]
  assert pair(e,e,Q)==L(qe)and pair(e,P,Q)==L()
  transverse=pair(e,[x-y for x,y in zip(b,v,strict=True)],Q);assert transverse
  assert all(ep[i]==L(int(i==anchor))for i in E)
  B=[[K(G[i][j])-(K(int(i==anchor))*int(j==anchor))/qe for j in E]for i in E]
  piv=pd(B)
  BI=[[K(I[i][j])+K(I[i][ai]*I[ai][j])/gamma for j in range(20)]for i in range(20)]
  assert all(sum(B[i][k]*BI[k][j]for k in range(20))==K(int(i==j))for i in range(20)for j in range(20))
  # Exact H4 replacement; the other local finite components are unchanged.
  neighbor=next(j for j in cc if G[anchor][j]<0)
  cos5=K(1)/(2*(K(2)-1/qe));assert cos5==K(F(3,8),F(1,8))
  witness=recipe['witness'];assert witness['kind']=='obtuse';release=witness['release'];k=E.index(release);j=witness['new_wall']
  D=[sum(L(BI[k][l])*W[i][d]for l,i in enumerate(E))-L(BI[k][ai]/qe)*e[d]for d in range(22)]
  assert D==[L(sum(I[k][l]*W[i][d]for l,i in enumerate(E)))-L(K(I[k][ai])/gamma)*L(0,z[d])for d in range(22)]
  dp=[dot(D,co)for co in C];qd=pair(D,D,Q)
  assert qd==L(BI[k][k])and qd.sign()>0 and pair(P,D,Q)==L()and pair(e,D,Q)==L()
  assert all(dp[i]==L(int(i==release))for i in E)
  assert dp[j].sign()<0 and pp[j]<0
  step=L(pp[j])/dp[j];assert step.sign()>0
  x=[L(a)-step*d for a,d in zip(P,D,strict=True)];xp=[L(a)-step*d for a,d in zip(pp,dp,strict=True)]
  assert [dot(x,co)for co in C]==xp and pair(x,e,Q)==L()
  qx=L(qp)+step*step*qd;assert pair(x,x,Q)==qx and qx.sign()<0 and pair(x,h,Q).sign()<0
  active=[i for i,y in enumerate(xp)if not y]
  assert active==witness['active']==sorted([i for i in E if i!=release]+[j])and all(y.sign()<=0 for y in xp)
  # The retained19 Gram is a PD principal minor of B. The new root has
  # positive Schur complement because qx/qp>0; hence all20 are independent.
  schur=dp[j]*dp[j]/qd+L(pp[j]*pp[j]/qp)
  assert schur==dp[j]*dp[j]/qd*qx/qp and schur.sign()>0
  a,bb=witness['bad_pair'];assert a in active and bb in active
  qa=L(G[a][a])-ep[a]*ep[a]/qe;qb=L(G[bb][bb])-ep[bb]*ep[bb]/qe
  p=L(G[a][bb])-ep[a]*ep[bb]/qe
  assert qa.sign()>0 and qb.sign()>0 and p.sign()>0 and (qa*qb-p*p).sign()>0
  rows.append({'scout_case':caseid,'source_corner':c['id'],'source_edge_index':seed['edge_index'],'endpoint_orbits':seed['endpoint_orbits'],'edge_walls':E,'A4_component':cc,'A4_endpoint':anchor,'position':'1/2','radical_sign':1,'radicand':RAD.enc(),'radicand_field_norm':str(RAD.a*RAD.a-5*RAD.b*RAD.b),'section_normal':enc(e),'section_norm':qe.enc(),'initial_point':enc(P),'initial_norm':str(qp),'edge_transverse_pair':transverse.enc(),'old_PD_pivots':oldpiv,'projected_PD_pivots':piv,'new_H4_squared_cosine':cos5.enc(),'release':release,'new_wall':j,'edge_step':step.enc(),'actual_vertex':enc(x),'actual_vertex_norm':qx.enc(),'actual_vertex_future_pair':pair(x,h,Q).enc(),'actual_vertex_incident_walls':active,'all_210_vertex_pairings':enc(xp),'new_root_schur_complement':schur.enc(),'actual_obtuse_ridge':[a,bb],'ridge_projected_Gram':[[qa.enc(),p.enc()],[p.enc(),qb.enc()]]})
  print('EXACT',si,'actual obtuse ridge',a,bb,flush=True)
 assert len(rows)==9 and len({tuple(r['endpoint_orbits'])for r in rows})==9
 return {'status':'PASS_NINE_EXACT_A4_TO_H4_EDGE_SECTIONS_WITH_ACTUAL_OBTUSE_RIDGES','sources':{p:sha(p)for p in[RAW,WALK,SCOUT,f'scripts/{TAG}.py']},'field_convention':'K=Q(sqrt5), a+b*sqrt5; each row L=K(sqrt(radicand)>0), coordinate encodes [a,b] with a,b in K','cases':rows,'scope':'Nine exact literal midpoint/+branch sections chosen from108 numerical proposals. Every section has an actual finite H4-containing corner and finite volume inherited from the source, but also an actual obtuse ridge at another finite vertex. No all-parameter or complete edge-orbit exclusion.'}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');aa=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if aa.verify:assert json.loads(p.read_text())==out
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'])
