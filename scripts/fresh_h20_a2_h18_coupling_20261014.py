#!/usr/bin/env python3
"""All-label direct A2 coupling obstruction for the fixed 40-wall H18 source.

Only 40 source pairings and three scalar margins are needed to certify the
common strict escape for the entire 138-normal relaxed profile bank.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_b2_h18_coupling_20261014_review as A

K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_a2_h18_coupling_20261014'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
PREVIOUS='results/fresh_h20_b2_h18_coupling_20261014.json'
HELPER='scripts/fresh_h20_b2_h18_coupling_20261014_review.py'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def p5(a=0,b=0):return(F(a),F(b))
def plus5(x,y):return(x[0]+y[0],x[1]+y[1])
def scale5(x,s):return(x[0]*s,x[1]*s)
def mul5(x,y):return(x[0]*y[0]+5*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def sign5(x):
 a,b=x
 if not b:return(a>0)-(a<0)
 if not a:return(b>0)-(b<0)
 if(a>0)==(b>0):return(a>0)-(a<0)
 return((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))
def enc5(x):return dict(a=str(x[0]),b=str(x[1]),basis='1,sqrt5')
def encode(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,K):return x.enc()
 if isinstance(x,dict):return{k:encode(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[encode(v)for v in x]
 return x

def run():
 d=json.loads((ROOT/SOURCE).read_text());oldcert=json.loads((ROOT/PREVIOUS).read_text())
 for f,h in oldcert['sources'].items():assert sha(f)==h,f
 G=[[A.dec(x)for x in row]for row in d['gram']];assert len(G)==40
 basis=oldcert['source_basis'];assert len(set(basis))==19
 Q=[[G[i][j]for j in basis]for i in basis];J=A.inverse(Q)
 old=[A.mv(J,[G[j][i]for j in basis])for i in range(40)]
 cov=[A.mv(Q,n)for n in old]
 assert all(A.dot(old[i],cov[j])==G[i][j]for i in range(40)for j in range(40))
 norms=[G[i][i]for i in range(40)];assert all(q in(K(1),K(2))for q in norms)
 # Complete positive-projection contact list on a norm2 A2 support.
 support_gram=[[K(2),K(-1)],[K(-1),K(2)]]
 support_inverse=[[K(F(2,3)),K(F(1,3))],[K(F(1,3)),K(F(2,3))]]
 assert A.inverse(support_gram)==support_inverse
 assert sum((z for row in support_inverse for z in row),K())==2
 contact_squares={3:p5(1),4:p5(2),5:p5(F(3,2),F(1,2))}
 residual_squares={m:plus5(p5(2),scale5(a,F(-2,3)))for m,a in contact_squares.items()}
 assert all(sign5(h)>0 for h in residual_squares.values())
 assert F(2,3)*3==2 # m=6 is null, every larger finite label is nonpositive.
 assert 2*(F(2,3)+F(1,3))==2 # two minimal pi3 contacts are already null.
 # Complete source-pi3 transition table. Ordinary/special has product
 # -sqrt(h/2) in (-1,0). For two specials compare positive terms by squares.
 table=[]
 for ma in(0,3,4,5):
  for mb in(0,3,4,5):
   colors=['irrelevant']if not ma or not mb else['same','opposite']
   for color in colors:
    if ma==mb==0:kind='pi3'
    elif not ma or not mb:kind='forbidden_between_minus1_and0'
    else:
     a,b=contact_squares[ma],contact_squares[mb]
     diff=plus5(plus5(a,b),p5(-3))
     if color=='same':diff=plus5(diff,mul5(a,b))
     s=sign5(diff)
     if s==0:kind='right'
     elif s>0:kind='forbidden_positive'
     else:
      assert ma==mb==3 and color=='opposite'
      kind='forbidden_minus_one_third'
    table.append(dict(first=ma,second=mb,color=color,relation=kind))
 allowed=[(r['first'],r['second'],r['color'])for r in table if r['relation']in('right','pi3')]
 assert allowed==[(0,0,'irrelevant'),(3,3,'same'),(3,4,'opposite'),(4,3,'opposite')]
 # Ordinary/nonordinary status propagates across every source-pi3 edge.
 edges=[(i,j)for i in range(40)for j in range(i)if G[i][j].sign()<0 and 4*G[i][j]*G[i][j]==norms[i]*norms[j]]
 forced=sorted(set(i for e in edges for i in e));assert len(edges)==27 and len(forced)==24
 assert G[0][1]==G[1][2]==-1 and G[0][2]==0
 seen={0};tree=[]
 while len(seen)<len(forced):
  e=next((e for e in edges if(e[0]in seen)!=(e[1]in seen)),None);assert e is not None
  u,v=e if e[0]in seen else e[::-1];tree.append([u,v]);seen.add(v)
 assert sorted(seen)==forced
 # If one node were special, all24 would be special. In particular0,2
 # would then have strictly positive pairing, whatever their colors.
 eligible=sorted(set(range(40))-set(forced))
 assert eligible==[17]+list(range(23,36))+[38,39]
 # Reuse the literal exact source point, not a numerical optimizer result.
 x=[K(a,b)/23 for a,b in[(142,-56),(282,-84),(422,-112),(351,-84),
   (280,-56),(233,-42),(186,-28),(139,-14),(92,0),(45,14),(-2,28),
   (-49,42),(-96,56),(-119,56),(-142,56),(-165,56),(-188,56),(-211,56),(211,-56)]]
 assert x==list(map(A.dec,oldcert['source_point']))
 h=[sum((r[j]for r in old),K())for j in range(19)]
 xp=[A.dot(x,c)for c in cov];hp=[A.dot(h,c)for c in cov]
 qx=A.dot(x,A.mv(Q,x));qh=A.dot(h,A.mv(Q,h));xh=A.dot(x,A.mv(Q,h))
 assert qx==K(F(3857,529),F(-2800,529))and qx.sign()<0
 assert qh==K(-2952,-1280)and xh==K(-16,-40)
 assert all(v.sign()<=0 for v in xp)and all(v.sign()<0 for v in hp)
 rt=K(0,1)
 for i in eligible:assert(xp[i]+(K(1)if norms[i]==1 else rt)).sign()<=0
 s=F(3207,1000);epsilon=F(1,10**7)
 margins={m:plus5(scale5(residual_squares[m],s*s),scale5(contact_squares[m],-1))for m in(3,4,5)}
 assert all(sign5(z)>0 for z in margins.values())
 q_unperturbed=s*s*qx+2;assert q_unperturbed.sign()>0
 z=[s*a+epsilon*b for a,b in zip(x,h)]
 zp=[A.dot(z,c)for c in cov]
 assert all(v.sign()<0 for v in zp)
 for i in eligible:assert(zp[i]+s*(K(1)if norms[i]==1 else rt)).sign()<0
 E=z+[K(-1),K(-1)];qE=A.dot(z,A.mv(Q,z))+2
 assert qE.sign()>0
 center_source=[100*a for a in h];C=center_source+[K(-1),K(-1)]
 cp=[A.dot(center_source,c)for c in cov]
 for i in eligible:assert(cp[i]+s*(K(1)if norms[i]==1 else rt)).sign()<0
 assert all(v.sign()<0 for v in cp)
 qC=10000*qh+2;EC=100*A.dot(z,A.mv(Q,h))+2
 assert qC.sign()<0 and EC.sign()<0
 # The two A2 support products of E and C are both -1. For a special
 # profile, sqrt(h/q)*(z,n)+a < -s*sqrt(h)+a < 0 by the squared margins.
 profiles=[dict(source_wall=i,labels=[0]if i in forced else[0,3,4,5],
    colors=[]if i in forced else[0,1])for i in range(40)]
 assert 2+len(forced)+7*len(eligible)==138
 paths=[f'scripts/{TAG}.py',HELPER,SOURCE,PREVIOUS]
 return encode(dict(status='PASS_ALL_DIRECT_A2_H18_COUPLINGS_HAVE_COMMON_STRICT_ESCAPE',
  sources={p:sha(p)for p in paths},source_dimension=18,target_dimension=20,
  source_basis=basis,source_gram=G,source_form=Q,source_normals=old,
  support_gram=support_gram,support_inverse=support_inverse,
  complete_contact_profiles=dict(ordinary=dict(contact=0,residual_square=2),
    single_contact=[dict(label=m,contact_square=enc5(contact_squares[m]),residual_square=enc5(residual_squares[m]),colors=[0,1])for m in(3,4,5)]),
  pair_rule='sqrt(h_i*h_j/(q_i*q_j))*G_ij + a_i*a_j*(2/3 if same_color else 1/3); the second term is zero for an ordinary profile.',
  pi3_transition_table=table,pi3_edges=edges,forcing_tree=tree,orthogonal_witness=[0,2],forced_ordinary=forced,eligible_special=eligible,
  profiles=profiles,relaxed_bank_count=138,source_point=x,source_point_norm=qx,
  source_point_pairings=xp,source_center=h,source_center_norm=qh,source_center_pairings=hp,
  rational_scale=s,strict_perturbation=epsilon,squared_profile_margins={str(m):enc5(v)for m,v in margins.items()},
  unperturbed_escape_norm=q_unperturbed,strict_escape=E,escape_norm=qE,
  escape_source_pairings=zp,support_escape_pairings=[-1,-1],strict_center=C,
  center_norm=qC,center_source_pairings=cp,escape_center_pairing=EC,
  exact_checks=dict(source_gram_entries=1600,pi3_component_vertices=24,pi3_edges=27,source_point_gaps=16,positive_profile_margins=3,bank_strict_signs=138),
  theorem='Every coherent direct A2 coupling of the fixed40-facet H18 source has infinite volume, uniformly over all finite inducing-contact labels.',
  scope='Only the two A2 support walls and walls inducing the fixed40 source facets. Additional noninducing walls, other source metrics and generalH20 constructions are outside this gate. No assignment enumeration, full chamber census or finite-volume example claimed.'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved A2 all-label coupling replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status']);print(d['escape_norm'])
