#!/usr/bin/env python3
"""All-label direct A1xA1 coupling obstruction for the fixed H18 source.

The single-contact label tail is closed analytically by a forced source
anchor. No finite-label cutoff or assignment enumeration is assumed.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_h20_b2_h18_coupling_20261014_review as A

K=A.K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_a1a1_h18_coupling_20261015'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
PREVIOUS='results/fresh_h20_b2_h18_coupling_20261014.json'
HELPER='scripts/fresh_h20_b2_h18_coupling_20261014_review.py'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
 if isinstance(x,F):return str(x)
 if isinstance(x,K):return x.enc()
 if isinstance(x,dict):return{k:enc(v)for k,v in x.items()}
 if isinstance(x,(list,tuple)):return[enc(v)for v in x]
 return x
def sign5(a,b):
 a,b=F(a),F(b)
 if not b:return(a>0)-(a<0)
 if not a:return(b>0)-(b<0)
 if(a>0)==(b>0):return(a>0)-(a<0)
 return((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))

def run():
 d=json.loads((ROOT/SOURCE).read_text());oldcert=json.loads((ROOT/PREVIOUS).read_text())
 for p,h in oldcert['sources'].items():assert sha(p)==h,p
 G=[[A.dec(x)for x in row]for row in d['gram']];assert len(G)==40
 basis=oldcert['source_basis'];assert len(set(basis))==19
 Q=[[G[i][j]for j in basis]for i in basis];J=A.inverse(Q)
 old=[A.mv(J,[G[j][i]for j in basis])for i in range(40)]
 cov=[A.mv(Q,n)for n in old]
 assert all(A.dot(old[i],cov[j])==G[i][j]for i in range(40)for j in range(40))
 norms=[G[i][i]for i in range(40)];assert all(q in(K(1),K(2))for q in norms)
 edges=[(i,j)for i in range(40)for j in range(i)if G[i][j].sign()<0 and 4*G[i][j]*G[i][j]==norms[i]*norms[j]]
 forced=sorted(set(i for e in edges for i in e));assert len(edges)==27 and len(forced)==24
 assert all(norms[i]==2 for i in forced)
 seen={0};tree=[]
 while len(seen)<len(forced):
  e=next((e for e in edges if(e[0]in seen)!=(e[1]in seen)),None);assert e is not None
  u,v=e if e[0]in seen else e[::-1];tree.append([u,v]);seen.add(v)
 assert sorted(seen)==forced
 clique=[0,2,4]
 assert all(G[i][j]==0 for i in clique for j in clique if i!=j)
 # Special profiles have nonempty contact support in two colors. Three
 # pairwise orthogonal source walls would require three pairwise-disjoint
 # nonempty subsets. Check this finite set statement literally.
 subsets=[{0},{1},{0,1}]
 assert not any(not(a&b or a&c or b&c)for a in subsets for b in subsets for c in subsets)
 eligible=sorted(set(range(40))-set(forced));core=[17,23,24]
 assert eligible==[17]+list(range(23,36))+[38,39]
 anchors=[]
 for i in eligible:
  g=K(-1)if i in core else K(-2)
  j=next(j for j in forced if G[i][j]==g)
  assert norms[i]==1 if i in core else norms[i]in(K(1),K(2))
  anchors.append(dict(source_wall=i,ordinary_anchor=j,source_pairing=g,source_norm=norms[i]))
 # Positive projection initially allows a single arbitrary finite label,
 # or exactly (3,3),(3,4),(4,3),(3,5),(5,3). For double contacts m,n>=4
 # already cost>=2. For (3,n), positivity requires n<6.
 initial_doubles=[[3,3],[3,4],[4,3],[3,5],[5,3]]
 # Exact scalar endpoints for the all-label anchor proof. For m>=9,
 # c=cos(pi/9) is the root of f(c)=4c^3-3c-1/2 on (sqrt(3)/2,1).
 # At u=sqrt(7/8), f(u)=(u-1)/2<0 and f'>0 there. Thus h9<1/4.
 u2=F(7,8)
 assert 4*u2-3==F(1,2) and F(3,4)<u2<1 and 12*u2-3>0
 # h8=1-sqrt2/2 is strictly between1/4 and1/2. Hence m7,8 give
 # 1<4h<2, and m>=9 give0<4h<1 for the norm-one -2 anchor.
 h8=K(1,F(-1,2));assert(h8-F(1,4)).sign()>0 and(h8-F(1,2)).sign()<0
 # m5: h5=(5-sqrt5)/4. A norm-one -2 anchor has squared product
 # 5-sqrt5 strictly between phi^2 and3. A norm-two -2 anchor has
 # squared product(5-sqrt5)/2 strictly between1 and2.
 assert sign5(F(7,2),F(-3,2))>0 and sign5(-2,1)>0
 assert sign5(F(3,2),F(-1,2))>0 and sign5(F(-1,2),F(1,2))>0
 # Double(3,5) has h=(3-sqrt5)/4, so even the largest anchor square4h<1.
 assert sign5(3,-1)>0 and sign5(-2,1)>0
 # Every surviving profile has a contact sum and residual square in a
 # small exact alphabet. Contact0 denotes right; individual a_m=2cos(pi/m).
 profiles={
  'O':dict(labels=[0,0],h=F(2),contact_sum_square=K()),
  '30':dict(labels=[3,0],h=F(3,2),contact_sum_square=K(1)),
  '03':dict(labels=[0,3],h=F(3,2),contact_sum_square=K(1)),
  '40':dict(labels=[4,0],h=F(1),contact_sum_square=K(2)),
  '04':dict(labels=[0,4],h=F(1),contact_sum_square=K(2)),
  '60':dict(labels=[6,0],h=F(1,2),contact_sum_square=K(3)),
  '06':dict(labels=[0,6],h=F(1,2),contact_sum_square=K(3)),
  '33':dict(labels=[3,3],h=F(1),contact_sum_square=K(4)),
  '34':dict(labels=[3,4],h=F(1,2),contact_sum_square=K(3,2)),
  '43':dict(labels=[4,3],h=F(1,2),contact_sum_square=K(3,2)),
 }
 core_profiles=['O','40','04','33'];other_profiles=list(profiles)
 assignments=[dict(source_wall=i,profiles=['O']if i in forced else core_profiles if i in core else other_profiles)for i in range(40)]
 assert 2+sum(len(r['profiles'])for r in assignments)==168
 anchor_gates=[]
 for row in anchors:
  i=row['source_wall'];g=row['source_pairing'];q=norms[i]
  for name in assignments[i]['profiles']:
   h=profiles[name]['h'];g2=g*g*h/q
   if g2==1:relation='pi3'
   elif g2==2:relation='pi4'
   elif g2==3:relation='pi6'
   elif g2==4:relation='parallel'
   else:assert(g2-4).sign()>0;relation='dotted'
   anchor_gates.append(dict(source_wall=i,profile=name,squared_product=g2,relation=relation))
 x=[K(a,b)/23 for a,b in[(142,-56),(282,-84),(422,-112),(351,-84),
   (280,-56),(233,-42),(186,-28),(139,-14),(92,0),(45,14),(-2,28),
   (-49,42),(-96,56),(-119,56),(-142,56),(-165,56),(-188,56),(-211,56),(211,-56)]]
 assert x==list(map(A.dec,oldcert['source_point']))
 h=[sum((n[j]for n in old),K())for j in range(19)]
 xp=[A.dot(x,c)for c in cov];hp=[A.dot(h,c)for c in cov]
 qx=A.dot(x,A.mv(Q,x));qh=A.dot(h,A.mv(Q,h));xh=A.dot(x,A.mv(Q,h))
 assert qx==K(F(3857,529),F(-2800,529))and qx.sign()<0
 assert qh==K(-2952,-1280)and xh==K(-16,-40)
 assert all(v.sign()<=0 for v in xp)and all(v.sign()<0 for v in hp)
 rt=K(0,1)
 for i in eligible:assert(xp[i]+(K(1)if norms[i]==1 else rt)).sign()<=0
 s=F(2);eps=F(1,100000)
 margins={name:s*s*p['h']-p['contact_sum_square']/4 for name,p in profiles.items()}
 assert all(v.sign()>0 for v in margins.values())
 z=[s*a+eps*b for a,b in zip(x,h)];zp=[A.dot(z,c)for c in cov]
 assert all(v.sign()<0 for v in zp)
 for i in eligible:assert(zp[i]+s*(K(1)if norms[i]==1 else rt)).sign()<0
 E=z+[K(F(-1,2)),K(F(-1,2))];qE=A.dot(z,A.mv(Q,z))+1
 assert (s*s*qx+1).sign()>0 and qE.sign()>0
 C=[100*a for a in h]+[K(F(-1,2)),K(F(-1,2))]
 cp=[100*v for v in hp]
 for i in eligible:assert(cp[i]+s*(K(1)if norms[i]==1 else rt)).sign()<0
 qC=10000*qh+1;EC=100*A.dot(z,A.mv(Q,h))+1
 assert qC.sign()<0 and EC.sign()<0
 return enc(dict(status='PASS_ALL_DIRECT_A1A1_H18_COUPLINGS_HAVE_COMMON_STRICT_ESCAPE',
  sources={p:sha(p)for p in[f'scripts/{TAG}.py',HELPER,SOURCE,PREVIOUS]},
  arithmetic='Fraction Q(sqrt2), rational comparisons in Q(sqrt5), and explicit monotonic trigonometric tail proof; no numerical acceptance',
  source_basis=basis,source_form=Q,source_normals=old,source_gram=G,
  support_gram=[[2,0],[0,2]],support_inverse=[[F(1,2),0],[0,F(1,2)]],
  pi3_edges=edges,forcing_tree=tree,orthogonal_triple=clique,forced_ordinary=forced,
  eligible_special=eligible,core=core,ordinary_anchors=anchors,
  initial_positive_double_profiles=initial_doubles,
  all_label_tail=dict(cos_pi9_comparison_square=u2,cubic_remainder_coefficient=F(1,2),h8=h8,
    proof='f(c)=4c^3-3c-1/2 increases on c>1/2. f(sqrt(7/8))=(sqrt(7/8)-1)/2<0, so h9=2sin^2(pi/9)<1/4; h decreases with m.'),
  allowed_profiles=profiles,per_source_profiles=assignments,complete_anchor_relations=anchor_gates,
  relaxed_bank_count=168,
  pair_rule='sqrt(h_i*h_j/(q_i*q_j))*G_ij+(a_i0*a_j0+a_i1*a_j1)/2; a_im=2cos(pi/m), with a_i0=0 for a right contact.',
  source_point=x,source_point_norm=qx,source_point_pairings=xp,source_center=h,
  source_center_norm=qh,source_center_pairings=hp,scale=s,strict_perturbation=eps,
  squared_profile_margins=margins,strict_escape=E,escape_norm=qE,
  unperturbed_escape_norm=s*s*qx+1,escape_source_pairings=zp,support_escape_pairings=[-1,-1],
  strict_center=C,center_norm=qC,center_source_pairings=cp,escape_center_pairing=EC,
  exact_checks=dict(source_gram_entries=1600,forcing_component_vertices=24,source_gap_constraints=16,ordinary_anchor_checks=len(anchor_gates),bank_strict_signs=168),
  scope='All coherent direct A1xA1 couplings with two support walls and inducing walls over the fixed40 source facets, allowing all finite contact labels. Additional noninducing walls and changed source metrics remain outside. No finite label cutoff, assignment enumeration or whole chamber census is claimed.'))

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
 if a.verify:assert d==json.loads(out.read_text());print('PASS saved all-label A1A1 coupling replay')
 else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status']);print(d['escape_norm'])
