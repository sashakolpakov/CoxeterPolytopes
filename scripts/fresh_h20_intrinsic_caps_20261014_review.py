#!/usr/bin/env python3
"""Independent accepted-profile and intrinsic-geometry check.

Uses SymPy radical expressions and rational matrix operations. It checks every
saved accepted profile and transport, but does not repeat the producer's full
rejected-profile enumeration. Completeness is reviewed separately in the note.
"""
import argparse,json,hashlib,math
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import sympy as S
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_intrinsic_caps_20261014_review'
SOURCE='results/fresh_h20_intrinsic_caps_20261014_endpoint.json'
RAD=[1,2,3,6,5,10,15,30]
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def exact(a):return S.expand(sum(S.Rational(x)*S.sqrt(n)for x,n in zip(a,RAD,strict=True)))
def sign(a):
 a=S.expand(a)
 if a==0:return 0
 lo=hi=F(0);scale=10**90
 # Match the squarefree radicals independently, using decimal enclosures.
 residual=a
 for n in RAD:
  coeff=a.coeff(S.sqrt(n))if n!=1 else a.as_coeff_Add()[0]
  assert coeff.is_Rational
  residual-=coeff*S.sqrt(n);v=F(coeff)
  b=math.isqrt(n*scale*scale);l=F(b,scale);h=l if b*b==n*scale*scale else F(b+1,scale)
  lo+=v*(l if v>=0 else h);hi+=v*(h if v>=0 else l)
 assert S.expand(residual)==0
 assert lo>0 or hi<0
 return 1 if lo>0 else -1

def run():
 p=json.loads((ROOT/SOURCE).read_text())
 for f,h in p['sources'].items():assert sha(f)==h
 assert p['seed_type']=='affine A11 x affine A2 x affine E6'
 ids=p['seed_order'];G=S.Matrix(p['normalized_seed_Gram']);marks=list(map(S.Rational,p['positive_affine_marks']));idx={w:i for i,w in enumerate(ids)}
 assert G.shape==(22,22)and G.rank()==19 and len(set(ids))==22
 assert all(G[i,i]==2 for i in range(22)) and all(G[i,j]in(-1,0)for i in range(22)for j in range(i))
 seen=set();components=[]
 for i in range(22):
  if i in seen:continue
  todo=[i];seen.add(i)
  for a in todo:
   for b in range(22):
    if b not in seen and G[a,b]:seen.add(b);todo.append(b)
  components.append(sorted(todo))
 assert sorted(map(len,components))==[3,7,12]
 for comp in components:
  M=G.extract(comp,comp);m=S.Matrix([marks[i]for i in comp]);assert M*m==S.zeros(len(comp),1)and all(x>0 for x in m)
 counts=Counter();summaries=[];mag={3:S.Integer(1),4:S.sqrt(2),5:(1+S.sqrt(5))/2,6:S.sqrt(3)}
 for typ in p['types']:
  e=typ['edge_support'];ii=list(map(idx.get,e));B=G.extract(ii,ii);I=B.inv()
  assert B==S.Matrix(typ['PD19_Gram'])and I==S.Matrix([[S.Rational(x)for x in r]for r in typ['inverse_Gram']])
  assert all(B[:j,:j].det()>0 for j in range(1,20))and all(x>=0 for x in I)
  assert min(I[i,i]for i in range(19))>=S.Rational(1,2)
  local=Counter();five=[]
  for r in typ['profiles']:
   contact=dict(r['contacts']);assert len(contact)==len(r['contacts'])and set(contact)<=set(e)
   c=S.Matrix([mag[contact[w]]if w in contact else 0 for w in e])
   cost=S.expand((c.T*I*c)[0]);residual=S.expand(2-cost)
   assert cost==exact(r['norm_two_projection_square'])and residual==exact(r['norm_two_schur_residual'])
   assert cost/2==exact(r['normalized_spatial_cost'])
   assert [S.expand(x/2)for x in c]==list(map(exact,r['cosine_vector']))
   sg=sign(residual);assert sg>=0 and r['endpoint']==('finite'if sg>0 else'ideal_PSD')
   local[r['endpoint']]+=1;counts['accepted_profiles']+=1
   for ratio in r['omitted_contact_ratios']:
    omit=ratio['omitted'];comp=next(v for v in components if idx[omit]in v)
    a=sum(marks[j]for j in comp)/marks[idx[omit]]
    b=S.expand(sum(marks[idx[w]]*c[j]/2 for j,w in enumerate(e)if idx[w]in comp)/marks[idx[omit]])
    assert S.Rational(ratio['numerator_over_radius'])==a and exact(ratio['subtract'])==b
    assert ratio['numerator_has_triangle_inradius_s']==(len(comp)==3)
   if any(m==5 for m in contact.values()):
    # Find the connected endpoint block containing the new node.
    conn={19};todo=[19]
    for v in todo:
     near=[j for j,w in enumerate(e)if c[j]]if v==19 else [j for j in range(19)if j!=v and B[v,j]]+([19]if c[v]else[])
     for j in near:
      if j not in conn:conn.add(j);todo.append(j)
    name={2:'I2(5)',3:'H3',4:'H4'}[len(conn)]
    if name=='H3':assert cost==(3+S.sqrt(5))/3
    if name=='H4':assert cost==(17+3*S.sqrt(5))/12
    five.append({'profile':r['id'],'contacts':r['contacts'],'block':name,'residual':str(residual)})
  assert dict(local)==typ['counts']
  tails=typ['dihedral_tails']
  for tail in tails:
   i=e.index(tail['wall']);assert B[i,i]==2 and all(B[i,j]==0 for j in range(19)if j!=i)
   assert I[i,i]==S.Rational(1,2)and len(tails)==1
  summaries.append({'type':typ['E6_deletion_type'],'multiplicity':typ['literal_edge_multiplicity'],'counts':dict(local),'label5_blocks':five,'isolated_I2_tails':len(tails)})
 transports=p['literal_edge_transports'];assert len(transports)==252 and len({tuple(x['omitted'])for x in transports})==252
 for tr in transports:
  perm=tr['seed_permutation'];assert set(perm)==set(ids)
  pp=[idx[x]for x in perm];assert G.extract(pp,pp)==G and [marks[i]for i in pp]==marks
  base=p['types'][tr['type']]
  assert [perm[idx[x]]for x in base['omitted']]==tr['omitted']
  assert [perm[idx[x]]for x in base['edge_support']]==tr['edge_support_in_profile_order']
 assert Counter(t['type']for t in transports)=={0:108,1:108,2:36}
 # Formal Lorentz Gram identities for an arbitrary finite/ideal profile.
 h,rho,cost=S.symbols('h rho cost',positive=True)
 # Coordinates in (p,V,D), with p²=cost, V²=D²=pV=pD=0, VD=-1.
 H=S.Matrix([[cost,0,0],[0,0,-1],[0,-1,0]])
 r=S.Matrix([1,h,-rho/(2*h)]);x=S.Matrix([0,1,rho/(2*h*h)])
 assert S.expand((r.T*H*r)[0])==cost+rho
 assert S.expand((r.T*H*x)[0])==0 and S.expand((x.T*H*x)[0])==-rho/(h*h)
 assert S.expand((r.T*H*S.Matrix([0,0,1]))[0])==-h
 a,b,d=S.symbols('R_i R_j distance_square',positive=True)
 pair=(a*a+b*b-d)/(a*b)
 assert S.expand(pair.subs(d,a*a+b*b))==0
 assert S.expand(pair.subs(d,(a+b)**2))==-2
 paths=[SOURCE,'scripts/fresh_h20_intrinsic_caps_20261014_endpoint.py',str(Path(__file__).resolve().relative_to(ROOT))]
 return {'status':'PASS_INTRINSIC_ENDPOINT_ACCEPTED_PROFILE_AND_GEOMETRY_REVIEW','sources':{f:sha(f)for f in paths},'counts':dict(counts),'edge_types':summaries,'transports':252,'finite_profiles':sum(r['counts']['finite']for r in summaries),'ideal_PSD_profiles':sum(r['counts']['ideal_PSD']for r in summaries),'formal_geometry':'r=p+hV-rho2*D/(2h), p²+rho2=2; endpoint x=V+rho2*D/(2h²), x²=-rho2/h²; cap-ball pair=(Ri²+Rj²-distance²)/(Ri Rj).','scope':'Independent arithmetic replay of all260 accepted profiles, all omitted-wall scalar formulas and252 transports; source/proof review of complete support/label enumeration. Does not independently re-enumerate rejected profiles or certify a globally finite-volume cap selection.'}

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run();p=ROOT/'results'/f'{TAG}.json'
 if args.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status']);print(out['finite_profiles'],out['ideal_PSD_profiles'])
