#!/usr/bin/env python3
"""Independent new-domain, elementary-sign and nine-polynomial review."""
import argparse,ast,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
from math import comb
from pathlib import Path
import fresh_h9_star_single_third_opposite_20261015_review as A
R,P=A.R,A.P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_star_third_four_20261015';PRE='fresh_h9_star_single_third_opposite_20261015'
S=A.S;CORE=A.CORE;Z,O=R.Z,R.O;k=A.k;pa=A.pa;ps=A.ps;pm=A.pm;pp=A.pp

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(v):return list(map(str,v))
def dec(v):return tuple(map(F,v))
def run():
 seen={}
 def read(p):
  if p in seen:return seen[p]
  d=json.loads((ROOT/p).read_text());seen[p]=d
  for f,h in d.get('sources',{}).items():
   assert sha(f)==h,f
   if f.endswith('.json'):read(f)
  return d
 d=read(f'results/{TAG}.json');prior=read(f'results/{PRE}.json');pr=read(f'results/{PRE}_review.json')
 assert d['status']=='PASS_SINGLE_PI4_THIRD_SHARED_OR_OPPOSITE_ANCHORS'
 assert pr['status']=='PASS_INDEPENDENT_OPPOSITE_SINGLE_THIRD_DOMAINS_AND_POLYNOMIALS'and pr['sources'][f'results/{PRE}.json']==sha(f'results/{PRE}.json')
 geom=read('results/fresh_h9_tree_cube_wedges_20261014.json')['representatives'][0]['geometry'];V=[set(v['incident'])for v in geom['vertices']];MN=set(map(frozenset,geom['minimal_nonfaces']))
 assert len(V)==378 and len(MN)==23 and d['rank_core']==CORE
 assert any(set(S)<=v for v in V)and all(any(set(p)<=v for v in V)for p in combinations(CORE,2))
 faces={frozenset(v&set(CORE))for v in V};maximal={v for v in faces if not any(v<w for w in faces)}
 assert len(maximal)==9 and maximal==set(map(frozenset,d['complete_maximal_core_faces']))
 guardcounts=Counter();domain_summaries=[];canonical_by_mode={}
 for mode in('opposite','shared'):
  own={11:13,17:12 if mode=='opposite'else 13}
  base={(3,4):5,(6,7):7,(9,10):7,(12,13):7,(4,14):4}
  for a,r in own.items():base[tuple(sorted((a,r)))]=3
  saved=[r for r in d['positive_kernel_guards']if r['anchor_mode']==mode]
  assert len(saved)==len(prior['positive_kernel_guards'])==54
  for old,row in zip(prior['positive_kernel_guards'],saved):
   nodes=list(old['ordered_support'])
   if mode=='shared'and old['kind']in('mutual_and_third_contact_forbidden','mutual_and_owner_foreign_forbidden','nonright_owner_foreign_forces_other_third_right'):nodes=[13 if x==12 else x for x in nodes]
   assert nodes==row['ordered_support']and row['kind']==old['kind']and row['edges']==old['edges']and row['conditional']==old['conditional']
   assert len(nodes)==len(set(nodes))and set(nodes)==set(row['support'])<=V[row['actual_vertex']]
   assert set(nodes)<=set(CORE)
   facts=dict(base)
   for name,value in row['conditional'].items():
    if name=='trial_owner_internal_at_least':facts[3,4]=value
    elif name=='trial_mutual_label_at_least':facts[11,17]=value
    elif name in('foreign_contacts','third_contacts'):
     for pair in value:facts[tuple(sorted(pair))]=3
    else:
     assert name in('trial_nonright','trial_label_at_least_four','owner_foreign','third_contact','trial_mutual_nonright','foreign_contact','trial_third_nonright')
     facts[tuple(sorted(value))]=4 if name=='trial_label_at_least_four'else 3
   assert sorted([list(p)+[m]for p,m in facts.items()])==row['available_edge_lower_bounds']
   assert all(facts.get(tuple(sorted((nodes[i],nodes[j]))),0)>=m for i,j,m in row['edges'])
   rad=3 if any(m==6 for _,_,m in row['edges'])else 2
   def qa(x,y):return x[0]+y[0],x[1]+y[1]
   def qm(x,y):return x[0]*y[0]+rad*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
   def qs(x,c):return x[0]*c,x[1]*c
   def parse(text):
    def f(t):
     if isinstance(t,ast.Constant):return F(t.value),F(0)
     if isinstance(t,ast.UnaryOp):assert isinstance(t.op,ast.USub);return qs(f(t.operand),-1)
     if isinstance(t,ast.Call):assert t.func.id=='sqrt'and t.args[0].value==rad;return F(0),F(1)
     if isinstance(t.op,ast.Add):return qa(f(t.left),f(t.right))
     if isinstance(t.op,ast.Sub):return qa(f(t.left),qs(f(t.right),-1))
     if isinstance(t.op,ast.Mult):return qm(f(t.left),f(t.right))
     assert isinstance(t.op,ast.Div)and isinstance(t.right,ast.Constant);return qs(f(t.left),F(1,t.right.value))
    return f(ast.parse(text,mode='eval').body)
   cos={3:(F(1,2),F(0)),4:(F(0),F(1,2)),6:(F(0),F(1,2))}
   G=[[(F(1),F(0))if i==j else(F(0),F(0))for j in nodes]for i in nodes]
   for i,j,m in row['edges']:G[i][j]=G[j][i]=qs(cos[m],-1)
   w=list(map(parse,row['positive_kernel']));assert all(R.sign_quad(*z,rad)>0 for z in w)
   for rr in G:
    total=(F(0),F(0))
    for x,y in zip(rr,w):total=qa(total,qm(x,y))
    assert total==(0,0)
   guardcounts[mode+':'+row['kind']]+=1
  budgets=[r for r in d['projection_budgets']if r['anchor_mode']==mode];assert len(budgets)==3
  for row in budgets:
   a=row['leaf'];drop,keep=(3,4)if a==14 else(25-own[a],own[a])
   assert(row['removed_own_root'],row['retained_own_root'])==(drop,keep)
   assert set(row['support'])==(set(S)-{drop})|{a}and len(row['support'])==8 and set(row['support'])<=V[row['actual_vertex']]
   assert F(row['own_cost'])==(F(1,2)if a==14 else F(1,4))
   if a==14:assert F(row['foreign_pi3_cost_lower'])==F(1,2)and F(row['sum_lower'])==1
  assert F(1,2)+2*F(1,4)==1
  literal=[]
  for a,b in product((None,3,4),repeat=2):
   for dx,dy,mu in product((0,1),(0,1),(2,3,4,5)):
    if a is not None and b is not None:continue
    if dx and dy:continue
    if mu!=2 and(a is not None or b is not None or dx or dy):continue
    if(a==4 or b==4)and(dx or dy):continue
    literal.append([a,b,dx+2,dy+2,mu,None if a is None and b is None else 5])
  assert len(literal)==14
  rows=next(x for x in d['contact_domains']if x['anchor_mode']==mode)
  assert sorted(literal,key=str)==sorted(rows['literal_patterns'],key=str)
  canon=set()
  for a,b,dx,dy,mu,m in literal:
   if b is not None:a,b,dx,dy=b,a,dy,dx
   canon.add((None if a is None else a-3,dx-2,dy-2,mu,m))
  assert len(canon)==10 and canon==set(map(tuple,rows['canonical_patterns']))
  assert sum(p[0]is None for p in canon)==6 and all(not(p[1]and p[2])for p in canon)
  canonical_by_mode[mode]=canon;domain_summaries.append(dict(mode=mode,literal=14,canonical=10,no_foreign=6,foreign=4))
 # Independent permutation expansion of the elementary identity, over six variables.
 one6={(0,)*6:O};variables=[]
 for i in range(6):e=[0]*6;e[i]=1;variables.append({tuple(e):O})
 D,mu,h,dx,dy,v=variables
 aa=pa(one6,ps(D,k(F(-1,4))));off=pa(ps(mu,k(-1)),ps(pm(h,D),k(F(-1,4))))
 G=[[aa,off,ps(dx,k(-1))],[off,aa,ps(dy,k(-1))],[ps(dx,k(-1)),ps(dy,k(-1)),v]]
 identity=pa(pa(pm(v,pa(pm(aa,aa),ps(pm(off,off),k(-1)))),ps(pm(aa,pa(pm(dx,dx),pm(dy,dy))),k(-1))),ps(pm(pm(off,dx),dy),k(2)))
 assert P.detpoly(G)==identity
 # Both modes' anchor minors are negative at every permitted t,mu.
 assert F(100,19)>4 and 1-F(100,19)/2<0 and 1-F(7,16)*F(100,19)<0
 c5=R.COS[5];beta=R.COS[4];q=R.sub(O,R.mul(c5,c5))
 assert R.sign(q)>0 and R.sign(R.sub(q,k(F(1,2))))<0
 # The strict third diagonal v=1-1/(2q)<0 persists for every owner>=5.
 assert R.sign(R.sub(O,R.mul(k(F(1,2)),R.inverse(q))))<0
 one={(0,0,0):O};t={(0,0,1):O};den=pa(one,ps(pp(t,2),k(-1)))
 # Check the anchor-minor identities after multiplying by (1-t²)^2.
 for mode in('shared','opposite'):
  aa=pa(den,ps(one,k(F(-1,4))));off=ps(one if mode=='shared'else t,k(F(-1,4)))
  # mu dependence reconstructed as a second independent variable.
  m={(1,0,0):O};off=pa(off,ps(pm(m,den),k(-1)))
  minor=pa(pm(aa,aa),ps(pm(off,off),k(-1)))
  expected=pa(pa(pm(den,den),ps(den,k(F(-1,2)if mode=='shared'else F(-7,16)))),ps(pm(pm(m,one if mode=='shared'else t),den),k(F(-1,2))))
  expected=pa(expected,ps(pm(pp(m,2),pp(den,2)),k(-1)))
  assert minor==expected
 # The already independently proved cubic gives the full discrete tail split;
 # verify its rational bracket again rather than relying on floating cosines.
 C=[F(-1),F(-4),F(4),F(8)];lower=F(31,50);upper=2*F(901,1000)**2-1
 assert sum(a*lower**i for i,a in enumerate(C))<0<sum(a*upper**i for i,a in enumerate(C))
 assert 24*F(1,2)**2+8*F(1,2)-4>0 and 48*F(1,2)+8>0
 assert (1+lower)/2==F(81,100)and 1/(1-F(81,100))==F(100,19)
 assert R.sign(R.sub(R.scale(R.add(k(2),(F(0),F(1),F(0),F(0))),F(1,4)),k(F(23,25)**2)))>0
 checks=[];keys=set();counts=Counter();foreign_patterns={(0,0,0),(0,0,1),(0,1,0),(1,0,0)}
 for row in d['foreign_cases']:
  mode=row['anchor_mode'];foreign=row['foreign_owner_root'];x,y=map(F,row['third_anchor_cosines']);branch=row['branch']
  key=(mode,foreign,int(2*x),int(2*y),branch);assert key not in keys;keys.add(key)
  assert(foreign,int(2*x),int(2*y))in foreign_patterns and row['owner_internal']==5
  ff=R.mul(k(F(1,4)),R.inverse(q));T=R.mul(R.mul(beta,c5 if foreign==0 else O),R.mul(k(F(1,2)),R.inverse(q)))
  aa=pa(den,ps(one,k(F(-1,4))));ab=ps(one if mode=='shared'else t,k(F(-1,4)))
  at=ps(den,R.neg(R.add(k(x),T)));bt=ps(den,k(-y));cc=ps(den,R.sub(O,R.mul(R.mul(beta,beta),R.inverse(q))))
  M=[[pa(aa,ps(den,R.neg(ff))),ab,at],[ab,aa,bt],[at,bt,cc]]
  poly=ps(A.divide_square(P.detpoly(M),2),q)
  assert poly==A.parse(row['scaled_determinant'],{'t':t})
  exception=mode=='opposite'and(foreign,int(2*x),int(2*y))==(0,0,1)
  assert branch in({'last7','last_ge8'}if exception else{'all_last_ge7'})
  lo,hi,sign=F(9,10),F(1),1
  if branch=='last7':hi,sign=F(901,1000),-1
  if branch=='last_ge8':lo=F(23,25)
  assert list(map(F,row['interval']))==[lo,hi]and row['target_sign']==sign
  if exception:
   old=next(z for z in prior['continuous_cases']if z['branch']=='third5_last7')
   assert poly==A.compose(A.parse(old['scaled_determinant'],{'beta':{(1,0,0):O},'c':{(0,1,0):O},'t':t}),[{(0,0,0):c5},{(0,0,0):beta},t])
  result=A.compose(poly,[one,one,pa({(0,0,0):k(lo)},ps(t,k(hi-lo)))])
  assert all(a==b==0 for a,b,j in result)
  values=[R.total([R.scale(result.get((0,0,j),Z),F(comb(i,j),comb(2,j)))for j in range(i+1)])for i in range(3)]
  assert values==list(map(dec,row['bernstein_coefficients']))and all(R.sign(a)==sign for a in values)
  counts[sign]+=3;checks.append(dict(mode=mode,foreign=foreign,third_cosines=list(map(str,(x,y))),branch=branch,sign=sign,
    polynomial=[dict(exponents=list(e),coefficient=enc(v))for e,v in sorted(poly.items())],bernstein=list(map(enc,values))))
 expected={(m,f,x,y,b)for m in('shared','opposite')for f,x,y in foreign_patterns for b in(['last7','last_ge8']if(m,f,x,y)==('opposite',0,0,1)else['all_last_ge7'])}
 assert keys==expected and len(keys)==9 and counts=={1:24,-1:3}
 assert len(d['families'])==len(d['literal_orientations'])==8
 for i,(mode,orev,lrev)in enumerate(product(('opposite','shared'),(False,True),(False,True))):
  perm=list(range(18))
  if orev:perm[3],perm[4]=4,3
  if lrev:perm[12],perm[13]=13,12
  assert d['literal_orientations'][i]==dict(anchor_mode=mode,permutation=perm)
  assert{frozenset(perm[x]for x in z)for z in V}==set(map(frozenset,V))and{frozenset(perm[x]for x in z)for z in MN}==MN
  fixed=[[3,14,2],[4,14,4],[11,12,2],[11,13,3],[12,17,2 if mode=='shared'else 3],[13,17,3 if mode=='shared'else 2]]
  bounds=[[3,4,5],[6,7,7],[9,10,7],[12,13,7]]
  for name,rows in(('fixed',fixed),('bounds',bounds)):assert sorted([sorted((perm[a],perm[b]))+[m]for a,b,m in rows])==d['families'][i][name]
 batch=read(f'results/{PRE}_search_batch_014.json');labels={tuple(r[:2]):r[2]for r in batch['solution']['pair_labels']}
 assert batch['cumulative_additional_nodes']==d['covered_saved_cursor']==60015
 assert any(all(labels[a,b]==m for a,b,m in f['fixed'])and all(labels[a,b]>=m for a,b,m in f['bounds'])for f in d['families'])
 paths=[f'scripts/{TAG}.py',f'results/{TAG}.json',f'scripts/{TAG}_review.py',f'scripts/{PRE}_review.py',
   'scripts/fresh_h9_star_single_third_20261014_review.py','scripts/fresh_h9_star_lastblock_20261014_review.py']
 return dict(status='PASS_INDEPENDENT_SINGLE_PI4_THIRD_BOTH_ANCHOR_MODES',sources={p:sha(p)for p in paths},
  arithmetic='Separate Fraction Qsqrt2sqrt5 and sparse permutation determinant arithmetic; exact quadratic-field positive kernels. No producer or SymPy imports.',
  recursive_manifests=len(seen),literal_guards=108,projection_budgets=6,guard_counts=dict(guardcounts),domains=domain_summaries,
  no_foreign_identity=True,no_foreign_regions=12,foreign_polynomial_regions=9,positive_coefficients=24,negative_coefficients=3,
  independent_polynomials=checks,literal_orientations=8,covered_cursor=60015,
  scope='Owner internal>=5, otherthree>=7, singlepi4third, sharedoropposite singlepi3anchors, allnoncorefree. Fullnewdomain andpolynomialreview; nowholeincidenceexclusion.')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}_review.json'
 if a.verify:assert d==json.loads(out.read_text())
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['literal_guards'],d['foreign_polynomial_regions'],d['positive_coefficients'],d['negative_coefficients'])
