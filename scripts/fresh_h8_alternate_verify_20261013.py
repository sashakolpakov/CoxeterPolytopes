#!/usr/bin/env python3
"""Independent Fraction-only replay of alternate released-H8 cap certificates.

No project imports, NumPy, SciPy or symbolic-algebra package. The cap
proposal searches are not rerun; every accepted normal, field identity,
conflict and strict escape used in the finite-bank obstruction is replayed.
"""
from fractions import Fraction as F
from itertools import combinations, product
from functools import lru_cache
from pathlib import Path
import argparse
import json

Z=(F(0),F(0));O=(F(1),F(0));P=(F(0),F(1))
def parse(x):return F(x['a']),F(x['b'])
def enc(x):return {'a':str(x[0]),'b':str(x[1])}
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def sub(x,y):return add(x,neg(y))
def scale(x,k):return x[0]*k,x[1]*k
def mul(x,y):return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def inv(x):
    n=x[0]*x[0]+x[0]*x[1]-x[1]*x[1]
    return (x[0]+x[1])/n,-x[1]/n
def div(x,y):return mul(x,inv(y))
def conj(x):return x[0]+x[1],-x[1]
def sign(x):
    a=2*x[0]+x[1];b=x[1]
    if not b:return (a>0)-(a<0)
    if not a or (a>0)==(b>0):return (b>0)-(b<0)
    return ((a>0)-(a<0))*((a*a>5*b*b)-(a*a<5*b*b))
def integral(x):return all(c.denominator==1 for c in x)
def total(values):
    out=Z
    for x in values:out=add(out,x)
    return out
def dot(x,y):return total(mul(a,b) for a,b in zip(x,y,strict=True))
def mv(a,v):return tuple(dot(row,v) for row in a)
def tr(a):return tuple(zip(*a))
def mm(a,b):return tuple(tuple(dot(row,col) for col in tr(b)) for row in a)
def pair(x,a,y):return dot(x,mv(a,y))
def matrix(x):return tuple(tuple(map(parse,row)) for row in x)
def eye(n):return tuple(tuple(O if i==j else Z for j in range(n)) for i in range(n))

def rref(a):
    a=[list(row) for row in a];cols=len(a[0]);pivots=[]
    for j in range(cols):
        k=len(pivots)
        pivot=next((i for i in range(k,len(a)) if a[i][j]!=Z),None)
        if pivot is None:continue
        a[k],a[pivot]=a[pivot],a[k]
        q=a[k][j];a[k]=[div(x,q) for x in a[k]]
        for i in range(len(a)):
            if i==k:continue
            q=a[i][j];a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[k])]
        pivots.append(j)
        if len(pivots)==len(a):break
    return a,pivots

def inverse(a):
    n=len(a);rows,pivots=rref(tuple(tuple(row)+ident for row,ident in zip(a,eye(n))))
    assert pivots==list(range(n))
    return tuple(tuple(row[n:]) for row in rows)

def determinant(a):
    a=[list(row) for row in a];out=O
    for j in range(len(a)):
        k=next((i for i in range(j,len(a)) if a[i][j]!=Z),None)
        if k is None:return Z
        if k!=j:a[k],a[j]=a[j],a[k];out=neg(out)
        pivot=a[j][j];out=mul(out,pivot)
        for i in range(j+1,len(a)):
            q=div(a[i][j],pivot)
            a[i]=[sub(x,mul(q,y)) for x,y in zip(a[i],a[j])]
    return out

def kernel_one(a):
    rows,pivots=rref(a);free=[i for i in range(len(a[0])) if i not in pivots]
    if len(free)!=1:return None
    out=[Z]*len(a[0]);out[free[0]]=O
    for row,pivot in zip(rows,pivots):out[pivot]=neg(row[free[0]])
    return tuple(out)

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h8_alternate_verify_20261013.json'
def ve(x):return tuple(map(parse,x))
def ea(x,y):return add(x[0],y[0]),add(x[1],y[1])
def en(x):return neg(x[0]),neg(x[1])
def es(x,y):return ea(x,en(y))
def em(x,y,d):return add(mul(x[0],y[0]),mul(d,mul(x[1],y[1]))),add(mul(x[0],y[1]),mul(x[1],y[0]))
def sg(x,d):
 a,b=x;sa,sb=sign(a),sign(b)
 if not sb:return sa
 if not sa or sa==sb:return sb
 return sa*sign(sub(mul(a,a),mul(d,mul(b,b))))
def pe(x):return parse(x['rational']),parse(x['radical'])
def angle(q):return sign(sub(q,scale(O,4)))>=0 or q in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P))
def bsg(x,d,e):
 a=x[:2];b=x[2:];sa,sb=sg(a,d),sg(b,d)
 if not sb:return sa
 if not sa or sa==sb:return sb
 q=em(b,b,d);return sa*sg(es(em(a,a,d),(mul(e,q[0]),mul(e,q[1]))),d)
def bm(x,y,d,e):
 out=[Z]*4
 for i,a in enumerate(x):
  if a==Z:continue
  for j,b in enumerate(y):
   if b==Z:continue
   v=mul(a,b)
   if i&j&1:v=mul(v,d)
   if i&j&2:v=mul(v,e)
   out[i^j]=add(out[i^j],v)
 return tuple(out)

def check_escape(row,I,center,old,caps=()):
 assert row['status']=='EXACT_STRICT_SPACELIKE_ESCAPE'
 y=ve(row['y']);q=pair(y,I,y);cp=pair(y,I,center)
 assert q==parse(row['norm']) and sign(q)>0
 assert cp==parse(row['center_pairing']) and sign(cp)<0
 assert all(sign(x)>0 for x in y) and all(sign(dot(a,y))>0 for a in old)
 assert all(sg((dot(a,y),dot(b,y)),d)>0 for d,a,b in caps)
 return y

def verify_base(label,deleted,source):
 raw=json.loads((ROOT/f'data/enumerations/fresh_h8_alternate_{label}_{deleted}_20261013.json').read_text())
 oldQ=matrix(source['lattice_gram']);original=matrix(source['simple_roots'])
 extra=scale(P,2)if label==3 else sub(P,O)
 Q=tuple(tuple((*r,Z))for r in oldQ)+((*(Z for _ in range(9)),extra),)
 Rall=tuple(tuple((*r,O if i==0 else Z))for i,r in enumerate(original));Rall=(tuple(div(x,P)for x in Rall[0]),*Rall[1:])
 ids=[i for i in range(11)if i!=deleted];R=tuple(Rall[i]for i in ids);G=tuple(tuple(pair(r,Q,s)for s in R)for r in R);I=inverse(G)
 assert matrix(raw['host_gram'])==Q and matrix(raw['seed_normals'])==R and matrix(raw['seed_gram'])==G and matrix(raw['seed_gram_inverse'])==I
 assert len(rref(R)[1])==raw['seed_rank']==10 and raw['retained_original_walls']==ids
 assert all(sign(G[i][j])<=0 and angle(div(scale(mul(G[i][j],G[i][j]),4),mul(G[i][i],G[j][j])))for i in range(10)for j in range(i))
 h=(*ve(source['geometry']['interior_point']),Z);center=tuple(neg(pair(h,Q,r))for r in R)
 assert center==ve(raw['strict_center_y']) and sign(pair(center,I,center))<0 and all(sign(x)>0 for x in center)
 assert ve(raw['dual_vertex_norms'])==tuple(I[i][i]for i in range(10))
 old=[];norms=[]
 for cap in raw['caps']:
  p=ve(cap['seed_pairings']);c=mv(I,p);r=tuple(total(mul(a,w[j])for a,w in zip(c,R))for j in range(10));q=pair(r,Q,r)
  assert r==ve(cap['normal']) and c==ve(cap['seed_coefficients']) and q==parse(cap['norm']) in (O,scale(O,2))
  assert tuple(pair(r,Q,w)for w in R)==p and dot(p,c)==q and sign(dot(c,center))>0
  assert all(sign(v)<=0 and angle(div(scale(mul(v,v),4),mul(q,G[j][j])))for j,v in enumerate(p))
  old.append(c);norms.append(q)
 check_escape(raw['volume_gate'],I,center,old)
 polarraw=json.loads((ROOT/f'data/enumerations/fresh_h8_alternate_{label}_{deleted}_polar_20261013.json').read_text());expected=[];polars=[]
 for i in range(10):
  if sign(I[i][i])>0 and angle(div(scale(O,4),mul(I[i][i],G[i][i]))):expected.append(i)
 assert [r['seed_index']for r in polarraw['polar_caps']]==expected
 for row in polarraw['polar_caps']:
  i=row['seed_index'];c=tuple(neg(x)for x in I[i]);q=I[i][i]
  assert c==ve(row['seed_coefficients']) and q==parse(row['norm']) and pair(c,G,c)==q
  assert row['original_wall']==ids[i] and parse(row['four_cosine_squared'])==div(scale(O,4),mul(q,G[i][i]))
  assert sign(dot(c,center))>0 and mv(G,c)==tuple(neg(O)if j==i else Z for j in range(10))
  polars.append(c)
 check_escape(polarraw['volume_gate'],I,center,old+polars)
 return raw,Q,R,G,I,center,old,norms,polars

def verify_coupled(path,base,include_polar):
 raw=json.loads(path.read_text());_,Q,R,G,I,center,old,norms,polars=base
 old=list(old);norms=list(norms)
 if include_polar:
  old+=polars;norms+=[pair(c,G,c)for c in polars]
 caps=[];target_index=0
 for rd in raw['rounds']:
  target=ve(rd['target_escape_y']);count=rd['exact_new_caps']
  assert len(caps)+count==rd['total_new_caps']
  for row in raw['caps'][target_index:target_index+count]:
   d=parse(row['discriminant']);assert sign(d)>0
   pa,pb=ve(row['pairings_rational']),ve(row['pairings_radical']);ca,cb=mv(I,pa),mv(I,pb)
   assert ca==ve(row['coefficients_rational'])and cb==ve(row['coefficients_radical'])
   norm=(add(dot(pa,ca),mul(d,dot(pb,cb))),scale(dot(pa,cb),2));assert sg(es(norm,(scale(O,2),Z)),d)==0
   assert all(sg((a,b),d)<=0 for a,b in zip(pa,pb))
   u,v=row['free_seed_indices'];assert u<v and sg(es(en((pa[u],pb[u])),(scale(O,2),Z)),d)>0 and sg(es(en((pa[v],pb[v])),(scale(O,2),Z)),d)>0
   f=ve(row['fixed_pairing_magnitudes']);assert f[u]==f[v]==Z
   for j in range(10):
    if j not in (u,v):assert pa[j]==neg(f[j])and pb[j]==Z
    squared=em((pa[j],pb[j]),(pa[j],pb[j]),d);ratio=em(squared,(div(scale(O,2),G[j][j]),Z),d)
    assert sg(es(ratio,(scale(O,4),Z)),d)>0 or any(sg(es(ratio,(q,Z)),d)==0 for q in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)))
   a=row['anchor_base_index'];M=parse(row['anchor_pairing_magnitude']);assert sg(ea((dot(pa,old[a]),dot(pb,old[a])),(M,Z)),d)==0
   assert angle(div(scale(mul(M,M),2),norms[a]))
   A,B,C=ve(row['norm_polynomial_descending']);assert d==sub(mul(B,B),scale(mul(A,C),4));x=pe(row['free_x'])
   L=div(sub(M,dot(old[a],f)),old[a][v]);T=div(neg(old[a][u]),old[a][v]);w0=list(f);w0[v]=L;w0=tuple(w0)
   w1=tuple(O if j==u else T if j==v else Z for j in range(10))
   assert A==pair(w1,I,w1) and B==scale(pair(w0,I,w1),2) and C==sub(pair(w0,I,w0),scale(O,2))
   assert sg(ea(ea(em((A,Z),em(x,x,d),d),em((B,Z),x,d)),(C,Z)),d)==0 and sg(ea(x,(pa[u],pb[u])),d)==0
   assert pe(row['target_escape_pairing'])==(dot(ca,target),dot(cb,target)) and sg((dot(ca,target),dot(cb,target)),d)<0
   assert sg((dot(ca,center),dot(cb,center)),d)>0
   # Reconstruct actual Lorentz normals in both field coefficients.
   ra=tuple(total(mul(c,r[j])for c,r in zip(ca,R))for j in range(10));rb=tuple(total(mul(c,r[j])for c,r in zip(cb,R))for j in range(10))
   assert tuple(pair(ra,Q,r)for r in R)==pa and tuple(pair(rb,Q,r)for r in R)==pb
   assert sg(es((add(pair(ra,Q,ra),mul(d,pair(rb,Q,rb))),scale(pair(ra,Q,rb),2)),(scale(O,2),Z)),d)==0
   caps.append((d,ca,cb))
  target_index+=count
  if rd['volume_gate']['status']=='EXACT_STRICT_SPACELIKE_ESCAPE':check_escape(rd['volume_gate'],I,center,old,caps)
  else:assert rd['volume_gate']['status']=='UNDECIDED_NO_NUMERIC_ESCAPE'
 assert target_index==len(raw['caps'])
 certified=[r['volume_gate']for r in raw['rounds']if r['volume_gate']['status']=='EXACT_STRICT_SPACELIKE_ESCAPE']
 forged=dict(certified[-1]);forged['y']=[enc(x)for x in center]
 try:check_escape(forged,I,center,old,caps)
 except AssertionError:pass
 else:raise AssertionError('timelike forged escape accepted')
 return {'file':str(path.relative_to(ROOT)),'exact_new_caps':len(caps),'verified_rounds':len(raw['rounds']),'certified_escape_rounds':len(certified),'all_archived_escape_certificates_pass':True,'last_all_cap_gate':raw['rounds'][-1]['volume_gate']['status'],'forged_timelike_escape_rejected':True}

def check_selection_and_sat(base):
 data,Q,R,G,I,center,old,norms,polars=base
 raw=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json').read_text())
 selection=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_selection_20261013.json').read_text())
 sat=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_sat_20261013.json').read_text())
 graphsat=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_graph_sat_20261013.json').read_text())
 nodes=[]
 for i,(c,q)in enumerate(zip(old,norms)):nodes.append((q,Z,c,(Z,)*10,mv(G,c),(Z,)*10,'base',i))
 for i,c in enumerate(polars):nodes.append((pair(c,G,c),Z,c,(Z,)*10,mv(G,c),(Z,)*10,'polar',i))
 for i,row in enumerate(raw['caps']):nodes.append((scale(O,2),parse(row['discriminant']),ve(row['coefficients_rational']),ve(row['coefficients_radical']),ve(row['pairings_rational']),ve(row['pairings_radical']),'coupled',i))
 assert selection['nodes']==[{'source':n[6],'index':n[7]}for n in nodes]
 def npair(i,j):
  a,b=nodes[i],nodes[j]
  return (dot(a[4],b[2]),dot(a[5],b[2]),dot(a[4],b[3]),dot(a[5],b[3]))
 seenpairs=set()
 def case(row,angles):
  chosen=row.get('selected_node_indices',row.get('selected_indices'));assert len(chosen)==len(set(chosen))
  for i,j in combinations(chosen,2):
   if (i,j,angles)in seenpairs:continue
   seenpairs.add((i,j,angles));a,b=nodes[i],nodes[j];p=npair(i,j);assert bsg(p,a[1],b[1])<=0
   if angles:
    factor=div(scale(O,4),mul(a[0],b[0]));rat=tuple(mul(factor,x)for x in bm(p,p,a[1],b[1]))
    assert bsg((sub(rat[0],scale(O,4)),*rat[1:]),a[1],b[1])>0 or any(bsg((sub(rat[0],q),*rat[1:]),a[1],b[1])==0 for q in (Z,O,scale(O,2),scale(O,3),add(O,P),add(scale(O,2),P)))
  chosenold=[nodes[i][2]for i in chosen if nodes[i][1]==Z];chosencaps=[(nodes[i][1],nodes[i][2],nodes[i][3])for i in chosen if nodes[i][1]!=Z]
  check_escape(row['volume_gate'],I,center,chosenold,chosencaps)
 for row in selection['cases']:case(row,True)
 for row in sat['checks']:case(row,False)
 for row in graphsat['checks']:case(row,True)
 # Exact escape-cover clauses: each cutter list must be COMPLETE over all493.
 clauses=sat['escape_clauses']+graphsat['new_escape_clauses']
 for row in clauses:
  y=check_escape(row['certificate'],I,center,[])
  assert row['cutters']==[i for i,n in enumerate(nodes)if sg((dot(n[2],y),dot(n[3],y)),n[1])<0]
 for row in sat['positive_pair_conflicts']:
  i,j=row['indices'];p=npair(i,j);assert p==ve(row['pairing']) and bsg(p,nodes[i][1],nodes[j][1])>0
 # A separate standard-library DPLL verifies abstract accepted-edge UNSAT.
 # The search graph's geometric completeness is an independent question.
 N=len(nodes);neighbors=[0]*N
 for key in selection['edges']:
  i,j=map(int,key.split(','));assert 0<=i<j<N;neighbors[i]|=1<<j;neighbors[j]|=1<<i
 masks=[sum(1<<i for i in row['cutters'])for row in clauses];hit=[0]*N
 for k,mask in enumerate(masks):
  for i in range(N):
   if mask>>i&1:hit[i]|=1<<k
 visits=0
 @lru_cache(None)
 def solve(available,unsatisfied):
  nonlocal visits
  visits+=1
  if not unsatisfied:return True
  choices=None;bits=unsatisfied
  while bits:
   lsb=bits&-bits;k=lsb.bit_length()-1;bits-=lsb;candidate=masks[k]&available
   if not candidate:return False
   if choices is None or candidate.bit_count()<choices.bit_count():choices=candidate
  while choices:
   bit=choices&-choices;i=bit.bit_length()-1;choices-=bit
   if solve(available&neighbors[i],unsatisfied&~hit[i]):return True
   available&=~bit
  return False
 assert graphsat['status']=='UNSAT_FOR_ACCEPTED_EDGE_GRAPH'
 assert not solve((1<<N)-1,(1<<len(masks))-1)
 all_clause_states=visits
 core=json.loads((ROOT/'data/enumerations/fresh_h8_alternate_graph_core_20261013.json').read_text());ci=core['clause_indices']
 assert core['sources']==[clauses[i]['source']for i in ci] and core['cutter_counts']==[len(clauses[i]['cutters'])for i in ci]
 union=sorted(set(j for i in ci for j in clauses[i]['cutters']));assert union==core['union_nodes']
 assert not solve(sum(1<<i for i in union),sum(1<<i for i in ci))
 compact_visits=visits
 for i,j in ((20,114),(21,114)):
  a,b=nodes[i],nodes[j];p=npair(i,j);assert bsg(p,a[1],b[1])<0
  rat=tuple(mul(div(scale(O,4),mul(a[0],b[0])),x)for x in bm(p,p,a[1],b[1]))
  assert bsg((sub(rat[0],scale(O,4)),*rat[1:]),a[1],b[1])==0
  neighbors[i]|=1<<j;neighbors[j]|=1<<i
 solve.cache_clear();visits=0
 assert not solve(sum(1<<i for i in union),sum(1<<i for i in ci))
 return {'dictionary_nodes':N,'greedy_subset_escapes':len(selection['cases']),'nonpositive_SAT_subset_escapes':len(sat['checks']),'accepted_graph_SAT_subset_escapes':len(graphsat['checks']),'exact_complete_escape_cover_clauses':len(clauses),'exact_positive_pair_conflicts':len(sat['positive_pair_conflicts']),'standard_library_accepted_graph_UNSAT':True,'DPLL_full_clause_states':all_clause_states,'DPLL_total_states_including_compact_core':compact_visits,'independent_five_clause_core_UNSAT':True,'core_nodes':len(union),'parallel_edges_added':[[20,114],[21,114]],'core_with_parallel_edges_UNSAT':True,'DPLL_core_with_parallel_states':visits,'scope':'Core UNSAT includes both parallel edges. Completeness of its angle graph is supplied by the separate independent54-node graph audit.'}

def main():
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args()
 source=json.loads((ROOT/'data/known_examples/bugaenko_compact_h8.json').read_text());oldQ=matrix(source['lattice_gram'])
 assert oldQ[0][0]==scale(P,-2) and all(oldQ[0][i]==oldQ[i][0]==Z for i in range(1,9))
 assert all(sign(determinant(tuple(tuple(oldQ[i][j]for j in range(1,n+1))for i in range(1,n+1))))>0 for n in range(1,9))
 bases={(l,k):verify_base(l,k,source)for l,k in ((3,5),(3,6),(4,4))}
 rows=[verify_coupled(ROOT/'data/enumerations/fresh_h8_alternate_coupled_20261013.json',bases[(4,4)],False)]
 extra=ROOT/'data/enumerations/fresh_h8_alternate_polar_coupled_20261013.json'
 if extra.exists():rows.append(verify_coupled(extra,bases[(3,6)],True))
 out={'status':'PASS','standard_library_only':True,'seed_comparisons':[{'label':l,'omitted_wall':k,'caps':len(v[6]),'coxeter_polar_caps':len(v[8]),'hyperideal_dual_rays':sum(sign(v[4][i][i])>0 for i in range(10))}for(l,k),v in bases.items()],'coupled_experiments':rows,'selection_audit':check_selection_and_sat(bases[(3,6)]),'scope':'Exact finite-dictionary certificates only; numerical proposal scans are not exhaustive and no unrestricted completion obstruction is asserted.'}
 if args.verify:assert out==json.loads(OUTPUT.read_text())
 else:OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps(out,indent=2))
if __name__=='__main__':main()
