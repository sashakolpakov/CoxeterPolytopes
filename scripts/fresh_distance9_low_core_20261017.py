#!/usr/bin/env python3
"""Complete small-core domain for selective distance9; no dimension claim.

All labels are bounded analytically before exhaustive graph enumeration.
PSD recognition uses exact tree Schur complements and the affine-cycle
kernel criterion; joins use exact Q(sqrt5) leaf pivots plus at most five
remaining symbolic cycle nodes. No floating-point arithmetic is used.
"""
import argparse,json,sys
from pathlib import Path
from itertools import combinations,product,permutations
from functools import lru_cache
from collections import Counter,deque
from hashlib import sha256
import sympy as s
from fresh_selective_distance10_dimension116_20261017_review import K,F
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance9_low_core_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
BASE='fresh_selective_distance10_dimension116_20261017'
REVIEW=BASE+'_review'

def digest(p):
 h=sha256()
 with p.open('rb')as f:
  for b in iter(lambda:f.read(1<<20),b''):h.update(b)
 return h.hexdigest()

BETA={0:K(0),3:K(1),4:K(2),5:K(F(3,2),F(1,2)),6:K(3)}

def topology(n,code):
 es=list(combinations(range(n),2));adj=[{}for _ in range(n)]
 for(i,j),m in zip(es,code):
  if m:adj[i][j]=adj[j][i]=m
 return adj

def components(adj):
 rem=set(range(len(adj)));out=[]
 while rem:
  start=min(rem);rem.remove(start);part=[start]
  for i in part:
   for j in adj[i]:
    if j in rem:rem.remove(j);part.append(j)
  out.append(part)
 return out

@lru_cache(None)
def verdict(n,code):
 adj=topology(n,code);strict=True
 for part in components(adj):
  nedges=sum(len(adj[i])for i in part)//2
  if nedges>=len(part):
   # A PSD acute Coxeter cycle is exactly an isolated simple affine A cycle.
   if nedges==len(part)and all(len(adj[i])==2 and all(m==3 for m in adj[i].values())for i in part):strict=False;continue
   return False,False
  assert nedges==len(part)-1
  diagonal={i:K(2)for i in part};a={i:dict(adj[i])for i in part}
  while a:
   i=next(i for i in a if len(a[i])<=1);q=diagonal[i]
   if q<0:return False,False
   if not q:
    if a[i]:return False,False
    strict=False
   if a[i]:
    j=next(iter(a[i]));diagonal[j]-=BETA[a[i][j]]/q;del a[j][i]
   del a[i];del diagonal[i]
 return True,strict


def subcode(n,code,ids):
 adj=topology(n,code)
 return tuple(adj[i].get(j,0)for i,j in combinations(ids,2))

def determinant_tree(n,code,ids):
 adj=topology(n,code)
 @lru_cache(None)
 def rec(nodes):
  if not nodes:return K(1)
  i=nodes[0];rest=nodes[1:];v=2*rec(rest)
  for j in rest:
   if j in adj[i]:v-=BETA[adj[i][j]]*rec(tuple(k for k in rest if k!=j))
  return v
 # Called only on a spherical deletion: a forest.
 h=topology(len(ids),subcode(n,code,ids));assert all(sum(len(h[i])for i in p)//2==len(p)-1 for p in components(h))
 return rec(tuple(ids))

def canonical(n,code):
 adj=topology(n,code)
 return min(tuple(adj[i].get(j,0)for i,j in combinations(perm,2))for perm in permutations(range(n)))

def distance(adj,marks):
 d={i:0 for i in marks};q=deque(marks)
 while q:
  i=q.popleft()
  for j in adj[i]:
   if j not in d:d[j]=d[i]+1;q.append(j)
 assert len(d)==len(adj)
 return [d[i]for i in range(len(adj))]


def census_order(n):
 alphabet=(0,3,4,5,6)if n==4 else(0,3,4,5);es=list(combinations(range(n),2))
 ids=[tuple(i for i in range(n)if i!=b)for b in range(n)]
 # Complete preceding-order PSD table, not a selected diagram catalogue.
 table={c:verdict(n-1,c)for c in product(alphabet,repeat=(n-1)*(n-2)//2)}
 enc=[[(es.index((i,j)))for i,j in combinations(part,2)]for part in ids]
 found={};labelled=0
 for code in product(alphabet,repeat=len(es)):
  flags=[]
  for ix in enc:
   c=tuple(code[k]for k in ix);v=table[c]
   if not v[0]:break
   flags.append(v[1])
  if len(flags)<n:continue
  adj=topology(n,code)
  if len(components(adj))!=1 or verdict(n,code)[0]:continue
  labelled+=1;cc=canonical(n,code)
  found[cc]=found.get(cc,0)+1
 print('ORDER',n,'labelled',labelled,'canonical',len(found),flush=True)
 rows=[]
 for index,(code,count)in enumerate(sorted(found.items())):
  finite=[b for b in range(n)if verdict(n-1,subcode(n,code,ids[b]))[1]]
  contexts=[];adj=topology(n,code)
  for black in combinations(finite,2):
   dd=distance(adj,black);ports=[]
   for i,d in enumerate(dd):
    if d<2:continue
    vals=[]
    for b in black:
     keep=ids[b];num=determinant_tree(n,code,tuple(j for j in keep if j!=i));den=determinant_tree(n,code,keep)
     assert den>0;vals.append(num/den)
    ports.append(dict(node=i,distance=d,inverse_coefficients=[v.save()for v in vals],allowed3=max(vals)<2,allowed4=max(vals)<1))
   if ports:contexts.append(dict(black=list(black),distances=dd,ports=ports))
  rows.append(dict(index=index,order=n,code=list(code),edges=[[i,j,m]for(i,j),m in zip(es,code)if m],finite_deletions=finite,labelled_multiplicity=count,contexts=contexts))
 print('CONTEXTS',n,sum(len(r['contexts'])for r in rows),'PORTS',sum(len(c['ports'])for r in rows for c in r['contexts']),'ALLOWED3',sum(p['allowed3']for r in rows for c in r['contexts']for p in c['ports']),flush=True)
 return dict(order=n,all_labelled_minimal_indefinite=labelled,canonical_types=len(rows),diagrams=rows)

def fromcoeff(v):return K(F(v[0]),F(v[1]))
def fullcode(n,edges):
 lookup={tuple(sorted((i,j))):m for i,j,m in edges}
 return tuple(lookup.get((i,j),0)for i,j in combinations(range(n),2))

def sym(q):return s.Rational(q.a.numerator,q.a.denominator)+s.Rational(q.b.numerator,q.b.denominator)*s.sqrt(5)

def inertia(edges,n):
 # Eliminate all available tree leaves exactly over squared-bond Q(sqrt5).
 a={i:{}for i in range(n)};diag={i:K(2)for i in range(n)};pos=neg=zero=0;steps=[]
 for i,j,m in edges:a[i][j]=a[j][i]=m
 while a:
  i=next((i for i in a if len(a[i])<=1),None)
  if i is None:break
  q=diag[i]
  if not a[i]:
   pos+=int(q>0);neg+=int(q<0);zero+=int(q==0);steps.append(dict(nodes=[i],pivot=q.save()));del a[i];del diag[i]
  elif q:
   j=next(iter(a[i]));w=BETA[a[i][j]];diag[j]-=w/q
   pos+=int(q>0);neg+=int(q<0);steps.append(dict(nodes=[i],pivot=q.save(),neighbor=j,squared_bond=w.save()));del a[j][i];del a[i];del diag[i]
  else:
   j=next(iter(a[i]));pos+=1;neg+=1;steps.append(dict(zero_leaf_pair=[i,j]))
   for k in list(a[j]):
    if k!=i:del a[k][j]
   del a[i];del a[j];del diag[i];del diag[j]
 # Any cyclic remainder has at most5 nodes, evaluated symbolically exactly.
 ids=list(a);h=s.zeros(len(ids));theta={3:s.Integer(1),4:s.sqrt(2),5:(1+s.sqrt(5))/2,6:s.sqrt(3)}
 for i,x in enumerate(ids):
  h[i,i]=sym(diag[x])
  for j,y in enumerate(ids):
   if x!=y and y in a[x]:h[i,j]=-theta[a[x][y]]
 while h.rows:
  i=next((i for i in range(h.rows)if h[i,i]!=0),None)
  if i is not None:
   q=s.simplify(h[i,i]);assert q.is_positive or q.is_negative
   pos+=int(bool(q>0));neg+=int(bool(q<0));rest=[j for j in range(h.rows)if j!=i]
   steps.append(dict(cyclic_pivot=str(q)))
   h=(h.extract(rest,rest)-h.extract(rest,[i])*h.extract([i],rest)/q).applyfunc(s.simplify)
  else:
   pair=next(((i,j)for i in range(h.rows)for j in range(i+1,h.rows)if h[i,j]!=0),None)
   if pair is None:zero+=h.rows;break
   i,j=pair;rest=[k for k in range(h.rows)if k not in pair];block=h.extract(pair,pair)
   assert block.det()<0;pos+=1;neg+=1;steps.append(dict(cyclic_two_pivot=[str(x)for x in block]))
   h=(h.extract(rest,rest)-h.extract(rest,pair)*block.inv()*h.extract(pair,rest)).applyfunc(s.simplify)
 return dict(signature=[pos,neg,zero],steps=steps)

def join_and_overlap(low,big,high):
    assert all(p['diagram']==25 for p in high)and len(high)==7
    hcode=fullcode(10,big['edges'])
    rows=[];types={};good={};counts=Counter()
    for order in('4','5'):
     for d in low[order]['diagrams']:
      n=d['order'];code=tuple(d['code'])
      for context in d['contexts']:
       u=context['black'];ww=tuple(i for i in range(n)if i not in u)
       for p in context['ports']:
        if not p['allowed3']:continue
        x=p['node'];lv=max(fromcoeff(v)for v in p['inverse_coefficients'])
        lx=determinant_tree(n,code,tuple(i for i in ww if i!=x))/determinant_tree(n,code,ww)
        for hp in high:
         if p['distance']+hp['distance']<9:continue
         v=hp['black'];y=hp['node'];hww=tuple(i for i in range(10)if i not in v)
         hy=determinant_tree(10,hcode,tuple(i for i in hww if i!=y))/determinant_tree(10,hcode,hww)
         hv=K(max(map(F,hp['inverse_diagonals'])))
         for m,bs in((3,K(1)),(4,K(2))):
          if bs*lv>=2 or bs*hv>=2:continue
          row=dict(low=[n,d['index'],u,x],high=[25,v,y],cross_label=m,white_inverse_product=(bs*lx*hy).save())
          if bs*lx*hy>=1:row['kind']='white_PD_rejection'
          else:
           key=(n,d['index'],x,y,m)
           if key not in types:
            ed=d['edges']+[[i+n,j+n,mm]for i,j,mm in big['edges']]+[[x,y+n,m]]
            types[key]=inertia(ed,n+10)
            print('TYPE',key,types[key]['signature'],flush=True)
           row['kind']='two_negative'if types[key]['signature'][1]>=2 else'survive'
           row['inertia_type']=list(key)
           ck=(n,d['index'],tuple(u),tuple(v));good.setdefault(ck,[]).append((x,y,m))
          counts[n,m,row['kind']]+=1;rows.append(row)
    print('JOINS',dict(counts),'types',len(types),'multiple',[(k,v)for k,v in good.items()if len(v)>1],flush=True)
    # Low5 radius3 can meet only high radius7 in an overlap, with one shared root.
    ov=[]
    for d in low['5']['diagrams']:
     for c in d['contexts']:
      for p in c['ports']:
       if p['distance']!=3:continue
       for hp in high:
        if hp['distance']!=7:continue
        row=dict(low=[5,d['index'],c['black'],p['node']],high=[25,hp['black'],hp['node']])
        if not p['allowed3']:
         row['kind']='inverse_rejection';ov.append(row);continue
        mp={hp['node']:p['node']};nextnode=5
        for i in range(10):
         if i not in mp:mp[i]=nextnode;nextnode+=1
        ed=d['edges']+[[mp[i],mp[j],m]for i,j,m in big['edges']]
        code=fullcode(nextnode,ed);black=c['black']+[mp[i]for i in hp['black']];white=[i for i in range(nextnode)if i not in black]
        bad=None
        for b in black:
         ids=tuple(sorted(white+[b]));sub=subcode(nextnode,code,ids)
         if not verdict(len(ids),sub)[1]:bad=ids;break
        if bad:row.update(kind='actual_edge_PD_rejection',bad_nodes=bad)
        else:row.update(kind='survive',edges=ed,black=black,inertia=inertia(ed,nextnode))
        ov.append(row)
    print('OVERLAPS',Counter(r['kind']for r in ov),flush=True)
    assert len(rows)==104 and len(types)==23
    assert counts=={(4,3,'two_negative'):33,(4,4,'white_PD_rejection'):2,(5,3,'two_negative'):60,(5,3,'white_PD_rejection'):2,(5,4,'white_PD_rejection'):7}
    multiple=[]
    for key,es in good.items():
     if len(es)<2:continue
     n,index,u,v=key;assert n==5 and len(es)==2
     d=next(d for d in low['5']['diagrams']if d['index']==index)
     edges=d['edges']+[[i+n,j+n,m]for i,j,m in big['edges']]+[[i,j+n,m]for i,j,m in es]
     code=fullcode(n+10,edges);adj=topology(n+10,code);white=[i for i in range(n+10)if i not in list(u)+[n+x for x in v]]
     witness=None
     for size in(3,4):
      for ids in combinations(white,size):
       h=topology(size,subcode(n+10,code,ids))
       if all(len(h[i])==2 and all(m==3 for m in h[i].values())for i in range(size)):
        witness=list(ids);break
      if witness:break
     assert witness
     multiple.append(dict(low=[n,index,list(u)],high=[25,list(v)],edges=list(map(list,es)),white_affine_cycle=witness,positive_kernel=[1]*len(witness)))
    assert len(multiple)==6
    assert Counter(r['kind']for r in ov)==dict(inverse_rejection=24,survive=6)
    return dict(joins=rows,join_counts={str(k):v for k,v in counts.items()},inertia_types={str(k):v for k,v in types.items()},overlaps=ov,multiple_edge_obstructions=multiple,
     low_first_simple_joins=95,low_first_possible_double_joins=9,low_first_two_negative=93,low_first_white_rejections=11,
     low_first_overlap_candidates=30,low_first_surviving_overlaps=6,ordered_surviving_overlaps=12)
TH={3:s.Integer(1),4:s.sqrt(2),5:(1+s.sqrt(5))/2,6:s.sqrt(3)}

def dist(g,marks):
 d={i:0 for i in marks};q=deque(marks)
 while q:
  i=q.popleft()
  for j in range(g.rows):
   if i!=j and g[i,j]and j not in d:d[j]=d[i]+1;q.append(j)
 return[d.get(i)for i in range(g.rows)]

def parts(g):
 rem=set(range(g.rows));out=[]
 while rem:
  p=[min(rem)];rem.remove(p[0])
  for i in p:
   for j in list(rem):
    if g[i,j]:p.append(j);rem.remove(j)
  out.append(sorted(p))
 return out

def pdtest(g):
 a=g.copy()
 for k in range(a.rows):
  q=s.simplify(a[k,k])
  if q<=0:return False
  for i in range(k+1,a.rows):
   for j in range(k+1,a.rows):a[i,j]=s.simplify(a[i,j]-a[i,k]*a[k,j]/q)
 return True
def extensions(joined):
    
    out=[]
    for r in joined['overlaps']:
     if r['kind']!='survive':continue
     n=max(max(i,j)for i,j,m in r['edges'])+1;g=2*s.eye(n)
     for i,j,m in r['edges']:g[i,j]=g[j,i]=-TH[m]
     u=r['black'][:2];v=r['black'][2:];white=[i for i in range(n)if i not in u+v]
     z=g.nullspace();iv=None if z else g.inv().applyfunc(s.simplify);du,dv=dist(g,u),dist(g,v)
     cases=[]
     for a,b in product(u,v):
      ids=white+[a,b];ids.sort();h=g.extract(ids,ids);allowed=sorted(set(u+v)-{a,b});blocks=[];marked=None
      for part in parts(h):
       actual=[ids[i]for i in part];hh=h.extract(part,part);ker=hh.nullspace()
       if ker:
        assert len(ker)==1;zz=ker[0]
        if all(x<0 for x in zz):zz=-zz
        assert all(x>0 for x in zz);blocks.append(dict(nodes=actual,affine=True))
       else:
        assert pdtest(hh);inv=hh.inv().applyfunc(s.simplify)
        allowed+= [actual[i]for i in range(len(part))if inv[i,i]<=2];blocks.append(dict(nodes=actual,inverse_diag=[str(inv[i,i])for i in range(len(part))]))
       if a in actual and b in actual:marked=actual
      if not marked:continue
      c=dict(preferred=[a,b],allowed=sorted(allowed),marked=marked,blocks=blocks)
      if z:
       c['radical']=[str(x)for x in z[0]]
       c['zero_allowed']=[i for i in allowed if z[0][i]==0]
       if not c['zero_allowed']:
        pos=[i for i in allowed if z[0][i]>0];neg=[i for i in allowed if z[0][i]<0]
        c['shortcut']=min(max(du[i]for i in pos)+max(dv[j]for j in neg)+2,max(dv[i]for i in pos)+max(du[j]for j in neg)+2)
      else:
       c['inverse_block']=[[str(iv[i,j])for j in sorted(allowed)]for i in sorted(allowed)]
       c['candidates']=[dict(node=i,diagonal=str(iv[i,i]),du=du[i],dv=dv[i],in_marked=i in marked)for i in sorted(allowed)if iv[i,i]<=2]
      cases.append(c)
     # Every extra cross edge must shorten or fail an actual-edge PD test.
     lowshared=r['low'][3];left=[i for i in range(5)if i!=lowshared];right=[i for i in range(5,n)]
     short=[];bad=[];unresolved=[]
     for i,j in product(left,right):
      if g[i,j]:continue
      gg=g.copy();gg[i,j]=gg[j,i]=-1
      if min(dist(gg,u)[b]for b in v)<=9:short.append([i,j]);continue
      fails=[]
      for b in u+v:
       ids=sorted(white+[b])
       if not pdtest(gg.extract(ids,ids)):fails.append(b)
      if fails:bad.append(dict(edge=[i,j],failed_marks=fails))
      else:unresolved.append([i,j])
     assert not unresolved and not bad and len(short)==36
     choices=[]
     for c in cases:
      if z:
       if not c['zero_allowed']and c['shortcut']<=9:
        c['method']='actual_radical_relation';c['component_order_bound']=len(c['marked']);choices.append(c)
      else:
       left=[i for i in c['allowed']if i<5];right=[i for i in c['allowed']if i>=5]
       if not left or not right:continue
       if not all(iv[i,j]>0 for block in(left,right)for i in block for j in block):continue
       if not all(iv[i,i]>2 for i in c['allowed']):continue
       bound=max(du[i]for i in left)+max(dv[i]for i in right)+2
       if bound>9:continue
       c.update(method='positive_inverse_end_blocks',end_sets=[left,right],shortcut=bound,component_order_bound=len(c['marked']))
       choices.append(c)
     assert choices
     selected=min(choices,key=lambda c:(len(c['allowed']),len(c['marked']),c['preferred']))
     assert selected['component_order_bound']<=12
     rr=dict(selected_corner=selected,low=r['low'],high=r['high'],corners=cases,extra_short=short,extra_pd_bad=bad,extra_unresolved=unresolved)
     out.append(rr)
     print('CASE',r['low'],r['high'],'extras',len(short),len(bad),unresolved,flush=True)
     for c in cases:print('CORNER',c['preferred'],'component',len(c['marked']),'allowed',c['allowed'],'zero',c.get('zero_allowed'),'shortcut',c.get('shortcut'),'candidates',c.get('candidates'),flush=True)
    assert len(out)==6 and max(r['selected_corner']['component_order_bound']for r in out)==12
    return dict(cases=out,low_first_extra_cross_shortcuts=sum(len(r['extra_short'])for r in out),component_order_cap=12)

def run():
 low={str(n):census_order(n)for n in(4,5)}
 assert [(low[str(n)]['all_labelled_minimal_indefinite'],low[str(n)]['canonical_types'])for n in(4,5)]==[(392,32),(1050,14)]
 assert all(sum(d>=3 for d in c['distances'])<=1 for r in low['5']['diagrams']for c in r['contexts'])
 cat=json.loads((ROOT/f'results/{CENSUS}.json').read_text());base=json.loads((ROOT/f'results/{BASE}.json').read_text())
 for package in(cat,base):
  for path,h in package['sources'].items():assert digest(ROOT/path)==h,(path,h)
 big=next(d for d in cat['diagrams']if d['index']==25)
 high=[p for p in base['census']['ports']if p['distance']>=6 and p['allowed']]
 joined=join_and_overlap(low,big,high);ext=extensions(joined)
 paths=[Path(__file__),ROOT/f'scripts/{REVIEW}.py']+[ROOT/f'{folder}/{stem}.{suffix}'for stem in(CENSUS,BASE)for folder,suffix in(('scripts','py'),('results','json'),('results','md'))]
 return dict(status='PASS_COMPLETE_DISTANCE9_SMALL_CORE_DOMAIN',sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
  classification=low,high_ports=high,disjoint_and_overlap=joined,extensions=ext,
  theorem='Under the actual quadrilateral-star and no marked white path<=9 hypotheses, all configurations containing a minimal indefinite core of order<=5 are excluded or have a preferred marked component isolated with order<=12.',
  scope='A complete small-core lemma, not the remaining high-rank census or an ambient dimension theorem. Label bounds and arbitrary-supporting-white isolation have written proofs; no numerical cutoffs, construction or novelty claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
 out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
 if args.verify:assert out==json.loads(path.read_text())
 else:path.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'32+14types/104joins/6overlaps',flush=True)
