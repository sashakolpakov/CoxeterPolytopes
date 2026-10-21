#!/usr/bin/env python3
"""Distance8 small-core research certificate; no dimension claim.

The order<=3 exclusion is analytic and all-label. The finite order4/5
catalogue is inherited intact; its new marked domains are recomputed.
"""
import argparse,json
from pathlib import Path
from itertools import combinations,product,permutations
from functools import lru_cache
from collections import Counter
import sympy as s
import fresh_distance9_low_core_20261017 as L
from fresh_selective_distance10_dimension116_20261017_review import F,K
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_distance8_low_core_20261017'
CENSUS='fresh_high_rank_quasilanner_census_20261017'
LOW='fresh_distance9_low_core_20261017'
HIGH='fresh_selective_distance9_dimension108_20261017'
TH={3:s.Integer(1),4:s.sqrt(2),5:(1+s.sqrt(5))/2,6:s.sqrt(3)}


def gram(n,edges):
 g=2*s.eye(n)
 for i,j,m in edges:g[i,j]=g[j,i]=-TH[m]
 return g

def encode(g):return [[str(g[i,j])for j in range(g.cols)]for i in range(g.rows)]

def inertia_small(g):
 h=g.copy();sig=[0,0,0];pivots=[]
 while h.rows:
  i=next((i for i in range(h.rows)if h[i,i]!=0),None)
  if i is not None:
   q=s.simplify(h[i,i]);assert q.is_positive or q.is_negative,q
   sig[0 if q>0 else 1]+=1;rest=[j for j in range(h.rows)if j!=i];pivots.append(str(q))
   h=(h.extract(rest,rest)-h.extract(rest,[i])*h.extract([i],rest)/q).applyfunc(s.simplify)
  else:
   pair=next(((i,j)for i in range(h.rows)for j in range(i+1,h.rows)if h[i,j]!=0),None)
   if pair is None:sig[2]+=h.rows;break
   rest=[k for k in range(h.rows)if k not in pair];block=h.extract(pair,pair);assert block.det()<0;sig[0]+=1;sig[1]+=1
   h=(h.extract(rest,rest)-h.extract(rest,pair)*block.inv()*h.extract(pair,rest)).applyfunc(s.simplify);pivots.append('hyperbolic_pair')
 return dict(signature=sig,pivots=pivots)


def union(d,e,match):
 mp={j:i for i,j in match};n=d['order']
 for j in range(e['order']):
  if j not in mp:mp[j]=n;n+=1
 edges={tuple(sorted((i,j))):m for i,j,m in d['edges']}
 for i,j,m in e['edges']:
  k=tuple(sorted((mp[i],mp[j])))
  if k in edges:assert edges[k]==m
  edges[k]=m
 return n,[[i,j,m]for(i,j),m in sorted(edges.items())],[mp[j]for j in range(e['order'])]


def edge_pd(n,edges,black):
 code=L.fullcode(n,edges);white=[i for i in range(n)if i not in black]
 for b in black:
  ids=tuple(sorted(white+[b]));cc=L.subcode(n,code,ids)
  if not L.verdict(len(ids),cc)[1]:return dict(mark=b,nodes=list(ids),code=list(cc),certificate=positive_witness(len(ids),cc))
 return None


def small_contexts(low):
 out=[]
 for order in('4','5'):
  for d in low['classification'][order]['diagrams']:
   n=d['order'];code=tuple(d['code']);adj=L.topology(n,code)
   for bb in combinations(d['finite_deletions'],2):
    dd=L.distance(adj,bb);ports=[]
    for i,dist in enumerate(dd):
     if dist==0:continue
     vals=[]
     for b in bb:
      keep=tuple(j for j in range(n)if j!=b)
      vals.append(L.determinant_tree(n,code,tuple(j for j in keep if j!=i))/L.determinant_tree(n,code,keep))
     ports.append(dict(node=i,distance=dist,inverse_coefficients=[q.save()for q in vals],allowed3=max(vals)<2,allowed4=max(vals)<1))
    out.append(dict(order=n,index=d['index'],black=list(bb),distances=dd,ports=ports))
 return out


def triangle_gate(high,diagrams):
 ports=[p for p in high['census']['ports']if p['distance']==7 and p['allowed']]
 assert [(p['diagram'],p['black'],p['node'])for p in ports]==[(25,[1,2],9),(25,[1,3],9),(25,[2,3],9)]
 g=gram(10,diagrams[25]['edges']);iv=g.inv();assert iv[9,9]==0
 a,b,x,t=s.symbols('a b x t',real=True);low=s.Matrix([[2,-x,-a],[-x,2,-b],[-a,-b,2]])
 coupling=s.zeros(10,3);coupling[9,2]=-t
 correction=coupling.T*iv*coupling;assert correction==s.zeros(3)
 assert s.expand(low.det())==8-2*(x*x+a*a+b*b)-2*a*b*x
 return dict(high_ports=ports,high_inverse_diagonal='0',low_gram=encode(low),schur_correction=encode(correction),
  low_determinant=str(s.expand(low.det())),scope='Every minimal indefinite3core, everyreal black-black parameter and everyreal cross magnitude; the high Lorentzian core contributes one negative direction and its zero Schur correction preserves the low negative direction. Order<=2 cannot have an unpaid join/overlap by marked-radius0.')


def disjoint(contexts,low_diagrams,big,high):
 hp=[p for p in high['census']['ports']if p['allowed']and p['distance']>=5]
 assert len(hp)==18 and min(max(map(F,p['inverse_diagonals']))for p in hp)==F(6,7)
 inv={i:gram(d['order'],d['edges']).inv().applyfunc(s.simplify)for i,d in big.items()if i in{p['diagram']for p in hp}}
 rows=[];types={};groups={};counts=Counter()
 for c in contexts:
  d=low_diagrams[c['order'],c['index']];n=d['order'];bb=c['black'];ww=tuple(i for i in range(n)if i not in bb);code=tuple(d['code'])
  for p in c['ports']:
   if not p['allowed3']:continue
   x=p['node'];lv=max(L.fromcoeff(v)for v in p['inverse_coefficients'])
   wx=L.determinant_tree(n,code,tuple(i for i in ww if i!=x))/L.determinant_tree(n,code,ww)
   for q in hp:
    if p['distance']+q['distance']<8:continue
    e=big[q['diagram']];y=q['node'];vb=q['black'];hv=K(max(map(F,q['inverse_diagonals'])));hww=tuple(i for i in range(e['order'])if i not in vb);hcode=L.fullcode(e['order'],e['edges'])
    wy=L.determinant_tree(e['order'],hcode,tuple(i for i in hww if i!=y))/L.determinant_tree(e['order'],hcode,hww)
    for m,bs in((3,K(1)),(4,K(2))):
     if bs*lv>=2 or bs*hv>=2:continue
     row=dict(low=[n,d['index'],bb,x],high=[q['diagram'],vb,y],cross_label=m,white_inverse_product=(bs*wx*wy).save())
     if bs*wx*wy>=1:row['kind']='white_PD_rejection'
     else:
      key=(n,d['index'],x,q['diagram'],y,m)
      if key not in types:
       h=gram(n,d['edges']);h[x,x]-=L.sym(bs)*inv[q['diagram']][y,y]
       cert=inertia_small(h);cert['high_inverse_diagonal']=str(inv[q['diagram']][y,y]);cert['schur_gram']=encode(h);types[key]=cert
      row['kind']='two_negative'if types[key]['signature'][1]else'survive';row['inertia_type']=list(key)
      key2=(n,d['index'],tuple(bb),q['diagram'],tuple(vb));groups.setdefault(key2,{}).setdefault((x,y),[]).append(m)
     counts[n,m,row['kind']]+=1;rows.append(row)
 print('LOW_DISJOINT',dict(counts),'types',len(types),'groups',len(groups),'multiple',Counter(len(g)for g in groups.values()),flush=True)
 multiple=[]
 for key,choices in sorted(groups.items()):
  if len(choices)<2:continue
  n,index,bb,di,vb=key;d=low_diagrams[n,index];e=big[di];nn,base,mp=union(d,e,());black=list(bb)+[mp[i]for i in vb]
  # Check every nonempty >=2edge set with the least permitted label3.
  # Any larger native product preserves a nonnegative Z obstruction.
  pairs=sorted(choices)
  for size in range(2,len(pairs)+1):
   for chosen in combinations(pairs,size):
    edges=base+[[i,mp[j],3]for i,j in chosen];bad=edge_pd(nn,edges,black)
    multiple.append(dict(low=[n,index,list(bb)],high=[di,list(vb)],cross_edges=list(map(list,chosen)),kind='actual_edge_PD_rejection'if bad else'survive',obstruction=bad))
 assert len(rows)==726 and len(types)==131
 assert counts=={(4,3,'two_negative'):294,(4,3,'white_PD_rejection'):11,(4,4,'two_negative'):17,(4,4,'white_PD_rejection'):6,(5,3,'two_negative'):342,(5,4,'two_negative'):15,(5,3,'white_PD_rejection'):22,(5,4,'white_PD_rejection'):19}
 assert Counter(len(g)for g in groups.values())=={1:249,2:117,3:51}
 assert len(multiple)==321 and all(r['kind']=='actual_edge_PD_rejection'for r in multiple)
 print('LOW_MULTIPLE',len(multiple),Counter(r['kind']for r in multiple),flush=True)
 return dict(high_ports=hp,rows=rows,counts={str(k):v for k,v in counts.items()},types={str(k):v for k,v in types.items()},multiple=multiple)


def overlaps(contexts,low_diagrams,big,high):
 hcontexts={}
 for r in high['census']['ports']:
  if r['distance']>=6:hcontexts.setdefault((r['diagram'],tuple(r['black'])),[]).append(r)
 records=[];count=0
 for c in contexts:
  d=low_diagrams[c['order'],c['index']];ga=L.topology(d['order'],tuple(d['code']));cp={p['node']:p for p in c['ports']}
  for (di,vb),hp in hcontexts.items():
   e=big[di];hcode=L.fullcode(e['order'],e['edges']);gb=L.topology(e['order'],hcode);hq={q['node']:q for q in hp}
   eligible=[(p['node'],q['node'])for p in c['ports']for q in hp if p['distance']+q['distance']>=9]
   for size in range(1,min(len({i for i,j in eligible}),len({j for i,j in eligible}))+1):
    for match in combinations(eligible,size):
     count+=1;ls={i for i,j in match};rs={j for i,j in match}
     if len(ls)!=size or len(rs)!=size:continue
     if any(ga[i].get(k,0)!=gb[j].get(l,0)for i,j in match for k,l in match):continue
     row=dict(low=[d['order'],d['index'],c['black']],high=[di,list(vb)],match=list(map(list,match)))
     bad=[]
     for i,j in match:
      other=[k for k in gb[j]if k not in rs]
      if other and max(L.fromcoeff(v)for v in cp[i]['inverse_coefficients'])>=2:bad.append(dict(side='low',node=i,external=other))
      other=[k for k in ga[i]if k not in ls]
      if other and max(map(F,hq[j]['inverse_diagonals']))>=2:bad.append(dict(side='high',node=j,external=other))
     if bad:row.update(kind='inverse_rejection',rejections=bad)
     else:
      n,edges,mp=union(d,e,match);black=c['black']+[mp[i]for i in vb];bad=edge_pd(n,edges,black)
      if bad:row.update(kind='actual_edge_PD_rejection',obstruction=bad)
      else:row.update(kind='survive',order=n,edges=edges,black=black,right_map=mp)
     records.append(row)
 assert count==630 and len(records)==523 and Counter(r['kind']for r in records)=={'inverse_rejection':422,'survive':66,'actual_edge_PD_rejection':35}
 print('LOW_OVERLAPS',count,len(records),Counter(r['kind']for r in records),flush=True)
 print('SURVIVING_CORES',Counter((r['low'][0],r['low'][1],r['high'][0],tuple(map(tuple,r['match'])))for r in records if r['kind']=='survive'),flush=True)
 return dict(candidate_matches=count,consistent_matches=len(records),records=records)


def extension_scout(overlap,big,ds):
 out=[];corecache={};summary=Counter()
 for record in overlap['records']:
  if record['kind']!='survive':continue
  n=record['order'];edges=record['edges'];code=L.fullcode(n,edges);g=gram(n,edges);u=record['black'][:2];v=record['black'][2:];white=[i for i in range(n)if i not in u+v]
  du,dv=L.dist(g,u),L.dist(g,v);key=str((record['low'][:2],record['high'][0],record['match']))
  if key not in corecache:
   sig=L.inertia(edges,n);ks=g.nullspace();iv=None if ks else g.inv().applyfunc(s.simplify)
   corecache[key]=dict(signature=sig['signature'],gram=encode(g),radical=[list(map(str,k))for k in ks],inverse=encode(iv)if iv is not None else None)
  core=corecache[key];ks=[s.Matrix([s.sympify(x)for x in z])for z in core['radical']];iv=s.Matrix([[s.sympify(x)for x in row]for row in core['inverse']])if core['inverse']is not None else None
  row={k:record[k]for k in('low','high','match','order','black')};row['core_key']=key;row['signature']=core['signature']
  corners=[];badcorner=None
  for p,q in product(u,v):
   omitted=sorted(set(u+v)-{p,q});ids=tuple(i for i in range(n)if i not in omitted);cc=L.subcode(n,code,ids)
   if not L.verdict(len(ids),cc)[0]:badcorner=dict(marks=[p,q],nodes=list(ids),code=list(cc));break
   allowed=omitted.copy();blocks=[];component=None
   h=g.extract(ids,ids)
   for part in L.parts(h):
    actual=[ids[i]for i in part];a=h.extract(part,part);z=a.nullspace()
    if z:
     assert len(z)==1;k=z[0]
     if all(x<0 for x in k):k=-k
     assert all(x>0 for x in k)
     blocks.append(dict(nodes=actual,kind='affine',positive_kernel=list(map(str,k))))
    else:
     assert L.pdtest(a);ji=a.inv().applyfunc(s.simplify);vals=[ji[i,i]for i in range(a.rows)]
     allowed += [actual[i]for i,x in enumerate(vals)if x<=2]
     blocks.append(dict(nodes=actual,kind='finite',inverse_diagonals=list(map(str,vals))))
    if p in actual and q in actual:component=actual
   c=dict(marks=[p,q],omitted=omitted,possible_contacts=sorted(allowed),blocks=blocks,component=component)
   if component is not None:
    if ks:
     z=ks[0];zeros=[i for i in allowed if z[i]==0];pos=[i for i in allowed if z[i]>0];neg=[i for i in allowed if z[i]<0]
     shortcuts=[max([du[i]for i in pos],default=0)+max([dv[i]for i in neg],default=0)+2,max([dv[i]for i in pos],default=0)+max([du[i]for i in neg],default=0)+2]
     c.update(zero_nodes=zeros,radical_shortcut=min(shortcuts),isolated=not zeros and min(shortcuts)<=8)
    else:
     ends=[]
     for radius in range(7):
      left=[i for i in allowed if du[i]<=radius];right=[i for i in allowed if dv[i]<=6-radius]
      if set(left+right)==set(allowed)and not(set(left)&set(right))and all(iv[i,j]>0 for a in(left,right)for i in a for j in a):ends.append([left,right])
     cand=[i for i in allowed if iv[i,i]<=2];c.update(end_partitions=ends,candidates=[dict(node=i,projection=str(iv[i,i]),relevant=i in component)for i in cand])
     if ends:c['minimum_multicontact']=str(min([iv[i,i]+iv[j,j]+2*iv[i,j]for a in ends[0]for i,j in combinations(a,2)],default=s.oo))
   corners.append(c)
  row['corners']=corners
  if badcorner:row.update(kind='adjacent_corner_indefinite',obstruction=badcorner)
  elif core['signature'][1]>1:row['kind']='two_negative'
  else:
   d=ds[record['low'][0],record['low'][1]];e=big[record['high'][0]];match=record['match'];mp=record['right_map'];left=[i for i in range(d['order'])if i not in{a for a,b in match}];right=[mp[i]for i in range(e['order'])if i not in{b for a,b in match}];extras=[];short=0
   for x,y in product(left,right):
    if g[x,y]:continue
    h=g.copy();h[x,y]=h[y,x]=-1
    if min(L.dist(h,u)[j]for j in v)<=8:short+=1;continue
    bad=edge_pd(n,edges+[[x,y,3]],u+v);extras.append(dict(edge=[x,y],kind='actual_edge_PD_rejection'if bad else'survive',obstruction=bad))
   row.update(kind='fixed_core',extra_shortcuts=short,extra_bonds=extras)
  summary[row['kind']]+=1;out.append(row)
  views=[]
  for c in corners:
   if c['component']is None:continue
   if ks:views.append((c['marks'],len(c['component']),c['zero_nodes'],c['radical_shortcut']))
   else:views.append((c['marks'],len(c['component']),bool(c['end_partitions']),c['candidates']))
  print('EXT',record['low'],record['high'],record['match'],core['signature'],row['kind'],views,flush=True)
 assert summary=={'fixed_core':66}and len(corecache)==14
 assert sum(r['extra_shortcuts']for r in out)==2184 and sum(len(r['extra_bonds'])for r in out)==12
 print('EXTENSION_SUMMARY',dict(summary),'cores',len(corecache),flush=True)
 return dict(cases=out,cores=corecache,counts=dict(summary))


def positive_witness(n,code):
 adj=L.topology(n,tuple(code));edges=[[i,j,m]for i in range(n)for j,m in adj[i].items()if i<j];g=gram(n,edges)
 # A shortest cycle has an all-positive norm<=0 vector, even with chords.
 cycles=[]
 for i,j,m in edges:
  todo=[i];parents={i:None}
  for a in todo:
   for b in adj[a]:
    if {a,b}=={i,j}or b in parents:continue
    parents[b]=a;todo.append(b)
   if j in parents:break
  if j in parents:
   p=[j]
   while p[-1]!=i:p.append(parents[p[-1]])
   cycles.append(p)
 if cycles:
  ids=min(cycles,key=lambda z:(len(z),z));v=s.Matrix([int(i in ids)for i in range(n)]);q=s.simplify((v.T*g*v)[0]);assert q<=0
  return dict(method='cycle',support=ids,vector=list(map(str,v)),square=str(q))
 diagonal={i:K(2)for i in range(n)};a={i:dict(adj[i])for i in range(n)};vectors={i:s.eye(n)[:,i]for i in range(n)};steps=[]
 while a:
  i=next(i for i in a if len(a[i])<=1);q=diagonal[i]
  if q<=0:
   v=vectors[i].applyfunc(s.simplify);assert all(x>=0 for x in v)
   assert s.simplify((v.T*g*v)[0]-L.sym(q))==0
   return dict(method='positive_leaf_schur',pivot=i,positive_eliminations=steps,vector=list(map(str,v)),square=str(L.sym(q)))
  if a[i]:
   j=next(iter(a[i]));m=a[i][j];diagonal[j]-=L.BETA[m]/q
   vectors[j]+=TH[m]*vectors[i]/L.sym(q);del a[j][i]
   steps.append(dict(leaf=i,neighbor=j,positive_pivot=q.save(),label=m))
  del a[i];del diagonal[i];del vectors[i]
 raise AssertionError('no nonpositive witness in rejected principal diagram')


def close_extensions(ext):
 rows=[]
 for r in ext['cases']:
  assert r['kind']=='fixed_core'and all(x['kind']=='actual_edge_PD_rejection'for x in r['extra_bonds'])
  core=ext['cores'][r['core_key']];iv=s.Matrix([[s.sympify(x)for x in row]for row in core['inverse']])if core['inverse']is not None else None
  choices=[]
  for c in r['corners']:
   if c['component']is None:continue
   bound=len(c['component'])
   if 'zero_nodes'in c:
    if c['isolated']:choices.append(dict(corner=c['marks'],method='radical_isolation',component_order_bound=bound,shortcut=c['radical_shortcut']))
    continue
   if not c['end_partitions']:continue
   rel=[x for x in c['candidates']if x['relevant']]
   if not rel:
    assert all(iv[i,i]>2 for i in set(c['possible_contacts'])&set(c['component']))
    choices.append(dict(corner=c['marks'],method='no_relevant_contact',component_order_bound=bound,end_sets=c['end_partitions'][0]))
   elif len(c['candidates'])==1:
    node=rel[0]['node'];p=iv[node,node];rho=s.simplify(2-p)
    assert p==(5+s.sqrt(5))/4 and rho==(3-s.sqrt(5))/4 and 0<rho<s.Rational(1,2)
    assert s.sympify(c['minimum_multicontact'])>2
    assert 2*p-2>0 and 2*p>2 and 2*rho<1
    choices.append(dict(corner=c['marks'],method='unique_sealed_label3_root',component_order_bound=bound+1,end_sets=c['end_partitions'][0],
       contact_node=node,projection_square=str(p),residual_square=str(rho),same_node_pair_lower=str(s.simplify(2*p-2)),
       core_orthogonal_contact_square_bound=str(s.simplify(2*rho)),minimum_multicontact=c['minimum_multicontact']))
  assert choices,(r['low'],r['high'],r['match'])
  chosen=min(choices,key=lambda c:(c['component_order_bound'],c['corner']))
  rows.append(dict(low=r['low'],high=r['high'],match=r['match'],core_key=r['core_key'],**chosen))
 assert len(rows)==66 and max(r['component_order_bound']for r in rows)==12
 assert Counter(r['method']for r in rows)=={'radical_isolation':39,'no_relevant_contact':24,'unique_sealed_label3_root':3}
 return dict(cases=rows,component_order_bound=12,method_counts=dict(Counter(r['method']for r in rows)))


def run():
 low=json.loads((ROOT/f'results/{LOW}.json').read_text());high=json.loads((ROOT/f'results/{HIGH}.json').read_text());cat=json.loads((ROOT/f'results/{CENSUS}.json').read_text());big={d['index']:d for d in cat['diagrams']}
 for obj in(low,high,cat):
  assert obj['status'].startswith('PASS')
  for path,h in obj.get('sources',{}).items():assert L.digest(ROOT/path)==h,path
 ds={(d['order'],d['index']):d for order in('4','5')for d in low['classification'][order]['diagrams']};contexts=small_contexts(low)
 tri=triangle_gate(high,big);print('TRIANGLE_ALL_REAL_PARAMETERS_PASS',flush=True)
 joins=disjoint(contexts,ds,big,high);ovs=overlaps(contexts,ds,big,high);ext=extension_scout(ovs,big,ds);closure=close_extensions(ext)
 paths=[Path(__file__),Path(L.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']+[ROOT/f'{folder}/{tag}.{ext}'for tag in(LOW,HIGH,CENSUS)for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
 return dict(status='PASS_COMPLETE_DISTANCE8_SMALL_CORE_DOMAIN',sources={str(p.relative_to(ROOT)):L.digest(p)for p in paths},triangle_all_label_gate=tri,small_contexts=contexts,disjoint=joins,overlaps=ovs,extension_geometry=ext,extension_closure=closure,scope='Complete all-label order<=5 disjoint/overlap domain; every surviving low-first overlap has a controlled actual marked component of order<=12 under the no-path<=8 hypothesis. No dimension claim.')

if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
 if a.verify:assert out==json.loads(p.read_text())
 else:p.write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],flush=True)
