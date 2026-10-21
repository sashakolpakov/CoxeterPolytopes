#!/usr/bin/env python3
"""Bounded mutual-label, complete fixed-profile three-endpoint belt diagnostic.

No new wall or global exclusion is claimed by numerical equation solving.
"""
import argparse,hashlib,itertools,json,math,time
from functools import lru_cache
import networkx as nx
import numpy as np
from scipy.optimize import least_squares
import fresh_h20_e8_partner_endpoints_20261017 as P
ROOT=P.ROOT;TAG='fresh_h20_e8_endpoint_belt_20261017';INPUT='results/'+P.TAG+'.json'
LABELS=(2,3,4,5,6,8,10)
def tree_key(g):
 def rec(v,p):return tuple(sorted((g[v][w]['label'],rec(w,v))for w in g[v]if w!=p))
 return min(rec(v,None)for v in nx.center(g))
MODELS={}
for name,g in P.models():MODELS.setdefault(len(g),set()).add(tree_key(g))
@lru_cache(None)
def finite(n,edges):
 g=P.graph(n,edges)
 for ids in nx.connected_components(g):
  h=g.subgraph(ids)
  if len(h)==1:continue
  if len(h.edges)!=len(h)-1:return False
  if len(h)==2:continue
  if tree_key(h)not in MODELS.get(len(h),set()):return False
 return True

def old_edges():
 out=[]
 for path in[list(range(1,11)),[12,13,14,15],[61,87]]:
  out.extend((i,j,3)for i,j in zip(path,path[1:]))
 out.extend([(1,159,3),(10,160,3),(159,161,4),(160,161,3)])
 return out

def corner(T,old,rows,bonds):
 ids=T+old+list(range(200,200+len(rows)));ix={v:i for i,v in enumerate(ids)}
 ee=[(ix[i],ix[j],m)for i,j,m in old_edges()if i in ix and j in ix]
 for k,row in enumerate(rows):
  ee.extend((ix[j],ix[200+k],m)for j,m in row['contacts']if j in ix)
 ee.extend((ix[200+i],ix[200+j],m)for i,j,m in bonds if m>2)
 return finite(len(ids),tuple(sorted((min(i,j),max(i,j),m)for i,j,m in ee)))

def run(solve=True):
 D=P.load(INPUT);cc=D['endpoint_censuses'];T=[1,2,3,4,5,6,7,8,9,10,12,13,14,15,61,87,151]
 allowed={};counts={}
 for i,j in[(0,1),(0,2),(1,2)]:
  shared=list(set(cc[i]['edge'])&set(cc[j]['edge']));one={}
  for a in cc[i]['profiles']:
   for b in cc[j]['profiles']:
    labs=[m for m in LABELS if corner(T,shared,[a,b],[(0,1,m)])]
    if labs:one[a['id'],b['id']]=labs
  allowed[i,j]=one;counts[f'{i},{j}']=dict(profile_pairs=len(one),labelled_pairs=sum(map(len,one.values())))
  print('PAIR',i,j,counts[f'{i},{j}'],flush=True)
 patterns=[]
 for a,b in allowed[0,1]:
  for c in range(len(cc[2]['profiles'])):
   if(a,c)not in allowed[0,2]or(b,c)not in allowed[1,2]:continue
   rows=[cc[0]['profiles'][a],cc[1]['profiles'][b],cc[2]['profiles'][c]]
   for lm in itertools.product(allowed[0,1][a,b],allowed[0,2][a,c],allowed[1,2][b,c]):
    if corner(T,[],rows,[(0,1,lm[0]),(0,2,lm[1]),(1,2,lm[2])]):patterns.append(dict(profiles=[a,b,c],mutual_labels=list(lm)))
 print('FINITE_BELT_PATTERNS',len(patterns),flush=True)
 src=[f'scripts/{TAG}.py',INPUT,f'scripts/{P.TAG}.py']
 out=dict(status='EXACT_FINITE_PROFILE_BELT_PATTERNS_WITH_BOUNDED_MUTUAL_LABELS',sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in src},mutual_label_domain=list(LABELS),fixed_profile_counts=[len(c['profiles'])for c in cc],pair_counts=counts,finite_patterns=patterns,scope='Complete combination of all159 fixed endpoint profiles with mutual labels2,3,4,5,6,8,10. Each shared-old corner and the three-new-root corner is literally spherical by finite-tree classification. The unbounded singleton151 endpoint families and other mutual I2 labels are not exhausted. This is a proposed seven-facet truncated-tetrahedron incidence pattern, not a complete class of H3 closures.')
 if not solve:return out
 G,Y,CT,J,q,H,live,verts=P.geometry();m=G['m'];families=[]
 for ci,c in enumerate(cc):
  z={}
  for row in c['profiles']:
   good=[br for br in c['families'][row['id']]['branches']if br['all_unprescribed_dotted_intervals']]
   assert len(good)==1
   br=good[0];assert len(br['all_unprescribed_dotted_intervals'])==1
   lo,hi=br['all_unprescribed_dotted_intervals'][0];assert hi is None
   z[row['id']]=dict(lo=lo,p=np.array(br['particular']),a=np.array(br['positive_k_coefficient']),b=np.array(br['inverse_k_coefficient']),dt=np.array(br['T_pairings']))
  families.append(z)
 solutions=[];attempts=0;best=None
 for pi,pat in enumerate(patterns):
  ff=[families[i][v]for i,v in enumerate(pat['profiles'])];off=np.array([[x['dt']@J@y['dt']for y in ff]for x in ff]);target=np.array([-2*np.cos(np.pi/l)if l>2 else 0 for l in pat['mutual_labels']])
  def rs(x):return np.array([f['p']+f['a']*(f['lo']*np.exp(t))+f['b']/(f['lo']*np.exp(t))for f,t in zip(ff,x)])
  def fun(x):
   rr=rs(x);gg=rr@m@rr.T+off
   return np.array([gg[0,1],gg[0,2],gg[1,2]])-target
  # Numerical discovery only: bounded log parameter and finitely many starts.
  for start in[(0,0,0),(.2,.2,.2),(1,1,1),(2,2,2),(.1,1,2),(2,.1,1),(1,2,.1)]:
   fit=least_squares(fun,np.array(start,float),bounds=(np.zeros(3),np.full(3,5.)),max_nfev=120,ftol=1e-10,xtol=1e-10,gtol=1e-10);attempts+=1
   err=float(max(abs(fun(fit.x))))
   if best is None or err<best['error']:best=dict(pattern=pi,error=err,log_parameters=fit.x.tolist())
   if err<1e-6:
    rr=rs(fit.x);record=dict(pattern=pi,profiles=pat['profiles'],mutual_labels=pat['mutual_labels'],log_parameters=fit.x.tolist(),parameters=[float(f['lo']*np.exp(t))for f,t in zip(ff,fit.x)],residual4=rr.tolist(),T_pairings=[f['dt'].tolist()for f in ff],error=err)
    if not any(s['pattern']==pi and max(abs(np.array(s['log_parameters'])-fit.x))<1e-5 for s in solutions):solutions.append(record);print('SOLUTION',record,flush=True)
    break
  if pi%100==0:print('SOLVED',pi,'of',len(patterns),'best',best,'solutions',len(solutions),flush=True)
 out.update(status='BOUNDED_NUMERICAL_THREE_ENDPOINT_BELT_SOLVER',numerical_solver=dict(log_parameter_bounds=[0,5],starts=7,max_nfev=120,attempts=attempts,best=best),solutions=solutions)
 return out
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');ap.add_argument('--profiles-only',action='store_true');a=ap.parse_args();d=run(not a.profiles_only);out=ROOT/f'results/{TAG}{"_profiles"if a.profiles_only else ""}.json'
 if a.verify:
  from fresh_h20_triangle_adjacent_h3_20261016 import close
  assert close(d,json.loads(out.read_text()))
 else:out.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],flush=True)
