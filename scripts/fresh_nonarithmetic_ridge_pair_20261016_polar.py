#!/usr/bin/env python3
"""Numerical discovery on every positive-E extreme direction of the 52-wall cone."""
import argparse,json
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import mpmath as mp
import fresh_nonarithmetic_escape49_20261016_polar as P
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_nonarithmetic_ridge_pair_20261016_polar'
CL='results/fresh_nonarithmetic_affine_filtered_20261016_cluster.json'
CO='results/fresh_nonarithmetic_affine_filtered_20261016_cone.json'
DOM='results/fresh_nonarithmetic_escape_cutter_domain_20261016.json'
END='results/fresh_nonarithmetic_ridge_pair_20261016_endpoints.json'
mp.mp.dps=80
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def mid(a):
 q=(Fraction(a[0])+Fraction(a[1]))/2;return mp.mpf(q.numerator)/q.denominator
def geometry():
 old=P.B.read(P.BATCH);cl=P.B.read(CL);co=P.B.read(CO)
 Q=mp.matrix([[P.scalar(P.D.encfield(P.D.lift(P.A.B.db(x))))for x in row]for row in old['form']])
 N=[mp.matrix([P.scalar(x)for x in row])for row in old['normals']]+[mp.matrix([mid(x)for x in c['normal_interval']])for c in cl['caps']]
 N=[n*mp.sqrt(2/(n.T*Q*n)[0])for n in N]
 H=mp.matrix([P.scalar(x)for x in old['center']]);E=mp.matrix([P.scalar(x)for x in co['strict_escape']['point']]);E/=-(E.T*Q*H)[0]
 H/=mp.sqrt(-(H.T*Q*H)[0]);return Q,N,H,E
def relation(g):
 if abs(g)<mp.mpf('1e-30'):return dict(kind='finite',label=2)
 if g<-2-mp.mpf('1e-30'):return dict(kind='dotted')
 if abs(g+2)<mp.mpf('1e-30'):return dict(kind='unresolved_parallel')
 if g>mp.mpf('1e-30'):return dict(kind='obtuse')
 m=mp.pi/mp.acos(-g/2);n=int(mp.nint(m))
 if n>=2 and abs(m-n)<mp.mpf('1e-27'):return dict(kind='finite',label=n)
 return dict(kind='finite_gap_proposal',effective_label=mp.nstr(m,45))
def run():
 Q,N,H,E=geometry();dom=P.B.read(DOM);ends=P.B.read(END)
 cov=[Q*n for n in N];ecov=Q*E;hcov=Q*H;stats=Counter();records=[];proposals=[]
 endpoints=[]
 for i,x in enumerate(ends['candidates']):
  u=mp.matrix([mp.mpf(z)for z in x['point']])
  if not any(max(abs(a-b)for a,b in zip(u,v))<mp.mpf('1e-25')for _,v in endpoints):endpoints.append((i,u))
 for saved in dom['records']:
  if saved['on_cutting_boundary']:continue
  ids=saved['active'];v=mp.matrix(P.null([list(cov[j])for j in ids]));v/=-(v.T*hcov)[0]
  q=(v.T*Q*v)[0];assert q>0;v*=mp.sqrt(2/q);ep=(v.T*ecov)[0];assert ep>0
  vals=[(v.T*c)[0]for c in cov];assert max(vals)<mp.mpf('1e-29')
  rels=[dict(wall=j,pairing=mp.nstr(g,45),**relation(g))for j,g in enumerate(vals)]
  bad=[r for r in rels if r['kind']not in['finite','dotted']]
  row=dict(active=ids,escape_pairing=mp.nstr(ep,45),center_pairing=mp.nstr((v.T*hcov)[0],45),relations=rels,bad=bad)
  stats['cutting_polars']+=1;stats['bad_contacts']+=len(bad)
  if not bad:
   row['point']=[mp.nstr(x,65)for x in v];row['endpoint_compatibility']=[]
   for i,u in endpoints:
    g=(v.T*Q*u)[0];r=relation(g)
    if r['kind']in['finite','dotted']:row['endpoint_compatibility'].append(dict(endpoint_index=i,pairing=mp.nstr(g,45),**r))
   proposals.append(row);print('PROPOSAL',ids,'E',row['escape_pairing'],'compatible_endpoints',len(row['endpoint_compatibility']),flush=True)
  records.append(row)
  if len(records)%20==0:print('screened',len(records),'proposals',len(proposals),flush=True)
 assert len(records)==146
 paths=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',P.BATCH,CL,CO,DOM,END]
 return dict(status='NUMERICAL_ALL_146_ESCAPE_CUTTING_POLARS_OF_FILTERED52',sources={p:sha(p)for p in paths},stats=dict(stats),endpoint_count=len(endpoints),endpoint_representatives=[i for i,_ in endpoints],proposals=proposals,records=records,scope='80-digit discovery from certified input enclosures, not exact polar rejection or acceptance. All146 extreme rays with positive renewed-E pairing are screened. Other candidate normals and boundary combinations remain outside this finite polar screen. The noncutting endpoint bank is numerical and bounded; all candidate pairs require exact certification.')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
 if a.verify:assert d==json.loads(p.read_text());print('PASS saved numerical cutting-polar replay')
 else:p.write_text(json.dumps(d,indent=2)+'\n')
 print(d['status'],d['stats'],'proposals',len(d['proposals']))
