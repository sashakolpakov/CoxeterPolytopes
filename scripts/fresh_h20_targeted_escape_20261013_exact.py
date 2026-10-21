#!/usr/bin/env python3
"""Exact new86th H20 facet that cuts the saved strict85-wall escape."""
import argparse
import copy
from collections import Counter
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label,rank

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
EXT='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'

def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()

def run():
 b,p,e=read(BASE),read(POLAR),read(EXT)
 assert p['source']==BASE and p['source_sha256']==sha(BASE)
 assert e['base_source']==POLAR and e['base_sha256']==sha(POLAR)
 Q=tuple(map(vec,b['form']));old=list(map(vec,b['normals']))+[vec(c['vector'])for c in p['cuts']]
 cov=[mv(Q,r)for r in old];ns=[dot(r,c)for r,c in zip(old,cov,strict=True)]
 h=vec(b['strict_center']);hc=mv(Q,h);w=vec(p['strict_escape']['vector']);wc=mv(Q,w)
 v=vec(p['rejected_polar_proposals'][27]['vector'])
 assert v[-1]==-O and dot(v,mv(Q,v))==O
 f=Extension(B,E(F(10)));r=tuple(f.lift(x)for x in v[:-1])+((E(F(-4)),O),)
 Qf=tuple(tuple(f.lift(x)for x in row)for row in Q);rc=fmv(f,Qf,r)
 q=fdot(f,r,rc);assert same(f,q,(E(F(76)),E(F(-24))))and f.sign(q)>0
 center=fdot(f,r,tuple(map(f.lift,hc)));cut=fdot(f,r,tuple(map(f.lift,wc)))
 assert same(f,center,f.rational(-110))and same(f,cut,(E(F(1225,128)),E(F(-3))))
 assert f.sign(center)<0 and f.sign(cut)>0
 pairs=[]
 for i,(c,n)in enumerate(zip(cov,ns,strict=True)):
  pp=fdot(f,r,tuple(map(f.lift,c)));kind=label(f,pp,q,f.lift(n))
  pairs.append({'old':i,'field':'Q(sqrt6,sqrt10)','pairing':f.encode(pp),'kind':kind})
 for i,c in enumerate(e['caps'],82):
  assert c['field']['kind']=='quadratic'and c['field']['base']=={'kind':'sqrt6'}
  g=Extension(f,f.lift(parse(c['field']['radicand'])))
  rr=tuple((f.lift(parse(x[0])),f.lift(parse(x[1])))for x in c['vector'])
  pp=fdot(g,tuple(map(g.lift,rc)),rr)
  kind=label(g,pp,g.lift(q),g.lift(f.lift(parse(c['norm']))))
  assert kind=='dotted'
  pairs.append({'old':i,'field':{'base':'Q(sqrt6,sqrt10)','radicand':c['field']['radicand']},
    'pairing':g.encode(pp),'kind':kind})
 assert len(pairs)==85
 right=[a['old']for a in pairs if a['kind']=='2'];assert len(right)==19
 assert rank([cov[i]for i in right])==19
 assert [a['old']for a in pairs if a['kind']=='4']==[26,34,44]
 t=(E(F(4)),-O)
 factor=f.mul(f.rational(3),t)
 vertex=tuple(f.mul(factor,f.lift(x))for x in v[:-1])+(f.rational(-2),)
 vertex_cov=fmv(f,Qf,vertex);vertex_q=fdot(f,vertex,vertex_cov)
 assert same(f,vertex_q,f.mul(f.rational(-6),q))and f.sign(vertex_q)<0
 assert f.sign(fdot(f,vertex,tuple(map(f.lift,hc))))<0
 vp=[]
 for i,c in enumerate(cov):
  pp=fdot(f,vertex,tuple(map(f.lift,c)))
  assert f.sign(pp)==(0 if i in right else -1)
  vp.append(f.encode(pp))
 for c in e['caps']:
  g=Extension(f,f.lift(parse(c['field']['radicand'])))
  rr=tuple((f.lift(parse(x[0])),f.lift(parse(x[1])))for x in c['vector'])
  pp=fdot(g,tuple(map(g.lift,vertex_cov)),rr);assert g.sign(pp)<0;vp.append(g.encode(pp))
 assert f.sign(fdot(f,r,vertex_cov))==0
 cusps=[]
 for c in b['cusp_rows']:
  pp=fdot(f,r,tuple(map(f.lift,mv(Q,vec(c['ray'])))))
  assert f.sign(pp)<0;cusps.append(f.encode(pp))
 return {'status':'EXACT_86_WALL_H20_TARGETED_CUT_VOLUME_PENDING',
  'source_sha256':{path:sha(path)for path in(BASE,POLAR,EXT)},
  'construction':'Copy first20 coordinates of rejected polar27; replace last coordinate−1 by−4+sqrt10.',
  'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':['10','0']},
  'normal':list(map(f.encode,r)),'norm':f.encode(q),'strict_center_pairing':f.encode(center),
  'old_escape_pairing':f.encode(cut),'old_escape_norm':p['strict_escape']['norm'],
  'old_pair_checks':pairs,'pair_counts':dict(sorted(Counter(a['kind']for a in pairs).items())),
  'actual_finite_vertex':{'vector':list(map(f.encode,vertex)),'norm':f.encode(vertex_q),
   'old_pairings':vp,'incident_walls':right+[85],
   'proof':'The19 old short normals have independent positive-definite span because their orthogonal complement inside the oldH19 form contains timelike v_old. The new normal is orthogonal to all19 and spacelike, giving rank20 positive definite; the displayed point is timelike and strict at every other wall.'},
  'inherited_cusp_strict_pairings':cusps,
  'scope':'Exact positive new norm, all85 old pair gates, strict common center, positive cut of archived strict escape, and all107 inherited cusp rays strictly retained. All86 walls are actual facets by acute projection. Volume must be checked separately.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run()
 if args.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k]for k in ['status','norm','strict_center_pairing','old_escape_pairing','pair_counts']}))
