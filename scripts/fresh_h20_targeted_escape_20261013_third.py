#!/usr/bin/env python3
"""Exact18-short two-angle cap cutting the common numerical vertexbank gap."""
import argparse,json,hashlib
from pathlib import Path
from fractions import Fraction as F
from collections import Counter
from fresh_dimension_bridge_20261013_verify import E,Z,O,parse,vec,dot,mv,inverse
from fresh_dimension_bridge_two_angle_20261013_independent import B,Extension,fdot,fmv,same,label
from fresh_h20_targeted_escape_20261013_vertex_bank import positive_solve

ROOT=Path(__file__).resolve().parents[1]
BASE='data/enumerations/fresh_dimension_bridge_20261013.json'
POLAR='data/enumerations/fresh_dimension_bridge_polar_20261013.json'
EXT='data/enumerations/fresh_dimension_bridge_two_angle_20261013_exact.json'
FIRST='data/enumerations/fresh_h20_targeted_escape_20261013_exact.json'
SECOND='data/enumerations/fresh_h20_targeted_escape_20261013_second.json'
SCOUT='data/enumerations/fresh_h20_targeted_escape_20261013_edge_directions.json'
TARGET='data/enumerations/fresh_h20_targeted_escape_20261013_bank_geometry.json'
OUT='data/enumerations/fresh_h20_targeted_escape_20261013_third.json'
def read(p):return json.loads((ROOT/p).read_text())
def sha(p):return hashlib.sha256((ROOT/p).read_bytes()).hexdigest()
def line(A,b):
 a=[list(r)+[x]for r,x in zip(A,b,strict=True)];piv=[];row=0
 for j in range(21):
  k=next((k for k in range(row,20)if a[k][j].sign()),None)
  if k is None:continue
  a[row],a[k]=a[k],a[row];d=a[row][j];a[row]=[x/d for x in a[row]]
  for k in range(20):
   if k==row:continue
   d=a[k][j]
   if d.sign():a[k]=[x-d*y for x,y in zip(a[k],a[row],strict=True)]
  piv.append(j);row+=1
  if row==20:break
 assert row==20;free,=set(range(21))-set(piv);u=[Z]*21;v=[Z]*21;v[free]=O
 for i,j in enumerate(piv):u[j]=a[i][-1];v[j]=-a[i][free]
 return tuple(u),tuple(v)

def run():
 base,polar,ext,first,second,scout,geometry=map(read,(BASE,POLAR,EXT,FIRST,SECOND,SCOUT,TARGET))
 source=scout['caps'][1];assert source['anchors']==[8,42]and source['labels']==[3,3]
 S=source['short_orthogonal_support'];assert len(S)==18 and 3 in S and 8 not in S
 Q=tuple(map(vec,base['form']));h=vec(base['strict_center']);hc=mv(Q,h)
 old=list(map(vec,base['normals']))+[vec(c['vector'])for c in polar['cuts']]
 cov=[mv(Q,r)for r in old];ns=[dot(r,c)for r,c in zip(old,cov,strict=True)]
 extras=[{'field':c['field'],'normal':c['vector'],'norm':[c['norm'],['0','0']]}for c in ext['caps']]+[
  {'field':first['field'],'normal':first['normal'],'norm':first['norm']},
  {'field':second['field'],'normal':second['normal'],'norm':second['norm']}]
 ids=S+[8,42];rhs=[Z]*18+[E(F(-1)),E(F(-3))]
 u,v=line([cov[i]for i in ids],rhs);aa=dot(v,mv(Q,v));bb=E(F(2))*dot(u,mv(Q,v));cc=dot(u,mv(Q,u))-E(F(2))
 dd=bb*bb-E(F(4))*aa*cc;assert aa.sign()and dd.sign()>0
 f=Extension(B,dd);Qf=tuple(tuple(f.lift(x)for x in row)for row in Q)
 raw=tuple(E(F(str(x)))for x in geometry['whole_bank_geometry']['point'])
 accepted=[]
 for branch in(-1,1):
  t=(-bb/(E(F(2))*aa),E(F(branch))/(E(F(2))*aa))
  r=tuple(f.add(f.lift(x),f.mul(t,f.lift(y)))for x,y in zip(u,v,strict=True));rc=fmv(f,Qf,r);q=fdot(f,r,rc)
  assert same(f,q,f.rational(2))
  center=fdot(f,r,tuple(map(f.lift,hc)))
  if f.sign(center)>=0:continue
  pairs=[];okay=True
  try:
   for i,(cv,n)in enumerate(zip(cov,ns,strict=True)):
    p=fdot(f,r,tuple(map(f.lift,cv)));pairs.append({'old':i,'pairing':f.encode(p),'kind':label(f,p,q,f.lift(n))})
   for i,c in enumerate(extras,82):
    g=Extension(f,f.lift(parse(c['field']['radicand'])))
    other=tuple((f.lift(parse(x[0])),f.lift(parse(x[1])))for x in c['normal'])
    oq=tuple(f.lift(parse(x))for x in c['norm'])
    p=fdot(g,tuple(map(g.lift,rc)),other);kind=label(g,p,g.lift(q),oq)
    pairs.append({'old':i,'pairing':g.encode(p),'kind':kind})
  except AssertionError:continue
  target=None
  for power in range(12,61):
   y=tuple(x+E(F(1,2**power))*z for x,z in zip(raw,h,strict=True));yc=mv(Q,y);yq=dot(y,yc);yh=dot(y,hc)
   if yq.sign()<=0 or yh.sign()>=0:continue
   ps=[dot(y,c)for c in cov]
   if any(p.sign()>=0 for p in ps):continue
   ep=[]
   for c in extras:
    g=Extension(B,parse(c['field']['radicand']));rr=tuple(map(g.parse,c['normal']))
    p=fdot(g,rr,tuple(map(g.lift,yc)))
    if g.sign(p)>=0:break
    ep.append(g.encode(p))
   if len(ep)!=5:continue
   cut=fdot(f,r,tuple(map(f.lift,yc)))
   if f.sign(cut)<=0:continue
   target={'vector':[x.enc()for x in y],'norm':yq.enc(),'center_pairing':yh.enc(),
    'old_wall_pairings':[p.enc()for p in ps]+ep,'new_cap_pairing':f.encode(cut),'perturbation_power':power};break
  if target:accepted.append((branch,r,rc,q,center,pairs,target))
 assert len(accepted)==1
 branch,r,rc,q,center,pairs,target=accepted[0]
 # Actual finite vertex at18short walls S, oldwall42, andnewcap.
 J=S+[42];M=[[dot(old[i],cov[j]).a for j in J]+[F(-3)if i==42 else F(0)]for i in J]
 M.append([F(-3)if i==42 else F(0)for i in J]+[F(2)])
 assert positive_solve(M,[F(0)]*20)is not None
 inv=inverse(M);incident=[tuple(map(f.lift,old[i]))for i in J]+[r]
 b=[fdot(f,x,tuple(map(f.lift,hc)))for x in incident]
 lam=[fdot(f,tuple(f.rational(a)for a in row),b)for row in inv]
 x=tuple(f.sub(f.lift(h[i]),fdot(f,lam,[a[i]for a in incident]))for i in range(21));xc=fmv(f,Qf,x);xq=fdot(f,x,xc)
 assert f.sign(xq)<0 and f.sign(fdot(f,x,tuple(map(f.lift,hc))))<0
 vertex_pairs=[]
 for i,cv in enumerate(cov):
  p=fdot(f,x,tuple(map(f.lift,cv)));assert f.sign(p)==(0 if i in J else -1);vertex_pairs.append(f.encode(p))
 for c in extras:
  g=Extension(f,f.lift(parse(c['field']['radicand'])))
  other=tuple((f.lift(parse(a[0])),f.lift(parse(a[1])))for a in c['normal'])
  p=fdot(g,tuple(map(g.lift,xc)),other);assert g.sign(p)<0;vertex_pairs.append(g.encode(p))
 assert f.sign(fdot(f,r,xc))==0
 return {'status':'EXACT_88_WALL_H20_18_SHORT_TARGET_CUT_VOLUME_PENDING',
  'source_sha256':{p:sha(p)for p in(BASE,POLAR,EXT,FIRST,SECOND,SCOUT,TARGET)},
  'scout_index':1,'short_orthogonal_support':S,'anchor_indices':[8,42],'anchor_labels':[3,3],
  'field':{'kind':'quadratic','base':{'kind':'sqrt6'},'radicand':dd.enc()},
  'normal':list(map(f.encode,r)),'norm':[['2','0'],['0','0']],'quadratic_branch':branch,
  'strict_center_pairing':f.encode(center),'old_pair_checks':pairs,
  'pair_counts':dict(sorted(Counter(p['kind']for p in pairs).items())),
  'strict_target':target,'actual_finite_vertex':{'vector':list(map(f.encode,x)),'norm':f.encode(xq),
   'incident_walls':J+[87],'old_pairings':vertex_pairs,'gram':[[str(a)for a in row]for row in M]},
  'scope':'All87pairrelations exact, norm2 andstrictcenter, actual20wallfinitevertex, andcutof anexactstrict87wallpointnear thecommon1593bank numericalescape. The1593bank itself remainsnumeric; thiscertificate doesnot claimits exactuniversalescape. Volume88mustbecheckedseparately.'}

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();out=run()
 if a.verify:assert out==read(OUT)
 else:(ROOT/OUT).write_text(json.dumps(out,indent=2)+'\n')
 print(json.dumps({k:out[k]for k in ['status','field','pair_counts']}))
