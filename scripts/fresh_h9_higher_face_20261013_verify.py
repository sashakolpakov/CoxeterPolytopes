#!/usr/bin/env python3
"""Exact orthogonal identities and rigorous signs for a linked seven-cap H5 roof."""
import argparse,json,sys
from fractions import Fraction as Q
from itertools import combinations
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from scripts import fresh_h9_remote_frontier_20261013_second_verify as S
V=S.V;K=S.K;N=S.N;C=S.C
PREFIX='data/enumerations/fresh_h9_higher_face_20261013'
COMMON=['seed1','seed2','compact1','compact2','compact3'];FREE='seed1'
T=['seed6','crosspair0','nextH4','compact0'];BASES=['wall17','facePrism']
SPECS=[{'corner':1,'profile':[2]*8,'anchor':'wall17','label':16},
 {'corner':2,'profile':[2]*8,'anchor':'roof0','label':12},
 {'corner':3,'profile':[2]*8,'anchor':'roof0','label':2},
 {'corner':5,'profile':[2]*8,'anchor':'roof0','label':3},
 {'corner':4,'profile':[2]*8,'anchor':'roof1','label':3},
 {'corner':7,'profile':[3]+[2]*7,'anchor':'roof0','label':5},
 {'corner':6,'profile':[3]+[2]*7,'anchor':'roof5','label':4}]
def key(a,b):return '|'.join(sorted((a,b)))
def edge(ci):return[n for i,n in enumerate(T)if i!=ci//2]+[BASES[ci%2]]+[n for n in COMMON if n!=FREE]
def cos2(m):
 if m==2:return V.ZERO
 if m==3:return V.rat(1)
 if m==4:return V.sqrt(V.rat(2))
 if m==5:return V.golden(K.parse({'a':'0','b':'1'}))
 if m==8:return V.sqrt(V.add(V.rat(2),V.sqrt(V.rat(2))))
 if m==12:return V.sqrt(V.add(V.rat(2),V.sqrt(V.rat(3))))
 if m==16:return V.sqrt(V.add(V.rat(2),V.sqrt(V.add(V.rat(2),V.sqrt(V.rat(2))))))
 raise ValueError(m)
def matinv(A):
 n=len(A);M=[list(r)+[V.rat(int(i==j))for j in range(n)]for i,r in enumerate(A)]
 for j in range(n):
  good=[i for i in range(j,n)if M[i][j][0]>0 or M[i][j][1]<0];assert good,('uncertified invertibility',j)
  p=max(good,key=lambda i:abs(sum(M[i][j])));M[j],M[p]=M[p],M[j];z=M[j][j];M[j]=[V.div(x,z)for x in M[j]]
  for i in range(n):
   if i==j:continue
   z=M[i][j];M[i]=[V.sub(x,V.mul(z,y))for x,y in zip(M[i],M[j],strict=True)]
 return[r[n:]for r in M]
def gram(ss,pairs):return[[V.rat(2)if i==j else V.neg(cos2(pairs[key(a,b)]['label']))for j,b in enumerate(ss)]for i,a in enumerate(ss)]
def inertia(ss,pairs):
 A=gram(ss,pairs);piv=[]
 for j in range(len(A)):
  z=A[j][j];piv.append(V.enc(z))
  if z[1]<0:return 'indefinite',piv
  assert z[0]>0,('possible affine or uncertified pivot',ss,j,V.enc(z))
  for a in range(j+1,len(A)):
   for b in range(j+1,len(A)):A[a][b]=V.sub(A[a][b],V.div(V.mul(A[a][j],A[j][b]),z))
 return 'positive_definite',piv

def old_state():
 d=V.read(S.PREFIX+'_caps.json');prior=S.check(d);assert prior==V.read(S.PREFIX+'_verify.json')
 G,I,h,g,roots,names,center,pairs,_=S.old_state();r,_=S.make(d['branch'],g,roots,names,center,pairs)
 for n in names:pairs[key(S.NEW,n)]={'kind':'finite','label':S.PROFILE[n]}if n in S.PROFILE else{'kind':'dotted'}
 return G,I,h,g,roots+[r],names+[S.NEW],center,pairs,prior

def make(spec,g,roots,names,center,pairs):
 ee=edge(spec['corner']);E=[roots[names.index(n)]for n in ee];f=roots[names.index(FREE)];M=C.gram(ee,pairs);Mi=K.inverse(M)
 assert N.T.P.T.psd(M)=='positive_definite'
 assert all(pairs[key(FREE,n)]=={'kind':'finite','label':2}for n in ee)
 b=[K.scale(K.O,-1)if m==3 else K.Z for m in spec['profile']];assert all(m in(2,3)for m in spec['profile'])
 coeff=K.mv(Mi,b);part=N.combine(list(map(V.golden,coeff)),E);qp=K.dot(b,coeff);qr=K.sub(K.scale(K.O,2),qp);assert K.sign(qr)>0
 cp=[V.pair(center,g,r)for r in E];lam=V.mv([[V.golden(x)for x in row]for row in Mi],cp);proj=N.combine(lam,E)
 vf=V.div(V.pair(center,g,f),V.rat(2));v0=[V.sub(V.sub(x,y),V.mul(vf,z))for x,y,z in zip(center,proj,f,strict=True)];qv=V.pair(v0,g,v0);assert qv[1]<0
 v=[V.div(x,V.sqrt(V.neg(qv)))for x in v0];u=[V.div(x,V.sqrt(V.rat(2)))for x in f]
 anchor=roots[names.index(spec['anchor'])];aa=V.neg(V.pair(v,g,anchor))
 if spec['anchor']=='wall17':assert pairs[key(FREE,'wall17')]=={'kind':'finite','label':2};bb=V.ZERO
 else:bb=V.neg(V.pair(u,g,anchor));assert bb[0]>0
 delta=V.sub(V.mul(aa,aa),V.mul(bb,bb));assert aa[0]>0 and delta[0]>0
 z=V.div(V.add(V.pair(part,g,anchor),cos2(spec['label'])),V.sqrt(V.golden(qr)));assert z[0]>=0 or spec['label']==2 and part==[V.ZERO]*10
 if spec['label']==2 and part==[V.ZERO]*10:z=V.ZERO
 root=V.sqrt(V.add(delta,V.mul(z,z)));sh=V.div(V.add(V.mul(aa,z),V.mul(bb,root)),delta);ch=V.div(V.add(V.mul(bb,z),V.mul(aa,root)),delta)
 assert sh[0]>0 and ch[0]>0
 r=[V.add(p,V.mul(V.sqrt(V.golden(qr)),V.sub(V.mul(sh,x),V.mul(ch,y))))for p,x,y in zip(part,v,u,strict=True)]
 # Exact defining identities: part lies in E, v perpendicular E,f and v²=-1,
 # u=f/sqrt2 has u²=1. The displayed sh,ch obey ch²-sh²=1 and aa*sh-bb*ch=z.
 # Therefore r²=qp+qr=2, prescribed E pairings hold, and r.anchor=-2cos(pi/m).
 return r,{'edge':ee,'profile':spec['profile'],'anchor':spec['anchor'],'label':spec['label'],'part_norm':K.enc(qp),'residual_norm':K.enc(qr),'timelike_projection_norm':V.enc(qv),'aa':V.enc(aa),'bb':V.enc(bb),'delta':V.enc(delta),'z':V.enc(z),'sinh':V.enc(sh),'cosh':V.enc(ch)}

def vertex(ss,roots,names,center,g,h,pairs):
 typ,piv=inertia(ss,pairs);assert typ=='positive_definite';M=gram(ss,pairs);Mi=matinv(M)
 assert all(x[0]>=0 or x[0]<=0<=x[1]and max(abs(x[0]),abs(x[1]))<10**12 for row in Mi for x in row)
 # Positive definiteness plus nonpositive off-diagonal Gram gives a nonnegative
 # exact inverse. Interval checks below certify strict signs of the actual point.
 E=[roots[names.index(n)]for n in ss];lam=V.mv(Mi,[V.pair(center,g,r)for r in E]);assert all(x[1]<0 for x in lam)
 proj=N.combine(lam,E);p=[V.sub(x,y)for x,y in zip(center,proj,strict=True)];q=V.pair(p,g,p);assert q[1]<0 and V.dot(p,list(map(V.golden,h)))[0]>0
 checks=[{'wall':n,'pairing':V.enc(V.pair(p,g,r))}for n,r in zip(names,roots,strict=True)if n not in ss];assert all(int(x['pairing'][1])<0 for x in checks)
 return{'incident':ss,'pivots':piv,'point':list(map(V.enc,p)),'norm':V.enc(q),'nonincident':checks}

def build():
 G,I,h,g,roots,names,center,pairs,prior=old_state();recipes=[];changed=[];affine=[];allnewvertices=[]
 for i,spec in enumerate(SPECS):
  new=f'roof{i}';r,recipe=make(spec,g,roots,names,center,pairs);ff=dict(zip(edge(spec['corner']),spec['profile'],strict=True));ff[spec['anchor']]=spec['label'];checks=[]
  assert V.dot(r,list(map(V.golden,h)))[0]>0
  for n,a in zip(names,roots,strict=True):
   p=V.pair(r,g,a)
   if n in ff:
    want=V.neg(cos2(ff[n]));assert p[0]<=want[1]and want[0]<=p[1],(new,n,p,want)
    pairs[key(new,n)]={'kind':'finite','label':ff[n]}
   else:assert p[1]<-2*V.SCALE,(new,n,p);pairs[key(new,n)]={'kind':'dotted'}
   checks.append({'wall':n,'kind':pairs[key(new,n)],'pairing':V.enc(p)})
  roots.append(r);names.append(new);recipes.append(recipe);changed.append(checks);stat={'new':new,'subsets':2**len(ff),'dotted_excluded':0,'positive_definite':0,'indefinite':0,'affine':0};positive=[];nn=list(ff)
  for mask in range(2**len(nn)):
   ss=[new]+[n for j,n in enumerate(nn)if mask>>j&1]
   if any(pairs[key(a,b)]['kind']=='dotted'for a,b in combinations(ss,2)):stat['dotted_excluded']+=1;continue
   typ,_=inertia(ss,pairs);stat[typ]+=1
   if typ=='positive_definite'and len(ss)==9:positive.append(ss)
  affine.append(stat);allnewvertices.extend(positive);print(new,stat,flush=True)
 # Recheck all new finite vertices against the final family, and all eight old
 # vertical edge endpoints. No mere Gram-pair intersection inference is used.
 vertices=[vertex(ss,roots,names,center,g,h,pairs)for ss in allnewvertices]
 cover=[]
 for ci in range(8):
  choices=[ss for ss in allnewvertices if set(edge(ci))<=set(ss)];assert choices,ci;cover.append({'corner':ci,'endpoints':choices})
 assert prior['affine_census']['affine']==prior['inherited_affine_circuits']==0
 out={'status':'PASS_EXACT_SEVEN_CAP_LINKED_ROOF','sources':{S.PREFIX+'_caps.json':V.sha(S.PREFIX+'_caps.json'),S.PREFIX+'_verify.json':V.sha(S.PREFIX+'_verify.json'),'scripts/'+Path(__file__).name:V.sha('scripts/'+Path(__file__).name),'scripts/'+Path(S.__file__).name:V.sha('scripts/'+Path(S.__file__).name),'scripts/'+Path(V.__file__).name:V.sha('scripts/'+Path(V.__file__).name)},'facets':len(roots),'new_norms':[2]*7,'recipes':recipes,'new_pairs':changed,'new_affine_census':affine,'old_affine_circuits':0,'actual_new_vertices':vertices,'old_vertical_edge_coverage':cover,'fixed_H5_walls':[n for n in COMMON if n!=FREE],'normals':[list(map(V.enc,r))for r in roots[-7:]],'scope':'All seven caps retain the prior106 actual facets and close each of the eight original prism vertical edges. This is not by itself compactness or finite volume.'}
 return out,(G,I,h,g,roots,names,center,pairs)

def run(verify=False):
 out,_=build()
 if verify:assert out==V.read(PREFIX+'_verify.json')
 else:(V.ROOT/(PREFIX+'_verify.json')).write_text(json.dumps(out,indent=2)+'\n')
 print(out['status'],'facets',out['facets'],'vertices',len(out['actual_new_vertices']),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();run(a.verify)
