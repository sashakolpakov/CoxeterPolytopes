#!/usr/bin/env python3
"""Proposal scout for fork-leaf cusp sections; no numerical nonexistence claim."""
from pathlib import Path
from itertools import combinations
from fractions import Fraction as F
import json,sys
import numpy as np
import sympy as sp
from scipy.optimize import linprog
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
import fresh_h20_global_fraction_verify_20261013 as S
Q,W,G,h=S.source();Q=np.array(Q,dtype=np.int64);W=np.array(W,dtype=np.int64);G=np.array(G,dtype=np.int64);h=np.array(h,dtype=np.int64);C=W@Q
P='fresh_h20_cusp_sections_20261013';OUT=ROOT/f'results/{P}_scout.json'

def exact(v):return [sp.Rational(int(x))for x in v]
def vals(v):return list(map(str,v))
def point_at(e,zero):
    # Linear max-slack candidate, then use cusp drift where necessary.
    A=np.column_stack([C,np.ones(210)]);A[zero,-1]=0
    eq=np.vstack([np.r_[e@Q,0],np.r_[h@Q,0]]+[np.r_[C[i],0]for i in zero])
    b=np.array([0,-1]+[0]*len(zero))
    lp=linprog(np.r_[np.zeros(22),-1],A_ub=A,b_ub=np.zeros(210),A_eq=eq,b_eq=b,bounds=[(None,None)]*22+[(0,1)],method='highs')
    if not lp.success:return None
    x=lp.x[:22]
    return {'point':x.tolist(),'slack':float(lp.x[-1]),'norm':float(x@Q@x)}

def scan():
    walk=json.loads((ROOT/'results/fresh_h21_source_corner_walk_20261013.json').read_text());rows=[]
    for cusp in walk['corner_orbits']:
      if cusp['type']!='ideal':continue
      u=np.array(cusp['primitive_ray'],dtype=np.int64);ids=cusp['active_walls']
      for a,b in combinations(ids,2):
       if G[a,a]!=G[b,b] or G[a,b]!=0:continue
       if any(G[a,i]!=G[b,i]for i in ids if i not in[a,b]):continue
       for a,b in [(a,b),(b,a)]:
        e0=(W[a]-W[b])/2;qe=float(e0@Q@e0);assert qe>0 and e0@Q@u==0
        comp=next(c for c in cusp['affine_components']if a in c['walls'])
        marks=dict(zip(comp['walls'],map(F,comp['positive_kernel'])));m=marks[a];assert marks[b]==m
        # Removing b: projected positive kernel has coefficient 2m on a.
        T=F(comp['ray_multiple'])/(2*m)
        keep=[i for i in ids if i!=b];v=C@e0;H=G[np.ix_(keep,keep)]-np.outer(v[keep],v[keep])/qe
        ratios=[];ok=True
        for i,j in combinations(range(len(keep)),2):
          r=H[i,j]**2/(H[i,i]*H[j,j]);ratios.append(r)
          if H[i,j]>1e-8 or min(abs(r-z)for z in[0,.25,.5,.75,1])>1e-8:ok=False
        if not ok:continue
        row={'cusp_id':cusp['id'],'leaf_pair':[a,b],'cusp':u.tolist(),'base_normal':e0.tolist(),'norm':qe,'upper_parameter':str(T),'local_rank':int(np.linalg.matrix_rank(H)),'samples':[]}
        for f in[.05,.25,.5,.75,.95]:
          t=float(T)*f;e=e0+t*u;v=C@e;HH=G-np.outer(v,v)/qe
          center=h-(h@Q@e)/qe*e
          # Project center onto pairs; this is a sufficient actual-ridge test.
          bad=[];tested=0
          cand=[]
          for i,j in combinations(range(210),2):
            if b in[i,j]or min(HH[i,i],HH[j,j])<1e-8:continue
            det=HH[i,i]*HH[j,j]-HH[i,j]**2
            if det<1e-8:continue
            ratio=HH[i,j]**2/(HH[i,i]*HH[j,j]);co=-HH[i,j]/np.sqrt(HH[i,i]*HH[j,j])
            if HH[i,j]<=1e-8 and min(abs(co-np.cos(np.pi/m))for m in range(2,41))<1e-7:continue
            # Favor small indices and local/nonlocal incidence mix.
            cand.append((i,j,ratio,co))
          cp=C@center
          for i,j,ratio,co in cand:
            A=HH[np.ix_([i,j],[i,j])];coef=np.linalg.solve(A,cp[[i,j]])
            rp=W[[i,j]]-np.outer(v[[i,j]],e)/qe;x=center-coef@rp
            xp=C@x
            if x@Q@x<-1e-6 and np.max(xp)<1e-7:
              bad.append({'pair':[i,j],'cosine':float(co),'cosine_squared':float(ratio),'point':x.tolist()});break
          if not bad:
            for i,j,ratio,co in cand[:120]:
              pp=point_at(e,[i,j]);tested+=1
              if pp and pp['slack']>1e-8 and pp['norm']<0:
                bad.append({'pair':[i,j],'cosine':float(co),'cosine_squared':float(ratio),**pp});break
          row['samples'].append({'fraction':f,'parameter':t,'bad_ridges':bad,'lp_tests':tested})
        rows.append(row);print('pencil',cusp['id'],a,b,'T',T,'samples',[(x['fraction'],[z['pair']for z in x['bad_ridges']])for x in row['samples']],flush=True)
        OUT.write_text(json.dumps({'status':'NUMERICAL_FORK_CUSP_SECTION_SCOUT','pencils':rows},indent=2)+'\n')
    return rows
if __name__=='__main__':scan()
