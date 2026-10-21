#!/usr/bin/env python3
"""Continuous LP screen and finite-anchor scout after releasing one old wall.

Retain19old incident walls right; new_plus right; new_minus labelm.
All floating results are proposals, not exact exclusion/existence claims.
"""
import json,math
from collections import Counter
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import fresh_h22_replacement_edges_20261013 as P
from fresh_h22_nonright_bridge_20261013_scout import roots
R=P.R;TAG=P.TAG+'_target_scout';OUT='results/'+TAG+'.json'
NEGATIVE=(7,14,16,17,21,23,25,28,29,30,31,48,50,63)
LABELS=(3,4,5,6,8,10,12,16,24,32)
ANCHOR_LABELS=tuple(range(2,13))
def f(a):return sum(float(c)*math.sqrt(n)for n,c in a.items())
def run():
    src,M,wr,cov,h,v,inc,U=P.setup();ids=list(wr);ix={i:k for k,i in enumerate(ids)}
    W=np.array([[f(a)for a in wr[i]]for i in ids]);M=np.array(M,float);G=W@M@W.T
    h=np.array([f(a)for a in h]);point=np.array([f(a)for a in v]);v=point/math.sqrt(float(U))
    vv=W@M@v;vv[[ix[i]for i in inc]]=0
    AA=G[np.ix_([ix[i]for i in inc],[ix[i]for i in inc])]
    # The exact incidence is oldGI direct-sum2I2.
    AA[-2:,:]=0;AA[:,-2:]=0;AA[-2,-2]=AA[-1,-1]=2
    AI=np.linalg.inv(AA);esc=np.array(src['strict_escape']['coordinates'],float);ec=M@esc
    minus=W[ix['new_minus']];ss=-G[ix['new_minus']];ss[[ix[i]for i in inc]]=0;ss[ix['new_minus']]=-2
    anchors=np.array([ix[i]for i in range(210,511)]);qn=np.einsum('ij,ij->i',W@M,W)
    qn=np.round(qn);qn[ix['new_plus']]=qn[ix['new_minus']]=2
    ad=-np.sqrt(2*qn[anchors,None])*np.cos(np.pi/np.array(ANCHOR_LABELS))[None,:];ad[:,0]=0
    bounds=np.sqrt(2*qn);counts=Counter();rows=[];accepted=[];near=[]
    for old in NEGATIVE:
        k=inc.index(old);dual=AI[k]@W[[ix[i]for i in inc]];e=dual/math.sqrt(AI[k,k]);ee=W@M@e
        ee[[ix[i]for i in inc]]=0;ee[ix[old]]=1/math.sqrt(AI[k,k])
        A=np.column_stack((vv,-ee));rhs=-ss
        # Y>=X and Y-X<=sqrt6 are necessary for allm>=3.
        A=np.vstack((A,[1,-1],[-1,1]));rhs=np.r_[rhs,0,math.sqrt(6)]
        objective=np.array([v@ec,-e@ec]);res=linprog(-objective,A_ub=A,b_ub=rhs,bounds=[(0,None),(0,None)],method='highs')
        row={'released':old,'linear_screen_status':res.message,'linear_screen_success':bool(res.success),'accepted':0}
        if res.success:
            row.update(linear_maximum_scaled_escape=float(-minus@ec-res.fun),linear_point=res.x.tolist())
        rows.append(row);counts['linear_feasible'if res.success else'linear_infeasible']+=1
        if res.success and row['linear_maximum_scaled_escape']<=1e-7:
            counts['linear_no_target_cut']+=1;continue
        if not res.success:continue
        for m in LABELS:
            c=math.cos(math.pi/m);s=-c*minus;sp=c*ss;residual=2*(1-c*c);kk=math.sqrt(residual)/2
            rr,deg=roots(kk*(vv[anchors]-ee[anchors])[:,None],sp[anchors,None]-ad,-kk*(vv[anchors]+ee[anchors])[:,None]);counts['degenerate']+=deg
            for tt,valid in rr:
                valid&=np.isfinite(tt)&(tt>1+1e-9);ai,mi=np.where(valid)
                t=tt[ai,mi];x=kk*(t-1/t);y=-kk*(t+1/t);counts['outgoing']+=len(t)
                cut=s@ec+x*(v@ec)+y*(e@ec)>1e-7
                ai,mi,t,x,y=[z[cut]for z in(ai,mi,t,x,y)];counts['target_cuts']+=len(t)
                for j in range(len(ids)):
                    if not len(t):break
                    good=sp[j]+x*vv[j]+y*ee[j]<=1e-6
                    ai,mi,t,x,y=[z[good]for z in(ai,mi,t,x,y)]
                counts['acute_target_cuts']+=len(t)
                for aj,mj,tj,xj,yj in zip(ai,mi,t,x,y,strict=True):
                    pp=sp+xj*vv+yj*ee;pp[[ix[i]for i in inc if i!=old]]=0;pp[ix['new_minus']]=-2*c;pp[anchors[aj]]=ad[aj,mj]
                    cos=-pp/bounds;js=np.where(cos<1-1e-6)[0];mm=np.rint(np.pi/np.arccos(np.clip(cos[js],-1,1))).astype(int)
                    error=float(max(abs(cos[js]-np.cos(np.pi/np.maximum(mm,2))),default=0))
                    r=s+xj*v+yj*e;rec={'released':old,'new_minus_label':m,'parent_anchor':ids[anchors[aj]],'anchor_label':ANCHOR_LABELS[mj],
                        't':float(tj),'x':float(xj),'y':float(yj),'normal':r.tolist(),'norm':float(r@M@r),'height':float(r[-1]),
                        'center_pairing':float(r@M@h),'escape_pairing':float(r@ec),'maximum_angle_error':error}
                    if error<1e-6 and r@M@h<0:
                        accepted.append(rec);row['accepted']+=1;counts['accepted']+=1
                    else:
                        near.append(rec);near.sort(key=lambda z:z['maximum_angle_error']);near=near[:20]
        print(old,row,flush=True)
    out={'status':'NUMERICAL_NONRIGHT_ADDITION_SCREEN','sources':{p:R.sha(p)for p in(P.SOURCE,'scripts/'+P.TAG+'.py','scripts/'+TAG+'.py','scripts/fresh_h22_nonright_bridge_20261013_scout.py')},
         'edge_labels':list(LABELS),'anchor_labels':list(ANCHOR_LABELS),'counts':dict(counts),'rows':rows,'accepted':accepted,'closest_failures':near,
         'scope':'Only14specified oldwall releases, right19oldwalls/new_plus, finitecontact new_minus; numericalLP relaxation and boundedlabel quadratic scout. No exact exclusion or existence claim.'}
    (R.P/OUT).write_text(json.dumps(out,indent=2)+'\n');print(dict(counts),OUT,flush=True);return out
if __name__=='__main__':run()
