#!/usr/bin/env python3
"""Bounded nonright spherical endpoint profiles on four frozen H22 vertices.

Numerical proposals only: prescribed edge/anchor labels and profile support
are explicit. No firstness, arithmetic completeness, or universal failure.
"""
import argparse,json,math,time
from pathlib import Path
from itertools import combinations,product
from collections import Counter
import numpy as np
import fresh_h22_contact_bridge_20261013_scout as S

TAG='fresh_h22_nonright_bridge_20261013'
ROOT=Path(__file__).resolve().parents[1]

def roots(A,B,C):
    """All ordinary quadratic roots; degenerate cases counted by caller."""
    aa=np.broadcast_to(A,B.shape);cc=np.broadcast_to(C,B.shape)
    disc=B*B-4*aa*cc;ok=disc>=-1e-10
    dd=np.sqrt(np.maximum(disc,0));den=2*aa
    regular=abs(aa)>1e-11
    out=[]
    for sg in (-1,1):
        tt=np.divide(-B+sg*dd,den,out=np.full(B.shape,np.nan),where=regular)
        if sg==1:
            linear=(~regular)&(abs(B)>1e-11)
            tt[linear]=-cc[linear]/B[linear]
        out.append((tt,ok&(regular|(abs(B)>1e-11))))
    return out,int(np.sum((~regular)&(abs(B)<=1e-11)))

def run(links,labels,anchorlabels):
    start=time.monotonic();M,W,h,G,n,meta,doc,vertices=S.load();WM=W@M
    types=doc['selected_original_types']+[23,24]
    anchors=np.array(list(range(210,511))+[812+2*types.index(k)for k in (0,2,6,22)])
    ac=np.cos(np.pi/np.array(anchorlabels));ac[np.array(anchorlabels)==2]=0
    ad=-np.sqrt(2*n[anchors,None])*ac[None,:]
    counts=Counter();edge_rows=[];accepted=[];near=[];seen=set()
    bounds=np.sqrt(2*n);hc=M@h
    target=np.array(list(map(S.Q,doc['strict_escape']['vector'])),float);tc=M@target
    # An acute failure is often visible in the extra caps; check those first.
    order=np.array(list(range(812,860))+list(range(210,812))+list(range(210)))
    for vi,source in enumerate(doc['compact_H2_vertices']):
        inc=source['incident'];A=G[np.ix_(inc,inc)];AI=np.linalg.inv(A)
        v=vertices[vi]/math.sqrt(-source['norm']);va=WM@v;va[inc]=0.
        for k,released in enumerate(inc):
            edge=[i for i in inc if i!=released];ix=[i for i in range(22)if i!=k]
            EI=AI[np.ix_(ix,ix)]-np.outer(AI[ix,k],AI[k,ix])/AI[k,k]
            dual=AI[k]@W[inc];e=dual/math.sqrt(AI[k,k]);ea=WM@e
            ea[inc]=0.;ea[released]=1/math.sqrt(AI[k,k])
            DD=EI@G[edge];DV=EI@W[edge]
            row={'vertex':vi,'released':released,'edge':edge,'profiles':0,'candidates':0,'acute':0,'accepted':0}
            edge_rows.append(row)
            # Single-link feasibility is a necessary filter for every pair.
            choices=[(j,m,-math.sqrt(2*n[edge[j]])*math.cos(math.pi/m))
                     for j in range(21)for m in labels
                     if 2*n[edge[j]]*math.cos(math.pi/m)**2*EI[j,j]<2-1e-9]
            profiles=[z for z in combinations(choices,links)
                      if len({x[0]for x in z})==links]
            for profile in profiles:
                js=np.array([x[0]for x in profile]);p=np.array([x[2]for x in profile])
                q=float(p@EI[np.ix_(js,js)]@p);R=2-q
                if R<=1e-9:counts['nonpositive_schur_residual']+=1;continue
                counts['spherical_profiles']+=1;row['profiles']+=1
                s=p@DV[js];sp=p@DD[js];kk=math.sqrt(R)/2
                polyA=kk*(va[anchors]-ea[anchors])[:,None]
                polyB=sp[anchors,None]-ad
                polyC=-kk*(va[anchors]+ea[anchors])[:,None]
                rr,deg=roots(polyA,polyB,polyC);counts['degenerate_equations']+=deg
                for tt,ok in rr:
                    ok&=np.isfinite(tt)&(tt>1+1e-9)
                    ai,mi=np.where(ok)
                    if not len(ai):continue
                    t=tt[ai,mi];xx=kk*(t-1/t);yy=-kk*(t+1/t)
                    counts['outgoing_roots']+=len(t);row['candidates']+=len(t)
                    # Remaining gates are numerical and deliberately tolerant.
                    good=(s[-1]+xx*v[-1]+yy*e[-1]>1e-8)&(s@hc+xx*(v@hc)+yy*(e@hc)<-1e-8)
                    counts['height_or_center_rejected']+=len(t)-int(sum(good))
                    ai,mi,t,xx,yy=[z[good]for z in (ai,mi,t,xx,yy)]
                    for wall in order:
                        if not len(t):break
                        good=sp[wall]+xx*va[wall]+yy*ea[wall]<=2e-6
                        ai,mi,t,xx,yy=[z[good]for z in (ai,mi,t,xx,yy)]
                    if not len(t):continue
                    row['acute']+=len(t);counts['acute_roots']+=len(t)
                    for aj,mj,tj,xj,yj in zip(ai,mi,t,xx,yy,strict=True):
                        pp=sp+xj*va+yj*ea;pp[edge]=0
                        for j,m,val in profile:pp[edge[j]]=val
                        anchor=int(anchors[aj]);pp[anchor]=ad[aj,mj]
                        cos=-pp/bounds;finite=np.where(cos<1-2e-6)[0]
                        mm=np.rint(np.pi/np.arccos(np.clip(cos[finite],-1,1))).astype(int)
                        errors=abs(cos[finite]-np.cos(np.pi/np.maximum(mm,2)))
                        maxerr=float(max(errors,default=0));r=s+xj*v+yj*e
                        rec={'vertex':vi,'released':released,'edge_profile':[[int(edge[j]),int(m)]for j,m,_ in profile],
                             'anchor':anchor,'anchor_label':int(anchorlabels[mj]),'t':float(tj),'x':float(xj),'y':float(yj),
                             'schur_residual':R,'maximum_angle_error':maxerr,'escape_pair':float(r@tc),
                             'height':float(r[-1]),'center_pair':float(r@hc),'norm':float(r@M@r),
                             'normal':r.tolist(),'finite_contacts':{str(i):int(m)for i,m in zip(finite,mm,strict=True)}}
                        if maxerr>2e-6:
                            counts['forbidden_finite_angle']+=1
                            near.append(rec);near.sort(key=lambda z:z['maximum_angle_error']);near=near[:20]
                            continue
                        key=tuple(np.round(r,7))
                        if key in seen:counts['duplicate']+=1;continue
                        seen.add(key);accepted.append(rec);row['accepted']+=1;counts['accepted']+=1
            print('edge',vi,released,'links',links,'profiles',row['profiles'],'acute',row['acute'],'accepted',row['accepted'],flush=True)
    out={'status':'BOUNDED_NUMERICAL_NONRIGHT_SPHERICAL_PROFILE_SCOUT','links':links,'edge_labels':list(labels),
         'anchor_labels':list(anchorlabels),'anchors':list(map(int,anchors)),'counts':dict(counts),
         'edge_rows':edge_rows,'accepted':accepted,'closest_angle_failures':near,
         'elapsed_seconds':time.monotonic()-start,
         'sources':{p:S.V.sha(p)for p in (S.P.OUT,'scripts/fresh_h22_contact_bridge_20261013_scout.py','scripts/'+Path(__file__).name)},
         'scope':'Numerical finite profile and anchor-label census on the specified88edge supports; no completeness beyond the stated labels/support size, and no exact existence claim from a survivor. Near-degenerate equations are counted rather than certified.'}
    path=ROOT/'results'/f'{TAG}_{links}link_scout.json';path.write_text(json.dumps(out,indent=2)+'\n')
    print(dict(counts),path,'seconds',out['elapsed_seconds'],flush=True)
    return out

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--links',type=int,choices=(1,2,3),default=1)
    p.add_argument('--labels',default='3,4,5,6');p.add_argument('--anchor-labels',default='2,3,4,5,6,8,10,12')
    a=p.parse_args();run(a.links,tuple(map(int,a.labels.split(','))),tuple(map(int,a.anchor_labels.split(','))))
