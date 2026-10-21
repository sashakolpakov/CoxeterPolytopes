#!/usr/bin/env python3
"""Continuous endpoint hyperbolas with strictly dotted nonincident pairs."""
import json,sys,itertools
from fractions import Fraction
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from scripts import fresh_h9_prism_completion_20261013_scout as S
V=S.V

def negative_intervals(a,b,c):
    if abs(a)<1e-11:
        if abs(b)<1e-11:return [(0.,np.inf)]if c< -1e-10 else []
        r=-c/b
        return [(0.,r)]if b>0 and r>0 else [(max(0.,r),np.inf)]if b<0 else []
    disc=b*b-4*a*c
    if disc<=0:return [(0.,np.inf)]if a<0 else []
    rr=sorted(((-b-np.sqrt(disc))/(2*a),(-b+np.sqrt(disc))/(2*a)))
    if a>0:return [(max(0.,rr[0]),rr[1])]if rr[1]>max(0.,rr[0])else []
    out=[]
    if rr[0]>0:out.append((0.,rr[0]))
    out.append((max(0.,rr[1]),np.inf));return out

def intersect(segments,other):
    return [(max(a,c),min(b,d))for a,b in segments for c,d in other if min(b,d)>max(a,c)+1e-9]

def run():
    d=S.setup();G=np.array(d['gram']);old=np.array(d['normals']);perm=[0,2,1,3,4,6,5,7,8,9]
    R=np.vstack([old,old[90:,perm]]);names=d['names']+['mirror'+n[0].upper()+n[1:]for n in d['names'][90:]]
    source=S.PREFIX+'_caps.json';recipe=V.read(source);prior=V.read(S.PREFIX+'_verify.json');cov=R@G
    y=np.array([float(Fraction(x['a']))for x in recipe['escape']]);w=-np.linalg.solve(G,y);h=np.array(d['h']);center=np.array(d['center'])
    edges=[e for e in prior['new_component_frontier']['edges']if e['status']=='open_to_light_cone'];rows=[];stats=[]
    for ei,edge in enumerate(edges):
        ids=[names.index(n)for n in edge['common_facets']];A=R[ids];M=A@G@A.T;Mi=np.linalg.inv(M)
        project=lambda v:v-(A@G@v)@Mi@A
        t=project(center);t/=np.sqrt(-t@G@t);e=project(w);e+=(e@G@t)*t;e/=np.sqrt(e@G@e)
        profiles=[np.zeros(8)]
        for i in range(8):
            for m in(3,5):
                p=np.zeros(8);p[i]=-2*np.cos(np.pi/m)
                if p@Mi@p<2-1e-8:profiles.append(p)
        for i,j in itertools.combinations(range(8),2):
            for m,n in itertools.product((3,5),repeat=2):
                p=np.zeros(8);p[i]=-2*np.cos(np.pi/m);p[j]=-2*np.cos(np.pi/n)
                if p@Mi@p<2-1e-8:profiles.append(p)
        stat={'edge':ei,'profiles':len(profiles),'target_feasible_branches':0,'all_dotted_feasible_branches':0}
        for pi,p in enumerate(profiles):
            base=p@Mi@A;rho=np.sqrt(2-p@Mi@p)
            for sigma in(-1,1):
                aa=rho*(t+sigma*e)/2;bb=rho*(-t+sigma*e)/2
                segments=intersect(negative_intervals(-aa@h,-base@h,-bb@h),negative_intervals(aa@y,base@y,bb@y))
                if not segments:continue
                stat['target_feasible_branches']+=1
                for j in range(len(R)):
                    if j in ids:continue
                    segments=intersect(segments,negative_intervals(cov[j]@aa,cov[j]@base+2,cov[j]@bb))
                    if not segments:break
                if not segments:continue
                stat['all_dotted_feasible_branches']+=1
                for lo,hi in segments:
                    x=2*lo+1 if np.isinf(hi)else (lo+hi)/2
                    r=base+aa*x+bb/x;pp=cov@r
                    assert r@y<0 and r@h>0 and max(pp[j]for j in range(len(R))if j not in ids)<-2+1e-6
                    rows.append({'edge_index':ei,'edge':edge,'profile_labels':[2 if abs(z)<1e-8 else 3 if abs(z+1)<1e-8 else 5 for z in p],
                                 'branch':sigma,'parameter_interval':[lo,None if np.isinf(hi)else hi],'parameter':x,'normal':r.tolist(),'target_cut':float(-r@y)})
        stats.append(stat);print(stat,flush=True)
    out={'status':'NUMERICAL_CONTINUOUS_ENDPOINT_FAMILY_SCREEN','source':source,'source_sha256':V.sha(source),'stats':stats,'candidates':rows,
         'scope':'Each new23 edge, spherical profiles with at most two nonright labels3/5. Both branches of the full norm-two endpoint hyperbola are intersected with target/center signs and all remaining wall pairings strictly below -2. Numerical interval endpoints; no exact universal obstruction.'}
    (V.ROOT/(S.PREFIX+'_continuous_scout.json')).write_text(json.dumps(out,indent=2)+'\n');print('DONE',len(rows))

if __name__=='__main__':run()
