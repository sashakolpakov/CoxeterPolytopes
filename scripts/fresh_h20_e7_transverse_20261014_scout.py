#!/usr/bin/env python3
"""Numerical scout only: 44 vertex-ray x 7 single-weight transverse caps."""
import json, math
from pathlib import Path
from fractions import Fraction
from collections import Counter
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_transverse_20261014_scout'
BASE='results/fresh_h20_e7_intrinsic_caps_20261014_cusps.json'
VERTICES='results/fresh_h20_e7_intrinsic_caps_20261014.json'

def load():
    d,b=[json.loads((ROOT/p).read_text())for p in [BASE,VERTICES]]
    radicals=np.sqrt([1,2,79,158,11,22,869,1738])
    def scalar(c):return sum(float(Fraction(a))*r for a,r in zip(c,radicals))
    def vector(v):return np.array([scalar(c)for c in v])
    Q=np.array([[float(Fraction(a))for a in row]for row in d['form']])
    N=np.array([vector(v)for v in d['normals']]);H=vector(d['strict_center']);E=vector(d['strict_escape']['point'])
    S=d['E7_support'];G=N[S]@Q@N[S].T;J=np.linalg.inv(G);W=J@N[S]
    V=np.array([vector(row['ambient_point'])for row in b['face_finite_vertices']])
    q=np.array([-float(Fraction(row['norm']))for row in b['face_finite_vertices']])
    U=V/np.sqrt(q[:,None]);norms=np.array(list(map(float,d['norms'])))
    assert np.max(np.abs(np.diag(U@Q@U.T)+1))<1e-7
    return d,b,Q,N,H,E,S,J,W,U,norms

def run():
    d,b,Q,N,H,E,S,J,W,U,norms=load()
    A=-W@Q@N.T;B=-U@Q@N.T;AE=-W@Q@E;BE=-U@Q@E
    AH=-W@Q@H;BH=-U@Q@H
    rows=[];candidates=[];counts=Counter();phi=(1+math.sqrt(5))/2
    def gate(c,vi,j):
        dd=J[j,j]
        if dd*c*c<=2+1e-10:return None
        r=-c*W[j]+math.sqrt(dd*c*c-2)*U[vi]
        pair=r@Q@N.T;cos=-pair/np.sqrt(2*norms)
        if r@Q@E<=1e-8:return None
        if max(pair)>1e-7 or r@Q@H>=-1e-8:return None
        # Numerical recognition only; no finite-label completeness claim.
        defects=[];labels=[]
        for i,x in enumerate(cos):
            if x>=1-1e-8:labels.append('dotted_or_parallel');continue
            if abs(x)<1e-7:labels.append(2);continue
            m=math.pi/math.acos(max(-1,min(1,x)));mi=round(m)
            if mi<3 or abs(m-mi)>1e-6:defects.append([i,float(x),float(m)])
            labels.append(mi)
        return dict(vertex=vi,E7_position=j,E7_wall=S[j],c=float(c),escape_pairing=float(r@Q@E),defects=defects,labels=labels)
    for vi in range(44):
        for j in range(7):
            dd=J[j,j];lower=max(0,math.sqrt(max(0,dd-2)));upper=min(math.sqrt(dd),AE[j]/BE[vi])
            constraints=[('Coxeter_c_at_least_one',lower)]
            impossible=[]
            for i,(a,bb)in enumerate(zip(A[j],B[vi])):
                if abs(bb)<1e-7:
                    if a>1e-7:impossible.append(i)
                elif bb>0:
                    if a/bb>lower:lower=a/bb;constraints.append((i,lower))
                else:assert bb>-1e-7,(vi,j,i,bb)
            if AH[j]/BH[vi]>lower:lower=AH[j]/BH[vi];constraints.append(('center',lower))
            rec=dict(vertex=vi,E7_position=j,E7_wall=S[j],weight_norm=float(dd),u_lower=float(lower),u_upper_cut=float(upper),last_lower_source=constraints[-1],zero_slope_blockers=impossible)
            if impossible or lower>=upper-1e-10:
                rec['status']='no_continuous_acute_cutter';counts[rec['status']]+=1
            else:
                clo=math.sqrt(2/(dd-lower*lower));chi=math.sqrt(2/(dd-upper*upper)) if upper<math.sqrt(dd)-1e-12 else None
                rec.update(status='continuous_interval',c_lower=clo,c_upper=chi)
                counts[rec['status']]+=1
                cs=[1,math.sqrt(2),phi,math.sqrt(3),2]
                # Contact equations to any retained wall: (a*c-b*sqrt(d*c²-2))=-sqrt(2n)*cos(pi/m).
                for i,(a,bb)in enumerate(zip(A[j],B[vi])):
                    if abs(bb)<1e-8:continue
                    for m in [2,3,4,5,6,8,10,12,16]:
                        k=math.sqrt(2*norms[i])*math.cos(math.pi/m)
                        aa=dd*bb*bb-a*a;ab=-2*a*k;ac=-2*bb*bb-k*k
                        rr=np.roots([aa,ab,ac]) if abs(aa)>1e-12 else [-ac/ab] if abs(ab)>1e-12 else []
                        for z in rr:
                            if abs(complex(z).imag)>1e-8:continue
                            c=float(complex(z).real)
                            if c>0 and dd*c*c>2 and abs(a*c-bb*math.sqrt(dd*c*c-2)+k)<1e-6:cs.append(c)
                checked=[]
                for c in sorted(cs):
                    if c<clo-1e-7 or (chi is not None and c>=chi+1e-7):continue
                    if checked and abs(c-checked[-1])<1e-8:continue
                    checked.append(c);ans=gate(c,vi,j)
                    if ans is not None:
                        if not ans['defects']:candidates.append(ans)
                rec['contact_candidates_tested']=len(checked)
                rec['discrete_c_at_least_one_gate']=[gate(c,vi,j)for c in [1,math.sqrt(2),phi,math.sqrt(3),2] if c>=clo-1e-7 and (chi is None or c<chi+1e-7)]
            rows.append(rec)
    return dict(status='NUMERICAL_SCOUT_NOT_CERTIFICATE',counts=dict(counts),E7_walls=S,weight_norms=np.diag(J).tolist(),escape_A=AE.tolist(),records=rows,candidates=candidates)

if __name__=='__main__':
    out=run();(ROOT/f'results/{TAG}.json').write_text(json.dumps(out,indent=2)+'\n')
    print(out['counts'],'Coxeter cutters',len(out['candidates']))
    print('intervals',[(r['vertex'],r['E7_wall'],r.get('c_lower'),r.get('c_upper'))for r in out['records']if r['status']=='continuous_interval'])
