#!/usr/bin/env python3
"""Bounded noncentral proposals on affine/Lorentzian three-wall supports.

First impose (e,r_l)=0, repairing a family of right angles; solve a second
ridge angle exactly along the resulting t=a+b*s line. The old central line
is explicitly excluded. Each accepted algebraic proposal needs a strict
timelike center, two actual repaired ridges, and an actual global ridge gate.
Failed center scouts are unresolved proposals, never geometric exclusions.
"""
from __future__ import annotations
import argparse
from collections import Counter
from fractions import Fraction
import json
from types import SimpleNamespace
import numpy as np
import sympy as sp
import fresh_h20_triple_central_20261013 as source
import fresh_h20_triple_angle_targets_20261013 as old
from fresh_h20_global_section_parameters_20261013 import Field,coeffs

ROOT=source.ROOT
OUTPUT=ROOT/'data/enumerations/fresh_h20_bivariate_nonelliptic_20261013.json'
u=old.t;P=old.P

def section(Q,W,G,h,i,j,k,a,b):
    e0=W[i]-a*W[k];e1=W[j]-b*W[k]
    e=[P(e0[d]+u*e1[d]) for d in range(22)]
    q=P(e0@Q@e0+2*u*(e0@Q@e1)+u*u*(e1@Q@e1))
    v=[P(w@Q@e0+u*(w@Q@e1)) for w in W]
    def gr(a,b):return q*int(G[a,b])-v[a]*v[b]
    he=P(h@Q@e0+u*(h@Q@e1))
    center=[P(q*int(h[d])-he*e[d]) for d in range(22)]
    return (e,q,v,gr),center,he

def scout_center(Q,W,h,sd,value):
    e=np.array([float(p.eval(value)) for p in sd[0]]);Qf=np.array(Q,dtype=float);Wf=np.array(W,dtype=float);hf=np.array(h,dtype=float)
    q=float(e@Qf@e)
    if q<=1e-8:return False
    C=q*hf-(hf@Qf@e)*e
    return C@Qf@C<0 and C@Qf@hf<0 and max(Wf@Qf@C)<-1e-8

def selection(Q,W,G,h,limit):
    orbits=json.loads(source.ORBIT_OUTPUT.read_text())['orbits'];chosen=[]
    cp=[-int(w@Q@h) for w in W]
    for orbit in orbits:
        ids=[x-1 for x in orbit['representative']];M=sp.Matrix(G[ids,:][:,ids])
        if M[:2,:2].det()>0 and M.det()>0:continue
        proposals=[]
        for k in ids:
            i,j=[z for z in ids if z!=k];value=sp.Rational(cp[i],cp[j]);seen=set()
            for l in range(210):
                if l in ids or G[l,k]>=0:continue
                a=sp.Rational(int(G[l,i]),int(G[l,k]));b=sp.Rational(int(G[l,j]),int(G[l,k]))
                if (a,b) in seen or (a==0 and b==0):continue
                seen.add((a,b))
                if (a,b)==(sp.Rational(cp[i],cp[k]),sp.Rational(cp[j],cp[k])):continue
                sd,C,he=section(Q,W,G,h,i,j,k,a,b)
                if not scout_center(Q,W,h,sd,value):continue
                assert sd[2][l].is_zero and not he.is_zero
                parameter=SimpleNamespace(lo=value,hi=value)
                pair=None
                for aa,bb in old.candidates(Q,W,h,sd,C,parameter):
                    if sd[3](aa,bb).eval(value)>=0:continue
                    ratio=sp.cancel(sd[3](aa,bb)**2/(sd[3](aa,aa)*sd[3](bb,bb)))
                    if u not in ratio.free_symbols:continue
                    pair=(aa,bb);break
                if pair is None:continue
                distance=abs(float(a+b*value)-2*cp[i]/cp[k])
                proposals.append((distance,i,j,k,l,a,b,pair,sd,C,he,value))
        if not proposals:continue
        row=min(proposals,key=lambda p:p[:8]);chosen.append((orbit,row))
        print('selected',orbit['representative'],'negative',row[3]+1,'aux',row[4]+1,'line',str(row[5]),str(row[6]),flush=True)
        if len(chosen)>=limit:break
    return chosen

def actual_ridge(Q,W,h,sd,C,F,a,b,cp):
    e,q,v,gr=sd;aa,bb,ab=map(F.rem,(gr(a,a),gr(b,b),gr(a,b)));D=F.rem(aa*bb-ab*ab)
    if min(F.sign(aa),F.sign(bb),F.sign(D))<=0:return None
    A=F.rem(bb*cp[a]-ab*cp[b]);B=F.rem(aa*cp[b]-ab*cp[a])
    vals=[F.rem(D*cp[l]-A*gr(l,a)-B*gr(l,b)) for l in range(210)]
    signs=[F.sign(p) for p in vals]
    if max(signs)>0:return None
    zeros=[l for l,p in enumerate(signs) if p==0]
    if a not in zeros or b not in zeros:return None
    cone=[]
    for l in zeros:
        if l in (a,b):continue
        al=F.rem(bb*gr(l,a)-ab*gr(l,b));be=F.rem(aa*gr(l,b)-ab*gr(l,a))
        if min(F.sign(al),F.sign(be))<0:return None
        if not F.rem(D*gr(l,l)-al*gr(l,a)-be*gr(l,b)).is_zero:return None
        cone.append({'wall_one_based':l+1,'alpha_numerator':coeffs(al),'beta_numerator':coeffs(be),'denominator':coeffs(D)})
    p=[F.rem(D*C[d]-A*(q*int(W[a,d])-v[a]*e[d])-B*(q*int(W[b,d])-v[b]*e[d])) for d in range(22)]
    square=old.pdot(Q,p,p,F);hp=old.pdot(Q,p,[P(int(x)) for x in h],F)
    assert F.sign(square)<0 and F.sign(hp)<0 and old.pdot(Q,p,e,F).is_zero
    return {'ridge_walls_one_based':[a+1,b+1],'projected_gram':[coeffs(aa),coeffs(ab),coeffs(bb)],
            'point_coordinate_polynomials':[coeffs(x) for x in p],'point_square':coeffs(square),'point_h_pairing':coeffs(hp),
            'zero_walls_one_based':[l+1 for l in zeros],'redundant_zero_wall_cones':cone}

def run_pencil(Q,W,G,h,orbit,row):
    distance,i,j,k,l,a,b,target,sd,rawC,he,value=row;e,q,v,gr=sd
    aa,bb=target;N=gr(aa,bb)**2;den=gr(aa,aa)*gr(bb,bb)
    equations={3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}
    factors={}
    for m,eq in equations.items():
        assert not eq.is_zero
        for f,power in sp.factor_list(eq)[1]:factors[tuple(coeffs(f.monic()))]=f.monic()
    rows=[]
    for key,poly in sorted(factors.items()):
        for (lo,hi),mult in poly.intervals(eps=sp.Rational(1,10**18)):
            if hi<=0:continue
            if lo==hi:lo-=sp.Rational(1,10**20);hi+=sp.Rational(1,10**20)
            if lo<=0:lo,hi=sp.polys.polytools.refine_root(poly,lo,hi,eps=sp.Rational(1,10**30))
            assert lo>0
            F=Field(poly,lo,hi);record={'minimal_polynomial':list(key),'isolating_interval':[str(lo),str(hi)],'approximate_s':float((lo+hi)/2)}
            if F.sign(q)<=0:
                record['status']='NO_SPACELIKE_SECTION_NORMAL';rows.append(record);continue
            if min(F.sign(gr(aa,aa)),F.sign(gr(bb,bb)))<=0 or F.sign(gr(aa,bb))>0:
                record['status']='TARGET_NOT_ACUTE_ELLIPTIC_PAIR';rows.append(record);continue
            labels=[m for m,eq in equations.items() if F.sign(eq)==0 and (m!=5 or F.sign(2*N-den)>0)]
            if not labels:
                record['status']='CONJUGATE_ANGLE_BRANCH';rows.append(record);continue
            C=[F.rem(p) for p in rawC]
            cp=[old.pdot(Q,C,[P(int(x)) for x in w],F) for w in W]
            if F.sign(old.pdot(Q,C,C,F))>=0 or F.sign(old.pdot(Q,C,[P(int(x)) for x in h],F))>=0 or any(F.sign(p)>=0 for p in cp):
                record.update(status='UNRESOLVED_STRICT_CENTER_SCOUT',target_labels=labels);rows.append(record);continue
            target_proof=actual_ridge(Q,W,h,sd,C,F,aa,bb,cp)
            if target_proof is None:
                record.update(status='UNRESOLVED_TARGET_RIDGE_SCOUT',target_labels=labels);rows.append(record);continue
            right=None
            for other in range(210):
                if other==l or G[l,other]!=0:continue
                right=actual_ridge(Q,W,h,sd,C,F,l,other,cp)
                if right is not None:break
            if right is None:
                record.update(status='UNRESOLVED_RIGHT_RIDGE_SCOUT',target_labels=labels);rows.append(record);continue
            bad=None
            for x,y in old.candidates(Q,W,h,sd,C,F):
                bad=old.witness(Q,W,h,sd,C,F,x,y,cp)
                if bad is not None:break
            record.update(isolating_interval=[str(F.lo),str(F.hi)],status='NON_COXETER_ACTUAL_RIDGE' if bad else 'OPEN_COXETER_CANDIDATE',
                          target_labels=labels,strict_interior=[coeffs(p) for p in C],
                          repaired_target_ridge=target_proof,repaired_right_ridge=right,actual_forbidden_ridge=bad)
            rows.append(record)
            print('candidate',[i+1,j+1,k+1],labels,record['approximate_s'],record['status'],flush=True)
    return {'triple_one_based':orbit['representative'],'triple_gram':[[int(G[x,y]) for y in (i,j,k)] for x in (i,j,k)],
            'positive_walls_one_based':[i+1,j+1],'negative_wall_one_based':k+1,'auxiliary_orthogonal_wall_one_based':l+1,
            'line_t_coefficients':[str(a),str(b)],'section_center_pairing':coeffs(he),
            'scout_parameter':str(value),'target_pair_one_based':[aa+1,bb+1],'target_labels':[3,4,5,6],
            'target_equations':{str(m):coeffs(eq) for m,eq in equations.items()},'candidates':rows}

def build(limit):
    Q,W,G,h=source.base.load();selected=selection(Q,W,G,h,limit);rows=[]
    for orbit,row in selected:
        rows.append(run_pencil(Q,W,G,h,orbit,row))
        OUTPUT.write_text(json.dumps({'complete':False,'pencils':rows},indent=2)+'\n')
    counts=Counter(c['status'] for p in rows for c in p['candidates'])
    out={'schema':'h20-nonelliptic-two-angle-proposals-v1','complete':True,'source':str(source.base.SOURCE.relative_to(ROOT)),
         'selection_limit':limit,'pencil_count':len(rows),'pencils':rows,'statuses':dict(counts),
         'scope':'Bounded noncentral proposal search, first selected nonelliptic triple orbits under the deterministic declared scout. Right-angle orthogonality line plus labels3,4,5,6 on a second pair. Every accepted candidate has a strict interior and two actual repaired ridges. Unresolved scouts are not exclusions; no completeness over nonelliptic supports or coefficient ratios.'}
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n');print('DONE',dict(counts),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--limit',type=int,default=12);a=p.parse_args();build(a.limit)
