#!/usr/bin/env python3
"""Finite-volume-preserving inward Coxeter-angle cuts with exact pruning.

Replace r by v=r+t*s with t<0, retain s, and preserve a common strict
timelike point. Then r=v+(-t)*s is implied: the new region is contained
in the previous finite-volume region. Only exact nonnegative cone
certificates permit redundant inequalities to be removed.
"""
from __future__ import annotations
import argparse
from fractions import Fraction as F
from math import gcd,isqrt,lcm
import json
from pathlib import Path
import numpy as np
from scipy.optimize import linprog
import sympy as sp
import fresh_h20_mixed_b3_20261013 as mixed
import fresh_h20_mixed_dictionary_20261013 as selection

ROOT=Path(__file__).resolve().parents[1]
OUTPUT=ROOT/'data/enumerations/fresh_h20_mixed_inward_20261013.json'

def sgn(x):
    a,b=x;d=lcm(a.denominator,b.denominator);return mixed.sign(int(a*d),int(b*d))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[0],-x[1]
def mul(x,y):return x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def div(x,y):
    z=mul(x,(y[0],-y[1]));d=y[0]**2-2*y[1]**2;assert d
    return z[0]/d,z[1]/d
def rat_sqrt(x):
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return F(a,b) if a*a==x.numerator and b*b==x.denominator else None
def field_sqrt(x):
    a,b=x
    if b==0:
        r=rat_sqrt(a)
        if r is not None:return r,F(0)
        r=rat_sqrt(a/2)
        return None if r is None else (F(0),r)
    h=rat_sqrt(a*a-2*b*b)
    if h is None:return None
    for h in (h,-h):
        u=rat_sqrt((a+h)/2)
        if u:
            ans=(u,b/(2*u));assert mul(ans,ans)==x
            return ans if sgn(ans)>0 else neg(ans)
    return None
def inner(a,b,Q):return tuple(F(v) for v in mixed.pairing(a,b,Q))

def check_terms(target,terms,R):
    out=np.zeros((2,22),dtype=object)
    for index,a,b in terms:
        a,b=F(a),F(b);assert sgn((a,b))>0
        out[0]+=a*R[index,0]+2*b*R[index,1];out[1]+=a*R[index,1]+b*R[index,0]
    assert np.array_equal(out,target)

def cone_terms(index,active,R):
    others=[i for i in active if i!=index]
    V=np.array(R[others,0],dtype=float)+np.sqrt(2)*np.array(R[others,1],dtype=float)
    y=np.array(R[index,0],dtype=float)+np.sqrt(2)*np.array(R[index,1],dtype=float)
    norms=np.max(abs(V),axis=1);yscale=max(abs(y));V/=norms[:,None];y/=yscale
    prop=linprog(np.ones(len(others)),A_eq=V.T,b_eq=y,bounds=(0,None),method='highs')
    if not prop.success:return None
    ids=[others[i] for i in np.where(prop.x>1e-8)[0]]
    if not ids:return None
    A=sp.Matrix(R[ids,0].T);B=sp.Matrix(R[ids,1].T)
    M=A.row_join(2*B).col_join(B.row_join(A));target=sp.Matrix(list(R[index,0])+list(R[index,1]))
    try:coeff,parameters=M.gauss_jordan_solve(target)
    except ValueError:return None
    if parameters.rows:return None
    terms=[];n=len(ids)
    for j,i in enumerate(ids):
        a,b=F(coeff[j]),F(coeff[n+j])
        if a or b:
            if sgn((a,b))<=0:return None
            terms.append([i,str(a),str(b)])
    check_terms(R[index],terms,R);return terms

def proposed_rotation(i,j,m,R,Q,center):
    r,s=R[i],R[j];p=inner(r,s,Q);q=inner(s,s,Q);rr=inner(r,r,Q)
    delta=add(rr,neg(div(mul(p,p),q)))
    if sgn(delta)<=0:return None
    h=(F(0),F(0)) if m==2 else field_sqrt(mul(div(delta,q),(F(1,3) if m==3 else F(3) if m==6 else F(1),F(0))))
    if h is None:return None
    if m==8:h=mul(h,(F(1),F(1)))
    t=add(neg(div(p,q)),neg(h))
    if sgn(t)>=0:return None
    d=lcm(t[0].denominator,t[1].denominator)
    A=d*r[0]+int(d*t[0])*s[0]+2*int(d*t[1])*s[1]
    B=d*r[1]+int(d*t[0])*s[1]+int(d*t[1])*s[0]
    g=gcd(*map(int,A),*map(int,B));v=np.array([A//g,B//g],dtype=object)
    if sgn(inner(v,(center,np.zeros(22,dtype=object)),Q))>=0:return None
    assert not selection.conflict(mixed.pairing(v,s,Q),mixed.pairing(v,v,Q),mixed.pairing(s,s,Q))
    # r = (g/d)*v + (-t)*s, both coefficients strictly positive.
    return v,{'replaced_index':i,'retained_index':j,'target_angle_m':m,
              'new_coefficient':[str(F(g,d)),'0'],'retained_coefficient':[str(-t[0]),str(-t[1])]}

def conflicts(active,R,Q):
    V=R[active];A,B=V[:,0],V[:,1];Ga=A@Q@A.T+2*B@Q@B.T;Gb=A@Q@B.T+B@Q@A.T
    edges=[];degree={i:0 for i in active}
    for a,i in enumerate(active):
        for b in range(a):
            j=active[b]
            if selection.conflict((int(Ga[a,b]),int(Gb[a,b])),(int(Ga[a,a]),int(Gb[a,a])),(int(Ga[b,b]),int(Gb[b,b]))):
                edges.append((i,j));degree[i]+=1;degree[j]+=1
    return edges,degree

def surviving_ridge(active,R,Q,f):
    old=json.loads((ROOT/'data/enumerations/h20_gluing_independent_audit_20261013.json').read_text())
    p=np.array(old['six_tile_patch']['non_coxeter_ridge']['relative_interior_timelike_point'],dtype=object)
    assert p@Q@p<0 and p@Q@f==0
    values=[mixed.pairing(R[i],(p,np.zeros(22,dtype=object)),Q) for i in active]
    if any(mixed.sign(*v)>0 for v in values):return None
    zeros=[i for i,v in zip(active,values) if v==(0,0)]
    if len(zeros)!=2:return None
    gram=[[mixed.pairing(R[i],R[j],Q) for j in zeros] for i in zeros]
    if gram!=[[(6,0),(-2,0)],[(-2,0),(6,0)]]:return None
    return {'indices':zeros,'normal_gram':[[6,-2],[-2,6]],'cosine_squared':'1/9',
            'relative_interior_timelike_point':list(map(int,p)),'square':int(p@Q@p),
            'all_other_current_inequalities_strict':True}

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--steps',type=int,default=12);parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    Q,W,G=mixed.data();old=json.loads(selection.OLD.read_text())['six_tile']['primitive_outward_normals']
    zero=np.zeros(22,dtype=object);R=np.array([(r,zero) for r in old],dtype=object)
    center=np.array(json.loads(selection.SELECTED.read_text())['common_strict_timelike_center'],dtype=object)
    assert center@Q@center<0 and center@Q@W[0]==0
    assert all(sgn(inner(r,(center,zero),Q))<0 for r in R)
    if args.verify:
        out=json.loads(OUTPUT.read_text());active=list(range(580));authoritative=R.copy();R=np.array(out['normal_pairs'],dtype=object)
        assert np.array_equal(R[:580],authoritative)
        assert out['common_strict_timelike_center']==list(map(int,center))
        seen=set(active)
        for row in out['steps']:
            i,j,k=row['replaced_index'],row['retained_index'],row['new_index'];assert i in active and j in active
            assert k not in seen and k==len(seen);seen.add(k)
            check_terms(R[i],[[k,*row['new_coefficient']],[j,*row['retained_coefficient']]],R)
            assert sgn(inner(R[k],(center,zero),Q))<0 and inner(R[k],(W[0],zero),Q)==(0,0)
            assert sgn(inner(R[k],R[k],Q))>0
            assert not selection.conflict(mixed.pairing(R[k],R[j],Q),mixed.pairing(R[k],R[k],Q),mixed.pairing(R[j],R[j],Q))
            active.remove(i);active.append(k)
            for drop in row['redundancy_proofs']:
                assert drop['removed_index'] in active and all(t[0] in active and t[0]!=drop['removed_index'] for t in drop['terms'])
                check_terms(R[drop['removed_index']],drop['terms'],R);active.remove(drop['removed_index'])
        assert active==out['final_active_indices']
        edges,degree=conflicts(active,R,Q);assert len(edges)==out['final_conflict_count']
        assert out['surviving_non_coxeter_ridge']==surviving_ridge(active,R,Q,W[0])
        print('PASS_EXACT_FINITE_VOLUME_CONTAINMENT_AND_PRUNING',len(out['steps']),'cuts',len(active),'remaining inequalities',len(edges),'conflicts',flush=True);return
    active=list(range(580));out={'source_patch':str(selection.OLD.relative_to(ROOT)),'steps':[],
        'common_strict_timelike_center':list(map(int,center)),
        'scope':'Every recorded region is nonempty, full-dimensional and contained in the original finite-volume six-tile patch. A Coxeter assertion additionally requires complete compatibility after exact pruning.'}
    for step in range(args.steps):
        edges,degree=conflicts(active,R,Q);print('Step',step,'walls',len(active),'conflicts',len(edges),flush=True)
        if not edges:out['status']='FINITE_VOLUME_COXETER_CANDIDATE';break
        chosen=None
        if step==0:chosen=proposed_rotation(0,74,8,R,Q,center)
        else:
            pairs=sorted(((i,j) for a,b in edges for i,j in ((a,b),(b,a))),key=lambda p:(-degree[p[0]],degree[p[1]],p))
            for i,j in pairs:
                for m in (2,3,4,6,8):
                    chosen=proposed_rotation(i,j,m,R,Q,center)
                    if chosen is not None:break
                if chosen is not None:break
        if chosen is None:out['status']='NO_INWARD_FIELD_ROTATION_PRESERVING_CENTER';break
        v,row=chosen;k=len(R);R=np.concatenate((R,v[None,:,:]),axis=0);row['new_index']=k
        check_terms(R[row['replaced_index']],[[k,*row['new_coefficient']],[row['retained_index'],*row['retained_coefficient']]],R)
        active.remove(row['replaced_index']);active.append(k);row['redundancy_proofs']=[]
        for i in list(active):
            terms=cone_terms(i,active,R)
            if terms is not None:row['redundancy_proofs'].append({'removed_index':i,'terms':terms});active.remove(i)
        row['retained_inequality_count']=len(active);out['steps'].append(row)
        print('Cut',row['replaced_index'],'toward',row['retained_index'],'m',row['target_angle_m'],'removed redundant',len(row['redundancy_proofs']),flush=True)
        out['normal_pairs']=R.tolist();out['final_active_indices']=active
        OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    else:out['status']='BOUNDED_INWARD_STEP_LIMIT'
    edges,degree=conflicts(active,R,Q);out['final_conflict_count']=len(edges)
    out['normal_pairs']=R.tolist();out['final_active_indices']=active
    out['surviving_non_coxeter_ridge']=surviving_ridge(active,R,Q,W[0])
    out['remaining_conflict_examples']=[list(e) for e in edges[:20]]
    OUTPUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],'steps',len(out['steps']),'walls',len(active),'conflicts',len(edges),flush=True)

if __name__=='__main__':main()
