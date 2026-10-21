#!/usr/bin/env python3
"""Two independent angle repairs in a three-cusp normal plane.

Numerics select proposals only. Normals are E=7u0-x*u3-12y*u4;
their exact square is336(x+y-xy). Seventeen original mirrors stay
orthogonal to every normal. Two common-mirror angle constraints are
quadratic equations, solved by a linear square-root relation.
"""
import itertools,json
from collections import Counter
import numpy as np
from scipy.optimize import linprog
import fresh_h20_new_route_20261013 as O

ROOT=O.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_two_cusp_repairs_20261013_scout.json'
SAMPLES=((.5,.5),(3.5,7/6),(7,.5),(2,1),(1,1))
LABELS=(3,4,5,6)


def state(Q,W,G,h,B,x,y):
    e=np.array([1,x,y])@B;q=e@Q@e
    if q<=1e-7:return None
    v=W@Q@e;H=G-np.outer(v,v)/q
    lp=linprog(np.r_[np.zeros(22),-1.],A_ub=np.c_[W@Q,np.ones(210)],b_ub=np.zeros(210),
        A_eq=np.c_[np.array([h@Q,e@Q]),np.zeros(2)],b_eq=[-1.,0.],
        bounds=[(None,None)]*22+[(0,None)],method='highs')
    if not lp.success or lp.x[-1]<1e-7:return None
    C=lp.x[:-1]
    if C@Q@C>=-1e-9:return None
    return e,q,H,C,W@Q@C


def actual(H,cp,i,j):
    a,b,c=H[i,i],H[j,j],H[i,j];det=a*b-c*c
    if min(a,b,det)<=1e-8:return False
    A=b*cp[i]-c*cp[j];B=a*cp[j]-c*cp[i]
    values=det*cp-A*H[:,i]-B*H[:,j]
    return max(values)<=1e-8*max(1,max(abs(values)))


def intersection(A,da,B,db):
    # A,B are affine-linear coefficient triples; q=336(x+y-xy).
    if da<1e-10 and db<1e-10:
        M=np.array([A[1:],B[1:]],float)
        if abs(np.linalg.det(M))<1e-9:return
        x,y=np.linalg.solve(M,-np.array([A[0],B[0]],float))
        yield float(x),float(y),0
        return
    if da<1e-10:lines=[np.array(A,float)]
    elif db<1e-10:lines=[np.array(B,float)]
    else:lines=[np.array(A,float)-sgn*np.sqrt(da/db)*np.array(B,float)for sgn in(-1,1)]
    V=np.array(A if da>1e-10 else B,float);delta=da if da>1e-10 else db
    for ni,L in enumerate(lines):
        if abs(L[2])>1e-10:
            aa,bb=-L[0]/L[2],-L[1]/L[2]
            c,d=V[0]+V[2]*aa,V[1]+V[2]*bb
            roots=np.roots([d*d+336*delta*bb,2*c*d-336*delta*(1+bb-aa),c*c-336*delta*aa])
            for z in roots:
                if abs(z.imag)<1e-8:yield float(z.real),float(aa+bb*z.real),ni
        elif abs(L[1])>1e-10:
            x=-L[0]/L[1];c=V[0]+V[1]*x;d=V[2]
            for z in np.roots([d*d,2*c*d-336*delta*(1-x),c*c-336*delta*x]):
                if abs(z.imag)<1e-8:yield float(x),float(z.real),ni


def run():
    Q,W,G,h=O.base.load();Q=np.array(Q,float);W=np.array(W,float);G=np.array(G,float);h=np.array(h,float)
    raw=json.loads(O.CUSP_SOURCE.read_text())['cusps']
    U=np.array([raw[i]['vector']for i in(0,3,4)],float);B=np.array([7*U[0],-U[1],-12*U[2]])
    V=np.rint(W@Q@B.T).astype(int);common=np.flatnonzero(np.all(V==0,axis=1))
    profiles={}
    for x,y in SAMPLES:
        sd=state(Q,W,G,h,B,x,y);assert sd is not None
        _,_,H,C,cp=sd
        for i in common:
            for j in range(210):
                if i==j or not G[i,j]or not any(V[j]):continue
                if not actual(H,cp,i,j):continue
                key=tuple(map(int,V[j]))+(int(G[j,j]),float(G[i,j]**2/G[i,i]))
                profiles.setdefault(key,set()).add((int(i),j))
    profiles=[{'v':list(key[:3]),'norm':key[3],'rho':key[4],
               'pairs':[[i+1,j+1]for i,j in sorted(pairs)]}for key,pairs in sorted(profiles.items())]
    counts=Counter();survivors=[];seen=set()
    for ia,ib in itertools.combinations(range(len(profiles)),2):
        A,Bp=profiles[ia],profiles[ib]
        for ma,mb in itertools.product(LABELS,repeat=2):
            da=A['norm']-A['rho']/np.cos(np.pi/ma)**2
            db=Bp['norm']-Bp['rho']/np.cos(np.pi/mb)**2
            if min(da,db)<-1e-8:continue
            counts['two_equation_templates']+=1
            for x,y,branch in intersection(A['v'],da,Bp['v'],db):
                if not np.isfinite(x+y)or min(x,y)<=1e-7:continue
                key=(round(x,8),round(y,8))
                if key in seen:continue
                seen.add(key);counts['distinct_positive_numeric_roots']+=1
                sd=state(Q,W,G,h,B,x,y)
                if sd is None:continue
                e,q,H,C,cp=sd;counts['strict_numeric_sections']+=1
                pairs=[]
                for profile in(A,Bp):
                    pair=next((pair for pair in profile['pairs']if actual(H,cp,pair[0]-1,pair[1]-1)),None)
                    if pair is None:break
                    pairs.append(pair)
                if len(pairs)<2:continue
                counts['two_actual_target_proposals']+=1
                bad=[]
                for i,j in itertools.combinations(range(210),2):
                    aa,bb,ab=H[i,i],H[j,j],H[i,j];det=aa*bb-ab*ab
                    if min(aa,bb,det)<1e-8:continue
                    ratio=ab*ab/(aa*bb)
                    if abs(ab)<1e-8 or(ab<=0 and(ratio>=.75-1e-7 or min(abs(ratio-z)for z in(0,.25,.5,(3+5**.5)/8))<1e-7)):continue
                    if actual(H,cp,i,j):bad.append([i+1,j+1]);break
                row={'profile_indices':[ia,ib],'target_labels':[ma,mb],'linear_branch':branch,
                     'approximate_parameters':[x,y],'target_pairs':pairs,'bad_pair_proposal':bad,
                     'status':'NUMERICAL_PROPOSAL_REQUIRES_EXACT_REPLAY'}
                survivors.append(row)
                print('two-target proposal',pairs,[ma,mb],x,y,'bad',bad,flush=True)
    result={'status':'NUMERICAL_TWO_PARAMETER_PROPOSALS_ONLY','cusp_indices':[0,3,4],
            'normal_formula':'7u0-x*u3-12y*u4','normal_square':'336(x+y-x*y)',
            'common_wall_labels':[int(i+1)for i in common],'samples':SAMPLES,'labels':LABELS,
            'profiles':profiles,'counts':dict(counts),'proposals':survivors,
            'scope':'Finite analytically generated two-angle templates; numerical interiors and ridge tests select proposals only. Not a proof of completeness, nonexistence, or any accepted Coxeter polytope.'}
    OUT.write_text(json.dumps(result,indent=2)+'\n')
    print('DONE',len(profiles),dict(counts),flush=True)


if __name__=='__main__':run()
