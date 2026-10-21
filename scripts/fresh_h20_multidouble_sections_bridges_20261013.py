#!/usr/bin/env python3
"""New cusp quotient directions on the eight-chamber ambient union.

Exact cusp-kernel coordinates construct a section interior even when the
orthogonal projection of the original center fails. Numerical target
proposals are bounded and never certify completeness.
"""
import json
from itertools import combinations
import numpy as np
import sympy as sp
from scipy.optimize import linprog
import fresh_h20_multidouble_sections_20261013 as M

D=M.D;O=M.O;P=M.P;t=M.t;old=M.old;ROOT=M.ROOT
OUT=ROOT/'data/enumerations/fresh_h20_multidouble_sections_bridges_20261013.json'
PROPOSALS=ROOT/'data/enumerations/fresh_h20_multidouble_sections_bridge_proposals_20261013.json'


def section(Q,R,H,h,row):
    e0=np.array(list(map(sp.Rational,row['base_normal'])),dtype=object)
    u=np.array(row['cusp_vector'],dtype=object);q=P(e0@Q@e0)
    e=[P(a+t*b)for a,b in zip(e0,u)]
    v=[P(a+t*b)for a,b in zip(R@Q@e0,R@Q@u)]
    return (e,q,v,lambda a,b:q*int(H[a,b])-v[a]*v[b]),e0,u


def cusp_center(Q,W,R,h,pencil,sd,F):
    e,q,v,gr=sd;u=np.array(pencil['cusp_vector'],dtype=object)
    incident=[i for i,w in enumerate(W)if w@Q@u==0]
    assert not pencil['central_cusp_link']['inactive_incident_walls_one_based']
    e0=np.array(list(map(sp.Rational,pencil['base_normal'])),dtype=object)
    rhs={};multipliers=[]
    for comp in pencil['central_cusp_link']['components']:
        ids=[i-1 for i in comp['walls_one_based']];marks=comp['positive_kernel_marks']
        mu=P(q.as_expr()*sp.Rational(comp['null_sum_cusp_multiple'])-t*sum(m*(W[i]@Q@e0)for m,i in zip(marks,ids)))
        assert F.sign(mu)>0;multipliers.append(O.coeffs(F.rem(mu)))
        for i in ids:rhs[i]=sp.cancel(-mu.as_expr()/(q.as_expr()*sum(marks)))
    matrix=sp.Matrix([list(W[i]@Q)for i in incident]+[list(u@Q),list(np.array(h,dtype=object)@Q)])
    values=sp.Matrix([rhs[i]for i in incident]+[-1,0])
    solution,parameters=matrix.gauss_jordan_solve(values);assert len(parameters)==0
    C0=[F.rem(P(x))for x in solution]
    assert F.rem(old.pdot(Q,C0,e,F)).is_zero
    assert F.rem(old.pdot(Q,C0,[P(x)for x in u],F))==P(-1)
    shift=1
    while True:
        C=[F.rem(p+shift*int(x))for p,x in zip(C0,u)]
        if (F.sign(old.pdot(Q,C,C,F))<0 and
            F.sign(old.pdot(Q,C,[P(int(x))for x in h],F))<0 and
            all(F.sign(old.pdot(Q,C,[P(int(x))for x in w],F))<0 for w in R)):break
        shift*=2;assert shift<2**50
    return C,{'unshifted_point':[O.coeffs(x)for x in C0],'shift':shift,
               'cusp_component_null_multiples':multipliers}


def numeric_pairs(Q,R,H,e,C,value,actual_only=True):
    ef=np.array([float(p.eval(value))for p in e]);cf=np.array([float(p.eval(value))for p in C])
    Qf=np.asarray(Q,float);Rf=np.asarray(R,float);q=ef@Qf@ef;v=Rf@Qf@ef
    HG=np.asarray(H,float)-np.outer(v,v)/q;cp=Rf@Qf@cf
    for a,b in combinations(range(len(R)),2):
        aa,bb,ab=HG[a,a],HG[b,b],HG[a,b];det=aa*bb-ab*ab
        if min(aa,bb,det)<1e-8:continue
        ratio=ab*ab/(aa*bb)
        if ab<=0 and(ratio>=.749999 or min(abs(ratio-z)for z in(0,.25,.5,(3+5**.5)/8))<1e-8):continue
        if actual_only:
            A=bb*cp[a]-ab*cp[b];B=aa*cp[b]-ab*cp[a]
            vals=det*cp-A*HG[:,a]-B*HG[:,b]
            if max(vals)>1e-7*max(1.,max(abs(vals))):continue
        yield a,b


def ridge_lp(Q,R,h,sd,F,a,b,forbidden=False):
    e,q,v,gr=sd;value=float((F.lo+F.hi)/2)
    Qf=np.asarray(Q,float);Rf=np.asarray(R,float);hf=np.asarray(h,float)
    ef=np.array([float(p.eval(value))for p in e]);other=[i for i in range(len(R))if i not in(a,b)]
    eq=np.array([hf@Qf,ef@Qf,Rf[a]@Qf,Rf[b]@Qf])
    lp=linprog(np.r_[np.zeros(22),-1.],A_ub=np.c_[Rf[other]@Qf,np.ones(len(other))],b_ub=np.zeros(len(other)),
               A_eq=np.c_[eq,np.zeros(4)],b_eq=[-1,0,0,0],bounds=[(None,None)]*22+[(0,None)],method='highs')
    if not lp.success or lp.x[-1]<1e-8:return None
    for denominator in (10**4,10**6,10**9):
        raw=[P(sp.Rational(str(float(x))).limit_denominator(denominator))for x in lp.x[:-1]]
        ep=old.pdot(Q,raw,e,F);C=[F.rem(q*x-ep*y)for x,y in zip(raw,e)]
        proof=D.ridge(Q,R,h,sd,C,F,a,b,forbidden)
        if proof is not None:return proof
    return None


def run():
    Q,W,G,h=O.base.load();Q=np.asarray(Q,dtype=np.int64);W=np.asarray(W,dtype=np.int64);G=np.asarray(G,dtype=np.int64)
    ambient=json.loads(M.OUT.read_text());R=np.array(ambient['runs'][1]['boundary_normals'],dtype=np.int64);H=R@Q@R.T
    scout=json.loads(PROPOSALS.read_text());rows=[]
    for group in scout:
        pencil=group['pencil'];sd,e0,u=section(Q,R,H,h,pencil)
        assert e0@Q@W[0]!=0
        assert set(map(tuple,R[(R@Q@u)==0]))==set(map(tuple,W[(W@Q@u)==0]))
        selected=sorted(group['numeric_proposals'],key=lambda x:(abs(x['approx']),x['label'],x['pair']))[:2]
        candidates=[]
        for proposal in selected:
            factor=sp.Poly.from_list(list(map(sp.Rational,proposal['polynomial'])),t,domain=sp.QQ)
            F=O.Field(factor,*map(sp.Rational,proposal['interval']));a,b=[j-1 for j in proposal['pair']];label=proposal['label']
            e,q,v,gr=sd;aa,ab,bb=gr(a,a),gr(a,b),gr(b,b);N=ab*ab;den=aa*bb
            equations={2:ab,3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}
            assert F.sign(equations[label])==0 and min(F.sign(aa),F.sign(bb))>0 and F.sign(ab)<=0
            if label==5:assert F.sign(2*N-den)>0
            C,centerproof=cusp_center(Q,W,R,h,pencil,sd,F)
            target=D.ridge(Q,R,h,sd,C,F,a,b)
            if target is None:target=ridge_lp(Q,R,h,sd,F,a,b)
            bad=None;value=float((F.lo+F.hi)/2)
            for x,y in numeric_pairs(Q,R,H,e,C,value):
                bad=D.ridge(Q,R,h,sd,C,F,x,y,True)
                if bad is not None:break
            if bad is None:
                for count,(x,y)in enumerate(numeric_pairs(Q,R,H,e,C,value,False)):
                    if count>=20:break
                    bad=ridge_lp(Q,R,h,sd,F,x,y,True)
                    if bad is not None:break
            candidate={**proposal,'normal':[O.coeffs(F.rem(x))for x in e],
                       'strict_center':[O.coeffs(x)for x in C],'cusp_center':centerproof,
                       'target_ridge':target,'actual_forbidden_ridge':bad,
                       'status':'NON_COXETER_ACTUAL_RIDGE'if bad else'OPEN_COXETER_CANDIDATE'}
            candidates.append(candidate);print(group['option'],proposal['approx'],'target',target is not None,candidate['status'],flush=True)
        rows.append({'option_index_in_prior_census':group['option'],'pencil':pencil,
                     'central_new_mirror_pairs':group['pairs'],
                     'numeric_cusp_feasible_root_proposals':len(group['numeric_proposals']),
                     'selected_rule':'Two smallest absolute parameters among declared numerical proposals, label then pair as tie breakers.',
                     'candidates':candidates})
        OUT.write_text(json.dumps({'complete':False,'groups':rows},indent=2)+'\n')
    data={'schema':'new-cusp-directions-in-eight-chamber-union-v1','complete':True,
          'ambient_source':str(M.OUT.relative_to(ROOT)),'ambient_chambers':8,
          'proposal_source':str(PROPOSALS.relative_to(ROOT)),
          'cusp_source':str(O.CUSP_SOURCE.relative_to(ROOT)),'groups':rows,
          'scope':'Four previously unmaterialized cusp quotient directions with nonzero r1 pairing; '
                  'two selected algebraic parameters per direction from12 declared new-mirror angle pairs '
                  'and labels2..6. No exhaustion of sections or target roots. Exact accepted interiors, '
                  'whole cusp links and actual ridge claims only.'}
    OUT.write_text(json.dumps(data,indent=2)+'\n')


if __name__=='__main__':run()
