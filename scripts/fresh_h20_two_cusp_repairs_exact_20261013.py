#!/usr/bin/env python3
"""Exact geometric replay of selected simultaneous three-cusp angle repairs."""
import json,itertools,hashlib
import numpy as np
import sympy as sp
import fresh_h20_two_cusp_repairs_20261013 as scout
import fresh_h20_joint_cusp_irrational_20261013 as prior
import fresh_h20_reflected_sections_double_20261013 as D

O=scout.O;ROOT=O.ROOT;P=O.P;t=O.t
OUT=ROOT/'data/enumerations/fresh_h20_two_cusp_repairs_exact_20261013.json'


def exact_parameters(row,profiles):
    A,B=[profiles[i]for i in row['profile_indices']]
    da=sp.Rational(A['norm'])-sp.Rational(str(A['rho']))/sp.cos(sp.pi/row['target_labels'][0])**2
    db=sp.Rational(B['norm'])-sp.Rational(str(B['rho']))/sp.cos(sp.pi/row['target_labels'][1])**2
    da,db=sp.simplify(da),sp.simplify(db)
    x,y=sp.symbols('x y');av=sum(sp.Integer(z)*v for z,v in zip(A['v'],(1,x,y)))
    bv=sum(sp.Integer(z)*v for z,v in zip(B['v'],(1,x,y)))
    if da==db==0:
        sol=sp.solve((av,bv),(x,y));return sol[x],sol[y]
    lines=[av]if da==0 else[bv]if db==0 else[av-sgn*sp.sqrt(da/db)*bv for sgn in(-1,1)]
    candidates=[]
    for line in lines:
        if sp.diff(line,y)!=0:
            yy=sp.solve(line,y)[0]
            eq=sp.factor(((av**2-da*336*(x+y-x*y))if da!=0 else(bv**2-db*336*(x+y-x*y))).subs(y,yy))
            for xx in sp.solve(eq,x):candidates.append((sp.simplify(xx),sp.simplify(yy.subs(x,xx))))
        else:
            xx=sp.solve(line,x)[0]
            eq=sp.factor(((av**2-da*336*(x+y-x*y))if da!=0 else(bv**2-db*336*(x+y-x*y))).subs(x,xx))
            for yy in sp.solve(eq,y):candidates.append((sp.simplify(xx),sp.simplify(yy)))
    ox,oy=row['approximate_parameters'];ranked=[]
    for xx,yy in candidates:
        if xx.is_real is False or yy.is_real is False:continue
        try:dx,dy=float(sp.N(xx,20)),float(sp.N(yy,20))
        except(TypeError,ValueError):continue
        ranked.append((abs(dx-ox)+abs(dy-oy),xx,yy))
    error,xx,yy=min(ranked,key=lambda z:z[0]);assert error<1e-6
    assert sp.simplify(av.subs({x:xx,y:yy})**2-da*336*(xx+yy-xx*yy))==0
    assert sp.simplify(bv.subs({x:xx,y:yy})**2-db*336*(xx+yy-xx*yy))==0
    return xx,yy


def number_field(x,y):
    theta=sp.to_number_field([x,y]);poly=sp.Poly.from_list(theta.minpoly.all_coeffs(),t,domain=sp.QQ).monic()
    value=float(sp.N(theta.as_expr(),30))
    interval=next((a,b)for(a,b),m in poly.intervals(eps=sp.Rational(1,10**18))if float(a)-1e-12<=value<=float(b)+1e-12)
    lo,hi=interval
    if lo==hi:lo-=sp.Rational(1,10**20);hi+=sp.Rational(1,10**20)
    field=O.Field(poly,lo,hi)
    convert=lambda v:sp.Poly.from_list(sp.to_number_field(v,theta).coeffs(),t,domain=sp.QQ)
    return field,convert(x),convert(y)


def ridge_lp(Q,W,h,sd,field,a,b,forbidden=False):
    """Alternate ridge point; numerical proposal followed by exact projection."""
    from scipy.optimize import linprog
    e,q,v,gr=sd;value=float((field.lo+field.hi)/2)
    ef=np.array([float(p.eval(value))for p in e]);Qf=np.array(Q,float);Wf=np.array(W,float);hf=np.array(h,float)
    aa,bb,ab=map(field.rem,(gr(a,a),gr(b,b),gr(a,b)));det=field.rem(aa*bb-ab*ab)
    if min(field.sign(aa),field.sign(bb),field.sign(det))<=0:return None
    excluded={a,b}
    for j in range(210):
        al=field.rem(bb*gr(j,a)-ab*gr(j,b));be=field.rem(aa*gr(j,b)-ab*gr(j,a))
        if min(field.sign(al),field.sign(be))<0:continue
        if not field.rem(det*gr(j,j)-al*gr(j,a)-be*gr(j,b)).is_zero:continue
        if all(field.rem(det*(q*int(W[j,d])-v[j]*e[d])-al*(q*int(W[a,d])-v[a]*e[d])-be*(q*int(W[b,d])-v[b]*e[d])).is_zero for d in range(22)):
            excluded.add(j)
    active=[i for i in range(210)if i not in excluded]
    lp=linprog(np.r_[np.zeros(22),-1.],A_ub=np.c_[Wf[active]@Qf,np.ones(len(active))],b_ub=np.zeros(len(active)),
        A_eq=np.c_[np.array([hf@Qf,ef@Qf,Wf[a]@Qf,Wf[b]@Qf]),np.zeros(4)],b_eq=[-1.,0.,0.,0.],
        bounds=[(None,None)]*22+[(0,None)],method='highs')
    if not lp.success or lp.x[-1]<1e-8:return None
    for den in(1000,100000,10000000):
        raw=[P(sp.Rational(str(float(x))).limit_denominator(den))for x in lp.x[:-1]]
        pe=O.old.pdot(Q,raw,e,field);C=[field.rem(q*x-pe*y)for x,y in zip(raw,e)]
        proof=D.ridge(Q,W,h,sd,C,field,a,b,forbidden)
        if proof is not None:return proof
    return None


def star_diagnostic(Q,W,G,h):
    """The initially proposed finite star has an unavoidable obtuse cusp ridge."""
    import fresh_h20_joint_cusp_a2_20261013 as A2
    I=np.eye(22,dtype=object);gens=[I-np.outer(W[i],W[i]@Q)for i in(0,21)]
    group=[I];seen={tuple(I.flat)}
    for M in group:
        for S in gens:
            T=S@M;key=tuple(T.flat)
            if key not in seen:group.append(T);seen.add(key)
        assert len(group)<=6
    outside=[i for i in range(210)if i not in(0,21)]
    for M in group:
        assert np.array_equal(M.T@Q@M,Q)
        for i in outside:
            delta=M@W[i]-W[i];a=sp.Rational(2*(delta@Q@W[0])+delta@Q@W[21],3)
            b=sp.Rational(delta@Q@W[0]+2*(delta@Q@W[21]),3)
            assert a>=0 and b>=0 and np.array_equal(delta,a*W[0]+b*W[21])
    R=np.array(sorted({tuple(map(int,M@W[i]))for M in group for i in outside}),dtype=object)
    u=sum(W[[59,98,124]])*sp.Rational(1,2)
    assert W[0]@Q@u==0 and W[21]@Q@u==-1
    j=next(i for i in outside if G[0,i]==-1 and G[i,i]==2 and W[i]@Q@u==0
           and all(G[i,k]==0 for k in(59,98,124)))
    roots=[W[j],W[j]+W[0]];pair=[next(i for i,r in enumerate(R)if np.array_equal(r,v))for v in roots]
    assert [[r@Q@s for s in roots]for r in roots]==[[2,1],[1,2]]
    pp=A2.face(Q,W,G,h,sorted([0,j,59,98]));pm=A2.face(Q,W,G,h,sorted([0,j,124]))
    vp=R@Q@np.array(pp,dtype=object).reshape(-1);vm=R@Q@np.array(pm,dtype=object).reshape(-1)
    assert max(vp)<=0 and max(vm)<=0
    assert [i for i in range(len(R))if vp[i]==vm[i]==0]==sorted(pair)
    return {'seams':[1,22],'group_order':len(group),'boundary_count':len(R),
            'incident_seam':1,'nonincident_seam':22,'cusp_vector':list(map(int,u)),
            'neighbor_original_wall':j+1,'obtuse_boundary_normals':[list(map(int,r))for r in roots],
            'obtuse_gram':[[2,1],[1,2]],
            'positive_face_point':list(map(str,pp)),'negative_face_point':list(map(str,pm)),
            'scope':'For every s,t>0 in e=r60+s*r99-t*r125, positive interpolation of the two face points gives an actual 2pi/3 ridge in the six-star section. This closes this initial star, not other assemblies.'}


def run():
    Q,W,G,h=O.base.load();source=json.loads(scout.OUT.read_text());raw=json.loads(O.CUSP_SOURCE.read_text())['cusps']
    U=np.array([raw[i]['vector']for i in(0,3,4)],dtype=object)
    assert (U@Q@U.T).tolist()==[[0,-24,-2],[-24,0,-14],[-2,-14,0]]
    selected=[(i,r)for i,r in enumerate(source['proposals'])if any(m>=5 for m in r['target_labels'])]
    selected=[(0,source['proposals'][0])]+selected
    rows=[]
    for index,row in selected:
        x,y=exact_parameters(row,source['profiles']);field,xp,yp=number_field(x,y)
        e=[field.rem(7*int(U[0,d])-xp*int(U[1,d])-12*yp*int(U[2,d]))for d in range(22)]
        sd=prior.sd_from_e(Q,W,G,e,field);q=sd[1]
        assert field.rem(q-336*(xp+yp-xp*yp)).is_zero and field.sign(q)>0
        C=O.four.alternate_center(Q,W,h,sd,field);assert C is not None
        targets=[]
        for pair,label in zip(row['target_pairs'],row['target_labels']):
            a,b=[i-1 for i in pair];proof=D.ridge(Q,W,h,sd,C,field,a,b)
            if proof is None:proof=ridge_lp(Q,W,h,sd,field,a,b)
            aa,bb,ab=sd[3](a,a),sd[3](b,b),sd[3](a,b);N=ab*ab;den=aa*bb
            eq={3:4*N-den,4:2*N-den,5:16*N*N-12*N*den+den*den,6:4*N-3*den}[label]
            assert field.rem(eq).is_zero and field.sign(ab)<0
            if label==5:assert field.sign(2*N-den)>0
            targets.append({'label':label,'pair':pair,'actual_ridge':proof})
        a,b=[i-1 for i in row['bad_pair_proposal'][0]]
        bad=D.ridge(Q,W,h,sd,C,field,a,b,True)
        if bad is None:bad=ridge_lp(Q,W,h,sd,field,a,b,True)
        if bad is None:
            value=float((field.lo+field.hi)/2)
            ef=np.array([float(p.eval(value))for p in e]);cf=np.array([float(p.eval(value))for p in C])
            Qf=np.array(Q,float);Wf=np.array(W,float);vf=Wf@Qf@ef
            H=np.array(G,float)-np.outer(vf,vf)/(ef@Qf@ef);cp=Wf@Qf@cf
            for aa,bb in itertools.combinations(range(210),2):
                dd=H[aa,aa]*H[bb,bb];n=H[aa,bb]**2
                if min(H[aa,aa],H[bb,bb],dd-n)<=1e-8:continue
                if abs(H[aa,bb])<1e-8:continue
                if H[aa,bb]<0 and(n/dd>=.75-1e-7 or min(abs(n/dd-z)for z in(0,.25,.5,(3+5**.5)/8))<1e-7):continue
                if not scout.actual(H,cp,aa,bb):continue
                bad=D.ridge(Q,W,h,sd,C,field,aa,bb,True)
                if bad is not None:break
        assert bad is not None
        result={'scout_proposal':index,'parameter_expressions':[str(x),str(y)],
                'field_polynomial':O.coeffs(field.p),'field_interval':[str(field.lo),str(field.hi)],
                'parameters':[O.coeffs(xp),O.coeffs(yp)],'normal':[O.coeffs(p)for p in e],
                'normal_square':O.coeffs(q),'strict_center':[O.coeffs(p)for p in C],
                'repaired_actual_ridges':targets,'actual_forbidden_ridge':bad,
                'status':'FINITE_VOLUME_NON_COXETER_SECTION_WITH_TWO_ACTUAL_REPAIRS'if all(r['actual_ridge']for r in targets)else'FINITE_VOLUME_NON_COXETER_SECTION_TARGET_INCIDENCE_PARTLY_UNRESOLVED'}
        rows.append(result)
        print('EXACT',index,'degree',field.p.degree(),row['target_labels'],'actualtargets',[r['actual_ridge']is not None for r in targets],'bad',row['bad_pair_proposal'],flush=True)
    out={'status':'EXACT_TWO_PARAMETER_SECTION_BATCH_NO_COXETER_SURVIVOR',
         'ambient_source':str(O.base.SOURCE.relative_to(ROOT)),
         'cusp_source':str(O.CUSP_SOURCE.relative_to(ROOT)),
         'scout_source':str(scout.OUT.relative_to(ROOT)),
         'cusp_indices':[0,3,4],'cusp_vectors':[list(map(int,u))for u in U],
         'cusp_gram':[[0,-24,-2],[-24,0,-14],[-2,-14,0]],
         'common_wall_labels':source['common_wall_labels'],
         'selected_proposals':[i for i,r in selected],'cases':rows,
         'initial_star_diagnostic':star_diagnostic(Q,W,G,h),
         'scope':'Nine literal sections selected from the two-independent-parameter three-cusp plane. Every section has exact finite-volume inheritance and an actual forbidden ridge. Target angle equations are exact; actual incidence is claimed only where its certificate is present. No compactness claim, all-parameter exclusion, or exhaustive numerical-angle search claim.'}
    OUT.write_text(json.dumps(out,indent=2)+'\n');print('COMPLETE',len(rows),flush=True)


if __name__=='__main__':run()
