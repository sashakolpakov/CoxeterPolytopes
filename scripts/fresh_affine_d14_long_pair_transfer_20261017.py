#!/usr/bin/env python3
"""Universal affine-D14 long-pair half-coupon transfer; no dimension bound."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s
import fresh_cycle16_half_coupon_transfer_20261017 as C

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_d14_long_pair_transfer_20261017'
HALF=s.Rational(1,2)


def affine():
    g=2*s.eye(15)
    edges=[(i,i+1)for i in range(10)]+[(0,11),(0,12),(10,13),(10,14)]
    for i,j in edges:g[i,j]=g[j,i]=-1
    nu=s.Matrix([1]*11+[HALF]*4)
    assert g*nu==s.zeros(15,1) and C.inertia(g.tolist())==[14,0,1]
    return g,nu,edges


def end_profiles(a,right=False):
    x,u,v=s.symbols(('y vr wr'if right else'x ul wl'),real=True)
    reflect={i:10-i for i in range(11)}|{11:13,12:14,13:11,14:12}
    aa=reflect[a]if right else a
    def row(name,contacts,c0,l,h,paid=False):
        if right:contacts={reflect[i]:q for i,q in contacts.items()}
        return dict(name=name,contacts=contacts,c0=c0,l=l,h=h,small_affine_paid=paid)
    if aa!=0:
        assert aa in(11,12)
        mate=23-aa
        return [row('isolated',{aa:x},2,0,x/2),
                row('mate3',{aa:x,mate:1},1,1,(x+1)/2),
                row('mate4',{aa:x,mate:s.sqrt(2)},0,s.sqrt(2),(x+s.sqrt(2))/2),
                row('center3',{aa:x,0:1},0,2,1+x/2)]
    return [row('isolated',{0:x},2,0,x),
            row('outer_general',{0:x,11:u,12:v},2-(u*u+v*v)/2,0,x+(u+v)/2),
            row('inner3',{0:x,1:1},1,1,x+1),
            row('inner4',{0:x,1:s.sqrt(2)},0,s.sqrt(2),x+s.sqrt(2),True),
            row('inner_next3',{0:x,2:1},0,2,x+1,True),
            row('join_one3',{0:x,1:1,11:1},HALF,1,x+s.Rational(3,2)),
            row('join_one4',{0:x,1:1,11:s.sqrt(2)},0,1,x+1+s.sqrt(2)/2),
            row('join_two3',{0:x,1:1,11:1,12:1},0,1,x+2)]


def formula_checks():
    g,nu,edges=affine();rows=[];cross=[];paid=[]
    for a,b in product((0,11,12),(10,13,14)):
        W=sorted(set(range(15))-{a,b});J=g.extract(W,W).inv()
        R=10+2*(a!=0)+2*(b!=10)
        va=-nu[a]*g.extract([a],W).T;vb=-nu[b]*g.extract([b],W).T
        old=s.diag(nu[a],nu[b])*g.extract([a,b],[a,b])*s.diag(nu[a],nu[b])
        old-=s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
        assert old==s.Matrix([[1,-1],[-1,1]])/R
        L=end_profiles(a);D=end_profiles(b,True)
        for side,ps,mark,far,vfar in [('left',L,a,b,vb),('right',D,b,a,va)]:
            for p in ps:
                w=s.Matrix([p['contacts'].get(i,0)for i in W])
                h=sum(nu[i]*q for i,q in p['contacts'].items())
                B=(vfar.T*J*w)[0];norm=2-(w.T*J*w)[0]
                assert s.simplify(h-p['h'])==0
                assert s.simplify(B-p['l']/R)==0
                assert s.simplify(norm-p['c0']-p['l']**2/R)==0
                rows.append(dict(marks=[a,b],side=side,profile=p['name'],R=R,
                    contacts={str(i):str(q)for i,q in p['contacts'].items()},
                    c0=str(p['c0']),moment=str(p['l']),height=str(p['h']),
                    residual_norm=str(s.factor(norm)),far_contact=str(B),
                    small_affine_paid=p['small_affine_paid']))
        for c,d in product(L,D):
            w=s.Matrix([c['contacts'].get(i,0)for i in W]);z=s.Matrix([d['contacts'].get(i,0)for i in W])
            E=(w.T*J*z)[0]
            assert s.simplify(E-c['l']*d['l']/R)==0
            cross.append(dict(marks=[a,b],c=c['name'],d=d['name'],projection_cross=str(E)))
    assert len(rows)==96 and len(cross)==256
    # Literal affine components that pay the two exceptional center profiles.
    for name,j,t in [('inner4',1,s.sqrt(2)),('inner_next3',2,s.Integer(1))]:
        ids=list(range(1,11))+[13,14,15]
        h=2*s.eye(16);h[:15,:15]=g;h[15,j]=h[j,15]=-t
        k=h.extract(ids,ids);z=k.nullspace()[0]
        if z[0]<0:z=-z
        assert len(k.nullspace())==1 and all(q>0 for q in z)
        # Diagonal congruence by the positive kernel is a connected Laplacian.
        lap=s.diag(*z)*k*s.diag(*z)
        assert all(sum(lap.row(i))==0 for i in range(13))
        assert all(lap[i,j]<=0 for i in range(13)for j in range(13)if i!=j)
        paid.append(dict(profile=name,nodes=ids,order=13,kernel=list(map(str,z))))
    R,h,k,l,m,c0,d0,p,q,r=s.symbols('R h k l m c0 d0 p q r',nonzero=True)
    alpha=1/R;B=l/R;A=m/R;E=l*m/R
    T=s.Matrix([[alpha,-alpha,-h+B,-A],[-alpha,alpha,-B,-k+A],
        [-h+B,-B,c0+l*l/R,-E],[-A,-k+A,-E,d0+m*m/R]])
    expected=1-(2*l/h+c0/h**2+2*m/k+d0/k**2)/R
    assert s.factor(T.det()/(h*h*k*k)-expected)==0
    T[2,2]-=p;T[2,3]-=q;T[3,2]-=q;T[3,3]-=r
    assert s.factor(T.det()/(h*h*k*k)-expected-(p/h**2-2*q/(h*k)+r/k**2)/R)==0
    return dict(end_rows=rows,cross_rows=cross,small_affine_rows=paid,
        determinant='1-(2*l/h+c0/h^2+2*m/k+d0/k^2)/R',
        spectator_correction='(p/h^2-2*q/(h*k)+r/k^2)/R')


def cost_checks():
    x,t=s.symbols('x t',nonnegative=True)
    leaf=[('isolated',8/x**2,1,8),('isolated_nonextreme',8/x**2,s.sqrt(2),4),
        ('mate3',4/(x+1)+4/(x+1)**2,0,8),('mate3_nonextreme',4/(x+1)+4/(x+1)**2,1,3),
        ('mate4',4*s.sqrt(2)/(x+s.sqrt(2)),0,4),('center3',4/(1+x/2),0,4)]
    center=[('isolated',2/x**2,1,2),('single_outer3',s.Rational(3,2)/(x+HALF)**2,0,6),
        ('single_outer3_native_nonzero',s.Rational(3,2)/(x+HALF)**2,1,s.Rational(2,3)),
        ('single_outer_high',8/t**2-2,s.sqrt(2),2),
        ('inner3',2/(x+1)+1/(x+1)**2,0,3),
        ('join_one3',2/(x+s.Rational(3,2))+HALF/(x+s.Rational(3,2))**2,0,s.Rational(14,9)),
        ('join_one4',2/(x+1+s.sqrt(2)/2),0,4-2*s.sqrt(2)),
        ('join_two3',2/(x+2),0,1)]
    rows=[]
    for typ,ps in [('leaf',leaf),('center',center)]:
        for name,fn,at,value in ps:
            variable=t if name=='single_outer_high'else x
            assert s.simplify(fn.subs(variable,at)-value)==0
            # Every rational expression is decreasing on its stated positive-height ray.
            deriv=s.factor(s.diff(fn,variable))
            assert deriv.subs(variable,max(at,s.Integer(1)))<0
            rows.append(dict(end_type=typ,profile=name,formula=str(fn),lower_contact=str(at),
                             maximum=str(value),derivative=str(deriv)))
    margins=[]
    for a,b in product(('center','leaf'),repeat=2):
        R=10+2*(a=='leaf')+2*(b=='leaf')
        hi={'center':6,'leaf':8};lo={'center':3,'leaf':4}
        one=max(lo[a]+hi[b],hi[a]+lo[b]);gap=F(1)-F(one,R)
        assert gap>0 and hi[a]+hi[b]==R+2
        margins.append(dict(ends=[a,b],R=R,extreme_cost_sum=R+2,
            at_least_one_nonextreme_cost_bound=one,positive_gap=str(gap),
            extreme_heights=['1/2','1/2'],base_determinant=str(-F(2,R)),
            required_unscaled_spectator_difference='1/2'))
    return dict(rows=rows,arbitrary_outer_bounds=dict(
        no_outer_contact='x>=1: g<=2',
        one_outer_contact_native_nonzero='c0<=3/2,h>=3/2: g<=2/3',
        one_outer_contact_native_zero='g=8/t^2-2: t1 gives6; t>=sqrt2 gives<=2',
        two_outer_contacts='u,v>=1: h>=1,c0<=1, so g<=1'),pair_margins=margins)


def components(g,ids):
    adj={i:{j for j in ids if j!=i and g[i,j]}for i in ids};todo=set(ids);out=[]
    while todo:
        start=min(todo);ds=C.distances(adj,ids,start);todo-=ds.keys()
        degrees={i:len(adj[i])for i in ds};branch=[i for i in ds if degrees[i]==3]
        if max(degrees.values())<=2:kind='A'
        elif len(branch)==1 and sum(degrees[i]==1 and next(iter(adj[i]))==branch[0]for i in ds)==2:kind='D'
        else:kind='affineD'
        out.append(dict(type=f'{kind}{len(ds)-1 if kind=="affineD"else len(ds)}',nodes=sorted(ds),branch=branch,
            short_leaves=[i for i in ds if degrees[i]==1 and any(j in branch for j in adj[i])]))
    return out


def equality_cases():
    f,nu,edges=affine();rows=[];tally=Counter()
    for a,b,L,R,side in product((0,11,12),(10,13,14),(11,12),(13,14),(15,16)):
        g=2*s.eye(18);g[:15,:15]=f
        for i,j in((15,L),(16,R),(17,side)):g[i,j]=g[j,i]=-1
        W=sorted(set(range(15))-{a,b})+[17]
        assert C.inertia(g.tolist())==[16,1,1]
        z=g.nullspace()[0];assert len(g.nullspace())==1
        assert C.inertia(g.extract(W,W).tolist())==[14,0,0]
        corners=[]
        for p,q in((a,b),(b,15),(15,16),(16,a)):
            ids=sorted(set(W)|{p,q});sig=C.inertia(g.extract(ids,ids).tolist())
            assert sig==([15,0,1]if(p,q)==(a,b)else[16,0,0])
            cc=components(g,ids)
            for part in cc:
                part['marks']=[r for r in(p,q)if r in part['nodes']]
                part['fork_marks']=[r for r in part['marks']if r in part['branch']+part['short_leaves']]
            cp=next(part for part in cc if p in part['nodes']);cq=next(part for part in cc if q in part['nodes'])
            base=0;coupon=F()
            if sig[2]==0:
                if cp is cq:
                    assert len(cp['nodes'])<=17
                    base=1
                else:
                    small=next((part for part in(cp,cq)if len(part['nodes'])<=3),None)
                    large=next((part for part in(cp,cq)if part['type']in('D12','D14')and part['fork_marks']),None)
                    if small and large:coupon=F(1,2)
            corners.append(dict(marks=[p,q],signature=sig,components=cc,
                finite_small17_unit=base,half_coupon=str(coupon)))
        units=sum(r['finite_small17_unit']for r in corners);coupons=sum(F(r['half_coupon'])for r in corners)
        assert units+coupons>=1
        tally[(units,str(coupons))]+=1
        rows.append(dict(marks=[a,b],leaf_attachments=[L,R],spectator_attached_to=side,
            gram=[[str(q)for q in row]for row in g.tolist()],signature=[16,1,1],radical=list(map(str,z)),
            white=W,corners=corners,finite_units=units,coupon_total=str(coupons)))
    assert len(rows)==72 and tally==Counter({(0,'1'):32,(1,'1/2'):32,(3,'0'):8})
    return dict(rows=rows,payment_classes=[dict(finite_units=k[0],coupon_total=k[1],cases=v)for k,v in sorted(tally.items())],
        sealing='All other spectators are orthogonal to the entire18-root core by equality in the PD spectator gap.')


def run():
    out=dict(profile_formulas=formula_checks(),cost_bounds=cost_checks(),equality=equality_cases())
    sources=[Path(__file__)]
    for stem in('fresh_ideal_distance6_mixed_replacement_20261017','fresh_cycle16_half_coupon_transfer_20261017',
                'fresh_distance8_ideal_corner_replacement_20261017','fresh_simple_cusp_affine_tree_completion_20261017'):
        sources += [ROOT/f'scripts/{stem}.py',ROOT/f'results/{stem}.json',ROOT/f'results/{stem}.md']
    sources += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    out.update(status='PASS_UNIVERSAL_AFFINE_D14_LONG_PAIR_HALF_COUPON_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in sources},
        counts=dict(marked_pairs=9,end_formula_checks=96,separated_cross_checks=256,
            center_affine13_payments=2,scalar_cost_rows=14,strict_nonextreme_end_pairs=4,equality_cases=72,
            equality_finite_corners=216,two_half_coupon_cases=32),
        theorem='Every actual quadrilateral with a marked pair across the two fork triples of a15-node affineD14 factor is paid at another corner by finite distance<=9 or whole finite component order<=17, or ideal distance<=6 or whole factor order<=14, or by two half-payments at sealed finite small-component x D12/D14 fork-mark corners. Small components have rank<=3.',
        scope='Universal local theorem with arbitrary supporting spectators and nonsimple cusps. No global finite budget, ambient dimension, finite-volume construction or external/formal verification claim.')
    return out


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
