#!/usr/bin/env python3
"""Universal nine-long-pair transfer for the14-node affine D13 factor."""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import sympy as s
import fresh_affine_d14_long_pair_transfer_20261017 as M

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_affine_d13_long_pair_transfer_20261017'
C=M.C

def affine():
    g=2*s.eye(14)
    edges=[(i,i+1)for i in range(9)]+[(0,10),(0,11),(9,12),(9,13)]
    for i,j in edges:g[i,j]=g[j,i]=-1
    nu=s.Matrix([1]*10+[s.Rational(1,2)]*4)
    assert g*nu==s.zeros(14,1) and C.inertia(g.tolist())==[13,0,1]
    return g,nu

def profiles(a,right=False):
    # Contract one unused middle spine edge in the frozen all-rank table.
    mapping={i:i for i in range(5)}|{i:i-1 for i in range(6,15)}
    inv={v:k for k,v in mapping.items()}
    out=[]
    for p in M.end_profiles(inv[a],right):
        assert all(i in mapping for i in p['contacts'])
        out.append(dict(p,contacts={mapping[i]:v for i,v in p['contacts'].items()}))
    return out

def formulas():
    g,nu=affine();rows=[];cross=[];paid=[]
    for a,b in product((0,10,11),(9,12,13)):
        W=sorted(set(range(14))-{a,b});J=g.extract(W,W).inv();R=9+2*(a!=0)+2*(b!=9)
        va=-nu[a]*g.extract([a],W).T;vb=-nu[b]*g.extract([b],W).T
        old=s.diag(nu[a],nu[b])*g.extract([a,b],[a,b])*s.diag(nu[a],nu[b])
        old-=s.Matrix.hstack(va,vb).T*J*s.Matrix.hstack(va,vb)
        assert old==s.Matrix([[1,-1],[-1,1]])/R
        L,D=profiles(a),profiles(b,True)
        for side,ps,vfar in(('left',L,vb),('right',D,va)):
            for p in ps:
                w=s.Matrix([p['contacts'].get(i,0)for i in W])
                height=sum(nu[i]*q for i,q in p['contacts'].items())
                far=(vfar.T*J*w)[0];norm=2-(w.T*J*w)[0]
                assert s.simplify(height-p['h'])==0
                assert s.simplify(far-p['l']/R)==0
                assert s.simplify(norm-p['c0']-p['l']**2/R)==0
                rows.append(dict(marks=[a,b],side=side,profile=p['name'],R=R,
                    contacts={str(i):str(v)for i,v in p['contacts'].items()},c0=str(p['c0']),
                    moment=str(p['l']),height=str(p['h']),residual_norm=str(s.factor(norm)),far_contact=str(far)))
        for c,d in product(L,D):
            w=s.Matrix([c['contacts'].get(i,0)for i in W]);z=s.Matrix([d['contacts'].get(i,0)for i in W])
            e=(w.T*J*z)[0];assert s.simplify(e-c['l']*d['l']/R)==0
            cross.append(dict(marks=[a,b],c=c['name'],d=d['name'],value=str(e)))
    for name,j,t in(('inner4',1,s.sqrt(2)),('inner_next3',2,s.Integer(1))):
        h=2*s.eye(15);h[:14,:14]=g;h[14,j]=h[j,14]=-t
        ids=list(range(1,10))+[12,13,14];k=h.extract(ids,ids);nu=k.nullspace()[0]
        if nu[0]<0:nu=-nu
        assert len(k.nullspace())==1 and all(q>0 for q in nu)
        lap=s.diag(*nu)*k*s.diag(*nu)
        assert all(sum(lap.row(i))==0 for i in range(12))
        assert all(lap[i,j]<=0 for i in range(12)for j in range(12)if i!=j)
        paid.append(dict(profile=name,order=12,nodes=ids,positive_kernel=list(map(str,nu))))
    assert len(rows)==96 and len(cross)==256
    return dict(end_rows=rows,cross_rows=cross,whole_affine12_payments=paid,
        determinant='1-(g_c+g_d)/R; spectator correction=(x/h-y/k)^T G_Z^-1(x/h-y/k)/R')

def bounds():
    costs=M.cost_checks()['rows']
    margins=[]
    for a,b in product(('center','leaf'),repeat=2):
        R=9+2*(a=='leaf')+2*(b=='leaf');hi={'center':6,'leaf':8};lo={'center':3,'leaf':4}
        non=max(lo[a]+hi[b],hi[a]+lo[b]);gap=1-F(non,R)
        assert gap>=0 and hi[a]+hi[b]==R+3
        margins.append(dict(ends=[a,b],R=R,extreme_sum=R+3,one_nonextreme_sum_bound=non,
            nonextreme_gap=str(gap),zero_case='a center inner3 at native opposite contact0, opposite end extreme'if not gap else None,
            extreme_heights=['1/2','1/2'],required_unscaled_D='3/4'))
    return dict(native_cost_rows=costs,pair_margins=margins,
        outer_contacts='Center outer-only u,v satisfy u²+v²<4 by actual edge PD. The frozen all-label bounds remain unchanged.',
        equality_classification='The sole nonextreme cost3 center profile is inner3 with native opposite contact0; extreme costs force height1/2.')

def corner_data(g,W,a,b,c,d):
    rows=[]
    for p,q in((a,b),(b,c),(c,d),(d,a)):
        ids=sorted(set(W)|{p,q});sig=C.inertia(g.extract(ids,ids).tolist());parts=M.components(g,ids)
        for t in parts:
            t['marks']=[r for r in(p,q)if r in t['nodes']]
            t['fork_marks']=[r for r in t['marks']if r in t['branch']+t['short_leaves']]
        cp=next(t for t in parts if p in t['nodes']);cq=next(t for t in parts if q in t['nodes'])
        unit=0;half=F()
        if(p,q)!=(a,b):
            assert sig==[len(ids),0,0]
            if cp is cq:
                assert len(cp['nodes'])<=17;unit=1
            else:
                small=next((t for t in(cp,cq)if len(t['nodes'])<=5),None)
                large=next((t for t in(cp,cq)if t['type']in('D11','D13')and t['fork_marks']),None)
                if small and large:half=F(1,2)
        else:assert sig==[len(ids)-1,0,1]
        rows.append(dict(marks=[p,q],signature=sig,components=parts,unit=unit,half_coupon=str(half)))
    return rows

def zero_cases():
    f,_=affine();rows=[]
    for side in('left','right'):
        for far,L in product((9,12,13),(12,13))if side=='left'else product((0,10,11),(10,11)):
            a,b=(0,far)if side=='left'else(far,9)
            g=2*s.eye(16);g[:14,:14]=f
            edges=((14,1),(15,L))if side=='left'else((14,L),(15,8))
            for i,j in edges:g[i,j]=g[j,i]=-1
            W=sorted(set(range(14))-{a,b})
            assert C.inertia(g.tolist())==[14,1,1]
            assert C.inertia(g.extract(W,W).tolist())==[12,0,0]
            corners=corner_data(g,W,a,b,14,15)
            chosen=corners[1]if side=='left'else corners[3]
            assert chosen['unit']==1
            assert any(t['type']=='D12'and len(t['marks'])==2 for t in chosen['components'])
            rows.append(dict(marks=[a,b],center_inner3_side=side,opposite_leaf_attachment=L,
                signature=[14,1,1],corners=corners,chosen_unit_corner=chosen['marks']))
    assert len(rows)==12
    return rows

def extreme_cases():
    f,_=affine();rows=[];tally=Counter();types=Counter()
    for a,b,L,R,side in product((0,10,11),(9,12,13),(10,11),(12,13),(14,15)):
        g=2*s.eye(19);g[:14,:14]=f
        for i,j in((14,L),(15,R),(16,side),(16,17),(17,18)):g[i,j]=g[j,i]=-1
        W=sorted(set(range(14))-{a,b})+[16,17,18]
        assert C.inertia(g.tolist())==[17,1,1]
        assert C.inertia(g.extract(W,W).tolist())==[15,0,0]
        corners=corner_data(g,W,a,b,14,15)
        units=sum(r['unit']for r in corners);halves=sum(F(r['half_coupon'])for r in corners)
        assert units+halves>=1;tally[(units,str(halves))]+=1
        for r in corners:
            if F(r['half_coupon']):
                marked=[t for t in r['components']if t['marks']]
                sm=next(t for t in marked if len(t['nodes'])<=5)
                lg=next(t for t in marked if len(t['nodes'])>5)
                types[(sm['type'],lg['type'])]+=1
        rows.append(dict(marks=[a,b],leaf_attachments=[L,R],A3_attached_to=side,
            signature=[17,1,1],white=W,corners=corners,units=units,coupon_total=str(halves)))
    assert len(rows)==72 and tally==Counter({(0,'1'):32,(1,'1/2'):32,(3,'0'):8})
    assert types==Counter({('A5','D11'):24,('A2','D11'):24,('A1','D13'):24,('A4','D13'):24})
    return dict(rows=rows,payment_classes=[dict(units=k[0],halves=k[1],cases=v)for k,v in sorted(tally.items())],
        coupon_types=[dict(small=k[0],large=k[1],occurrences=v)for k,v in sorted(types.items())])

def run():
    d=dict(formulas=formulas(),bounds=bounds(),nonextreme_zero_cases=zero_cases(),extreme_equality=extreme_cases())
    paths=[Path(__file__)]
    for stem in('fresh_affine_d14_long_pair_transfer_20261017','fresh_small_spectator_difference_gap_20261017',
                'fresh_finite8_a17_d11_coupon_budget_20261017','fresh_finite8_ideal5_local_replacement_20261017'):
        paths += [ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    d.update(status='PASS_UNIVERSAL_AFFINE_D13_LONG_PAIR_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(marked_pairs=9,end_formulas=96,cross_formulas=256,small_affine12_payments=2,
            native_cost_rows=14,nonextreme_zero_cases=12,extreme_A3_cases=72,extreme_finite_corners=216),
        theorem='Every actual quadrilateral with a pair across the two fork triples of affineD13(order14) transfers to finite8/small17, ideal5/small13, or two sealed half-coupons: whole A1/A2/A4/A5 against whole D11/D13 marked at a fork center/short leaf.',
        scope='Universal local theorem with all native labels and arbitrary supporting spectators. No complete finite8 overlap coverage, full ideal small-factor rule or new ambient dimension bound.')
    return d

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
