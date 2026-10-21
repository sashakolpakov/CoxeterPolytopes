#!/usr/bin/env python3
"""Exact Lorentz-dual escape bounds for 144 fixed E7 contact profiles.

The residual ranges over the complete intrinsic space, not an edge/ray bank.
No optimizer or floating arithmetic is used in the verification.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from pathlib import Path
import fresh_h20_e7_residual_gate_20261014_exact_data as X

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_e7_residual_gate_20261014_dual'
RECIPES='results/fresh_h20_e7_residual_gate_20261014_dual_recipes.json'
K=X.K
Z=K.Z+K.Z


def lift(a):return a+K.Z
def scalar(a):return lift(K.scalar(a))
def add(a,b):return tuple(x+y for x,y in zip(a,b))
def scale(a,q):return tuple(q*x for x in a)
def sub(a,b):return add(a,scale(b,-1))


def mul(a,b):
    p=K.multiply(a[:8],b[:8]);q=K.multiply(a[8:],b[8:])
    return K.add(p,K.scale(q,3))+K.add(K.multiply(a[:8],b[8:]),K.multiply(a[8:],b[:8]))


def sign(a):
    sa,sb=K.sign(a[:8]),K.sign(a[8:])
    if not sb:return sa
    if not sa or sa==sb:return sb
    q=K.sub(K.multiply(a[:8],a[:8]),K.scale(K.multiply(a[8:],a[8:]),3))
    sq=K.sign(q)
    assert sq  # sqrt3 is outside Q(sqrt2,sqrt79,sqrt11).
    return sa*sq


def total(items):
    out=Z
    for x in items:out=add(out,x)
    return out


def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(v)for v in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x


def run():
    data=X.load();recipes=json.loads((ROOT/RECIPES).read_text())
    for p,h in recipes['sources'].items():assert X.digest(p)==h,p
    JI=K.inverse(data['J'])
    S=data['S'];I=data['I'];C=data['coords'];tau=data['times'];P=data['wp']
    a=data['a'];B=data['B'];weightE=data['weightE']
    coeff={3:scalar(1),4:lift((F(0),F(1))+(F(0),)*6),6:K.Z+K.scalar(1)}
    rows=[];seen=set();margin=F(1,100000)
    for recipe in recipes['records']:
        support=recipe['support'];labels=recipe['labels']
        assert len(support)==len(labels)and 1<=len(support)<=3
        assert support==sorted(set(support))and all(0<=j<7 for j in support)
        assert all(m in [3,4,6]for m in labels)
        assert len(support)<3 or labels==[3]*3
        key=tuple(zip(support,labels));assert key not in seen;seen.add(key)
        c=[Z]*7
        for j,label in key:c[j]=coeff[label]
        mass=sub(total(scale(mul(c[i],c[j]),I[i][j])for i in support for j in support),scalar(2))
        assert sign(mass)>=0
        lam={int(i):F(v)for i,v in recipe['old_wall_multipliers'].items()}
        mu={int(i):F(v)for i,v in recipe['nonnegative_contact_multipliers'].items()}
        assert all(0<=i<87 and v>0 for i,v in lam.items())
        assert all(0<=i<13 and v>0 for i,v in mu.items())
        ell=[K.add(a[j],K.scalar(mu.get(j,F(0))))for j in range(13)]
        h=B
        for i,v in lam.items():
            h=K.add(h,K.scale(tau[i],v))
            for j in range(13):ell[j]=K.add(ell[j],K.scale(C[j][i],v))
        hsign=K.sign(h);assert hsign>0
        delta=K.scale(K.multiply(h,h),F(1,2))
        for i in range(13):
            for j in range(13):delta=K.sub(delta,K.scale(K.multiply(ell[i],ell[j]),JI[i][j]))
        assert K.sign(delta)>0
        constant=total(mul(lift(weightE[j]),c[j])for j in support)
        for i,v in lam.items():constant=add(constant,scale(total(mul(lift(P[j][i]),c[j])for j in support),v))
        # f <= constant - sqrt(mass*delta), by the reverse Lorentz
        # Cauchy inequality. Certify the stronger uniform f<-1/100000.
        shifted=add(constant,scalar(margin))
        if sign(shifted)<0:
            reason='strictly_negative_constant';gap=None
        else:
            gap=sub(mul(mass,lift(delta)),mul(shifted,shifted))
            assert sign(gap)>0,(support,labels)
            reason='positive_squared_Lorentz_margin'
        rows.append(dict(support=support,E7_walls=[S[j]for j in support],labels=labels,
                         method=recipe['method'],contact_coefficients=[c[j]for j in support],
                         mass_square=mass,old_wall_multipliers=lam,nonnegative_contact_multipliers=mu,
                         time_coefficient=h,spatial_coefficients=ell,Lorentz_dual_square=delta,
                         upper_bound_constant=constant,certified_upper_bound=-margin,
                         sign_certificate=reason,squared_margin=gap))
    assert len(rows)==144 and len(seen)==144
    methods=dict(Counter(r['method']for r in rows))
    assert methods==dict(direct_Lorentz=129,retained_wall_Lorentz_dual=15)
    assert sum(6 not in r['labels']for r in rows)==69
    # The one omitted negative-mass profile cannot have a causal residual.
    assert I[5][5]-2==F(-1,2)
    sources=[X.BASE,X.PARENT,X.FACE,RECIPES,f'scripts/{TAG}.py',
             'scripts/fresh_h20_e7_residual_gate_20261014_exact_data.py',
             'scripts/fresh_h20_e8_interior_caps_20261014_exact_review.py']
    return enc(dict(status='PASS_144_FIXED_E7_PROFILES_ALL_INTRINSIC_RESIDUALS_MISS_ESCAPE',
                    sources={p:X.digest(p)for p in sources},
                    field_radicands=X.RADICANDS+[3*r for r in X.RADICANDS],
                    E7_support=S,face_corner_walls=data['A'],face_inverse=data['J'],
                    profile_count=144,methods=methods,label3_4_only_count=69,
                    uniform_pairing_upper_bound=-margin,records=rows,
                    negative_mass_profile=dict(support=[5],E7_wall=S[5],labels=[3],weight_square=F(3,2),required_residual_square=F(1,2)),
                    scope='For each listed normalized E7 contact profile and every intrinsic residual acute with all87 walls, the relaxed future-cone norm condition already implies r.escape<-1/100000. This proves no escape-cutting norm2 Coxeter cap in these144 fixed profiles. The separate100 positive numerical relaxations are not certified Coxeter caps and are not excluded. No all-label or all-profile theorem, no finite-volume construction.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['methods'])
