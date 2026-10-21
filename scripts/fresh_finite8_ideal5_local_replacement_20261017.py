#!/usr/bin/env python3
"""Ideal-five local transfers with finite8/leaf9 and explicit sealed coupons."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

import fresh_finite8_cycle_replacement_20261017 as R

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_finite8_ideal5_local_replacement_20261017'
A,B,C=R.A,R.B,R.C


def cycle6():
    rows=[];special=[];regions=[]
    for Q in range(9,14):
        keep=[(p,R.classification(p,6,Q)[0])for p in B.cycle_profiles()]
        keep=[(p,k)for p,k in keep if k in('unrestricted','requires_spectator_contact')]
        strict=[]
        for(c,ck),(d,dk)in product(keep,repeat=2):
            side='c'if ck=='requires_spectator_contact'else'd'if dk=='requires_spectator_contact'else None
            cc,dd=dict(c),dict(d)
            if side=='c':cc['c0']-=A.Q(F(1,2))
            if side=='d':dd['c0']-=A.Q(F(1,2))
            assert cc['c0'].sign()>=0 and dd['c0'].sign()>=0
            vals=[]
            for ub,vb in product((0,1),repeat=2):
                u=c['umax']if ub else A.Q();v=d['umax']if vb else A.Q()
                val=A.value(A.Q(F(1,6)),A.Q(F(1,Q)),cc,dd,u,v)
                vals.append(val)
                rows.append(dict(Q_lower=Q,all_Q_at_least=Q==13,c=c['name'],d=d['name'],
                    fixed_cauchy_side=side,upper_height_reciprocals=[ub,vb],value=val.wire()))
            low=R.minimum(vals)
            if low.sign()>0:
                strict+=vals;continue
            assert side is None and c['name']in('isolated','end3_1')and d['name']in('isolated','end3_1')
            high=[]
            for which in('c','d'):
                # Isolated h is1 or>=sqrt2; endpoint3 h=x+1 is1 or>=2.
                uc=A.Q(1)/A.SQ if c['name']=='isolated'else A.Q(F(1,2))
                vd=A.Q(1)/A.SQ if d['name']=='isolated'else A.Q(F(1,2))
                values=[A.value(A.Q(F(1,6)),A.Q(F(1,Q)),cc,dd,
                    (uc if which=='c'else c['umax'])if ub else A.Q(),
                    (vd if which=='d'else d['umax'])if vb else A.Q())for ub,vb in product((0,1),repeat=2)]
                gap=R.minimum(values);assert gap.sign()>0
                high.append(dict(side=which,minimum_gap=gap.wire(),values=[v.wire()for v in values]))
            at_one=A.value(A.Q(F(1,6)),A.Q(F(1,Q)),c,d,A.Q(1),A.Q(1))
            assert at_one==low and low.b==0
            required=-low.a/(F(1,6)+F(1,Q))
            assert 0<=required<F(1,2)
            outcome='excluded_by_spectator_gap'if required else'sealed_A11_finite_component'if Q==10 else'sealed_A1_A17_two_half_coupons'
            special.append(dict(Q=Q,c=c['name'],d=d['name'],minimum=low.wire(),
                high_height_checks=high,forced_heights=['1','1'],required_unscaled_D=str(required),outcome=outcome))
        regions.append(dict(Q_lower=Q,all_Q_at_least=Q==13,
            profiles=[dict(name=p['name'],kind=k)for p,k in keep],strict_minimum=R.minimum(strict).wire()))
    assert len(rows)==400 and len(special)==7
    assert Counter(r['outcome']for r in special)==dict(excluded_by_spectator_gap=5,sealed_A11_finite_component=1,sealed_A1_A17_two_half_coupons=1)
    return dict(rectangles=rows,special=special,regions=regions)


def cycle_equality(m,attachments):
    n=m+2;a,b=0,6
    g=[[F(2 if i==j else 0)for j in range(n)]for i in range(n)]
    for i,j in[(i,(i+1)%m)for i in range(m)]+attachments:g[i][j]=g[j][i]=F(-1)
    assert C.inertia(g)==[m,1,1]
    W=[i for i in range(m)if i not in(a,b)]
    assert C.inertia([[g[i][j]for j in W]for i in W])==[m-2,0,0]
    rows=[]
    for p,q in((a,b),(b,m),(m,m+1),(m+1,a)):
        ids=sorted(set(W)|{p,q});sig=C.inertia([[g[i][j]for j in ids]for i in ids])
        adj={i:{j for j in ids if i!=j and g[i][j]}for i in ids}
        cp=C.distances(adj,ids,p);cq=C.distances(adj,ids,q)
        if(p,q)==(a,b):
            assert sig==[m-1,0,1]
            rows.append(dict(marks=[p,q],kind='ideal',unit=0,coupon='0'));continue
        assert sig==[m,0,0]
        unit=int(q in cp and(len(cp)<=17 or cp[q]<=8))
        parts=[];coupon=F()
        for mark,ds in((p,cp),(q,cq)):
            leaves=[i for i in ds if len(adj[i]&ds.keys())==1]
            path=all(len(adj[i]&ds.keys())<=2 for i in ds)
            enddist=sorted(ds[i]for i in leaves)
            parts.append(dict(nodes=sorted(ds),mark=mark,type=f'A{len(ds)}'if path else'other',
                endpoint_distances=enddist))
        if q not in cp and sorted((len(cp),len(cq)))==[1,17]:
            large=next(t for t in parts if len(t['nodes'])==17)
            assert large['endpoint_distances']==[5,11]
            coupon=F(1,2)
        rows.append(dict(marks=[p,q],kind='finite',signature=sig,marked_components=parts,
            marked_distance=cp.get(q),unit=unit,coupon=str(coupon)))
    assert sum(r['unit']for r in rows)+sum(F(r['coupon'])for r in rows)>=1
    if m==18:assert sum(F(r['coupon'])for r in rows)==1
    else:assert m==16 and rows[2]['marked_components'][0]['type']=='A11'
    return dict(order=m,gram=[[str(x)for x in row]for row in g],signature=[m,1,1],white=W,corners=rows,
        sealing='D=0 and absence of a common spectator neighbor imply all spectators orthogonal to c,d and F.',
        scope='Conditional literal equality core. No finite-volume realization is asserted or required.')


def outer_exception():
    # D8 (0..7) plus a simple prolongation of a short fork leaf is affine E8.
    edges=[(i,i+1)for i in range(5)]+[(5,6),(5,7),(6,8)]
    g=[[F(2 if i==j else 0)for j in range(9)]for i in range(9)]
    for i,j in edges:g[i][j]=g[j][i]=F(-1)
    nu=[1,2,3,4,5,6,4,3,2]
    assert all(sum(g[i][j]*nu[j]for j in range(9))==0 for i in range(9))
    assert C.inertia(g)==[8,0,1]
    assert all(C.inertia([[g[i][j]for j in range(9)if j!=k]for i in range(9)if i!=k])==[8,0,0]for k in range(9))
    return dict(edges=edges,positive_kernel=nu,signature=[8,0,1],
        exclusion='For an outer D8 component this whole affine block lies in actual edge W+c, which is PD.')


def trees():
    rows=[]
    for P in(6,7,8,9):
        for c,d in product(B.tree_profiles(),repeat=2):
            if P<=8 and any(p['name']=='inner_next3'for p in(c,d)):
                rows.append(dict(P=P,c=c['name'],d=d['name'],kind='immediate_payment',length=P-1,
                    finite_reason='distance<=8',ideal_reason='distance<=5'if P==6 else'persistent feature within1 of new mark'))
                continue
            costs=[2*p['moment']/p['minimum_height']+p['c0']/p['minimum_height']**2 for p in(c,d)]
            gap=A.Q(1)-sum(costs,A.Q())/P
            assert gap.sign()>=0
            if gap.sign()==0:
                assert P==6 and c['name']==d['name']=='inner_end3'
                rows.append(dict(P=P,c=c['name'],d=d['name'],kind='sealed_finite_A7',forced_heights=['1','1'],
                    identity='F=[(3-2u-u^2)+(3-2v-v^2)]/6; u,v<=1; equality onlyu=v=1'))
            else:rows.append(dict(P=P,c=c['name'],d=d['name'],kind='positive_gap',gap=gap.wire()))
    assert Counter(r['kind']for r in rows)==dict(immediate_payment=45,positive_gap=210,sealed_finite_A7=1)
    # After deleting a,b the middle A5 has c,d at its endpoints: whole A7.
    edges=[(i,i+1,1)for i in range(6)];piv=C.tree_pd(range(7),edges)
    return dict(rows=rows,minimum_outer_B_rank=7,minimum_outer_D_rank=8,
        added_outer_exception=outer_exception(),sealed_middle=dict(type='A7',marked_distance=6,positive_pivots=piv))


def run():
    data=dict(cycle6=cycle6(),cycle18_equality=cycle_equality(18,[(0,18),(6,19)]),
        cycle16_end3_equality=cycle_equality(16,[(15,16),(7,17)]),tree=trees())
    paths=[Path(__file__)]
    for stem in('fresh_finite8_cycle_replacement_20261017','fresh_ideal_distance6_mixed_replacement_20261017',
                'fresh_cycle16_half_coupon_transfer_20261017','fresh_distance8_ideal_corner_replacement_20261017'):
        paths += [ROOT/f'{folder}/{stem}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    paths += [ROOT/'results/fresh_cusp_factor_payment_framework_20261017.md']
    data.update(status='PASS_FINITE8_LEAF9_IDEAL5_LOCAL_REPLACEMENTS',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        counts=dict(cycle6_rectangle_values=400,cycle6_special_pairs=7,high_height_gates=14,
            literal_cycle_equality_cores=2,tree_profile_cases=256,outer_affine_exception_kernels=1),
        theorem='Actual ideal cycle distances6..9 at order>=15 and classical tree distances6..9 outside feature radius5 transfer to finite8/whole-small17/marked-leaf9, ideal5 or distances6..9 near features within5, or sealed A1 x A15-middle/A17-position5or11 half-coupons.',
        scope='Universal local face transfer only. No complete high-core overlap census, global finite budget, full small-factor payment rule, or ambient dimension improvement.')
    return data


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],flush=True)
