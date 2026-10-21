#!/usr/bin/env python3
"""Exact all-contact depth8 finite7 proof, with one sealed A24 equality.

No finite-forest optimization or ambient dimension comparison is performed.
"""
import argparse
from collections import deque
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sympy as S
import fresh_depth10_finite7_zero_transfer_20261018 as D
import fresh_order13_sparse_ambient60_frontier_20261018 as R

ROOT=D.ROOT
TAG='fresh_depth8_finite7_profile_transfer_20261018'
STEMS=(R.TAG,D.TAG,D.B.TAG)
N,A,B=24,8,15
ALLOW=(0,1,6,7,9,22,23)
SQ=S.sqrt(2)


def field(x):
    """Return exact rational coordinates and sign, without numerical tests."""
    x=S.expand(x)
    a,b=x.coeff(SQ,0),x.coeff(SQ)
    assert a.is_Rational and b.is_Rational and S.expand(x-a-b*SQ)==0,x
    a,b=F(int(a.p),int(a.q)),F(int(b.p),int(b.q))
    if not b:sg=(a>0)-(a<0)
    elif not a:sg=(b>0)-(b<0)
    elif a*b>0:sg=(a>0)-(a<0)
    else:
        delta=a*a-2*b*b
        assert delta!=0
        sg=((a>0)-(a<0))*((delta>0)-(delta<0))
    return dict(rational=str(a),sqrt2=str(b),sign=sg)


def rat(q):return S.Rational(q.numerator,q.denominator)
def inv(n,i,j):return rat(D.B.inverse(n,i,j))


def path(n):
    return S.Matrix(n,n,lambda i,j:2 if i==j else -1 if abs(i-j)==1 else 0)


def projection(x,ids):
    blocks=[]
    for i in sorted(ids):
        if blocks and i==blocks[-1][-1]+1:blocks[-1].append(i)
        else:blocks.append([i])
    return S.expand(sum(x.get(i,0)*x.get(j,0)*inv(len(z),ii,jj)
        for z in blocks for ii,i in enumerate(z) for jj,j in enumerate(z)))


def profiles():
    W=[i for i in range(N) if i not in (A,B)]
    corner=[i for i in range(N) if i!=A]
    records=[];kept=[]
    for digits in product((0,1,2),repeat=len(ALLOW)):
        if sum(bool(q) for q in digits)>2:continue
        x={i:S.sqrt(q) for i,q in zip(ALLOW,digits) if q}
        edge=projection(x,W);adj=projection(x,corner)
        accepted=field(2-edge)['sign']>0 and field(2-adj)['sign']>=0
        row=dict(index=len(records),contacts={str(i):str(q) for i,q in x.items()},
            strict_edge_projection=field(edge),adjacent_corner_projection=field(adj),accepted=accepted)
        records.append(row)
        if accepted:kept.append(x)
    assert len(records)==99 and len(kept)==17
    expected=[{}, {23:1},{23:SQ},{22:1}, {9:1},{9:1,23:1},{9:SQ},
        {7:1},{7:1,23:1},{7:1,9:1},{7:SQ},{6:1},{1:1},
        {0:1},{0:1,23:1},{0:1,9:1},{0:SQ}]
    assert kept==expected
    allowed_by_short=[i for i in W if not (i>B and i-B<=5) and not (i<B and B-i<=5 and i>A)]
    # c cannot contact within five white steps of adjacent black b.
    excluded_arm=[]
    for i in W:
        if i<A and field(2-inv(8,i,i))['sign']<=0:excluded_arm.append(i)
        if i>B and field(2-inv(8,i-B-1,i-B-1))['sign']<=0:excluded_arm.append(i)
    assert sorted(set(allowed_by_short)-set(excluded_arm))==list(ALLOW)
    return kept,dict(candidate_count=99,accepted_count=17,candidates=records,
        strict_arm_excluded_positions=excluded_arm,allowed_positions=list(ALLOW),
        inverse_checks=[D.B.path_check(n) for n in (6,8,15,24)],
        minimum_white_inverse_diagonal='6/7',three_contact_projection_lower_bound='18/7',
        label_at_least5_square_strict_lower_bound='5/2',large_label_projection_strict_lower_bound='15/7')


def coefficient_certificate(ps):
    t,u,T,U=S.symbols('t u T U',nonnegative=True)
    J=S.Matrix(N,N,lambda i,j:inv(N,i,j))
    rows=[];discarded=[];positive=[];equal=[]
    for ci,c in enumerate(ps):
        for di,r in enumerate(ps):
            d={N-1-i:q for i,q in r.items()}
            short=[(i,j) for i in c for j in d if abs(i-j)+2<=6
                and not any(min(i,j)<k<max(i,j) for k in (A,B))]
            if short:
                discarded.append(dict(c_profile=ci,d_profile=di,short_cd_path=short[0]));continue
            x=S.Matrix([c.get(i,0)+(t if i==A else 0) for i in range(N)])
            y=S.Matrix([d.get(i,0)+(u if i==B else 0) for i in range(N)])
            q=(x.T*J*x)[0];rr=(y.T*J*y)[0];s=(x.T*J*y)[0]
            determinant=S.expand((q-2)*(rr-2)-s*s)
            for tc,ud in product((0,1),repeat=2):
                if not tc and not set(c)&{6,7,9}:continue
                if not ud and not set(d)&{14,16,17}:continue
                sub={t:1+T if tc else 0,u:1+U if ud else 0}
                p=S.Poly(S.expand(determinant.subs(sub)),T,U)
                terms=[dict(powers=list(m),coefficient=field(v)) for m,v in p.terms()]
                assert all(z['coefficient']['sign']>=0 for z in terms)
                q0=S.expand(q.subs({t:tc}));r0=S.expand(rr.subs({u:ud}))
                assert field(q0-2)['sign']>0 and field(r0-2)['sign']>0
                constant=S.expand(p.eval({T:0,U:0}));fc=field(constant)
                row=dict(c_profile=ci,d_profile=di,t_branch='1+T' if tc else '0',
                    u_branch='1+U' if ud else '0',q_lower=field(q0),r_lower=field(r0),
                    coefficients=terms,constant=fc)
                rows.append(row)
                if fc['sign']>0:positive.append(constant)
                else:
                    assert fc['sign']==0 and p.is_zero and (ci,di,tc,ud)==(4,4,0,0)
                    equal.append(row)
    assert len(rows)==581 and len(positive)==580 and len(equal)==1
    least=positive[0]
    for value in positive[1:]:
        if field(value-least)['sign']<0:least=value
    return dict(branch_count=len(rows),strict_positive_branches=580,unique_equality=equal,
        complete_branches=rows,short_cd_profile_pairs=discarded,minimum_positive_constant=field(least),
        arbitrary_opposite_products='Each is0 or1+nonnegative parameter; every coefficient has an exact Q(sqrt2) sign.')


def boundary_and_equality():
    boundaries=[]
    for name,g in [('A9',path(9)),('A10',path(10))]:
        boundaries.append((name,g))
    g=path(10);g[0,1]=g[1,0]=-SQ;boundaries.append(('B10',g))
    g=S.diag(path(9),S.Matrix([[2]]));g[1,9]=g[9,1]=-1;boundaries.append(('D10',g))
    checked=[]
    for name,g in boundaries:
        det=[field(g[:i,:i].det()) for i in range(1,len(g.tolist())+1)]
        assert all(z['sign']>0 for z in det)
        checked.append(dict(type=name,gram=[[str(q) for q in row] for row in g.tolist()],leading_determinants=det))
    gram=S.diag(path(24),S.Matrix([[2,0],[0,2]]))
    gram[9,24]=gram[24,9]=gram[14,25]=gram[25,14]=-1
    radical=S.Matrix([inv(N,i,9)-inv(N,i,14) for i in range(N)]+[1,-1])
    assert gram*radical==S.zeros(26,1)
    assert (radical[0],radical[23])==(S.Rational(1,5),-S.Rational(1,5))
    q,r,s=inv(N,9,9),inv(N,14,14),inv(N,9,14)
    assert (q,r,s)==(6,6,4)
    schur=S.Matrix([[2-q,-s],[-s,2-r]])
    assert schur==S.Matrix([[-4,-4],[-4,-4]])
    eigen=schur.eigenvals();assert eigen=={-8:1,0:1}
    pairings=[]
    for x,y in product((0,1,SQ),repeat=2):
        val=field(y-x-S.Rational(1,5));assert val['sign']!=0
        pairings.append(dict(c_contact=str(x),d_contact=str(y),radical_pairing=val))
    double=path(9);double[0,1]=double[1,0]=-SQ
    assert double.inv()[0,0]==S.Rational(9,2)
    assert inv(9,0,0)==S.Rational(9,10) and inv(8,0,0)==S.Rational(8,9)
    assert inv(24,0,0)==S.Rational(24,25)
    # In the equality corner remove a,b. Its middle connected component is
    # exactly c,9,...,14,d; both marked nodes are leaves of a whole A8.
    ids=[i for i in range(26) if i not in (A,B)]
    reachable={24};todo=[24]
    while todo:
        i=todo.pop()
        for j in ids:
            if i!=j and gram[i,j]!=0 and j not in reachable:reachable.add(j);todo.append(j)
    order=[24,*range(9,15),25]
    assert set(order)==reachable and gram.extract(order,order)==path(8)
    return dict(minimal_core_order_bound=10,disconnected_deletion_extra_roots_at_most=1,
        central_contact_finite_boundaries=checked,external_extra_contact_edge_inverse='8/9',
        equality_gram=[[str(q) for q in row] for row in gram.tolist()],
        equality_inertia=[24,1,1],positive_A24_block=True,schur_eigenvalues=['-8','0'],
        radical=list(map(str,radical)),endpoint_values=['1/5','-1/5'],
        endpoint_extension_inverse='24/25',simple_white_extension_inverse='9/10',
        double_white_extension_inverse='9/2',simple_extension_all_contact_pairs=pairings,
        double_extension_only_contacts=[0,0],double_extension_radical_pairing=field(-SQ/5),
        sealed_whole_A8=dict(nodes=order,marks=[24,25],distance=7,marked_degrees=[1,1],retained_payment='1'),
        outside_white_rule='Whole A24 is orthogonal to every remaining white; radical forces equal c/d contacts, giving short2 unless both vanish.')


def counts():
    rows=[]
    for kind in 'ABD':
        for n in range(1 if kind=='A' else 2 if kind=='B' else 4,81):
            edges=[(i,i+1,3) for i in range(n-1)] if kind!='D' else [(i,i+1,3) for i in range(n-3)]+[(n-3,n-2,3),(n-3,n-1,3)]
            if kind=='B':edges[-1]=(n-2,n-1,4)
            adj={i:[] for i in range(n)}
            for i,j,m in edges:adj[i].append(j);adj[j].append(i)
            features={i for i in adj if len(adj[i])!=2}|{i for a,b,m in edges if m!=3 for i in (a,b)}
            dist=[]
            for i in range(n):
                dd={i:0};todo=deque([i])
                while todo:
                    k=todo.popleft()
                    for j in adj[k]:
                        if j not in dd:dd[j]=dd[k]+1;todo.append(j)
                assert len(dd)==n;dist.append(dd)
            depth=[min(dist[i][f] for f in features) for i in range(n)]
            pairs=[(i,j) for i,j in combinations(range(n),2) if dist[i][j]==7 and min(depth[i],depth[j])>=8]
            expected=max(0,n-{'A':23,'B':24,'D':25}[kind]);assert len(pairs)==expected
            rows.append(dict(kind=kind,rank=n,deleted_pairs=pairs,count=expected))
    return dict(formulas={k:'max(0,n-'+str(v)+')' for k,v in dict(A=23,B=24,D=25).items()},
        exact_graph_checks_through_rank=80,rows=rows,all_rank_reason='Both marks lie on the unique feature-free spine interval; counting starting positions gives these formulas.')


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for s,doc in docs.items():
        assert doc['status'].startswith('PASS'),s
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[R.TAG]['finite']['exact_maximum']=='819/2'
    protected=docs[D.B.TAG]['triangle_interface']['retained_distance7_elementary_payments']
    assert len(protected)==14 and all(len(r['feature_witness']['path'])-1<=5<8 for r in protected)
    special=docs[D.B.TAG]['triangle_interface']['exceptional_distance7_marked_leaf']
    assert special['fixed_simple_leaf'] in special['alternate_pair']
    ps,domain=profiles();certificate=coefficient_certificate(ps)
    paths={Path(__file__),Path(D.B.__file__)}
    paths.update(ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in (('scripts','py'),('results','json'),('results','md')))
    return dict(status='PASS_DEPTH8_FINITE7_PROFILE_LOCAL_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
        current_baseline=R.TAG,feature_depth=8,finite_distance=7,new_weight='0',
        contact_domain=domain,schur_certificate=certificate,equality_and_central_boundary=boundary_and_equality(),
        deleted_pair_counts=counts(),inherited_triangle_protected_payments=protected,
        inherited_exceptional_marked_leaf=special,triangle_minimum=2,quadrilateral_minimum=1,
        ideal_weights_unchanged=True,finite_compensation='0',all_other_weights_and_coupons_unchanged=True,
        new_forest_maximum=None,new_dimension_bound=None,new_finite_volume_example=False,
        scope='Complete99 native white profiles,17 accepted,581 exact unbounded opposite-product polynomial branches; central boundary and sealed equality proof. Inherited all-label minimal-core theorem and current short6/triangle geometry remain written inputs. No new finite DP, ideal convolution or ambient theorem.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    result=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert result==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],'profiles',result['contact_domain']['accepted_count'],
        'branches',result['schur_certificate']['branch_count'],'equalities',len(result['schur_certificate']['unique_equality']),flush=True)
