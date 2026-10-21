#!/usr/bin/env python3
"""An exact cap direction satisfying the scalar H8 attachment necessities.

No source placement or global coverage is claimed. An isolated quadratic root
over the frozen quartic field defines the cap without expanding its extension.
"""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_h9_mixed_d4_coupling_20261014_review as A
from fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate import I,ei,SCALE,BITS
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h9_h8_global_cap_gap_20261014'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
SOURCE='data/known_examples/bugaenko_compact_h8.json'
GAP='results/fresh_h9_h8_facet_cover_20261014_gap_gate.json'
ANCHOR='results/fresh_h9_h8_corner_anchors_20261014_pair_gate.json'
V=['3.346178','4.2047','3.038572','1.952451','.909731','.150188','-1.149919','-.592165','-.912382','-.750588']
TL=F('1.00000001124723266977347773466025')
TU=F('1.00000001124723266977347773466026')
SMALL=A.E(A.Q(F(-1,4),F(1,4)));BOUND=4*SMALL
DIAM=A.E(A.Q(17,8))
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def check_sources(path):
    d=json.loads((ROOT/path).read_text())
    for p,h in d.get('sources',{}).items():assert sha(p)==h,p
    return d
def hull(a,b):return I.raw(a.lo,b.hi)
def ipair(x,y,Q):return sum((x[i]*sum((ei(Q[i][j])*y[j]for j in range(len(y))),I(0))for i in range(len(x))),I(0))
def run():
    b=check_sources(BASE);gap=check_sources(GAP);anchor=check_sources(ANCHOR)
    assert gap['exceptional_minimum_norm']==SMALL.a.enc() and gap['ordinary_minimum_norm']==BOUND.a.enc()
    Q,N,H=A.dm(b['form']),A.dm(b['final_normals']),A.dv(b['center'])
    assert len(N)==12 and A.pair(H,H,Q).sign()<0
    v=[A.E(F(x))for x in V];d=[x-y for x,y in zip(v,H)]
    q0=A.pair(H,H,Q);q1=2*A.pair(H,d,Q);q2=A.pair(d,d,Q)
    g0=[A.pair(H,n,Q)for n in N];g1=[A.pair(d,n,Q)for n in N]
    nn=[A.pair(n,n,Q)for n in N]
    assert nn[8]==2
    # Exact target: the normalized old wall8 has projected square -SMALL.
    poly=[g0[8]*g0[8]-(2+SMALL)*q0,
          2*g0[8]*g1[8]-(2+SMALL)*q1,
          g1[8]*g1[8]-(2+SMALL)*q2]
    pe=lambda t:poly[0]+t*poly[1]+t*t*poly[2]
    assert pe(TL).sign()*pe(TU).sign()==-1
    deriv=[poly[1]+2*t*poly[2]for t in(TL,TU)]
    assert deriv[0].sign()==deriv[1].sign()!=0
    t=hull(I(TL),I(TU));z=[ei(h)+t*ei(x)for h,x in zip(H,d)]
    q=ei(q0)+t*ei(q1)+t.square()*ei(q2);assert q.lo>0
    g=[ei(a)+t*ei(c)for a,c in zip(g0,g1)]
    assert all(x.hi<0 for x in g)
    hz=ei(q0)+t*ei(q1)/2;assert hz.hi<0
    kapp=[2*x.square()/(q*ei(n))-2 for x,n in zip(g,nn)]
    assert kapp[8].overlap(ei(SMALL))
    assert all(x.lo>ei(BOUND).hi for i,x in enumerate(kapp)if i!=8)
    assert all((g[i].square()-q*ei(nn[i])).lo>0 for i in range(12))
    # Scale every old normal to square2. The low wall8 is already norm2.
    a=[-g[i]*(2/ei(nn[i])).sqrt()/q.sqrt()for i in range(12)]
    G=[[ei(A.pair(x,y,Q))*(4/(ei(nn[i])*ei(nn[j]))).sqrt()for j,y in enumerate(N)]for i,x in enumerate(N)]
    low=[];high=[];dots=[]
    rho_s=[(2/k).sqrt()for k in kapp];rho_c=[((k+2)/k).sqrt()for k in kapp]
    distances=[];unions=[]
    for i,j in combinations(range(12),2):
        dc=(a[i]*a[j]-G[i][j])/(kapp[i]*kapp[j]).sqrt()
        assert dc.lo>SCALE and dc.hi<ei(DIAM).lo
        ds=(dc.square()-1).sqrt()
        u=dc*(rho_c[i]*rho_c[j]+rho_s[i]*rho_s[j])+ds*(rho_s[i]*rho_c[j]+rho_c[i]*rho_s[j])
        distances.append(dict(pair=[i,j],cosh=dc.enc()))
        unions.append(dict(pair=[i,j],cosh=u.enc()))
        rel=A.relation(A.pair(N[i],N[j],Q),nn[i],nn[j])
        if rel[0]=='dotted':dots.append([i,j])
        else:assert rel[0]in('right','finite')and rel[1]<=5
    assert dots==[[4,11],[5,6],[8,11]]
    selected=next(r for r in unions if r['pair']==[7,8]);assert int(selected['cosh'][0])>ei(DIAM).hi
    # General finite-pair diameter bound whenever both projected norms>=BOUND.
    # Decreasing either kappa increases its radius and the center distance.
    B=ei(BOUND);phi=ei(A.PHI)
    dc=(B+2+phi)/B;sr=(2/B).sqrt();cr=((B+2)/B).sqrt()
    finite_upper=dc*(cr.square()+sr.square())+(dc.square()-1).sqrt()*2*cr*sr
    assert finite_upper.hi<ei(DIAM).lo
    # Reconstruct the exact source diameter, rather than accepting a decimal.
    src=check_sources(SOURCE)
    def sf(x):return A.E(A.Q(F(x['a'])+F(x['b'])/2,F(x['b'])/2))
    SQ=[[sf(x)for x in row]for row in src['lattice_gram']]
    SN=[[sf(x)for x in row]for row in src['simple_roots']]
    SH=[sf(x)for x in src['geometry']['interior_point']]
    SV=[]
    for ids in src['geometry']['vertices']:
        w=A.kernel_one([A.mv(SQ,SN[i])for i in ids]);assert w is not None
        if A.pair(w,SH,SQ).sign()>0:w=[-x for x in w]
        nw=A.pair(w,w,SQ);assert nw.sign()<0
        assert all(A.pair(w,n,SQ).sign()<=0 for n in SN)
        SV.append((w,nw))
    assert len(SV)==41
    diameter_pairs=[]
    for i,j in combinations(range(41),2):
        w,nw=SV[i];x,nx=SV[j];p=A.pair(w,x,SQ);assert p.sign()<0
        diff=DIAM*DIAM*nw*nx-p*p;assert diff.sign()>=0
        if not diff:diameter_pairs.append([i,j])
    assert diameter_pairs==[[32,33]]
    # The cap itself is an actual facet of a coherent13-wall chamber.
    Hc=[ei(x)for x in H];f=[h-hz*x/q for h,x in zip(Hc,z)]
    fq=ipair(f,f,Q);assert fq.hi<0
    fp=[ipair(f,[ei(x)for x in n],Q)for n in N];assert all(x.hi<0 for x in fp)
    assert ipair(f,z,Q).overlap(I(0))
    paths=[BASE,SOURCE,GAP,ANCHOR,'scripts/fresh_h9_mixed_d4_coupling_20261014_review.py',
        'scripts/fresh_h9_mixed_d4_endpoint_profiles_20261014_universal_gate.py',f'scripts/{TAG}.py']
    return dict(status='PASS_SCALAR_H8_ATTACHMENT_NECESSITIES_REMAIN_FEASIBLE',
        sources={p:sha(p)for p in paths},interval_bits=BITS,interval_scale=str(SCALE),
        cap_recipe=dict(center_source=BASE,rational_endpoint=V,parameter_interval=[str(TL),str(TU)],
            polynomial_ascending=A.enc(poly),endpoint_polynomial_values=A.enc([pe(TL),pe(TU)]),
            derivative_endpoint_values=A.enc(deriv),coordinates=[x.enc()for x in z],norm=q.enc(),center_pairing=hz.enc()),
        old_pairings=[x.enc()for x in g],normalized_projected_norm_squares=[x.enc()for x in kapp],
        unique_low_wall=8,exact_low_norm=A.enc(SMALL),high_norm_threshold=A.enc(BOUND),
        source_diameter_cosh=A.enc(DIAM),source_diameter_pairs=diameter_pairs,
        projected_center_pair_cosh=distances,ball_union_pair_cosh=unions,
        ball_union_exceeds_source_diameter_pair=[7,8],actual_cap_facet=dict(point=[x.enc()for x in f],norm=fq.enc(),old_pairings=[x.enc()for x in fp]),
        finite_old_pair_high_norm_union_bound=finite_upper.enc(),old_dotted_pairs=dots,
        general_necessary_condition='If source diameter endpoints are covered by two old balls whose centers both have kappa^2>=sqrt5-1, their old wall pair must be dotted. Otherwise at least one covering center has low projected norm.',
        scope='An exact all-dotted cap direction satisfying the source norm-gap count, source-center diameter and ball-union diameter necessities. No isometric H8 placement, all source-sidewall Coxeter relations, whole-ball coverage, finite-volume completion, or placement failure is established. It only shows that these scalar necessities do not themselves yield a universal exclusion.')
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();d=run()
    path=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(path.read_text());print('SAVED_REPLAY_PASS')
    else:path.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['unique_low_wall'],d['old_dotted_pairs'])
