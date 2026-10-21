#!/usr/bin/env python3
"""Degree-50 abelian-base upgrade and all-label consecutive-arc escape gate.

Standalone Fraction/Q(sqrt2) arithmetic. Inherits the separately reviewed
modular factorization and physical 134-wall certificate; no discovery bank.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_asymmetric_arc_escape_gate_20261016'
BASE='results/fresh_h20_s61_multicap_20261014_exact.json'
CAP='results/fresh_h20_dotted_arc_cap_20261016.json'
CAP_REVIEW='results/fresh_h20_dotted_arc_cap_20261016_review.json'
ABEL='results/fresh_h20_chain_abelian_gate_20261016.json'
ABEL_REVIEW='results/fresh_h20_chain_abelian_gate_20261016_review.json'
SYMMETRIC_REVIEW='results/fresh_h20_symmetric_arc_escape_gate_20261016_review.json'

def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def mul(x,y):return (x[0]*y[0]+2*x[1]*y[1],x[0]*y[1]+x[1]*y[0])
def scale(x,a):return (x[0]*a,x[1]*a)
def sign(x):
    a,b=x
    if not b:return (a>0)-(a<0)
    if not a:return (b>0)-(b<0)
    if a*b>0:return (a>0)-(a<0)
    d=a*a-2*b*b;assert d
    return ((a>0)-(a<0))*((d>0)-(d<0))
def enc(x):return list(map(str,x))
def padd(a,b):
    c=[F(0)]*max(len(a),len(b))
    for i,x in enumerate(a):c[i]+=x
    for i,x in enumerate(b):c[i]+=x
    while len(c)>1 and not c[-1]:c.pop()
    return c

def pmul(a,b):
    c=[F(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return padd(c,[0])
def pscale(a,t):return padd([x*t for x in a],[0])
def matmul(a,b):return [[sum(a[i][k]*b[k][j]for k in range(len(b)))for j in range(len(b[0]))]for i in range(len(a))]
def pd(a):
    a=[row[:]for row in a];pivots=[]
    for i in range(len(a)):
        p=a[i][i];assert p>0;pivots.append(str(p))
        for j in range(i+1,len(a)):
            for k in range(i+1,len(a)):a[j][k]-=a[j][i]*a[i][k]/p
    return pivots

def manifest(paths):
    seen=set();references=0
    def visit(p):
        nonlocal references
        if p in seen:return
        seen.add(p)
        if not p.endswith('.json'):return
        doc=json.loads((ROOT/p).read_text())
        for q,h in doc.get('sources',{}).items():
            references+=1;assert sha(q)==h,(p,q);visit(q)
    for p in paths:visit(p)
    return dict(files=len(seen),source_references=references)

def run():
    paths=[BASE,CAP,CAP_REVIEW,ABEL,ABEL_REVIEW,SYMMETRIC_REVIEW]
    base,cap,caprev,abel,abelrev,symrev=[json.loads((ROOT/p).read_text())for p in paths]
    manifests=manifest(paths)
    assert caprev['sources'][CAP]==sha(CAP)
    assert abelrev['sources'][ABEL]==sha(ABEL)
    assert abel['tower_degree']==64 and abel['irreducible_degree_over_tower']==50
    reductions={r['p']:r for r in abel['modular_reductions']}
    assert reductions[2017]['irreducibility']['degree']==50
    assert reductions[359]['factor_degrees']==[13,37]
    for p,r in reductions.items():
        assert r['simple_hensel_roots'] and r['squarefree']
        assert 0<r['polynomial_leading_unit']<p
    # Frobenius type(13)(37), raised to13, fixes13 roots and cycles37.
    sigma=list(range(50))
    for lo,length in [(0,13),(13,37)]:
        for j in range(length):sigma[lo+j]=lo+(j+1)%length
    tau=[]
    for j in range(50):
        k=j
        for _ in range(13):k=sigma[k]
        tau.append(k)
    lengths=[];seen=set()
    for j in range(50):
        if j in seen:continue
        k=j;l=0
        while k not in seen:seen.add(k);l+=1;k=tau[k]
        lengths.append(l)
    assert sorted(lengths)==[1]*13+[37]
    blocks=[b for b in range(2,50)if 50%b==0]
    assert blocks==[2,5,10,25]
    assert all(b<37 and 50//b<37 for b in blocks)
    group=dict(degree=50,Frobenius_cycle_lengths=[13,37],power=13,
        powered_cycle_lengths=sorted(lengths),proper_block_sizes=blocks,
        block_counts=[50//b for b in blocks],primitive=True,
        commutator_transitive=True,degree_over_K_Qab=50,
        no_nonzero_polynomial_over_Qab_of_degree_below=50,
        scope='Degree50 is over K Qab; only the lower bound50 is asserted over Qab itself.')
    # Polynomial identities for consecutive roots, ordered left/middle/right.
    g=[F(-2),F(-4)];four_minus_g2=padd([4],pscale(pmul(g,g),-1))
    det=pscale(four_minus_g2,2)
    assert det==[0,F(-32),F(-32)]
    # f^T adj(G) f: coefficients of A^2,B^2,D^2,BD.
    numerator={'A2':four_minus_g2,'B2':[F(4)],'D2':[F(4)],'BD':pscale(g,-4)}
    assert numerator=={'A2':[0,F(-16),F(-16)],'B2':[F(4)],'D2':[F(4)],'BD':[F(8),F(16)]}
    outer_sum_norm=padd([4],pscale(g,2));assert outer_sum_norm==[0,F(-8)]
    # If B=D=0,A>0 then either next available contact is (2C+1)A>0.
    centers=list(range(108,132));next_walls=[]
    for center in centers:
        available=center+2 if center+2<=132 else center-2
        assert 107<=available<=132 and available not in [center-1,center,center+1]
        next_walls.append([center,available])
    S=base['support'];J=[[F(q)for q in row]for row in base['support_inverse']]
    G=[[F(q)for q in row]for row in base['support_gram']]
    assert len(S)==18 and all(G[i][i]==2 for i in range(18))
    assert all(G[i][j]<=0 for i in range(18)for j in range(18)if i!=j)
    assert matmul(G,J)==[[F(i==j)for j in range(18)]for i in range(18)]
    support_pivots=pd(G);assert all(q>=0 for row in J for q in row)
    assert min(J[i][i]for i in range(18))==F(1,2)
    alo,ahi=map(F,base['parameter']['isolating_interval']);assert alo>0 and alo*alo>280
    inside=[S.index(i)for i in range(1,11)];outside=[i for i in range(18)if i not in inside]
    assert [[J[i][j]for j in inside]for i in inside]==[
        [F(min(i,j)*(11-max(i,j)),11)for j in range(1,11)]for i in range(1,11)]
    assert all(J[i][j]==0 for i in inside for j in outside)
    total=F(281,70);small=F(141,70)
    assert 4+4/(alo*alo)<total and total-2==small
    assert total-F(40,11)<F(1,2)
    assert 4*F(18,11)>total
    assert F(24,11)>small and F(10,11)*F(5,2)>small and 2*F(18,11)>small
    # E0=x-e/100-omega5 is weakly feasible for old134, by acuteness:
    # E0=Eeps+eps*sum_{j!=5}omega_j and omega_j.N<=0 offS.
    saved=cap['strict_escape'];eb=list(map(F,saved['dual_coefficients']))
    assert saved['support']==S and saved['frame_coefficients']==['1','-1/100','0']
    assert eb==[F(-1)if i==5 else F(-1,100000)for i in S]
    assert len(saved['all134_pair_bounds'])==134
    assert all(F(pair[1])<0 for pair in saved['all134_pair_bounds'])
    assert F(saved['strict_center_pair_bounds'][1])<0
    e0=[F(-1)if i==5 else F(0)for i in S]
    z=[J[i][S.index(5)]for i in range(18)]
    assert all(z[i]==0 for i in outside)
    K=F(9999,10000);norm=J[S.index(5)][S.index(5)]-K
    assert norm==F(190011,110000)>0
    old_bounds={str(i):str(-int(i==5))for i in S}
    endpoint=[]
    for wall in [1,10]:
        j=S.index(wall);u0=4*J[j][j]-2
        lead=J[j][j]*K-z[j]*z[j]
        margin=u0*lead-2*z[j]*z[j]
        assert lead>0 and margin>0
        endpoint.append(dict(wall=wall,J=str(J[j][j]),z=str(z[j]),u_lower=str(u0),
                             linear_coefficient=str(lead),strict_margin=str(margin)))
    zero=(F(0),F(0));one=(F(1),F(0));sqrt2=(F(0),F(1))
    slots=[1,2,9,10];domains=[[zero,one,sqrt2],[zero,one],[zero,one],[zero,one,sqrt2]]
    kept=[];rejected=[]
    for co in product(*domains):
        d=[zero]*18
        for i,x in zip(slots,co):d[S.index(i)]=x
        cost=zero;pay=zero
        for i in inside:
            pay=add(pay,scale(d[i],z[i]))
            for j in inside:cost=add(cost,scale(mul(d[i],d[j]),J[i][j]))
        key=[0 if x==zero else 3 if x==one else 4 for x in co]
        if sign(add(cost,(-small,F(0))))>=0:
            rejected.append(dict(labels=key,cost=enc(cost)));continue
        assert cost[1]==0
        margin=add((K*cost[0],F(0)),scale(mul(pay,pay),-1))
        if cost[0]:assert sign(margin)>0
        else:assert pay==zero and margin==zero
        kept.append(dict(labels=key,A10_cost=str(cost[0]),escape_A10_pair_magnitude=enc(pay),
            strict_squared_margin=enc(margin),
            proof='u>A10_cost>0 and K*A10_cost>z_A^2'if cost[0]else'u>0 and z_A=0; no quantitative gap required'))
    assert len(kept)==8 and len(rejected)==28
    expected=[[0,0,0,0],[3,0,0,0],[4,0,0,0],[0,3,0,0],[0,0,3,0],[0,0,0,3],[0,0,0,4],[3,0,0,3]]
    assert sorted(r['labels']for r in kept)==sorted(expected)
    paths.append(f'scripts/{TAG}.py')
    return dict(status='PASS_DEGREE50_AB_BASE_AND_ALL_CONSECUTIVE_ARCS_COMMON_WEAK_ESCAPE',
        sources={p:sha(p)for p in paths},recursive_manifest_audit=manifests,
        Galois_upgrade=group,centers=centers,
        triple=dict(Gram_outer_polynomial=list(map(str,g)),determinant=list(map(str,det)),
            adjugate_quadratic_coefficients={k:list(map(str,v))for k,v in numerator.items()},
            outer_sum_norm=list(map(str,outer_sum_norm)),
            norm_formula='A^2/2-((B+D)^2+4*C*B*D)/(8*C*(C+1))',
            abelian_polynomial='8*C*(C+1)*(A^2/2-2+Q)-(B+D)^2-4*C*B*D',
            nonzero_constant='-(B+D)^2',zero_outer_next_walls=next_walls,
            u_upper='2/C = 2+4/a^2'),
        support=S,support_Gram_LDL_pivots=support_pivots,
        all_label_budget=dict(total_support_cost_upper=str(total),
            A10_cost_upper_if_dotted_outside=str(small),minimum_inverse_diagonal='1/2'),
        weak_escape=dict(frame_coefficients=['1','-1/100','0'],dual_coefficients=list(map(str,e0)),
            norm=str(norm),frame_negative_norm=str(K),support_pairings=old_bounds,
            off_support_old134_strict=True,strict_center_pair_negative=True,
            derivation='E0=Eeps+eps*sum_{j!=5}omega_j; J>=0 and acuteness give every added term nonpositive offS and negative on a strict center.'),
        dotted_A10_endpoint_cases=endpoint,A10_slots=slots,A10_profiles=kept,A10_excluded_profiles=rejected,
        conclusion='Every coherent norm2 cap retaining old133 with finite/parallel contacts to any consecutive triple pairs strictly negatively with E0. Any finite coherent collection of these caps, including the134th cap, leaves infinite volume after a small strict-interior perturbation of E0.',
        scope='All24 consecutive triple centers; unequal outer contacts allowed; all integer finite labels and parallel contacts; arbitrary real coefficients and all18 support Coxeter contacts. No arithmeticity of the cap is assumed. Does not cover caps without such a consecutive triple or changes to retained oldwalls. Degree50 modular arithmetic and original physical chamber are inherited from independently reviewed frozen dependencies.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    d=run();p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],len(d['A10_profiles']),len(d['dotted_A10_endpoint_cases']),flush=True)
