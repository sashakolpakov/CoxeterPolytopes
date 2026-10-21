#!/usr/bin/env python3
"""Exact graph and scalar certificate for a refinement of Prokhorov's bound.

No geometry/search producer is imported. Published geometric hypotheses and
the all-rank proofs are stated in the accompanying report. Integer tree-arm
counts are independent of the reviewer's proposed graph BFS implementation.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from math import comb
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension_bound_prokhorov_20261015'

def path_pairs(n,k):
    return sum(n-j for j in range(1,min(k,n-1)+1))

def arm_pairs(arms,k):
    # Center/arm, same-arm, and different-arm unordered pairs, respectively.
    ans=sum(min(a,k) for a in arms)
    ans+=sum(path_pairs(a,k) for a in arms)
    ans+=sum(sum(i+j<=k for i in range(1,a+1) for j in range(1,b+1))
             for z,a in enumerate(arms) for b in arms[z+1:])
    return ans

def row(name,n,arms=None):
    f7=path_pairs(n,7) if arms is None else arm_pairs(arms,7)
    f15=path_pairs(n,15) if arms is None else arm_pairs(arms,15)
    w3=2*f7+f15
    return dict(type=name,rank=n,arms=arms,pairs_at_distance_at_most7=f7,
                pairs_at_distance_at_most15=f15,three_times_weight=w3,
                deficit=29*n-w3)

def nikulin(n,l,k):
    a,b=n//2,(n+1)//2
    return F(comb(n-l,n-k)*(comb(a,l)+comb(b,l)),comb(a,k)+comb(b,k))

def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,list):return [enc(v)for v in x]
    if isinstance(x,dict):return {k:enc(v)for k,v in x.items()}
    return x

def run():
    small=[]
    for n in range(1,17):
        small.append(row('path',n))
        if n>=4:small.append(row('D',n,[1,1,n-3]))
        if n in (6,7,8):small.append(row('E',n,[1,2,n-4]))
    for x in small:
        assert x['deficit']>0
        assert 17*x['deficit']>=173*x['rank']
    # Labels on paths include A/B/C/F/H and all finite I2(m), without cutoff.
    tail=[]
    for n in (17,18,29,990,991):
        a=row('path',n);d=row('D',n,[1,1,n-3])
        for k in (7,15):
            assert path_pairs(n,k)==k*n-k*(k+1)//2
            assert arm_pairs([1,1,n-3],k)==k*n-k*(k+1)//2+1
        assert a['deficit']==176 and d['deficit']==173
        tail.extend([a,d])
    # A second disconnected-graph proof only needs the complete-graph bound.
    fallback=[dict(rank=r,deficit_lower=F(r*(61-3*r),2))for r in range(1,17)]
    assert all(x['deficit_lower']>=F(13*x['rank'],2)for x in fallback)
    assert F(13*27,2)>173
    # Complete exact rank-partition optimization is a check, not the tail proof.
    best_component={}
    for n in range(1,101):
        rows=[row('path',n)]
        if n>=4:rows.append(row('D',n,[1,1,n-3]))
        if n in (6,7,8):rows.append(row('E',n,[1,2,n-4]))
        best_component[n]=max(x['three_times_weight']for x in rows)
    dp=[0]
    for n in range(1,101):
        dp.append(max(dp[n-r]+best_component[r]for r in range(1,n+1)))
        if n>=17:assert dp[n]==29*n-173
    scalar=[]
    for n in range(18,1002):
        W=F(29*(n-1)-173,3)
        C=W/(n-1)
        assert C>0
        assert 96*C+68==996-F(5536,n-1)
        P=n*n-997*n+6532
        assert P==(n-991)**2+985*(n-991)+586
        if n>=991:assert P>0 and n>=96*C+68
        if n in (989,990,991,992,995,996):
            scalar.append(dict(dimension=n,edge_rank=n-1,weight_upper=W,
                               C=C,strict_dimension_upper=96*C+68,
                               contradiction_polynomial=P,contradiction=P>=0))
    # Reconstruct the incidence calculation independently of the simplified lemma.
    parity=[]
    for n in range(18,1002):
        A=nikulin(n,1,3);B=nikulin(n,2,3)
        p=F(3*n+6,4*(n-3))if n%2==0 else F(3*n+9,4*(n-2))
        W=F(29*n-202,3)
        gap=7-B-p-F(2)*W*A/((n-1)*(n-2))
        if n%2==0:
            assert A==F(12*(n-1),n-4)and B==F(6*(n-2),n-4)
            num=n**3-1000*n*n+9604*n-19824
            den=4*(n-4)*(n-3)*(n-2)
        else:
            assert A==F(12*n,n-3)and B==F(6*(n-1),n-3)
            num=n**3-997*n*n+6679*n-147
            den=4*(n-3)*(n-2)*(n-1)
        assert gap==F(num,den)
        if n>=991:assert gap>0
        if 989<=n<=996:parity.append(dict(dimension=n,lower_minus_upper=gap))
    # Positive coefficients prove the entire parity tails beyond the finite checks.
    # Substitute n=992+t in the even numerator and n=991+t in the odd numerator.
    even_shift=[1634832,977796,1976,1]
    odd_shift=[726256,976868,1976,1]
    for t in (0,1,2,100):
        n=992+t
        assert sum(c*t**i for i,c in enumerate(even_shift))==n**3-1000*n*n+9604*n-19824
        n=991+t
        assert sum(c*t**i for i,c in enumerate(odd_shift))==n**3-997*n*n+6679*n-147
    src=f'scripts/{TAG}.py'
    return enc(dict(status='PASS_EXACT_PROKHOROV_WEIGHT_REFINEMENT',
        sources={src:sha256((ROOT/src).read_bytes()).hexdigest()},
        primary_reference='M.N. Prokhorov, Math. USSR-Izv.28(1987),401-411, Lemmas1-2 and Proposition5',
        primary_url='https://www.mathnet.ru/php/getFT.phtml?jrnid=im&option_lang=eng&paperid=1484&what=fullteng',
        scope='Finite-sided finite-volume hyperbolic Coxeter polytopes; all cusps and finite labels. Published geometric weighting lemma is inherited, not computationally re-proved.',
        weights='1 at graph distance1..7;1/3 at8..15;0 otherwise; unordered pairs',
        universal_graph_bound='3W(T)<=29N-173 for every spherical diagram rankN>=17',
        graph_coverage='Every irreducible finite Coxeter diagram is a path, D_n, or E6/E7/E8 as an unlabelled graph. All finite I2(m) share the rank2 path.',
        small_connected=small,large_rank_formula_checks=tail,
        large_rank_proof='For n>=17 and k=7,15: f_path=k*n-k(k+1)/2; f_D=f_path+1. Thus deficits176 and173. Small components have17*deficit>=173*rank; summing proves the disconnected case.',
        fallback_without_small_census=dict(small_rank_complete_graph_bounds=fallback,
            proof='If any component rank>=17 its deficit>=173. Otherwise total deficit>=13N/2>173 forN>=27. Every component has nonnegative deficit.'),
        partition_dp_maxima=dp,dimension_checks=scalar,
        universal_quadratic='n>=991 gives n²−997n+6532=(n−991)²+985(n−991)+586>0, contradicting published Lemma2 withC=W/(n−1)',
        direct_incidence_parity_checks=parity,positive_shift_coefficients=dict(even=even_shift,odd=odd_shift),
        derived_upper_dimension=990,
        novelty='No publication-priority claim. This is an exact refinement of a published proof, pending separate mathematical/source review before external theorem announcement.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved Prokhorov graph/scalar certificate')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
