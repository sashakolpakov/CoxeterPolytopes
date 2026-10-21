#!/usr/bin/env python3
"""Independent multiquadratic reconstruction of every allowed A2 profile."""
import argparse,json
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import sympy as s
ROOT=Path(__file__).resolve().parents[1];TAG='fresh_h20_a2_h18_coupling_20261014_review'
PRIMARY='results/fresh_h20_a2_h18_coupling_20261014.json'
SOURCE='data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'
B2='results/fresh_h20_b2_h18_coupling_20261014.json'
K=s.QQ.algebraic_field(s.sqrt(2),s.sqrt(3),s.sqrt(5));Z,O=K.zero,K.one
rt2,rt3,rt5=[K.from_sympy(s.sqrt(n))for n in (2,3,5)]
phi=(O+rt5)/2
RAD=[1,2,3,5,6,10,15,30]
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def dec(v):return K.from_sympy(s.Rational(v['a'])+s.Rational(v['b'])*s.sqrt(2))
def read(p):
    d=json.loads((ROOT/p).read_text())
    for f,h in d.get('sources',{}).items():assert sha(f)==h,f
    return d
@lru_cache(None)
def coefficients(v):
    expr=s.expand(K.to_sympy(v));out=[]
    for n in RAD[1:]:out.append(expr.coeff(s.sqrt(n)))
    a=s.expand(expr-sum(c*s.sqrt(n)for c,n in zip(out,RAD[1:])))
    assert a.is_Rational and all(c.is_Rational for c in out)
    return tuple(F(c)for c in [a]+out)
@lru_cache(None)
def interval(v):
    scale=10**40;lo=hi=F(0)
    for c,n in zip(coefficients(v),RAD):
        if n==1:a=b=F(1)
        else:
            a=F(isqrt(n*scale*scale),scale);b=a+F(1,scale)
            assert a*a<n<b*b
        lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
    return lo,hi
def sign(v):
    if not v:return 0
    lo,hi=interval(v);assert lo>0 or hi<0,(v,lo,hi)
    return 1 if lo>0 else -1
def dot(a,b):return sum((x*y for x,y in zip(a,b)),Z)
def mv(M,v):return [dot(r,v)for r in M]
def pair(a,b,Q):return dot(a,mv(Q,b))
def inverse(M):
    n=len(M);a=[r[:]+[O if i==j else Z for j in range(n)]for i,r in enumerate(M)]
    for j in range(n):
        i=next(i for i in range(j,n)if a[i][j]);a[j],a[i]=a[i],a[j]
        q=a[j][j];a[j]=[v/q for v in a[j]]
        for i in range(n):
            if i!=j:
                q=a[i][j];a[i]=[x-q*y for x,y in zip(a[i],a[j])]
    return [r[n:]for r in a]
def enc(v):return list(map(str,coefficients(v)))
def run():
    primary,raw,b2=map(read,(PRIMARY,SOURCE,B2));G=[[dec(x)for x in row]for row in raw['gram']]
    ids=list(range(19));Q0=[[G[i][j]for j in ids]for i in ids];J=inverse(Q0)
    roots=[mv(J,[G[j][i]for j in ids])for i in range(40)]
    assert all(pair(roots[i],roots[j],Q0)==G[i][j]for i in range(40)for j in range(40))
    edges=[(i,j)for i,j in combinations(range(40),2)if sign(G[i][j])<0 and 4*G[i][j]**2==G[i][i]*G[j][j]]
    seen={0}
    while True:
        more=seen|{j for i,j in edges if i in seen}|{i for i,j in edges if j in seen}
        if more==seen:break
        seen=more
    assert len(edges)==27 and len(seen)==24 and set(i for e in edges for i in e)==seen
    assert (0,1)in edges and (1,2)in edges and not G[0][2]
    eligible=sorted(set(range(40))-seen);assert eligible==b2['eligible_special']
    # Complete normalized norm2 inducing profiles: ordinary or single3/4/5.
    # The minimum two-contact cost is2; one label>=6 also costs>=2.
    contacts=[O,rt2,phi];residuals=[2*O-s.Rational(2,3)*t*t for t in contacts]
    factors=[2*rt3/3,rt2*rt3/3,(rt5-O)*rt2*rt3/6]
    assert all(a*a==p and sign(p)>0 for a,p in zip(factors,residuals))
    assert 2*(s.Rational(2,3)*O)+2*s.Rational(1,3)*O==2*O
    assert s.Rational(2,3)*rt3*rt3==2*O
    # An ordinary/special source pi3 pair gives -sqrt(p/2), strictly
    # between -1 and0. Every pair of nonzero A2 weights has positive product.
    assert all(sign(p)>0 and sign(2*O-p)>0 for p in residuals)
    transition=[]
    ts={0:Z,3:O,4:rt2,5:phi};fs={0:rt2,3:factors[0],4:factors[1],5:factors[2]}
    for first in (0,3,4,5):
        for second in (0,3,4,5):
            for color in (['irrelevant']if not first or not second else['same','opposite']):
                g=-fs[first]*fs[second]/2+ts[first]*ts[second]*(s.Rational(2,3)if color=='same'else s.Rational(1,3))
                if g==-O:kind='pi3'
                elif not g:kind='right'
                elif sign(g)>0:kind='forbidden_positive'
                elif g==-O/3:kind='forbidden_minus_one_third'
                else:
                    assert sign(g)<0 and sign(g+O)>0
                    kind='forbidden_between_minus1_and0'
                transition.append(dict(first=first,second=second,color=color,relation=kind))
    assert transition==primary['pi3_transition_table']
    A2=[[2*O,-O],[-O,2*O]];weights=[[2*O/3,O/3],[O/3,2*O/3]]
    assert all(pair(weights[i],weights[j],A2)==(2*O/3 if i==j else O/3)for i in range(2)for j in range(2))
    Q=[row+[Z,Z]for row in Q0]+[[Z]*19+row for row in A2]
    # Symmetric elimination with exact one- and two-dimensional pivots.
    M=[r[:]for r in Q];pos=neg=zero=0
    while M:
        diagonal=next((i for i in range(len(M))if M[i][i]),None)
        if diagonal is not None:
            order=[diagonal]+[i for i in range(len(M))if i!=diagonal]
            M=[[M[i][j]for j in order]for i in order];q=M[0][0]
            if sign(q)>0:pos+=1
            else:neg+=1
            M=[[M[i][j]-M[i][0]*M[0][j]/q for j in range(1,len(M))]for i in range(1,len(M))]
        else:
            edge=next(((i,j)for i in range(len(M))for j in range(i)if M[i][j]),None)
            if edge is None:zero+=len(M);break
            order=list(edge)+[i for i in range(len(M))if i not in edge]
            M=[[M[i][j]for j in order]for i in order];q=M[0][1];pos+=1;neg+=1
            M=[[M[i][j]-(M[i][0]*M[1][j]+M[i][1]*M[0][j])/q for j in range(2,len(M))]for i in range(2,len(M))]
    assert (pos,neg,zero)==(20,1,0)
    bank=[[Z]*19+[O,Z],[Z]*19+[Z,O]];profile_records=[]
    for i,r in enumerate(roots):
        assert G[i][i]in(O,2*O)
        norm_sqrt=O if G[i][i]==O else rt2
        ordinary=[rt2*x/norm_sqrt for x in r]+[Z,Z];bank.append(ordinary)
        if i in eligible:
            for label,t,factor in zip((3,4,5),contacts,factors):
                for color,w in enumerate(weights):
                    n=[factor*x/norm_sqrt for x in r]+[-t*x for x in w]
                    bank.append(n);profile_records.append(dict(source=i,label=label,color=color))
    assert len(bank)==138 and all(pair(n,n,Q)==2*O for n in bank)
    x=[dec(z)for z in b2['source_point']];h=[sum((r[j]for r in roots),Z)for j in range(19)]
    xn=pair(x,x,Q0);assert xn==K.from_sympy((s.Rational(3857)-2800*s.sqrt(2))/529)
    gaps=[pair(x,r,Q0)for r in roots]
    assert all(sign(g)<=0 for g in gaps)
    assert all(sign(gaps[i]+(O if G[i][i]==O else rt2))<=0 for i in eligible)
    assert pair(h,h,Q0)==-2952*O-1280*rt2 and pair(x,h,Q0)==-16*O-40*rt2
    assert all(sign(pair(h,r,Q0))<0 for r in roots)
    a=K.from_sympy(s.Rational(3207,1000));epsilon=K.from_sympy(s.Rational(1,10**9))
    margin=a*a*residuals[-1]-phi*phi;assert sign(margin)>0
    unperturbed=a*a*xn+2*O;assert sign(unperturbed)>0
    E=[a*u+epsilon*v for u,v in zip(x,h)]+[-O,-O]
    H=[100*u for u in h]+[-O,-O]
    enorm=pair(E,E,Q);hnorm=pair(H,H,Q);eh=pair(E,H,Q)
    assert sign(enorm)>0 and sign(hnorm)<0 and sign(eh)<0
    ep=[pair(E,n,Q)for n in bank];hp=[pair(H,n,Q)for n in bank]
    assert all(sign(v)<0 for v in ep+hp)
    paths=[PRIMARY,SOURCE,B2,f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_ALL_DIRECT_A2_H18_ESCAPE',sources={p:sha(p)for p in paths},
        arithmetic='Separate SymPy multiquadratic field; rigorous rational square-root intervals for every sign. No producer/project arithmetic imports.',
        scalar_basis=RAD,ambient_inertia=[pos,neg,zero],transition_cases=len(transition),source_gram_entries=1600,pi3_edges=[list(e)for e in edges],forced_ordinary=sorted(seen),
        induced_path=[0,1,2],eligible_special=eligible,bank_count=138,special_profiles=profile_records,
        escape_scale='3207/1000',center_perturbation='1/1000000000',worst_profile_margin=enc(margin),
        unperturbed_escape_norm=enc(unperturbed),escape_norm=enc(enorm),center_norm=enc(hnorm),escape_center_pairing=enc(eh),
        escape=[enc(v)for v in E],center=[enc(v)for v in H],
        bank_escape_pairings=[enc(v)for v in ep],bank_center_pairings=[enc(v)for v in hp],
        escape_norm_interval=list(map(str,interval(enorm))),
        scope='All seven normalized inducing profiles at every eligible fixed-source facet, all finite integer labels. The exact138-normal union has a common strict spacelike escape and strict timelike center; every coherent direct A2 assignment therefore has infinite volume. Extra noninducing walls and different sources remain outside.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['bank_count'],d['escape_norm_interval'])
