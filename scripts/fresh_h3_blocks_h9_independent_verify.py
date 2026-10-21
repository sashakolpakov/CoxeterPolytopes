#!/usr/bin/env python3
"""Fraction-only H3 triage and H3+A6 root, alphabet, and escape audit."""
import argparse
from functools import lru_cache
from fractions import Fraction as F
from itertools import product,permutations
import json
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,neg,sub,scale,mul,div,conj,sign,integral,total,dot,mv,tr,mm,pair,matrix,
    eye,inverse,determinant,bounded_integers,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read,valid_coxeter_value
from scripts.fresh_multibasepoint_h9_targeted_verify import escape
from scripts.fresh_rank3_h9_newplane_verify import integral_square_root

H3=tuple(tuple(scale(O,2) if i==j else neg(P) if {i,j}=={0,1} else neg(O) if abs(i-j)==1 else Z for j in range(3)) for i in range(3))
A6=tuple(tuple(scale(O,2) if i==j else neg(O) if abs(i-j)==1 else Z for j in range(6)) for i in range(6))

def diagonal_block(*blocks):
    n=sum(map(len,blocks));out=[[Z]*n for _ in range(n)];offset=0
    for B in blocks:
        for i,row in enumerate(B):
            for j,v in enumerate(row):out[offset+i][offset+j]=v
        offset+=len(B)
    return tuple(map(tuple,out))

def reflected(r,B):
    c=mv(B,r);q=dot(r,c)
    return tuple(tuple(sub(O if i==j else Z,div(scale(mul(r[i],c[j]),2),q)) for j in range(len(r))) for i in range(len(r)))

@lru_cache(None)
def brute_block(B,target,q):
    I=inverse(B);alphabet=tuple(d for d in bounded_integers(scale(target,2),scale(target,2)) if sign(d)>=0 and integral(div(scale(d,2),q)))
    rows=[]
    for weights in product(alphabet,repeat=len(B)):
        y=mv(I,weights)
        if not all(map(integral,y)):continue
        energy=dot(weights,y)
        if sign(sub(target,energy))>=0 and sign(conj(sub(target,energy)))>=0:rows.append((energy,y))
    return tuple(rows)

def zero_block(B,q):return brute_block(B,q,q)

def zero_census(blocks,pool):
    roots=set()
    for q in pool:
        for rows in product(*(zero_block(B,q) for B in blocks)):
            if total(e for e,_ in rows)==q:roots.add((q,(Z,*(x for _,r in rows for x in r))))
    return roots

def finite_roots(B,expected):
    roots=set(eye(len(B)));todo=list(roots)
    while todo:
        r=todo.pop();products=mv(B,r)
        for j in range(len(B)):
            s=tuple(sub(x,products[j]) if i==j else x for i,x in enumerate(r))
            assert pair(s,B,s)==scale(O,2)
            if s not in roots:roots.add(s);todo.append(s)
        assert len(roots)<=expected
    assert len(roots)==expected
    assert sum(all(sign(x)>=0 for x in mv(B,r)) for r in roots)==1
    return expected

RES=((0,0),(1,0),(0,1),(1,1))
def rmul(x,y):return ((x[0]*y[0]+x[1]*y[1])%2,(x[0]*y[1]+x[1]*y[0]+x[1]*y[1])%2)
def radd(x,y):return ((x[0]+y[0])%2,(x[1]+y[1])%2)
def line(v):
    f=next(x for x in v if x!=(0,0));iv=next(x for x in RES if rmul(f,x)==(1,0))
    return tuple(rmul(iv,x) for x in v)
def rd(v,w):
    x=(0,0)
    for a,b in zip(v,w):x=radd(x,rmul(a,b))
    return x

def structure():
    d=read('results/fresh_h3_blocks_h9_structure_20261013.json')
    assert matrix(d['H3_gram'])==H3 and matrix(d['H3_inverse'])==inverse(H3)
    A=matrix(d['H3_to_2I3']);assert mm(tr(A),mm(tuple(tuple(scale(O,2) if i==j else Z for j in range(3)) for i in range(3)),A))==H3
    dual=tuple(map(parse,d['H3_discriminant_generator']));assert pair(dual,H3,dual)==scale(O,F(3,2))
    assert all(integral(x) for x in mv(H3,dual)) and any(not integral(x) for x in dual) and all(integral(scale(x,2)) for x in dual)
    lines={line(v) for v in product(RES,repeat=3) if any(x!=(0,0) for x in v) and rd(v,v)==(0,0)}
    assert sorted(lines)==[tuple(map(tuple,v)) for v in d['spatial_isotropic_discriminant_lines']] and len(lines)==5
    assert all(rd(v,w)!=(0,0) for v in lines for w in lines if v!=w)
    orbits={tuple(sorted({line(tuple(v[i] for i in p)) for p in permutations(range(3))})) for v in lines}
    assert sorted(map(len,orbits))==[2,3]
    lines4={line(v) for v in product(RES,repeat=4) if any(x!=(0,0) for x in v) and rd(v,v)==(0,0)}
    planes=set()
    for v,w in product(lines4,repeat=2):
        if v==w or rd(v,w)!=(0,0):continue
        vectors={tuple(radd(rmul(a,x),rmul(b,y)) for x,y in zip(v,w)) for a,b in product(RES,repeat=2)}
        planes.add(tuple(sorted({line(z) for z in vectors if any(x!=(0,0) for x in z)})))
    assert len(lines4)==21 and len(planes)==5
    assert sorted(lines4)==[tuple(map(tuple,v)) for v in d['even_host_isotropic_lines']]
    assert sorted(planes)==[tuple(tuple(map(tuple,v)) for v in pl) for pl in d['even_host_maximal_isotropic_planes']]
    assert set.intersection(*(set(p) for p in planes))=={((1,0),)*4}
    for row in d['hosts']:
        alpha=parse(row['alpha']);B=diagonal_block(((neg(alpha),),),H3,H3,H3)
        assert B==matrix(row['gram']) and determinant(B)==parse(row['determinant'])
        pool=tuple(map(parse,row['norm_pool']));assert pool==((O,scale(O,2),scale(O,4)) if alpha==P else (scale(O,2),scale(O,4)))
        assert all(integral(scale(x,2)) for line_ in inverse(B) for x in line_)
        expected=zero_census((H3,H3,H3),pool)
        assert expected=={(parse(r['norm']),tuple(map(parse,r['root']))) for r in row['zero_height_dominant_census']} and len(expected)==3
        targetalpha=parse(row['comparison_target_alpha']);diag=diagonal_block(((neg(targetalpha),),),eye(9))
        target=matrix(row['target_E8_gram'])
        for key,Btarget,comparison in (('to_diagonal',diag,'diagonal'),('to_E8_plus_scalar',target,'E8_plus_scalar')):
            T=matrix(row[key]);assert tuple(tuple(scale(x,2) for x in line_) for line_ in mm(tr(T),mm(Btarget,T)))==B
            record=next(r for r in row['explicit_group_distinctions'] if r['comparison']==comparison)
            j=record['source_simple_root_index'];root=tuple(O if i==j else Z for i in range(10))
            transformed=mv(T,root);assert transformed==tuple(map(parse,record['transported_root']))
            reflection=reflected(transformed,Btarget)
            assert all(integral(x) for line_ in reflected(root,B) for x in line_)
            for failure in record['nonintegral_reflection_entries']:
                x=reflection[failure['row']][failure['column']];assert x==parse(failure['value']) and not integral(x)
    return {'H3_finite_roots':finite_roots(H3,30),'triple_zero_census':3,'spatial_maximal_overclosures':5,
        'even_host_maximal_overclosures':5,'maximal_overclosure_block_permutation_orbits':[2,3],
        'field_similarities_verified':4,'explicit_reflection_group_distinctions_in_displayed_placements':4,
        'scope':'Does not rule out an unidentified conjugacy of the full groups. No positive-height H3-cubed search was run.'}

def valuation7(x):
    def rational_val(f):
        if not f:return 10000
        a,b=abs(f.numerator),f.denominator;v=0
        while a%7==0:a//=7;v+=1
        while b%7==0:b//=7;v-=1
        return v
    return min(rational_val(x[0]),rational_val(x[1]))

def h3_a6(name):
    d=read('results/fresh_h3_blocks_h9_a6_'+name+'_h30_20261013.json');alpha=parse(d['alpha'])
    gram=diagonal_block(((neg(alpha),),),H3,A6);assert gram==matrix(d['gram'])
    det=determinant(gram);assert det==parse(d['determinant'])==mul(scale(alpha,-14),sub(scale(O,2),P))
    assert valuation7(det)==1 and all((i*i-i-1)%7 for i in range(7))
    comparisons={'E8_phi':neg(P),'E8_2phi':scale(P,-2),'E8_sqrt5':sub(O,scale(P,2)),
        'E8_3phi':scale(P,-3),'E8_1plus4phi':sub(neg(O),scale(P,4)),'E8_minus1plus3phi':sub(O,scale(P,3)),
        'E8_even_natural':scale(P,-4),'E7_A2_phi':scale(P,-6),'E7_A2_2phi':scale(P,-12)}
    assert all(valuation7(v)==0 for v in comparisons.values())
    pool=tuple(map(parse,d['norm_pool']));expected=tuple(scale(O,q) for q in ((1,2,4,7,14,28) if alpha==P else (2,4,14,28)))
    assert pool==expected and all(integral(scale(x,14)) for row in inverse(gram) for x in row)
    # Inert prime divisors of28 give precisely these square-unit classes;
    # evenness removes odd classes for alpha=2phi.
    assert all((i*i-i-1)%2 for i in range(2))
    assert all(sign(determinant(tuple(row[:k] for row in B[:k])))>0 and sign(conj(determinant(tuple(row[:k] for row in B[:k]))))>0
        for B in (H3,A6) for k in range(1,len(B)+1))
    assert sign(alpha)>0 and sign(conj(alpha))<0
    zero=zero_census((H3,A6),pool)
    assert zero=={(parse(r['norm']),tuple(map(parse,r['root']))) for r in d['zero_height_dominant_census']} and len(zero)==2
    walls=matrix(d['all_walls']);inside=tuple(map(parse,d['strict_interior']));assert len(walls)==d['wall_count']
    assert d['completed'] and parse(d['completed_height'])==parse(d['height'])==scale(O,30)
    assert sign(inside[0])>0 and sign(pair(inside,gram,inside))<0
    assert walls[:9]==tuple(tuple(neg(O) if i==j else Z for i in range(10)) for j in range(1,10))
    cov=tuple(mv(gram,r) for r in walls);norms=tuple(dot(r,n) for r,n in zip(walls,cov));heights=[]
    for i,(r,n,q) in enumerate(zip(walls,cov,norms)):
        assert q in pool and all(map(integral,r)) and all(integral(div(scale(x,2),q)) for x in n)
        if q in (scale(O,4),scale(O,28)):assert not all(integral(scale(x,F(1,2))) for x in r)
        assert sign(dot(inside,n))<0
        height=div(mul(alpha,mul(r[0],r[0])),q);heights.append(height)
        assert sign(height)>=0 and sign(sub(scale(O,30),height))>=0
        for j in range(i):
            p=dot(r,cov[j]);assert sign(p)<=0
            assert valid_coxeter_value(div(scale(mul(p,p),4),mul(q,norms[j])))
    assert all(sign(sub(b,a))>=0 for a,b in zip(heights,heights[1:]))
    controls={(div(mul(alpha,mul(x,x)),q),q,x) for q in pool
        for x in bounded_integers(div(scale(q,30),alpha),neg(div(q,alpha)))
        if sign(x)>0 and integral(div(scale(mul(alpha,x),2),q))}
    assert controls=={(parse(row['height']),parse(row['norm']),parse(row['control'])) for row in d['levels']} and len(controls)==len(d['levels'])
    assert sum(row['roots'] for row in d['levels'])+9==len(walls)
    v=escape(d['hyperideal_witness'],walls,gram)
    reversal=(0,1,2,3,9,8,7,6,5,4)
    reverse=lambda r:tuple(r[i] for i in reversal)
    assert {reverse(r) for r in walls}==set(walls)
    assert tuple(tuple(gram[i][j] for j in reversal) for i in reversal)==gram
    average=tuple(scale(add(x,y),F(1,2)) for x,y in zip(v,reverse(v)))
    assert reverse(average)==average and sign(pair(average,gram,average))>0
    assert all(sign(pair(average,gram,r))<=0 for r in walls)
    h3weights=mv(H3,v[1:4]);a6weights=mv(A6,v[4:])
    assert all(x==Z for x in h3weights[1:]) and all(x==Z for i,x in enumerate(a6weights) if i!=2)
    return {'alpha':name,'walls':len(walls),'complete_zero_census':2,'controls':len(controls),'all_pairs_coxeter':True,
        'escape_norm':enc(pair(v,gram,v)),'prime7_discriminant_valuation':1,
        'A6_reversal_permutation':reversal,'original_escape_H3_fundamental_coefficients':list(map(enc,h3weights)),
        'original_escape_A6_fundamental_coefficients':list(map(enc,a6weights)),
        'A6_reversal_fixed_escape':list(map(enc,average)),'A6_reversal_fixed_escape_norm':enc(pair(average,gram,average)),
        'compared_prior_determinants_with_valuation_zero':{k:enc(v) for k,v in comparisons.items()}}

def verify():
    result={'status':'INDEPENDENT_H3_TRIAGE_AND_PRIME7_H9_AUDIT_PASS','structure':structure(),'A6_finite_roots':finite_roots(A6,42),
        'hosts':[h3_a6(n) for n in ('phi','2phi')],
        'scope':'Independent zero-height and bounded time-control censuses, all roots/reflections/angles/escapes, structural similarities and determinant distinction. Positive-height spatial enumeration is source-reviewed, not independently rerun in full.'}
    assert integral_square_root(2,0) is None
    result['two_new_host_determinant_ratio_nonsquare']=2
    return result

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args();result=verify();print(json.dumps(result,sort_keys=True))
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
