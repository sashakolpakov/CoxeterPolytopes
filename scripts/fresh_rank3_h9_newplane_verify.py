#!/usr/bin/env python3
"""Independent exact replay of the new binary slice and 28-subspace screen.

Standard-library Fraction arithmetic only; the full bounded binary search is
repeated with direct integer square tests, without project Golden arithmetic.
"""
import argparse
import json
from itertools import combinations
from math import isqrt
from pathlib import Path
from scripts.fresh_h4_blocks_independent_verify import (
    Z,O,P,parse,enc,add,neg,sub,scale,mul,div,conj,sign,integral,total,dot,mv,tr,mm,
    pair,matrix,inverse,
)
from scripts.fresh_multibasepoint_h9_verify import ROOT,read
from scripts.fresh_multibasepoint_h9_targeted_verify import escape

def integer_sign(a,b):
    u=2*a+b
    if not b:return (u>0)-(u<0)
    if not u or (u>0)==(b>0):return (b>0)-(b<0)
    return ((u>0)-(u<0))*((u*u>5*b*b)-(u*u<5*b*b))

def integral_square_root(a,b):
    """Solve (c+d*phi)^2=a+b*phi by norm and trace integer squares."""
    if a<0:return None
    n=a*a+a*b-b*b
    if n<0:return None
    k=isqrt(n)
    if k*k!=n:return None
    trace=2*a+b
    for signed in set((k,-k)):
        d2=trace+2*signed;fiveb2=trace-2*signed
        if d2<0 or fiveb2<0 or fiveb2%5:continue
        d=isqrt(d2);v=isqrt(fiveb2//5)
        if d*d!=d2 or 5*v*v!=fiveb2:continue
        for signed_v in set((v,-v)):
            if d*signed_v!=b or (d-signed_v)%2:continue
            c=(d-signed_v)//2
            assert (c*c+signed_v*signed_v,2*c*signed_v+signed_v*signed_v)==(a,b)
            return c,signed_v
    return None

def independent_binary_census(limit,q):
    # beta'=82+107phi'>2, so |t'|<1. For t>0 this implies b>=1,
    # b<(limit+1)/sqrt5<(limit+1)/2, and a=floor(b/phi) or ceil(b/phi).
    total_t=1;parity=0;squares=0;solutions=[]
    for b in range(1,(limit+1)//2+2):
        f=(isqrt(5*b*b)-b)//2
        for a in (f,f+1):
            if integer_sign(limit-a,-b)<0:continue
            aa=a*a+b*b;bb=2*a*b+b*b
            c=82*aa+107*bb;d=107*aa+189*bb
            if integer_sign(q-c-d,d)<0:continue
            total_t+=1
            if a%2 or b%2!=(1 if q==1 else 0):continue
            parity+=1
            if (d-c+q)%2 or (c-q)%2:continue
            ya=(d-c+q)//2;yb=(c-q)//2
            if integer_sign(ya,yb)<=0:continue
            squares+=1;x=integral_square_root(ya,yb)
            if x is not None:solutions.append({'t':(a,b),'x':x})
    return {'norm':q,'alphabet_size_including_zero':total_t,'parity_admissible':parity,
        'integral_square_tests':squares,'integral_solutions':len(solutions),'solutions':solutions}

def check_screen(base,gram,walls,v):
    screen=read('results/fresh_rank4_h9_slice_screen_20261013.json')
    assert screen['source']=='results/fresh_rank3_h9_cones_20261013.json' and screen['source_case']==0
    assert screen['wall_count']==len(walls)==382 and tuple(map(parse,screen['target_escape']))==v
    ei=inverse(tuple(tuple(gram[i][j] for j in range(1,9)) for i in range(1,9)))
    omega=tuple((Z,*(ei[k][j] for k in range(8)),Z) for j in range(8))
    time=(O,*(Z for _ in range(9)));scalar=(*(Z for _ in range(9)),O)
    cov=tuple(mv(gram,r) for r in walls);target=mv(gram,v)
    seen=set();closed=0;geometric=[]
    for row in screen['slices']:
        i,j=row['fundamental_weight_indices'];assert 0<=i<j<8 and (i,j) not in seen;seen.add((i,j))
        basis=(time,omega[i],omega[j],scalar)
        if row['status']=='EXACT_SLICE_OBSTRUCTION':
            ids=row['wall_indices'];coeff=tuple(map(parse,row['nonnegative_coefficients']))
            assert len(ids)==len(coeff) and all(0<=k<len(walls) for k in ids) and all(sign(c)>=0 for c in coeff)
            for b in basis:assert dot(b,target)==total(mul(c,dot(b,cov[k])) for c,k in zip(coeff,ids))
            closed+=1
        else:
            assert row['status']=='EXACT_GEOMETRIC_CUTTER' and (i,j)==(4,7)
            r=tuple(map(parse,row['normal']));coeff=tuple(map(parse,row['slice_coefficients']))
            assert r==tuple(total(mul(c,b[k]) for c,b in zip(coeff,basis)) for k in range(10))
            q=pair(r,gram,r);assert q==parse(row['norm'])==sub(scale(O,16),scale(P,9)) and sign(q)>0
            assert pair(r,gram,v)==parse(row['escape_pairing']) and sign(pair(r,gram,v))>0
            assert all(sign(dot(r,n))<=0 for n in cov)
            inside=tuple(map(parse,base['cases'][0]['strict_interior']));assert sign(pair(r,gram,inside))<0
            assert all(map(integral,r)) and integral(div(O,r[-1])) # Primitive due to a unit scalar coordinate.
            assert any(not integral(div(scale(c,2),q)) for c in mv(gram,r))
            geometric.append({'fundamental_weight_indices':[i,j],'norm':enc(q),'reflective':False,'primitive':True})
    assert seen==set(combinations(range(8),2)) and closed==27 and len(geometric)==1
    return {'exact_closed_slices':closed,'geometric_cutters':geometric}

def rational_norm_witnesses():
    d=read('results/fresh_rank3_h9_newplane_rational_norms_20261013.json')
    beta=add(scale(O,82),scale(P,107));alpha=scale(P,2);delta=mul(alpha,beta)
    assert delta==parse(d['extension_radicand'])
    for row in d['solutions']:
        q=row['norm'];a=parse(row['relative_constant']);b=parse(row['relative_generator_coefficient'])
        t=parse(row['t']);x=parse(row['x'])
        assert sub(mul(a,a),mul(delta,mul(b,b)))==scale(beta,q)
        assert t==neg(div(a,beta)) and x==b and sign(t)>0 and sign(x)>0
        assert sub(mul(beta,mul(t,t)),mul(alpha,mul(x,x)))==scale(O,q)
        assert not integral(t) and not integral(x)
    assert {r['norm'] for r in d['solutions']}=={1,2}
    return {'exact_rational_solutions':2,'integral_solutions':0}

def verify():
    d=read('results/fresh_rank3_h9_newplane_20261013.json');base=read(d['source'])
    assert d['source_case']==0 and d['completed_parameter_box']
    gram=matrix(base['gram']);assert gram==matrix(d['gram'])
    case=base['cases'][0];walls=matrix(case['all_walls']);v=escape(case['hyperideal_witness'],walls,gram)
    alpha=scale(P,2);assert parse(d['alpha'])==alpha
    ei=inverse(tuple(tuple(gram[i][j] for j in range(1,9)) for i in range(1,9)))
    direction=tuple(add(scale(ei[i][4],2),mul(scale(P,2),ei[i][7])) for i in range(8))+(add(O,P),)
    assert direction==tuple(map(parse,d['direction'])) and integral(div(O,direction[-1]))
    z=(Z,*direction);beta=pair(z,gram,z)
    assert beta==parse(d['direction_norm'])==add(scale(O,82),scale(P,107))
    assert sign(sub(conj(beta),scale(O,2)))>0
    lower=parse(d['compatibility_lower_ratio']);upper=parse(d['escape_upper_ratio'])
    assert lower==add(scale(O,4),scale(P,3)) and upper==div(pair(z,gram,v),mul(alpha,v[0]))
    values=[div(pair(z,gram,r),mul(alpha,r[0])) for r in walls if r[0]!=Z]
    assert lower in values and all(sign(sub(lower,x))>=0 for x in values) and sign(sub(upper,lower))>0
    reports=[]
    for q,source in zip((1,2),d['norm_searches'],strict=True):
        row=independent_binary_census(d['positive_t_bound'],q)
        assert parse(source['norm'])==scale(O,q)
        for key in ('alphabet_size_including_zero','parity_admissible','integral_square_tests','integral_solutions'):
            assert row[key]==source[key],(key,row[key],source[key])
        assert row['integral_solutions']==0
        reports.append(row)
    assert d['cuts']==d['cases']==[]
    return {'status':'INDEPENDENT_NEW_BINARY_SLICE_AND_28_SUBSPACES_PASS','positive_t_bound':d['positive_t_bound'],
        'binary_searches':reports,'subspace_screen':check_screen(base,gram,walls,v),'rational_norm_witnesses':rational_norm_witnesses(),
        'scope':'Complete binary parameter bound independently rerun. All27 fixed-escape slice obstructions replayed. No unbounded conclusion follows from this bounded replay; the separate two-projection certificate excludes completion confined to S or U(S).'}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path);a=p.parse_args();result=verify()
    print(json.dumps(result,sort_keys=True))
    if a.output:a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
