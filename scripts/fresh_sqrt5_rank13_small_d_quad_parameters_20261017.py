#!/usr/bin/env python3
"""Finite-anchor arithmetic classification of the joint small-component D quad list."""
import argparse,hashlib,json,math
from collections import Counter
from fractions import Fraction as F
from itertools import combinations,permutations
from pathlib import Path
import fresh_sqrt5_adjacent_discount_witness_20261016 as A

ROOT=Path(__file__).resolve().parents[1];TAG='fresh_sqrt5_rank13_small_d_quad_parameters_20261017'
INPUT='results/fresh_sqrt5_rank13_small_d_joint_quads_20261017.json'
ANCHOR='results/fresh_sqrt5_compact_quad_finite_anchors_20261017.json'
REVIEW='results/fresh_sqrt5_compact_quad_finite_anchors_review_20261017.json'
THIRD='results/fresh_sqrt5_rank13_triangle_third_costs_20261017.json'
K=A.K;Z=A.Z;PHI=K(F(1,2),F(1,2))

def rsqrt(x):
    if x<0:return None
    a=math.isqrt(x.numerator);b=math.isqrt(x.denominator)
    return F(a,b)if a*a==x.numerator and b*b==x.denominator else None

def ksqrt(x):
    if not x:return Z
    if x.sign()<0:return None
    if not x.b:
        q=rsqrt(x.a)
        if q is not None:return K(q)
        q=rsqrt(x.a/5)
        return K(0,q)if q is not None else None
    n=rsqrt(x.a*x.a-5*x.b*x.b)
    if n is None:return None
    for sign in(-1,1):
        q=rsqrt((x.a+sign*n)/2)
        if q is None or not q:continue
        r=K(q,x.b/(2*q));assert r*r==x
        return r if r.sign()>0 else -r
    return None

def solve(coeff):
    assert max(coeff,default=0)<=2 and min(coeff,default=0)>=0
    c,b,a=[coeff.get(i,Z)for i in range(3)]
    if not a:
        if b:return[-c/b],False
        return[],not bool(c)
    d=ksqrt(b*b-4*a*c)
    if d is None:return[],False
    roots=[(-b+d)/(2*a),(-b-d)/(2*a)]
    for r in roots:assert a*r*r+b*r+c==Z
    return list(dict.fromkeys(roots)),False

def decode(p):return{tuple(m):F(c)for m,c in p}
def add(p,q):
    d=p.copy()
    for m,c in q.items():d[m]=d.get(m,F(0))+c
    return{m:c for m,c in d.items()if c}
def mul(p,q):
    d={}
    for(a,b),c in p.items():
        for(e,f),h in q.items():d[a+e,b+f]=d.get((a+e,b+f),F(0))+c*h
    return{m:c for m,c in d.items()if c}
def determinant(M):
    ans={};n=len(M)
    for p in permutations(range(n)):
        term={(0,0):F((-1)**sum(p[i]>p[j]for i in range(n)for j in range(i+1,n)))}
        for i,j in enumerate(p):term=mul(term,M[i][j])
        ans=add(ans,term)
    return ans

def exponent(x,n):
    if n>=0:
        out=K(1)
        for _ in range(n):out=out*x
        return out
    return K(1)/exponent(x,-n)

def valid_square(x):
    return x.integral()and x.sign()>0 and x.conj().sign()>0 and(4-x.conj()).sign()>0 and ((x-4).sign()>0 or x in(1+PHI,2+PHI))

def valid_direct(x):
    return x.integral()and x.sign()<0 and valid_square(x*x)

def coupling(row):
    n=row['nodes'];g={i:set()for i in range(n)}
    for(i,j),m in zip(combinations(range(n),2),row['labels']):
        if m==3:g[i].add(j);g[j].add(i)
    unseen=set(g);where={};count=0
    while unseen:
        seen={min(unseen)};queue=list(seen)
        for i in queue:
            for j in g[i]:
                if j not in seen:seen.add(j);queue.append(j)
        unseen-=seen
        for i in seen:where[i]=count
        count+=1
    return sorted((where[0],where[2]))==sorted((where[1],where[3]))

def parameters(p,fu,fv,coupled,direct,squares):
    result={};stats=Counter()
    for axis in(0,1):
        known=(fu,fv)[axis];other_known=(fu,fv)[1-axis]
        for anchor in direct if known else squares:
            coeff={}
            for powers,value in p.items():
                a,b=powers[axis],powers[1-axis]
                if known:
                    if not other_known:assert b%2==0
                    power=b if other_known else b//2
                    term=K(value)*exponent(anchor,a)
                elif other_known:
                    assert a%2==0;power=b;term=K(value)*exponent(anchor,a//2)
                elif coupled:
                    assert(a-b)%2==0;power=b;term=K(value)*exponent(anchor,(a-b)//2)
                else:
                    assert a%2==b%2==0;power=b//2;term=K(value)*exponent(anchor,a//2)
                coeff[power]=coeff.get(power,Z)+term
            roots,identity=solve(coeff);stats['anchor_tests']+=1
            if identity:
                stats['positive_null_identity_exclusions']+=1;continue
            for root in roots:
                values=[None,None];sq=[None,None];uv=None
                if known:values[axis]=anchor;sq[axis]=anchor*anchor
                else:sq[axis]=anchor
                if other_known:values[1-axis]=root;sq[1-axis]=root*root
                elif not known and coupled:
                    if not root.integral()or root.sign()<=0:continue
                    uv=root;sq[1-axis]=root*root/anchor
                else:sq[1-axis]=root
                if not all(valid_square(x)for x in sq):continue
                if any(x is not None and not valid_direct(x)for x in values):continue
                if fu and fv:uv=values[0]*values[1]
                key=tuple(z for x in sq for z in(x.a,x.b))
                result[key]=dict(values=values,squares=sq,uv=uv)
    stats['distinct_integral_parameter_pairs']=len(result)
    return list(result.values()),stats

def evaluate(p,param):
    u,v=param['values'];x,y=param['squares'];t=param['uv'];answer=Z
    for(a,b),c in p.items():
        if u is not None and v is not None:z=exponent(u,a)*exponent(v,b)
        elif u is not None:
            assert b%2==0;z=exponent(u,a)*exponent(y,b//2)
        elif v is not None:
            assert a%2==0;z=exponent(x,a//2)*exponent(v,b)
        elif a%2==b%2==0:z=exponent(x,a//2)*exponent(y,b//2)
        else:
            assert a%2==b%2==1 and t is not None
            z=t*exponent(x,(a-1)//2)*exponent(y,(b-1)//2)
        answer=answer+c*z
    return answer

def run():
    source=json.loads((ROOT/INPUT).read_text());domain=json.loads((ROOT/ANCHOR).read_text())
    review=json.loads((ROOT/REVIEW).read_text());third=json.loads((ROOT/THIRD).read_text())
    for doc in(source,domain,review,third):
        for p,h in doc['sources'].items():assert hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    direct=[K(*map(F,x['value']))for x in domain['direct_product_anchors']]
    squares=[K(*map(F,x['square']))for x in domain['squared_product_anchors']]
    records=[];counts=Counter();cache={}
    for rec in source['records']:
        symbolic=rec['symbolic'];p=decode(symbolic['determinant']);fu=symbolic['u_in_K'];fv=symbolic['v_in_K']
        coupled=not fu and not fv and coupling(rec['row'])
        key=(tuple(sorted(p.items())),fu,fv,coupled)
        if key not in cache:cache[key]=parameters(p,fu,fv,coupled,direct,squares)
        pars,stats=cache[key];counts.update(stats)
        M=[[decode(q)for q in row]for row in symbolic['residual']]
        opp=[determinant([[M[i][j]for j in pair]for i in pair])for pair in[(0,2),(1,3)]]
        edge=[determinant([[M[i][j]for j in ids]for i in ids])for ids in[(0,1,2),(0,1,3)]]
        accepted=[];rejected=Counter()
        for param in pars:
            assert evaluate(p,param)==Z
            ov=[evaluate(q,param)for q in opp]
            if any(x.sign()>=0 for x in ov):rejected['opposite_pair_not_hyperbolic']+=1;continue
            ev=[evaluate(q,param)for q in edge]
            if any(x.sign()!=-1 or x.conj().sign()!=1 for x in ev):rejected['physical_or_conjugate_edge_signature']+=1;continue
            accepted.append(dict(direct_products=[None if x is None else x.enc()for x in param['values']],
                squared_products=[x.enc()for x in param['squares']],mutual_product=None if param['uv']is None else param['uv'].enc(),
                opposite_determinants=[x.enc()for x in ov],edge_determinants=[x.enc()for x in ev]))
        records.append(dict(index=rec['index'],u_in_K=fu,v_in_K=fv,mutual_product_cycle=coupled,
            parameter_statistics=dict(stats),rejected=dict(rejected),accepted=accepted))
        counts['records_with_parameters']+=bool(accepted)
        if rec['index']%100==0:print('SMALL_D_QUAD_PARAMETERS',rec['index'],len(source['records']),counts['records_with_parameters'],flush=True)
    targets=[]
    for target in source['targets']:
        maps=[i for i,m in enumerate(target['maps'])if records[m['record']]['accepted']]
        pairs={tuple(target['maps'][i]['ports'])for i in maps}
        tri=next(t for t in third['targets']if t['index']==target['index'])
        # The third-cost certificate's refined corner set is intentionally
        # left external if its schema changes; this package reports quads.
        targets.append(dict(index=target['index'],forest=target['forest'],accepted_map_indices=maps,
            possible_quad_corners=[list(p)for p in sorted(pairs)],possible_quad_corner_count=len(pairs),
            refined_triangle_half_charge_corner_count=tri['half_charge_corner_count']))
    paths=[f'scripts/{TAG}.py',INPUT,ANCHOR,REVIEW,THIRD,'scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py']
    return dict(status='PASS_EXACT_FINITE_ANCHOR_JOINT_SMALL_D_QUADRILATERAL_PARAMETERS',
        sources={p:hashlib.sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},counts=dict(counts),records=records,targets=targets,
        completeness='Either member is in its appropriate8/59 anchor domain. Substitute both possibilities. A field-valued remaining product solves a quadratic overK; a square-only product solves a linear equation in its square. If both products have only their squares inK and close a common fixed-component cycle, solve instead for their integral mutual product. All exact K square roots are tested. Identity specializations are excluded by the frozen positive-null lemma.',
        opposite_pair_gate='Opposite residual pairs in a compact acute quadrilateral must be hyperbolic: a positive-definite pair would give another finite vertex inside all four halfspaces by projecting a common strict center, contradicting opposite incidence; a singular pair gives a feasible positive null ray. Both determinants are therefore strictly negative.',
        signature_gate='Every adjacent residual2-block is already rational positive definite. Each old-edge3-block has negative physical determinant and positive conjugate determinant, and the full4-block determinant vanishes. Thus its physical/conjugate inertias are(2,1,1)/(3,0,1).',
        scope='Complete necessary arithmetic parameter gate on the frozen D5+2D4,2D5+A3,3D4+A1-based joint quad census. Retained parameters are not a claim of global compact realization. This only controls quadrilaterals at these three vertex types; other selected and unselected types still require bounds. No dimension theorem is claimed.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args();data=run()
    out=ROOT/f'results/{TAG}.json'
    if a.verify:assert data==json.loads(out.read_text())
    else:out.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],[(t['forest'],t['possible_quad_corner_count'])for t in data['targets']],flush=True)
