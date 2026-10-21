#!/usr/bin/env python3
"""Exact finite budgets with one added quarter coupon on E6 and D11.

The local geometric justification is a separate theorem. This package proves
the conditional counting statement, with every last-item upper inequality and
an actual attaining forest, using the frozen prior finite payment rule.
"""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_a2_d11_full_coupon_forest_budget_20261017 as B

ROOT=B.ROOT
TAG='fresh_e6_d11_quarter_coupon_forest_budget_20261017'


def cross(a,b):
    return B.cross(a,b)+(2 if {a,b}=={('E',6),('D',11)} else 0)


def rep(t):
    return ('E',6) if t[1]==6 else B.representative(t)


def residual_type(r):
    return rep(('A' if r<=9 else 'D',r))


def costs(key):
    c=B.residual_costs(key)
    c[6]+=2*key[3]
    return c


def actual_cost(ts):
    return sum(map(B.base,ts))+sum(cross(a,b) for i,a in enumerate(ts) for b in ts[i+1:])


def dominance():
    types=([('A',r) for r in range(1,76)]+[('B',r) for r in range(2,76)]+
           [('D',r) for r in range(4,76)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    rows=[]
    for t in types:
        a=rep(t)
        assert a in types and a[1]==t[1] and B.base(a)>=B.base(t)
        gains=[]
        for other in types:
            gain=cross(a,other)-cross(t,other)
            assert gain>=0,(t,a,other)
            if gain:gains.append([B.label(other),gain])
        rows.append(dict(original=B.label(t),representative=B.label(a),
            doubled_base_gain=B.base(a)-B.base(t),positive_cross_gains=gains))
    return dict(types=len(types),partner_checks=len(types)**2,rows=rows)


def table(key):
    c=costs(key);v=[0]+[None]*B.LIMIT
    for r in B.RESIDUAL:
        for n in range(r,B.LIMIT+1):
            if v[n-r] is not None:
                q=v[n-r]+c[r]
                if v[n] is None or q>v[n]:v[n]=q
    assert all(x is not None for x in v)
    last=[0]+[max(r for r in B.RESIDUAL if r<=n and v[n]==v[n-r]+c[r]) for n in range(1,76)]
    return dict(interaction_key=list(key),doubled_values=v,last_rank=last)


def unpack(t,n):
    out=[]
    while n:
        r=t['last_rank'][n];assert 0<r<=n
        out.append(residual_type(r));n-=r
    return out


def large_types(row):
    return [('A' if r==17 else 'D',r) for r,c in zip(range(10,18),row['component_counts']) for _ in range(c)]


def run(verify=False):
    old=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    for p,h in old['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,p
    domain=B.count_domain();keys=sorted({tuple(r['interaction_key']) for r in domain})
    assert domain==old['certificate']['large_count_vectors']
    path=ROOT/f'results/{TAG}.json'
    saved=json.loads(path.read_text()) if verify else None
    if verify:
        tables={tuple(t['interaction_key']):t for t in saved['residual_tables']}
        assert len(tables)==len(saved['residual_tables']) and set(tables)==set(keys)
    else:tables={k:table(k) for k in keys}
    comparisons=bellman=attain=0
    for row in domain:
        large=large_types(row);key=tuple(row['interaction_key'])
        assert actual_cost(large)==row['doubled_base_cost']
        for r,w in costs(key).items():
            t=residual_type(r)
            assert w==B.base(t)+sum(cross(t,u) for u in large)
            comparisons+=1
    for key,t in tables.items():
        v=t['doubled_values'];last=t['last_rank'];c=costs(key)
        assert len(v)==len(last)==76 and v[0]==last[0]==0
        for n in range(1,76):
            assert isinstance(v[n],int)
            for r,w in c.items():
                if r<=n:assert v[n]>=v[n-r]+w;bellman+=1
            r=last[n];assert r in c and r<=n and v[n]==v[n-r]+c[r];attain+=1
    answers=[];global_checks=0
    for n in range(76):
        candidates=[]
        for row in domain:
            if row['used_rank']>n:continue
            t=tables[tuple(row['interaction_key'])];left=n-row['used_rank']
            value=row['doubled_base_cost']+t['doubled_values'][left]
            ts=B.canonical(large_types(row)+unpack(t,left))
            assert sum(r for _,r in ts)==n and actual_cost(ts)==value
            candidates.append((value,tuple(map(B.label,ts))))
        value,ts=max(candidates)
        assert all(value>=v for v,_ in candidates);global_checks+=len(candidates)
        assert value==old['certificate']['dimension_maxima'][n]['doubled_cost']
        answers.append(dict(dimension=n,doubled_cost=value,witness_components=ts,
            finite_curvature=str(B.curvature(n)),finite_slack=str(B.curvature(n)-F(value,2))))
    paths=[Path(__file__)]+[ROOT/f'{folder}/{B.TAG}.{ext}' for folder,ext in
        (('scripts','py'),('results','json'),('results','md'))]
    out=dict(status='PASS_E6_D11_QUARTER_COUPON_FINITE_BUDGET',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        conditional_coupon='A quarter between each of the two long-arm endpoint leaves of whole E6 and each of the two short fork leaves of separate whole D11; total added cost1 per component pair.',
        geometric_face_coverage_proved_here=False,type_domination=dominance(),
        large_count_vectors=domain,residual_tables=[tables[k] for k in keys],dimension_maxima=answers,
        counts=dict(large_vectors=len(domain),interaction_tables=len(keys),interaction_formula_checks=comparisons,
            bellman_upper_inequalities=bellman,attaining_predecessors=attain,global_alternatives=global_checks,typed_global_attainers=76),
        unchanged_maxima_through=75,finite_max65='487',curvature65='480',finite_max66='495',curvature66='495',
        scope='Complete exact conditional finite counting with all inherited coupons and the single E6/D11 addition. Rank6 is represented by E6, not A6; residual doubled cost gains2*d11. All typed domination comparisons and Bellman upper/attainment checks pass. Inherited geometric payment theorems and complete large-count domain are bound inputs. No ideal budget, face-coverage proof, new ambient theorem or example is asserted here.')
    out=json.loads(json.dumps(out))
    if verify:assert out==saved
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],'all maxima unchanged through75',flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();run(args.verify)
