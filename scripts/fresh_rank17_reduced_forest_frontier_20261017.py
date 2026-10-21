#!/usr/bin/env python3
"""Conditional exact finite maxima for explicit reduced rank17 base costs.

This counts an abstract rule; geometric coverage of its proposed deletions
is a separate obligation. All existing cross coupons, including A11, remain.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_ideal5_weighted_dimension68_20261017 as G

ROOT=G.ROOT;TAG='fresh_rank17_reduced_forest_frontier_20261017'
OLD='fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
A11='fresh_a11_cycle_coupon_domination_20261017'
MAXIMUM=75


def solve(a17,d17):
    residual_ranks=tuple(range(1,10))+tuple(range(17,MAXIMUM+1))
    def base(r):
        return 2*d17 if r==17 else r*(r-1) if r<=16 else 16*r-64
    def direct(cs):
        c=Counter(cs);a=sum(i*c[f'A{i}']for i in range(1,6));t=sum(c[f'D{i}']for i in range(11,17))
        cost=sum(2*a17 if s=='A17' else base(int(s[1:])) for s in cs)
        return cost+a*(3*t+2*c['A17'])+2*(c['A1']*(c['D11']+c['D12'])+c['A2']*c['B10']+c['A4']*c['D13'])
    tables={};best=[None]*(MAXIMUM+1);large=[]
    def residual(key):
        if key in tables:return tables[key]
        t,u,v,b,c,w=key
        weights={r:base(r)+(r*(3*t+2*u)+(2*(b+w)if r==1 else 2*v if r==2 else 2*c if r==4 else 0) if r<=5 else 0)for r in residual_ranks}
        values=[0];witness=[()];last=[0]
        for n in range(1,MAXIMUM+1):
            val,ws,r=max((values[n-r]+weights[r],tuple(sorted(witness[n-r]+(r,))),r)for r in residual_ranks if r<=n)
            values.append(val);witness.append(ws);last.append(r)
        row=dict(coupon_counts=list(key),doubled_optima=values,last_rank=last)
        tables[key]=(values,witness,row);return tables[key]
    def visit(r,left,counts):
        if r<=17:
            for c in range(left//r+1):visit(r+1,left-r*c,counts+[c])
            return
        used=MAXIMUM-left
        key=(sum(counts[1:7]),counts[7],counts[0],counts[2],counts[3],counts[1])
        values,witness,_=residual(key)
        basecost=sum(c*(2*a17 if j==17 else base(j))for j,c in zip(range(10,18),counts))
        types=tuple(('A17'if j==17 else 'B10'if j==10 else f'D{j}')for j,c in zip(range(10,18),counts)for _ in range(c))
        large.append(dict(component_counts=counts,used_rank=used,coupon_counts=list(key),doubled_base_cost=basecost))
        for rest in range(left+1):
            ws=tuple(sorted(types+tuple(('A'if j<=9 else 'D')+str(j)for j in witness[rest]),key=lambda s:(int(s[1:]),s)))
            candidate=(basecost+values[rest],ws)
            if best[used+rest]is None or candidate>best[used+rest]:best[used+rest]=candidate
    visit(10,MAXIMUM,[])
    assert len(large)==1608 and len(tables)==675
    maxima=[]
    for n,(cost,ws)in enumerate(best):
        assert sum(int(s[1:])for s in ws)==n and direct(ws)==cost
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        maxima.append(dict(dimension=n,doubled_cost=cost,witness_components=list(ws),
            finite_curvature_coefficient=str(beta),finite_slack=str(beta-F(cost,2))))
    return dict(assumed_A17_base=a17,assumed_D17_base=d17,
        assumed_B17_base_upper_bound=d17,
        residual_ranks=residual_ranks,large_count_vectors=large,
        residual_tables=[tables[k][2]for k in sorted(tables)],dimension_maxima=maxima,
        counts=dict(large_vectors=len(large),residual_tables=len(tables),dimensions=76))


def run():
    inputs={}
    for stem in(OLD,OLD+'_review',A11,G.TAG):
        d=json.loads((ROOT/f'results/{stem}.json').read_text())
        assert d['status'].startswith('PASS')
        for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
        inputs[stem]=d
    variants=[solve(131,131),solve(114,127)]
    profiles=[]
    for n in range(64,70):
        rows=[G.profile(rs)for rs in G.terminal_shapes(n-1)]
        minimum=min(rows,key=lambda z:F(z['margin']))
        profiles.append(dict(ambient_dimension=n,profiles=len(rows),minimum_margin=minimum['margin'],
            minimizing_ranks=minimum['ranks'],all_positive=all(F(r['margin'])>0 for r in rows)))
    assert [(r['ambient_dimension'],r['all_positive'])for r in profiles]==[(64,False),(65,False),(66,False),(67,True),(68,True),(69,True)]
    assert profiles[3]['minimum_margin']=='851/165888'
    assert [variants[0]['dimension_maxima'][n]['finite_slack']for n in(67,68,69)]==['-3/2','3','16']
    assert [variants[1]['dimension_maxima'][n]['finite_slack']for n in(67,68,69)]==['15/2','15','24']
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in inputs
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_CONDITIONAL_RANK17_FINITE_BUDGET_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        variants=variants,ideal_comparisons=profiles,new_dimension_bound=None,
        conditional_conclusion='If the simultaneous rank17 reductions preserve face coverage without extra coupons, A131/B<=131/D131 fits n68 and would give ceiling67. The stronger A114/B<=127/D127 rule fits n67 and n68 and would give ceiling66. The unchanged ideal rule is positive there; the completed68 theorem supplies n>=69.',
        inherited_domination='Ranks other than17 retain their frozen type domination and coupons. A17 is kept separately; B17 must cost no more than D17. A11 cycle coupons are dominated by unchanged D11 exactly as before.',
        missing_obligation='Complete geometric face coverage for simultaneous finite A17/B17/D17 reductions, and final independently reviewed global composition.',
        scope='Exact abstract finite counting under explicitly assumed rank17 costs, with all inherited coupons. Not a face theorem or new dimension claim. The finite-type reduction and current convex-order lemma are inherited; scalar arithmetic imports the frozen68 producer.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(v['assumed_A17_base'],v['assumed_D17_base'],[(n,v['dimension_maxima'][n]['doubled_cost'],v['dimension_maxima'][n]['finite_slack'])for n in(66,67,68,69)])for v in out['variants']],flush=True)
