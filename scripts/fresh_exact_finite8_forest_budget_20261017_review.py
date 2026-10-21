#!/usr/bin/env python3
"""No-producer certificate review of exact finite-forest payment maxima.

Checks type domination, the complete large-component domain, every Bellman
upper-bound alternative and an attaining predecessor, and every global
upper bound and attaining forest. No producer arithmetic is imported.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_exact_finite8_forest_budget_20261017'
TAG=PRIMARY+'_review'
INPUTS=('fresh_finite8_leaf9_half_coupon_budget_20261017',
        'fresh_finite8_a17_d11_coupon_budget_20261017',
        'fresh_finite8_a4_d13_neighbor_coupon_budget_20261017',
        'fresh_finite8_d12_b10_coupon_budget_20261017')
LIMIT=75

def digest(p):return sha256(p.read_bytes()).hexdigest()

def representative(r):
    if r<=9 or r==17:return('A',r)
    return('B',r)if r==10 else('D',r)

def base(t):
    kind,r=t
    if r<=17:return r*(r-1)
    assert kind in('A','B','D')
    return 16*r-(64 if kind=='D'else 68)

def cross(a,b):
    """Twice the weight on all unordered marked pairs between two components."""
    total=0
    for small,large in((a,b),(b,a)):
        if small[1]<=5:
            if large[0]=='D'and 11<=large[1]<=16:total+=3*small[1]
            elif large in(('A',13),('A',15)):total+=small[1]
            elif large==('A',17):total+=2*small[1]
        if(small,large)in((('A',1),('D',12)),(('A',2),('B',10)),(('A',4),('D',13))):
            total+=2
    return total

def forest_cost(types):
    return sum(map(base,types))+sum(cross(a,b)for i,a in enumerate(types)for b in types[i+1:])

def type_domination():
    # I2(*) denotes every I2(m),m>=4. Its rank/base/small-coupon profile is
    # independent of m; the exceptional m=3 is the separately present A2.
    types=([('A',r)for r in range(1,LIMIT+1)]+[('B',r)for r in range(2,LIMIT+1)]+
          [('D',r)for r in range(4,LIMIT+1)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    rows=[];tests=0
    for t in types:
        rep=representative(t[1]);assert rep in types and base(rep)>=base(t)
        deltas=[]
        for other in types:
            delta=cross(rep,other)-cross(t,other);assert delta>=0;tests+=1
            if delta:deltas.append(dict(other=list(other),doubled_increase=delta))
        rows.append(dict(original=list(t),representative=list(rep),doubled_base_increase=base(rep)-base(t),
            positive_cross_increases=deltas))
    return dict(types=len(types),pair_tests=tests,rows=rows,
        scope='The finite Coxeter classification and inherited same-component graph costs are assumed. All coupon coefficients are independently counted; I2(*) covers the infinite rank-two family.')

def verify(primary):
    assert primary['dimension_maximum']==LIMIT
    assert primary['large_coupon_ranks']==list(range(10,18))
    residual=list(range(1,10))+list(range(18,76))
    assert primary['residual_ranks']==residual
    for r in range(1,76):assert primary['representative_types'][str(r)]==''.join(map(str,representative(r)))
    # Enumerate by number of components and nondecreasing rank tuples,
    # independently of the producer's successive rank-count recursion.
    all_counts={}
    for number in range(8):
        for rs in combinations_with_replacement(range(10,18),number):
            used=sum(rs)
            if used>75:continue
            count=tuple(rs.count(r)for r in range(10,18));assert count not in all_counts
            key=(sum(11<=r<=16 for r in rs),rs.count(17),rs.count(10),rs.count(12),rs.count(13))
            all_counts[count]=(used,key,sum(base(representative(r))for r in rs))
    stored={tuple(r['component_counts']):r for r in primary['large_count_vectors']}
    assert len(stored)==len(primary['large_count_vectors'])==len(all_counts)==1608
    assert set(stored)==set(all_counts)
    for count,(used,key,cost)in all_counts.items():
        r=stored[count];assert r['used_rank']==used and tuple(r['coupon_counts'])==key and r['doubled_base_cost']==cost
    keys={x[1]for x in all_counts.values()};tables={tuple(r['coupon_counts']):r for r in primary['residual_tables']}
    assert len(tables)==len(primary['residual_tables'])==len(keys)==377 and set(tables)==keys
    table_checks=[];transitions=0;predecessors=0
    for key in sorted(keys):
        t,u,v,b,c=key;assert t>=b+c
        # Any D11..16 has the same broad small-component interaction. The
        # D12/D13 additions are retained separately in this witness multiset.
        large=[('D',11)]*(t-b-c)+[('D',12)]*b+[('D',13)]*c+[('A',17)]*u+[('B',10)]*v
        costs={r:base(representative(r))+sum(cross(representative(r),x)for x in large)for r in residual}
        row=tables[key];values=row['doubled_optima'];last=row['last_rank']
        assert len(values)==len(last)==76 and values[0]==last[0]==0
        checked=0
        for s in range(1,76):
            assert isinstance(values[s],int)
            for r in residual:
                if r>s:continue
                assert values[s]>=values[s-r]+costs[r],(key,s,r)
                checked+=1
            r=last[s];assert r in residual and 1<=r<=s
            assert values[s]==values[s-r]+costs[r],(key,s,'attainment')
            predecessors+=1
        transitions+=checked
        table_checks.append(dict(coupon_counts=list(key),doubled_optima=values,last_rank=last,
            verified_upper_bound_alternatives=checked,verified_attaining_predecessors=75))
    maxima=primary['dimension_maxima'];assert len(maxima)==76
    upper_tests=0;answers=[]
    for n,row in enumerate(maxima):
        assert row['dimension']==n
        opt=row['doubled_cost'];assert isinstance(opt,int)and 0<=opt<=16*n
        for used,key,cost in all_counts.values():
            if used<=n:
                assert opt>=cost+tables[key]['doubled_optima'][n-used]
                upper_tests+=1
        rs=row['witness_ranks'];assert sum(rs)==n and all(1<=r<=75 for r in rs)
        types=[representative(r)for r in rs]
        assert [''.join(map(str,t))for t in types]==row['witness_components']
        assert forest_cost(types)==opt
        count=tuple(rs.count(r)for r in range(10,18));assert count in all_counts
        used,key,cost=all_counts[count];rest=n-used
        assert opt==cost+tables[key]['doubled_optima'][rest]
        assert row['doubled_slack_below_8n']==16*n-opt
        beta2=None if n<3 else (n*(n-6)if n%2==0 else(n-1)*(n-5))//4
        assert row['doubled_curvature_slack']==(None if beta2 is None else beta2-opt)
        answers.append(dict(dimension=n,doubled_cost=opt,witness_components=row['witness_components'],
            doubled_curvature_slack=row['doubled_curvature_slack']))
    assert answers[69]['doubled_cost']==1096 and answers[69]['doubled_curvature_slack']==-8
    assert answers[70]['doubled_cost']==1106 and answers[70]['doubled_curvature_slack']==14
    return dict(type_domination=type_domination(),large_count_vectors=len(all_counts),
        counts_by_number_of_large_components=dict(Counter(sum(x)for x in all_counts)),
        bellman_tables=table_checks,dimension_maxima=answers,
        counts=dict(large_count_vectors=len(all_counts),residual_tables=len(tables),
            bellman_upper_bound_alternatives=transitions,attaining_predecessors=predecessors,
            global_upper_bound_alternatives=upper_tests,positive_rank_maxima=75),
        conclusion=dict(rank69_maximum='548',rank69_curvature='544',rank69_excess='4',
            rank70_maximum='553',rank70_curvature='560',rank70_slack='7'))

def main(compare=False,saved=False):
    primary=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert primary['status']=='PASS_EXACT_CURRENT_FINITE8_FOREST_MAXIMA'
    result=verify(primary)
    if compare:
        print('PASS_INDEPENDENT_BELLMAN_COMPARISON',result['counts'],result['conclusion'],flush=True);return
    paths=[Path(__file__)]
    for tag in(PRIMARY,*INPUTS):
        doc=json.loads((ROOT/f'results/{tag}.json').read_text());assert doc['status'].startswith('PASS')
        for p,h in doc['sources'].items():assert digest(ROOT/p)==h,(tag,p)
        paths.extend(ROOT/f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    result.update(status='PASS_INDEPENDENT_EXACT_FINITE8_BELLMAN_BUDGET_REVIEW',
        sources={str(p.relative_to(ROOT)):digest(p)for p in paths},
        scope='No producer imports. Independent count-vector domain, all Bellman upper bounds and attaining predecessors, all global bounds/attainment and coupon type domination. Inherited finite Coxeter classification and graph base costs are mathematically reviewed inputs; no physical realization or dimension theorem follows.')
    result=json.loads(json.dumps(result));path=ROOT/f'results/{TAG}.json'
    if saved:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],result['counts'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--compare-only',action='store_true');p.add_argument('--verify',action='store_true')
    a=p.parse_args();assert not(a.compare_only and a.verify);main(a.compare_only,a.verify)
