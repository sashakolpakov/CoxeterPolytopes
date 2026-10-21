#!/usr/bin/env python3
"""Independent Bellman certificates for conditional rank17 finite budgets.

Standard library only; no project/producer arithmetic imports. Geometric
coverage of the assumed simultaneous rank17 deletions remains a separate input.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PRIMARY='fresh_rank17_reduced_forest_frontier_20261017'
TAG=PRIMARY+'_review'
INPUTS=('fresh_a17_endpoint_d11_upgrade_forest_budget_20261017',
        'fresh_a17_endpoint_d11_upgrade_forest_budget_20261017_review',
        'fresh_a11_cycle_coupon_domination_20261017',
        'fresh_ideal5_weighted_dimension68_20261017')
LIMIT=75


def digest(path):return sha256(path.read_bytes()).hexdigest()
def parse(name):return name[0],int(name[1:])
def label(t):return t[0]+str(t[1])


def base(t,a17,d17):
    kind,r=t
    if r==17:return 2*(a17 if kind=='A'else d17)
    if r<=16:return r*(r-1)
    assert kind in('A','B','D')
    return 16*r-(64 if kind=='D'else 68)


def cross(a,b):
    """Doubled cross-coupon count, including the dominated A11 coupons."""
    out=0
    for small,large in((a,b),(b,a)):
        if small[1]<=5:
            if large[0]=='D'and 11<=large[1]<=16:out+=3*small[1]
            elif large in(('A',13),('A',15)):out+=small[1]
            elif large==('A',17):out+=2*small[1]
        if(small,large)in((('A',1),('D',11)),(('A',1),('D',12)),
                          (('A',2),('B',10)),(('A',4),('D',13)),
                          (('A',1),('A',11)),(('A',2),('A',11))):out+=2
    return out


def representative(t):
    kind,r=t
    if r<=9:return('A',r)
    if r==10:return('B',r)
    if r==17 and kind=='A':return t
    return('D',r)


def forest_cost(ts,a17,d17):
    return sum(base(t,a17,d17)for t in ts)+sum(cross(x,y)for i,x in enumerate(ts)for y in ts[i+1:])


def domination(a17,d17):
    types=([('A',r)for r in range(1,76)]+[('B',r)for r in range(2,76)]+
           [('D',r)for r in range(4,76)]+[('E',6),('E',7),('E',8),('F',4),('H',3),('H',4),('I2*',2)])
    rows=[]
    for t in types:
        rep=representative(t);assert rep in types
        inc=base(rep,a17,d17)-base(t,a17,d17);assert inc>=0
        deltas=[]
        for other in types:
            delta=cross(rep,other)-cross(t,other);assert delta>=0,(t,other,rep)
            if delta:deltas.append([label(other),delta])
        rows.append(dict(original=label(t),representative=label(rep),
            doubled_base_increase=inc,positive_cross_increases=deltas))
    return dict(types=len(types),cross_comparisons=len(types)**2,rows=rows,
        interpretation='Rank17 A and D are both retained; the asserted B17 upper cost is D17. All other same-component costs/classification are inherited. A11 unit/half coupons are explicitly dominated by D11. I2* denotes every I2(m),m>=4.')


def count_domain(a17,d17):
    rows={}
    # This enumerates component-number/rank multisets, not the producer's
    # consecutive-rank count recursion.
    for number in range(8):
        for rs in combinations_with_replacement(range(10,18),number):
            if sum(rs)>LIMIT:continue
            counts=tuple(rs.count(r)for r in range(10,18))
            assert counts not in rows
            types=[('A',r)if r==17 else('B',r)if r==10 else('D',r)for r in rs]
            key=(sum(11<=r<=16 for r in rs),rs.count(17),rs.count(10),rs.count(12),rs.count(13),rs.count(11))
            rows[counts]=(sum(rs),key,sum(base(t,a17,d17)for t in types))
    assert len(rows)==1608
    return rows


def verify_variant(p):
    a17,d17=p['assumed_A17_base'],p['assumed_D17_base']
    assert(a17,d17)in((131,131),(114,127))
    assert p['assumed_B17_base_upper_bound']==d17
    residual=tuple(range(1,10))+tuple(range(17,76))
    assert p['residual_ranks']==list(residual)
    domain=count_domain(a17,d17)
    stored={tuple(r['component_counts']):r for r in p['large_count_vectors']}
    assert len(stored)==len(p['large_count_vectors'])==1608 and set(stored)==set(domain)
    for count,(rank,key,cost)in domain.items():
        r=stored[count]
        assert(r['used_rank'],tuple(r['coupon_counts']),r['doubled_base_cost'])==(rank,key,cost)
    keys={key for rank,key,cost in domain.values()}
    tables={tuple(r['coupon_counts']):r for r in p['residual_tables']}
    assert len(tables)==len(p['residual_tables'])==len(keys)==675 and set(tables)==keys
    checked=[];alternatives=0
    for key in sorted(keys):
        t,u,v,b,c,w=key;assert t>=b+c+w
        # A dummy D14 has exactly the common broad small-component coupon
        # profile of the remaining D14/15/16 components. This is an interaction
        # multiset, not a claimed witness with the original total rank.
        large=([('D',14)]*(t-b-c-w)+[('D',12)]*b+[('D',13)]*c+
               [('D',11)]*w+[('A',17)]*u+[('B',10)]*v)
        types={r:('A',r)if r<=9 else('D',r)for r in residual}
        costs={r:base(types[r],a17,d17)+sum(cross(types[r],x)for x in large)for r in residual}
        q=tables[key];values=q['doubled_optima'];last=q['last_rank']
        assert len(values)==len(last)==76 and values[0]==last[0]==0
        count=0
        for n in range(1,76):
            assert isinstance(values[n],int)
            for r in residual:
                if r<=n:
                    assert values[n]>=values[n-r]+costs[r],(key,n,r)
                    count+=1
            r=last[n];assert r in residual and r<=n
            assert values[n]==values[n-r]+costs[r],(key,n,'attaining predecessor')
        alternatives+=count
        checked.append(dict(coupon_counts=list(key),doubled_optima=values,last_rank=last,
            all_alternative_bounds=count,attaining_predecessors=75))
    maxima=p['dimension_maxima'];assert len(maxima)==76
    upper=0;answers=[]
    for n,q in enumerate(maxima):
        assert q['dimension']==n
        value=q['doubled_cost'];assert isinstance(value,int)and value>=0
        for rank,key,cost in domain.values():
            if rank<=n:
                assert value>=cost+tables[key]['doubled_optima'][n-rank],(n,key)
                upper+=1
        ts=list(map(parse,q['witness_components']));assert sum(r for kind,r in ts)==n
        assert all(t==representative(t)for t in ts)
        assert forest_cost(ts,a17,d17)==value
        large=[r for kind,r in ts if 10<=r<=16 or(kind,r)==('A',17)]
        count=tuple(large.count(r)for r in range(10,18));rank,key,cost=domain[count]
        assert value==cost+tables[key]['doubled_optima'][n-rank]
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        assert F(q['finite_curvature_coefficient'])==beta and F(q['finite_slack'])==beta-F(value,2)
        answers.append(q)
    return dict(assumed_A17_base=a17,assumed_B17_base_upper_bound=d17,assumed_D17_base=d17,
        type_domination=domination(a17,d17),large_count_vectors=len(domain),
        counts_by_number_of_large_components=dict(Counter(sum(c)for c in domain)),
        bellman_tables=checked,dimension_maxima=answers,
        counts=dict(large_count_vectors=1608,residual_tables=675,
            bellman_upper_bound_alternatives=alternatives,attaining_predecessors=675*75,
            global_upper_bound_alternatives=upper,positive_rank_maxima=75))


def main(compare=False,verify=False):
    p=json.loads((ROOT/f'results/{PRIMARY}.json').read_text())
    assert p['status']=='PASS_CONDITIONAL_RANK17_FINITE_BUDGET_FRONTIER'and p['new_dimension_bound']is None
    rows=[verify_variant(x)for x in p['variants']]
    assert[(r['assumed_A17_base'],r['assumed_D17_base'])for r in rows]==[(131,131),(114,127)]
    assert[rows[1]['dimension_maxima'][n]['doubled_cost']for n in(67,68)]==[1008,1024]
    assert[rows[1]['dimension_maxima'][n]['witness_components']for n in(67,68)]==[['D67'],['D68']]
    assert[rows[0]['dimension_maxima'][n]['finite_slack']for n in(67,68)]==['-3/2','3']
    if compare:
        print('PASS_INDEPENDENT_CONDITIONAL_RANK17_BELLMAN',[(r['counts'],r['dimension_maxima'][67],r['dimension_maxima'][68])for r in rows],flush=True)
        return
    paths=[Path(__file__)]
    for tag in(PRIMARY,*INPUTS):
        doc=json.loads((ROOT/f'results/{tag}.json').read_text());assert doc['status'].startswith('PASS')
        for name,h in doc['sources'].items():assert digest(ROOT/name)==h,(tag,name)
        paths.extend(ROOT/f'{folder}/{tag}.{ext}'for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    result=dict(status='PASS_INDEPENDENT_CONDITIONAL_RANK17_BELLMAN_REVIEW',
        sources={str(path.relative_to(ROOT)):digest(path)for path in paths},variants=rows,
        new_dimension_bound=None,
        scope='No producer arithmetic imports. Complete independent large-count domain; each Bellman upper bound and attaining predecessor; every global upper bound and attaining typed forest; exact coupon type domination. Conditional on the proposed rank17 base costs and inherited finite graph/classification inputs. Does not independently reconstruct primary ideal scalar comparisons, prove geometric face coverage, or assert a dimension theorem.')
    result=json.loads(json.dumps(result));path=ROOT/f'results/{TAG}.json'
    if verify:assert result==json.loads(path.read_text())
    else:path.write_text(json.dumps(result,sort_keys=True,separators=(',',':'))+'\n')
    print(result['status'],[r['counts']for r in rows],flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--compare-only',action='store_true');ap.add_argument('--verify',action='store_true')
    args=ap.parse_args();assert not(args.compare_only and args.verify);main(args.compare_only,args.verify)
