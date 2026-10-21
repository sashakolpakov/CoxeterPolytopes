#!/usr/bin/env python3
"""Exact complete disjoint-base cone check at every exported coface slice.

No producer imports or optimizer. Reconstruct demands from actual-provenance
coefficients and common face masses, then test the complete saved domain.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from math import lcm
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_sqrt5_complete_coface_base_cones_20261017'
DOMAIN='fresh_simplex_base_packing_cone_20261017'
SLICES='fresh_sqrt5_labelled_coface_base_packing_20261017'


def digest(p):
    return sha256((ROOT/p).read_bytes()).hexdigest()


def run():
    domain=json.loads((ROOT/f'results/{DOMAIN}.json').read_text())
    data=json.loads((ROOT/f'results/{SLICES}.json').read_text())
    assert domain['status']=='PASS_EXACT_COMPLETE_SIMPLEX_BASE_PACKING_CONES'
    assert data['status']=='PASS_EXACT_LABELLED_COFACE_SUBSET_PACKING_DIAGNOSTIC'
    for doc in (domain,data):
        for p,h in doc['sources'].items():
            assert digest(p)==h,p
    cones={c['vertices']:c for c in domain['cones']}
    alpha={i:F(v) for i,v in data['face_multipliers']}
    assert all(x>0 for x in alpha.values())
    pf={}
    bases=data['actual_bases']
    for b in bases:
        p,f=b['provenance'],b['face']
        assert p not in pf or pf[p]==f
        pf[p]=f
    assert len(pf)==1858 and len(bases)==6016
    rows=[]
    violations=[]
    count=Counter()
    positive_by_dimension=Counter()
    max_utilization=F()
    max_row=None
    total_demand=F()
    for row in data['coface_slice_demands']:
        n=len(row['vertex_roots'])
        assert n==row['simplex_dimension']+1 and 2<=n<=5
        c=cones[n]
        coeffs=c['doubled_inequality_coefficients']
        d={}
        for mask,terms in row['block_coefficients']:
            assert mask in c['base_masks']
            assert len(dict(terms))==len(terms)
            assert all(F(v)>0 for _,v in terms)
            value=sum((F(v)*alpha.get(pf[p],F()) for p,v in terms),F())
            if value:
                d[mask]=value
        assert d=={m:F(v) for m,v in row['block_marginals']}
        mass=sum(d.values(),F())
        assert mass==F(row['total_block_mass'])
        assert mass==sum((alpha.get(bases[i]['face'],F()) for i in row['source_bases']),F())
        cap=row['coface_multiplicity']*alpha.get(row['simplex_face'],F())
        assert cap==F(row['capacity']) and cap>=0
        assert alpha.get(row['simplex_face'],F())==F(row['simplex_mass'])
        total_demand+=mass
        count['slices']+=1
        count['complete_inequalities']+=len(coeffs)
        if not d:
            count['zero_mass_slices']+=1
            # Every coefficient is nonnegative, and cap>=0. All loads are zero.
            continue
        count['positive_mass_slices']+=1
        count['evaluated_positive_slice_inequalities']+=len(coeffs)
        positive_by_dimension[n-1]+=1
        den=lcm(cap.denominator,*(x.denominator for x in d.values()))
        scaled=[int(d.get(m,F())*den) for m in c['base_masks']]
        loads=[sum(a*b for a,b in zip(w,scaled)) for w in coeffs]
        maximum=max(loads)
        optimum=F(maximum,2*den)
        maximizer=loads.index(maximum)
        slack=cap-optimum
        if cap:
            utilization=optimum/cap
            if utilization>max_utilization:
                max_utilization=utilization
                max_row=row['index']
        else:
            assert optimum==0,('positive demands at zero capacity',row['index'])
        entry=dict(slice=row['index'],simplex_dimension=n-1,
            simplex_face=row['simplex_face'],coface_label=row['coface_label'],
            block_marginals=row['block_marginals'],capacity=str(cap),
            exact_minimum_packing_cover=str(optimum),maximizing_inequality=maximizer,
            slack=str(slack),inequality_count=len(loads))
        rows.append(entry)
        if slack<0:
            violations.append(entry)
    assert total_demand==sum((alpha.get(b['face'],F()) for b in bases),F())
    assert count['slices']==34520 and count['positive_mass_slices']==870
    assert count['evaluated_positive_slice_inequalities']==65069
    assert dict(positive_by_dimension)=={1:382,2:101,3:172,4:215}
    assert not violations and max_utilization==1
    paths=[f'scripts/{TAG}.py',f'scripts/{DOMAIN}.py',f'results/{DOMAIN}.json',
           f'scripts/{SLICES}.py',f'results/{SLICES}.json']
    return dict(status='PASS_EXACT_COMPLETE_COFACE_BASE_PACKING_CONES',
        sources={p:digest(p) for p in paths},counts=dict(count),
        positive_slices_by_dimension=sorted(positive_by_dimension.items()),
        slice_checks=rows,violations=violations,
        maximum_capacity_utilization=str(max_utilization),maximizing_slice=max_row,
        input_bound=data['input_bound'],output_bound=data['input_bound'],
        alpha_changed=False,mu_changed=False,optimizer_invocations=0,
        scope='Every one of the34520 exported simultaneous-coface slices passes the complete nonnegative disjoint-base packing cone on at most five simplex vertices. Demands are reconstructed from actual-provenance coefficients and unchanged common alpha. This exhausts abstract disjointness capacities on these averaged slices, not additional geometric compatibility or coupling between different cofaces. No new ambient dimension bound or global optimum is claimed.')


if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--verify',action='store_true')
    args=p.parse_args()
    result=json.loads(json.dumps(run()))
    out=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'],'violations',len(result['violations']),flush=True)
