#!/usr/bin/env python3
"""Check independent algebra, inherited cases and the all-block domain split."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_t8_all_blocks_20261014_review'
PRIMARY = 'results/fresh_h9_t8_sparse_search_20261014_all_blocks.json'
ALGEBRA = 'results/fresh_h9_t8_all_blocks_20261014_review_algebra.json'
BASE = 'results/fresh_h9_t8_asymmetric_20261014.json'
SPARSE = 'results/fresh_h9_t8_sparse_frontier_20261014_review.json'
MIDDLE = 'results/fresh_h9_t8_middle_20261014_review.json'
NAMES = ('e','a','b','d','c','delta')

def digest(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def value(x): return tuple(map(F,x))+(F(0),)*(4-len(x))

def run():
    paths = [PRIMARY,ALGEBRA,BASE,SPARSE,MIDDLE]
    docs = {p:json.loads((ROOT/p).read_text()) for p in paths}
    bound = set()
    def bind(p):
        if p in bound: return
        bound.add(p)
        for q,expected in json.loads((ROOT/p).read_text()).get('sources',{}).items():
            assert digest(q) == expected, q
            if q.endswith('.json'): bind(q)
    for p in paths: bind(p)
    data, algebra, base, sparse, middle = [docs[p] for p in paths]
    assert data['status'] == 'PASS_T8_ALL_THREE_INTERNAL_LABELS_GE_4_SPARSE_FAMILY_EXCLUDED'
    assert algebra['status'] == 'PASS_INDEPENDENT_384_ALL_BLOCKS_QUARTIC_CASES'
    assert sparse['status'] == 'PASS_INDEPENDENT_SPARSE_T8_DEPENDENCY_REDUCTION'
    assert middle['status'] == 'PASS_INDEPENDENT_T8_UNRESTRICTED_MIDDLE_LABEL_EXTENSION'
    assert sorted(data['fixed_labels']) == [[1,14,2],[2,14,3],[4,14,3],[5,14,3]]
    assert sorted(data['finite_label_lower_bounds']) == [[1,2,7],[4,5,4],[7,8,4],[12,13,4]]
    original = {(b['id'],tuple(r['foreign_contacts']),tuple(r['third_contact_labels']))
                for b in base['blocks'] for r in b['cases']}
    assert len(original) == 408
    all_cases = {key+(second,) for key in original for second in (4,5)}
    changed = {key for key in all_cases if any(i in (4,5) for i in key[1])
               and (key[3] == 4 or any(i in (7,8) for i in key[1]))}
    inherited = all_cases-changed
    def primary_key(r): return r['block_id'],tuple(r['foreign_contacts']),tuple(r['third_contact_labels']),r['second_middle_label']
    assert {primary_key(r) for r in data['inherited_cases']} == inherited
    assert len(changed) == 384 and len(inherited) == 432
    separate = {(r['block'],tuple(r['foreign']),tuple(r['third_labels']),r['second_middle_label']):r
                for r in algebra['cases']}
    assert set(separate) == changed
    count = 0
    seen = set()
    for row in data['cases']:
        key = primary_key(row)
        assert key not in seen
        seen.add(key)
        own = separate[key]
        for label in ('foreign_costs','third_pairings'):
            assert list(map(value,row[label])) == list(map(value,own[label]))
        assert value(own['Q']) == value(data['changed_projection']['Q_at_first_t_lower'])
        variables = row['free_variables']
        indices = [NAMES.index(v) for v in variables]
        assert len(indices) == len(set(indices))
        for v,interval in zip(variables,row['variable_intervals']):
            assert list(map(value,interval)) == [value(own['lower_bounds'][v]),(F(1),F(0),F(0),F(0))]
        for name,saved in (('positive_residual_numerator',own['certificates'][0]),
                           ('negative_anchor_numerator',own['certificates'][1])):
            cert = row[name]
            assert all(saved['degree'][j] == 0 for j in range(6) if j not in indices)
            assert cert['degrees'] == [saved['degree'][j] for j in indices]
            for target,source in (('power_terms','unit_power_terms'),('positive_coefficients','positive_Bernstein_coefficients')):
                a = {tuple(r['index']):value(r['value']) for r in cert[target]}
                b = {tuple(r['index'][j] for j in indices):value(r['value']) for r in saved[source]}
                assert a == b and len(a) == len(cert[target])
                if target == 'positive_coefficients': count += len(a)
    assert seen == changed and count == 5472
    assert data['positive_coefficient_count'] == algebra['positive_coefficient_count'] == count
    minimum = F(261,2000),F(0),F(-31,1000),F(0)
    assert value(data['minimum_positive_coefficient']) == value(algebra['minimum_positive_coefficient']) == minimum
    assert minimum[0]**2 > 5*minimum[2]**2
    assert F(1,4)/(1-F(3,4)) == 1 and F(1,4)+1 > 1
    paths.append(f'scripts/{TAG}.py')
    return dict(status='PASS_INDEPENDENT_T8_ALL_INTERNAL_BLOCK_LABELS_EXTENSION',
        sources={p:digest(p) for p in paths}, recursively_bound_manifests=sorted(bound),
        changed_cases=384, inherited_rows=432, positive_coefficients=count,
        minimum_positive_coefficient=list(map(str,minimum)),
        domain_argument=[
            'All sparse guards use nonpositive symmetric Gram relaxations or nonnegative Schur-cost lower bounds. Increasing4-5 from4 preserves these, including the new affine stars.',
            'The first changed block p=4 is exactly the independently verified sparse theorem. At p=5,q=4, precisely256 cases with foreign4/5 contacts change.',
            'At p=q=5, precisely128 cases with one foreign contact in each block additionally change. A case without foreign4/5 uses the old q=5 bank with increasedQ; a case without foreign7/8 uses the p=5,q=4 bank literally.',
            'At either block label>=6, every foreign nonright contact there costs>=1 before the positive retained own cost, and is forbidden.',
            'Q=1/[4(1-t)]+1/[2(1-c_p)]-1 increases with the first high-block parameter t and the changed-block cosine c_p. With foreign4/5 absent, A,ac,bc do not depend on c_p. The inherited A<0 gives dN/dQ=-A>0.',
            'The exact coefficient comparisons rebuild both inverse blocks and third-wall projections independently. All nine own-support patterns, last-label ranges, and rank-eight residual contradiction retain their independently checked sparse domains.',
        ],
        scope='Only the literal T8 incidence with1-2>=7,4-5/7-8/12-13>=4, and fixed14-1=2,14-2=14-4=14-5=3. All other labels free within compact Coxeter requirements. No whole-incidence exclusion.')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify: assert result == json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['changed_cases'],result['positive_coefficients'])
