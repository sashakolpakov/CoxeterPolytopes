#!/usr/bin/env python3
"""Adding the cycle12 A11 coupons does not increase the revised finite maximum."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_a17_endpoint_d11_upgrade_forest_budget_20261017_review as R

ROOT=R.ROOT;TAG='fresh_a11_cycle_coupon_domination_20261017'
CYCLE='fresh_cycle12_ideal5_transfer_20261017'
D11='fresh_affine_d11_leaf_pair_transfer_20261017'

def added_twice(a,b):
    if b=='A11':a,b=b,a
    return 2 if a=='A11'and b in('A1','A2')else 0

def run():
    stems=(R.BASE,R.TAG,CYCLE,D11)
    data={stem:json.loads((ROOT/f'results/{stem}.json').read_text())for stem in stems}
    for stem,d in data.items():
        assert d['status'].startswith('PASS'),stem
        for name,h in d['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,(stem,name)
    assert data[CYCLE]['finite_coupons']==[
        dict(small='whole_A1',large='whole_A11',large_position=5,weight='1'),
        dict(small='whole_A2',large='whole_A11',large_position=5,weight='1/2')]
    assert data[D11]['old_cost']==66 and data[D11]['new_cost']==62
    types=['A'+str(n)for n in range(1,76)]+['B'+str(n)for n in range(2,76)]+['D'+str(n)for n in range(4,76)]
    # I2 denotes every non-A2 member of the arbitrary dihedral family.
    types+=['E6','E7','E8','F4','H3','H4','I2']
    assert len(types)==228
    assert R.internal('A11')==R.internal('D11')==110
    rows=[]
    for partner in types:
        old=R.cross('A11',partner,True);new=old+added_twice('A11',partner)
        dominant=R.cross('D11',partner,True)
        assert old==0 and dominant>=new
        rows.append(dict(partner=partner,old_A11_cross_twice=old,new_A11_cross_twice=new,
            D11_cross_twice=dominant,domination_slack_twice=dominant-new))
    primary=data[R.BASE];variant=primary['variants'][1];assert variant['D11_upgrade']
    maxima=[]
    for row in variant['dimension_maxima']:
        witness=row['witness_components']
        assert 'A11'not in witness
        cost=R.direct_forest(witness,True)
        cost+=sum(added_twice(a,b)for i,a in enumerate(witness)for b in witness[i+1:])
        assert cost==row['doubled_cost']
        maxima.append(dict(row))
    assert maxima[69]['doubled_cost']==1088 and maxima[69]['doubled_curvature_slack']==0
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_A11_CYCLE_COUPONS_WITH_UNCHANGED_EXACT_FINITE_MAXIMA',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        added_coupons=data[CYCLE]['finite_coupons'],partner_domination=rows,
        unchanged_same_component_cost='55',dominating_type='D11',dimension_maxima=maxima,
        counts=dict(partner_types=228,attaining_forests=76),
        finite_maximum_at_n69='544',finite_curvature_slack_at_n69='0',
        theorem='Replacing each whole finite A11 component by D11 preserves rank/base and weakly increases every total cross-coupon cost. The previous exact maxima are upper bounds for the enlarged rule. Their saved witnesses contain no A11 and still attain each maximum.',
        scope='Exact finite budget for A17 endpoint deletion, full A1/D11-short-leaf upgrades and both new cycle12 A11 coupons, including all inherited coupons. No ideal profile or ambient dimension assertion.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));p=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(p.read_text())
    else:p.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['counts'],flush=True)
