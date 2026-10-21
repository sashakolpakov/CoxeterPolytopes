#!/usr/bin/env python3
"""Conditional scalar frontier for an order12 affine cost cap65; no new bound."""
import argparse,json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_ideal5_weighted_dimension69_20261017_review as R

ROOT=R.ROOT;TAG='fresh_rank11_cost65_scalar_frontier_20261017'
FINITE='fresh_a17_endpoint_d11_upgrade_forest_budget_20261017'
def budget(r):return F(65,12)if r==11 else R.budget(r)

def run():
    finite=json.loads((ROOT/f'results/{FINITE}.json').read_text())
    assert finite['status']=='PASS_EXACT_A17_ENDPOINT_AND_CONDITIONAL_D11_FOREST_BUDGETS'
    for name,h in finite['sources'].items():assert sha256((ROOT/name).read_bytes()).hexdigest()==h,name
    gain=budget(11)+budget(13)-2*budget(12);assert gain==F(23,1092)>0
    groups=[]
    for d in(67,68):
        shapes=set()
        for high in(0,*range(18,d+1)):
            for low in range(10):
                for twelve in(0,1):
                    left=d-high-low-12*twelve
                    if left<0:continue
                    exceptional=([high]if high else[])+([low]if low else[])+[12]*twelve
                    for rest in R.anchors(left,(10,11,13,14,15,16,17)):
                        shapes.add(tuple(sorted(exceptional+rest)))
        rows=[]
        for ranks in sorted(shapes):
            row=R.scalar(ranks);cost=sum((budget(r)for r in ranks),F())
            assert F(row['normalized_cost'])-cost==F(ranks.count(11),12)
            row.update(normalized_cost=str(cost),margin=str(F(row['normalized_curvature'])-cost))
            assert F(row['margin'])>0
            rows.append(row)
        n=d+1;minimum=min(rows,key=lambda r:F(r['margin']))
        f=finite['variants'][1]['dimension_maxima'][n]
        beta=F(n*(n-6),8)if n%2==0 else F((n-1)*(n-5),8)
        groups.append(dict(ambient_dimension=n,profiles=rows,profile_count=len(rows),
            minimum_margin=minimum['margin'],minimizing_ranks=minimum['ranks'],
            current_revised_finite_maximum=str(F(f['doubled_cost'],2)),curvature_coefficient=str(beta),
            finite_slack=str(beta-F(f['doubled_cost'],2)),strict_global_constant=5*(n-1)))
    assert [(g['profile_count'],g['minimum_margin'],g['finite_slack'])for g in groups]==[
        (1906,'411499/3234816','-17'),(2060,'3641/13608','0')]
    stems=(R.PRIMARY,R.TAG,FINITE)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems
        for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_CONDITIONAL_RANK11_COST65_SCALAR_FRONTIER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        assumed_affine_order12_maximum=65,rank11_normalized_cost='65/12',groups=groups,
        concentration=dict(repeated_ranks=[10,11,13,14,15,16,17],optional_low_ranks=list(range(1,10)),
            optional_rank12=True,optional_high_rank_minimum=18,duplicate12_cost_increase=str(gain),
            low_rank_cap=10,termination='Every nontrivial concentration increases the bounded integer squared-rank sum.'),
        new_dimension_bound=None,
        missing_obligations=['Complete A/B/C order12 transfers with triangle and quadrilateral composition.',
            'Any new finite coupons from those transfers must be included in a complete finite budget.'],
        scope='Conditional exact ideal comparisons with a newly justified terminal domain retaining repeatedrank10. Existing changed finite rule fits ambient69 but not68. This is neither a promoted ideal envelope nor a new ambient dimension theorem.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],[(g['ambient_dimension'],g['profile_count'],g['minimum_margin'],g['finite_slack'])for g in out['groups']],flush=True)
