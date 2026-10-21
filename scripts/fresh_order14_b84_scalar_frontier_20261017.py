#!/usr/bin/env python3
"""Complete old/new concentration comparisons after the order14 B reduction.

Scalar diagnostic only. Existing negative ideal coefficients remain; no
ambient theorem or geometric realization of a profile is asserted.
"""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_order14_b84_scalar_frontier_20261017'
LOCAL='fresh_order14_affine_b_fork_projection_transfer_20261017'
BASE='fresh_ideal5_weighted_dimension66_20261017'


def budget(r,new):
    special={11:F(65,12),12:F(75,13),13:F(6) if new else F(43,7),
             14:F(32,5),15:F(107,16),16:F(120,17)}
    return F(r,2) if r<=10 else special[r] if r<=16 else F(163,23)


def domain(total,restricted=None):
    def visit(left,minimum,prefix,low,high,middle):
        if left==0:
            yield prefix;return
        for r in range(minimum,left+1):
            if (r<10 and low) or (r>17 and high) or (r==restricted and middle):
                continue
            yield from visit(left-r,r,prefix+(r,),low or r<10,high or r>17,
                             middle or r==restricted)
    return list(visit(total,1,(),False,False,False))


def curvature(ranks):
    h=[1]
    for r in ranks:
        out=[0]*(len(h)+r)
        for j,a in enumerate(h):
            for k in range(r+1):out[j+k]+=a
        h=out
    d=sum(ranks);v=prod(r+1 for r in ranks)
    assert len(h)==d+1 and h==h[::-1] and sum(h)==v
    g=sum(a*min(j,d-j) for j,a in enumerate(h))
    return F(5*g-d*v,2*v)


def comparison(ranks):
    k=curvature(ranks);old=sum((budget(r,False) for r in ranks),F());new=sum((budget(r,True) for r in ranks),F())
    assert old-new==F(ranks.count(13),7)
    return dict(ranks=list(ranks),curvature=str(k),old_cost=str(old),new_cost=str(new),
                old_margin=str(k-old),new_margin=str(k-new))


def summary(rows,key):
    m=min(rows,key=lambda r:F(r[key]))
    return dict(profiles=len(rows),minimum_margin=m[key],minimizing_ranks=m['ranks'],
                negative_profiles=sum(F(r[key])<0 for r in rows),zero_profiles=sum(F(r[key])==0 for r in rows))


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in (LOCAL,BASE)}
    for s,d in docs.items():
        assert d['status'].startswith('PASS')
        for p,h in d['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h,(s,p)
    assert docs[LOCAL]['new_order14_costs']==dict(A=70,B=84,C=84,D=82)
    gain_old12=budget(11,False)+budget(13,False)-2*budget(12,False)
    gain_new12=budget(11,True)+budget(13,True)-2*budget(12,True)
    gain_new13=budget(12,True)+budget(14,True)-2*budget(13,True)
    assert (gain_old12,gain_new12,gain_new13)==(F(23,1092),-F(19,156),F(11,65))
    old_domain=domain(65,12);new_domain=domain(65,13);common=domain(65)
    assert len(common)==len(set(common))==1796
    assert set(old_domain)<=set(common) and set(new_domain)<=set(common)
    rows=[comparison(p) for p in common];by={tuple(r['ranks']):r for r in rows}
    old_rows=[by[p] for p in old_domain];new_rows=[by[p] for p in new_domain]
    comparisons=dict(old_reduced_domain=summary(old_rows,'old_margin'),
        new_reduced_domain=summary(new_rows,'new_margin'),
        common_expanded_domain_old=summary(rows,'old_margin'),
        common_expanded_domain_new=summary(rows,'new_margin'))
    assert comparisons['old_reduced_domain']['minimum_margin']=='-1565/10976'
    assert comparisons['new_reduced_domain']['minimum_margin']=='-31895/248832'
    assert comparisons['common_expanded_domain_new']['negative_profiles']==10
    assert comparisons['common_expanded_domain_new']['minimum_margin']=='-31895/248832'
    negatives=[]
    for r in rows:
        if F(r['new_margin'])>=0:continue
        old_image=r['ranks'][:]
        while old_image.count(12)>=2:
            old_image.remove(12);old_image.remove(12);old_image.extend((11,13))
        old_image.sort();image=by[tuple(old_image)]
        assert tuple(old_image) in set(old_domain)
        assert F(image['old_margin'])<=F(r['old_margin'])<0
        negatives.append(dict(**r,already_in_old_terminal_domain=tuple(r['ranks']) in set(old_domain),
            old_concentration_image=old_image,old_image_margin=image['old_margin']))
    negatives.sort(key=lambda r:F(r['new_margin']))
    old_obstruction=by[(13,)*5]
    assert old_obstruction['new_margin']=='6275/10976'
    assert F(1565,10976)*F(14,5)==F(313,784)
    assert F(23,1092)*14==F(23,78)
    next_obstruction=by[(10,11,11,11,11,11)]
    assert next_obstruction['old_margin']==next_obstruction['new_margin']=='-31895/248832'
    assert F(31895,248832)*F(12,5)==F(6379,20736)
    # A separately unproved order12 max65->64 improvement would change only
    # c(11), leaving the corrected concentration moves valid. This measures
    # the next geometric obligation; it does not supply that obligation.
    assert all(r['ranks'].count(11)>0 for r in negatives)
    required=max((-F(r['new_margin'])*12/r['ranks'].count(11),r['ranks']) for r in negatives)
    assert required==(F(27960,41327),[10,11,12,16,16])
    conditional=[(F(r['new_margin'])+F(r['ranks'].count(11),12),r['ranks']) for r in new_rows]
    assert all(m>0 for m,p in conditional)
    conditional_minimum=min(conditional)
    paths=[Path(__file__)]+[ROOT/f'{folder}/{s}.{ext}' for s in docs
        for folder,ext in (('scripts','py'),('results','json'),('results','md'))]
    paths.append(ROOT/'results/fresh_distance14_unrestricted_dimension156_20261017.md')
    return dict(status='PASS_ORDER14_B84_COMPLETE_SCALAR_FRONTIER_DIAGNOSTIC',ambient_dimension=66,
        new_dimension_bound=None,comparisons=comparisons,common_domain_rows=rows,
        old_terminal_ranks=[list(p) for p in old_domain],new_terminal_ranks=[list(p) for p in new_domain],
        remaining_negative_profiles=negatives,
        preexisting_old_terminal_negatives=sum(r['already_in_old_terminal_domain'] for r in negatives),
        newly_admitted_repeated12_negatives=sum(not r['already_in_old_terminal_domain'] for r in negatives),
        previous_worst_profile=old_obstruction,next_worst_profile=next_obstruction,
        concentration=dict(old_duplicate12_gain=str(gain_old12),new_duplicate12_gain=str(gain_new12),
            new_duplicate13_gain=str(gain_new13),common_repeated_ranks=list(range(10,18)),
            old_at_most_one_middle=12,new_at_most_one_middle=13,
            low_cap=10,at_most_one_low=list(range(1,10)),high_constant_starts=17,at_most_one_high_minimum=18),
        single_profile_cost_thresholds=dict(order14_drop_strictly_above='313/784',
            duplicate12_rule_breaks_above_order14_drop='23/78',
            order12_drop_for_next_profile_strictly_above='6379/20736',
            warning='These three thresholds concern the named individual profiles or move only; none proves all new margins positive.'),
        conditional_next_order12=dict(geometric_proof_supplied=False,
            current_order12_costs=dict(A=60,B=65,C=65,D=62),
            required_maximum_cost_drop_strictly_above=str(required[0]),controlling_profile=required[1],
            hypothetical_new_order12_maximum=64,complete_profiles=len(conditional),
            all_positive=True,minimum_margin=str(conditional_minimum[0]),minimizing_ranks=conditional_minimum[1],
            meaning='Both B and C costs must be improved; this is a conditional ideal-only comparison, with no new local transfer or finite bound.'),
        finite_obstruction_in_completed66='D66 still costs496 against finite curvature495 before any separately proved finite refinement.',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        scope='Complete exact comparison in the old, corrected new and common expanded concentration domains. Cost decrease exposes existing secondary obstructions; it cannot create a negative actual profile. No ambient theorem or cusp realization is asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['comparisons'],
          'old/newly admitted remaining negatives',out['preexisting_old_terminal_negatives'],out['newly_admitted_repeated12_negatives'],flush=True)
