#!/usr/bin/env python3
"""Integrity and saved-headline audit, not a geometric/computational replay."""
import argparse,ast,json,re
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension59_complete_and_current_search_20261018'
STEMS=(
 'fresh_small_d16_coupon_and_ideal17_rebalance_20261018',
 'fresh_depth6_finite7_adaptive_transfer_20261018',
 'fresh_depth6_sparse_ambient60_forest_frontier_20261018',
 'fresh_depth6_dimension59_20261018',
 'fresh_raw76_extra_vertex_complete_family_20261018',
 'fresh_moving33_compact10_bridge_20261018',
 'fresh_compact10_three_chain_gate_20261018')

def run():
    hashes={}
    def digest(p):
        p=Path(p)
        if p not in hashes:
            h=sha256()
            with p.open('rb') as f:
                for b in iter(lambda:f.read(1<<20),b''):h.update(b)
            hashes[p]=h.hexdigest()
        return hashes[p]
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    bindings=0
    for stem,d in docs.items():
        assert d['status'].startswith(('PASS','BOUNDED')),stem
        for p,h in d['sources'].items():assert digest(ROOT/p)==h,(stem,p);bindings+=1
    local=docs[STEMS[1]];count=docs[STEMS[2]]
    assert local['triangle_minimum']==2 and local['quadrilateral_minimum']==1
    assert local['schur']['branch_count']==766 and local['schur']['positive_branches']==765
    assert local['endpoint_seal']['norm2_projection_lower']=='5/2'
    assert local['counts']['formulas']==count['conditional_all_rank_deleted_pair_formulas']
    finite=count['finite'];ideal=count['ideal']
    assert finite['exact_maximum']=='404' and finite['curvature']=='405' and finite['actual_attaining_cost']=='404'
    assert finite['complete_fixed_vectors']==1044 and finite['interaction_keys']==1042 and finite['above_curvature']==0
    assert ideal['complete_profiles']==2113 and ideal['summary']['absolute_minimum']=='54'
    assert all(F(r['normalized_margin'])>0 and F(r['absolute_margin'])>=54 for r in ideal['rows'])
    assert count['strict_global_constant']=='4500'
    proof=docs[STEMS[3]]
    assert proof['dimension_bound']==59 and proof['excluded_dimension']==60
    assert proof['conditional_count_hypothesis_discharged'] and proof['constant_counted_once']
    assert proof['finite_maximum']=='404' and proof['actual_ideal_absolute_surplus']=='54'
    assert proof['contradiction']=='0 >= 1*v_f + 54*v_i + 4500 > 0'
    raw=docs[STEMS[4]]
    assert raw['counts']['classes']==235 and len(raw['exported']['columns'])==242 and raw['optimizer_invocations']==0
    assert F(raw['exported']['repair']['bound'])==F(102699323852930280,8367588524424317)>12
    moving=docs[STEMS[5]];gate=docs[STEMS[6]]
    assert moving['summary']['second_branches']==560 and moving['summary']['native_pairs']==0
    assert gate['complete_low_signs']==92 and gate['high_families']==5 and gate['old44_control']
    sources=[ROOT/f'scripts/{s}.py' for s in STEMS]
    sources.extend(ROOT/f'scripts/{s}.py' for s in (
        'check_raw76_extra_vertex_complete_family_saved_20261018',
        'check_moving33_compact10_bridge_saved_20261018',
        'check_compact10_three_chain_gate_saved_20261018'))
    sources.append(Path(__file__))
    for p in sources:ast.parse(p.read_text());digest(p)
    reports=[ROOT/f'results/{s}.md' for s in STEMS]
    reports += [ROOT/f'results/{STEMS[2]}_domain_clarification.md',ROOT/f'results/{TAG}.md']
    checked_links=0
    def links(p,text):
        nonlocal checked_links
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if target.startswith(('http:','https:','#')):continue
            dest=(p.parent/target.split('#')[0]).resolve()
            if dest==ROOT/f'results/{TAG}_audit.json':continue
            assert dest.exists(),(p,target)
            checked_links+=1
    for p in reports:links(p,p.read_text());digest(p)
    for name,stop in (('README.md','The [preceding bound61 checkpoint]'),
                      ('RESEARCH_PLAN.md','The [preceding bound62 checkpoint]'),
                      ('GOAL.md','## Objective')):
        p=ROOT/name;opening=p.read_text().split(stop,1)[0]
        assert '59' in opening and '21' in opening
        links(p,opening);digest(p)
    return dict(status='PASS_DIMENSION59_CHECKPOINT_INTEGRITY_AND_HEADLINES',packets=len(STEMS),
        direct_bindings=bindings,files_hashed=len(hashes),source_AST_checks=len(sources),
        reports_and_openings=len(reports)+3,local_links=checked_links,
        sources={str(p.relative_to(ROOT)):h for p,h in sorted(hashes.items())},
        saved_arithmetic=dict(finite_maximum='404',curvature='405',finite_surplus='1',ideal_surplus='54',constant='4500'),
        scope='File integrity, syntax, local links and saved headline arithmetic only. Does not rerun geometry, catalogue generation, profile/pair/clique enumeration, optimization, dynamic programming, convolutions, global full-matrix replay or formal proof verification.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run();p=ROOT/f'results/{TAG}_audit.json'
    if args.verify:assert result==json.loads(p.read_text())
    else:assert not p.exists();p.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],{k:v for k,v in result.items() if k not in('sources','scope')})
