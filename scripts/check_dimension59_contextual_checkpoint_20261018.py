#!/usr/bin/env python3
"""Metadata/headline audit only; no mathematical or numerical search replay."""
import argparse,ast,json,re
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
TAG='dimension59_contextual_payments_and_moving64_20261018'
STEMS=(
 'fresh_depth6_dimension59_20261018',
 'fresh_rank13_contextual_coupon_rebalance_20261018',
 'fresh_order14_b_basic_distance9_transfer_20261018',
 'fresh_order14_c_basic_distance9_transfer_20261018',
 'fresh_order14_d_near_fork_half_transfer_20261018',
 'fresh_order14_contextual_ambient59_frontier_20261018',
 'fresh_raw56_raw57_extra_vertex_family_20261018',
 'fresh_joint64_native_compact_closure_20261018',
 'fresh_finite_rank14_existing_payment_transfer_20261018')
HELPERS=('check_raw56_raw57_extra_vertex_saved_20261018.py',
 'check_joint64_native_compact_saved_20261018.py')


def digest(p):
    h=sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def run():
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    checked={};bindings=0
    for stem,d in docs.items():
        for p,h in d['sources'].items():
            if p not in checked:checked[p]=digest(ROOT/p)
            assert checked[p]==h,(stem,p);bindings+=1
    old,root,b,c,d,new,ar,con,finite14=[docs[s] for s in STEMS]
    assert old['dimension_bound']==new['completed_unrestricted_dimension_bound']==59
    assert root['new_dimension_bound'] is new['new_dimension_bound'] is None
    assert b['new_cost']==c['new_cost']==d['new_cost']=='80'
    assert new['current_order14_family_costs']==dict(A=70,B=80,C=80,D=80)
    assert new['finite']['source_preflight_geometric_hypothesis_discharged']
    f=root['conditional_neighbor_preflight'];nf=new['finite']
    assert f['exact_maximum']==nf['exact_maximum']=='795/2'
    assert F(f['exact_maximum'])-F(f['curvature'])==6
    assert f['fixed_vectors']==len(f['alternatives'])==961 and f['interaction_keys']==960
    assert len(f['tables'])==1920 and len(f['above_curvature'])==14
    assert f['Bellman_upper_inequalities']==2613120 and f['attaining_predecessors']==113280
    assert all(len(t['values'])==len(t['last'])==60 for t in f['tables'])
    assert max(r['sixteenths_upper'] for r in f['alternatives'])==max(r['sixteenths_actual'] for r in f['alternatives'])==6360
    assert f['max_with_rank13']=='1563/4' and f['max_with_rank10']=='386'
    assert [w['new'] for w in f['small_obstructions_removed']]==['1561/4','1561/4']
    assert f['actual_D15_cost']=='201/2' and f['A1_D13_component_cost']=='2'
    assert len(root['combined_literal_interface']['rows'])==new['combined_literal_interface']['configurations']==216
    literal_key=lambda r:(r['source'],json.dumps(r['case'],sort_keys=True))
    predecessor={literal_key(r):r for r in finite14['combined_literal_interface']['rows']}
    assert len(predecessor)==216
    assert {literal_key(r) for r in root['combined_literal_interface']['rows']}==set(predecessor)
    assert Counter(len(r['corners']) for r in root['combined_literal_interface']['rows'])=={4:184,3:32}
    for row in root['combined_literal_interface']['rows']:
        before=predecessor[literal_key(row)]
        assert len(row['corners'])==len(before['corners'])
        assert [c['marks'] for c in row['corners']]==[c['marks'] for c in before['corners']]
    ideal=new['ideal'];assert ideal['complete_profiles']==len(ideal['rows'])==1890
    assert len({tuple(r['ranks']) for r in ideal['rows']})==1890
    assert sum(F(r['normalized_margin'])<0 for r in ideal['rows'])==ideal['summary']['negative_count']==7
    assert ideal['summary']['zero_count']==0 and ideal['old_negative_count']==10
    assert min(F(r['normalized_margin']) for r in ideal['rows'])==-F(5241,57460)
    assert ideal['summary']['normalized_minimizers']==[[9,9,12,12,16]]
    assert F(next(r['normalized_margin'] for r in ideal['rows'] if r['ranks']==[9,9,9,9,9,13]))>0
    assert not ideal['fresh_convolution_performed'] and ideal['negative_absolute_margin_not_extended_to_arbitrary_profiles']
    assert ar['counts']['family_classes']==len(ar['family'])==1408
    assert ar['counts']['marker_orientations']==21584
    assert ar['counts']['target_forests']==len(ar['forest_records'])==1212
    assert ar['counts']['exported_columns']==len(ar['exported']['columns'])==1212
    assert len({r['vertex'] for r in ar['exported']['columns']})==1212
    assert ar['counts']['anchor_tests']==4115 and ar['counts']['retained_cross_pairs']==13104
    assert ar['counts']['strictly_violated_columns']==60 and ar['optimizer_invocations']==0
    repair=ar['exported']['repair'];ratio=F(repair['objective_numerator'])/F(repair['normalizer'])
    assert ratio==F(repair['bound'])==F(14671331978990040,1220681360116301)>12
    assert len(ar['exported']['preserved_raw76_columns'])==242
    assert len(ar['exported']['preserved_multi_owner_costs'])==14
    assert con['summary']==dict(selected_prefixes=9,M_recipes=152,M_branches=304,M_native=63,R_banks=63,
        R_recipes=1323,R_branches=2646,native_five_caps=0,clean_five_caps=0,completed_disks=0,
        fewer_than_ten_ends=0,unresolved=0)
    caps=[m for r in con['records'] for m in r['bank']['caps']];assert len(caps)==63
    assert all(m['score'] is not None and m['finite_rays']==19 and m['open_ends']==13 for m in caps)
    assert con['exact_deletion_debit']['lost_finite_count']==19
    assert con['exact_deletion_debit']['compact_faces_by_dimension']=={'2':12,'3':1}
    assert con['replacement64_id']==71 and not con['accepted'] and not con['unresolved']
    reports=[ROOT/f'results/{s}.md' for s in STEMS]+[ROOT/f'results/{TAG}.md']
    scripts=[ROOT/f'scripts/{s}.py' for s in STEMS]+[ROOT/'scripts'/s for s in HELPERS]+[Path(__file__)]
    for p in scripts:ast.parse(p.read_text())
    links=0
    for report in reports:
        for link in re.findall(r'\]\(([^)]+)\)',report.read_text()):
            if '://' in link or link.startswith('#'):continue
            path=(report.parent/link.split('#')[0]).resolve()
            # The result being produced is referenced before its first write.
            if path==ROOT/f'results/{TAG}_audit.json':continue
            assert path.is_file(),(report,link);links+=1
    for name in ('README.md','RESEARCH_PLAN.md','GOAL.md'):
        opening=(ROOT/name).read_text()[:6000]
        assert f'results/{TAG}.md' in opening,name
    paths=reports+scripts+[ROOT/f'results/{s}.json' for s in STEMS]
    return dict(status='PASS_DIMENSION59_CONTEXTUAL_CHECKPOINT_METADATA',
        direct_bindings=bindings,unique_bound_files=len(checked),syntax_files=len(scripts),report_links=links,
        sources={str(p.relative_to(ROOT)):digest(p) for p in sorted(set(paths))},
        completed_bound=59,new_dimension_bound=None,finite_maximum='795/2',finite_curvature='783/2',
        finite_representative_obstructions=14,ideal_negative_profiles=7,order14_envelope='40/7',
        arithmetic_fixed_mass_ratio=str(ratio),new_finite_volume_example=False,
        scope='Metadata, local links and saved exact headline/cardinality audit only. No geometric census, finite Bellman proof, convolution, native field arithmetic, numerical cone or matrix optimization replay.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}_audit.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['direct_bindings'],out['unique_bound_files'],out['syntax_files'],out['report_links'])
