#!/usr/bin/env python3
"""Bounded D15 fork-coupon provenance, not a changed payment theorem.

Read frozen literal data and output inventories. No native embedding search,
inertia, spectator classification or forest dynamic program is run.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path
import fresh_finite_d16_d17_long_pair_half_transfer_20261017 as T
import fresh_d12_type_split_ambient58_forest_frontier_20261018 as C

ROOT=C.ROOT
TAG='fresh_small_d15_fork_coupon_provenance_20261018'
COMPOSE='fresh_d12_half_face_composition_20261018'
SENS='fresh_ambient58_obstruction_payment_sensitivities_20261018'
LITERAL='fresh_finite_d12_long_leaf_support_scout_20261018'
ROUTES=[
 ('counting_enlargements',[
  'fresh_broad_finite_half_coupons_20261017','fresh_finite8_a17_d11_coupon_budget_20261017'],
  'Budget enlargement to all small marks against D fork triples; not itself a geometric output.'),
 ('complete_original_routing',['fresh_finite8_ideal5_combined_payments_20261017'],
  'Low/disjoint/high, protected/adaptive and ordinary selected outputs are same-component; coupon replacements have the separate universal or sealed proofs below.'),
 ('cycle_tree_short_paths',[
  'fresh_finite8_cycle_replacement_20261017','fresh_finite8_ideal5_local_replacement_20261017',
  'fresh_cycle12_ideal5_transfer_20261017','fresh_cycle13_ideal5_transfer_20261017','fresh_cycle14_ideal5_transfer_20261017'],
  'Whole A1/A2 against path A11/A13/A15/A17 roles; short paths and same-component alternatives. A15 is not D15.'),
 ('universal_fork_triples',[
  'fresh_affine_d13_long_pair_transfer_20261017','fresh_affine_d14_long_pair_transfer_20261017'],
  'A1/A2/A4/A5 against D11/D13, or small rank at most3 against D12/D14; no D15 cross output.'),
 ('original_order16_17',[
  'fresh_order16_finite_cross_payments_20261017','fresh_ideal7_finite_cross_payments_20261017'],
  'D15 fork coupons have whole A1 only; other original outputs are A2/D13, A1/A15, A1/D16, A2/D14 or same-component units.'),
 ('selected_ideal_partial_transfers',[
  'fresh_d16_partial_weight_transfer_20261017','fresh_order16_d_distance9_half_transfer_20261018',
  'fresh_order16_bd_sparse_unit_transfer_20261018','fresh_order17_d5_d12_unit_coupon_zero_transfer_20261018',
  'fresh_order17_e6_half_coupon_transfer_20261018','fresh_order17_e6_unit_coupon_zero_transfer_20261018'],
  'Whole A4/D13, A4/D12, A5/D11, D5/D11, D5/D12 or E6/D11; complete supporting-tail exclusions precede terminal corners.'),
 ('order11_12_spectators',[
  'fresh_cycle11_half_pair_transfer_20261018','fresh_order11_bcd_half_pair_transfer_20261018',
  'fresh_order11_bc_zero_transfer_20261018','fresh_order12_b_center_terminal_transfer_20261018',
  'fresh_order12_bc_near_terminal_half_transfer_20261018'],
  'Small A components against B9/B10/B11/D10/D11; other outputs are same-component or independent affine/short paths.'),
 ('order13_spectators',[
  'fresh_order13_d_center_pair_zero_transfer_20261018','fresh_order13_d_center_leaf_half_transfer_20261018',
  'fresh_order13_bc_sparse_terminal_transfer_20261018','fresh_order13_bc_inner_terminal_transfer_20261018'],
  'A2/A3-D10, A1/D12, A4/B9, A1/B12; whole D12/D13 and B11/B13 outputs are same-component.'),
 ('order14_15_universal',[
  'fresh_order14_d_near_fork_half_transfer_20261018','fresh_order15_bd_terminal_coupon_transfer_20261018',
  'fresh_order15_bd_basic_fractional_transfer_20261018'],
  'A1/D13 neighbor, A4/D10, A4/A5-D11, D5/D10, A1/D14 neighbor; whole D14/D15 alternatives are same-component.'),
 ('current_order16_sixteenth',['fresh_order16_cd_sparse_sixteenth_transfer_20261018'],
  'A4/A12 unit, A4/D12 short-leaf half, and whole A1/D15 node11 sixteenth. Node11 is not center12 or leaf13/14.'),
 ('contextual_rules',[
  'fresh_rank13_contextual_coupon_rebalance_20261018','fresh_order17_center_contextual_zero_transfer_20261018'],
  'A1/A13 and A2/A13 middle, A2/D13 center, contextual A2/A2 with whole A10, or A2/D5 with whole A1/A9. D15 is a same-component half.'),
 ('later_finite_outputs',[COMPOSE],
  'The22 fully scoped finite transfer packets are separately inventoried by their564 frozen terminal corners; no cross-component D15 occurs.'),
 ('additional_finite_interfaces',[
  'fresh_finite_a15_endpoint_coupon_zero_20261018','fresh_rank15_sparse_payment_exchange_20261018',
  'fresh_finite_a17_endpoint_payment_support_20261017','fresh_finite_d17_long_pair_support_20261017',
  'fresh_finite_d16_d17_long_pair_half_transfer_20261017','fresh_large_finite_d_leaf_half_transfer_20261017'],
  'A2/A13, A2/D13, A1/A2-D14 and small/D10, retained same-component A16/A17/D15/D16/D17, or independent short8; no rank3–5/D15 fork coupon.'),
 ('small_coupon_rebalances',[
  'fresh_small_d16_coupon_and_ideal17_rebalance_20261018','fresh_d5_d12_alternate_coupon_rebalance_20261018'],
  'Prior broad-versus-geometric inventory; retained A1/D16 and new D5/A12 or D5/A1, not D15.'),
 ('deep_short_path_priorities',[
  'fresh_depth6_finite6_zero_transfer_20261018','fresh_depth6_finite7_adaptive_transfer_20261018',
  'fresh_depth6_finite8_zero_transfer_20261018'],
  'Independent shorter paths, fixed leaves, sealed whole A8/A13 or A1/A12 middle unit and A1/A15 middle halves.'),
 ('order10_priorities',[
  'fresh_order10_c_endpoint_zero_transfer_20261018','fresh_order10_b_terminal_half_transfer_20261018',
  'fresh_order10_cycle_half_transfer_20261018','fresh_order10_d_opposite_fork_fractional_transfer_20261018'],
  'Small whole components against B9/A9/D9, restored same-component A17 endpoint, finite joined/short/whole alternatives. No small/D15 cross output.'),
 ('new_D12_rule',[
  'fresh_d12_radical_saturated_half_replacements_20261018','fresh_d12_remaining_selected_half_replacements_20261018',
  'fresh_d12_low_protected_adaptive_interfaces_20261018','fresh_d12_universal_ideal_output_repairs_20261018'],
  'New whole A1/A2-A11 endpoint halves or retained same-component outputs; no small/D15 fork coupon.'),
]


def digest(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {str(k):enc(v)for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [enc(v)for v in x]
    return x


def d15_role(g,nodes,mark):
    """Recognize this exact simple tree; no classification or inertia search."""
    if len(nodes)!=15:return None
    adj={i:[j for j in nodes if i!=j and g[i][j]]for i in nodes}
    if sum(map(len,adj.values()))!=28:return None
    centers=[i for i in nodes if len(adj[i])==3]
    if len(centers)!=1 or max(map(len,adj.values()))!=3:return None
    center=centers[0];leaves=[i for i in adj[center]if len(adj[i])==1]
    if len(leaves)!=2:return None
    assert all(g[i][i]==2 and all(g[i][j]==-1 for j in adj[i])for i in nodes)
    end=[i for i in nodes if len(adj[i])==1 and i not in leaves]
    assert len(end)==1
    dist={end[0]:0};queue=end[:]
    for i in queue:
        for j in adj[i]:
            if j not in dist:dist[j]=dist[i]+1;queue.append(j)
    assert len(dist)==15 and dist[center]==12 and all(dist[i]==13 for i in leaves)
    images=[13,14]if mark in leaves else[dist[mark]]
    return dict(center=center,short_leaves=sorted(leaves),canonical_mark_images=images,
        role='fork_center'if mark==center else'short_fork_leaf'if mark in leaves else'long_arm')


def literal_inventory(docs):
    count=corner_count=0;hits=[];keys=[]
    for stem,case,g,white,cs in T.literal_inputs(docs):
        count+=1;keys.append([stem,enc(case)]);corner_count+=len(cs)
        for marks,*_ in cs:
            parts=T.A.components(g,sorted(set(white)|set(marks)))
            aa,bb=[next(p for p in parts if m in p)for m in marks]
            if aa==bb:continue
            for small,large,mark in ((aa,bb,marks[1]),(bb,aa,marks[0])):
                role=d15_role(g,large,mark)
                if role is None:continue
                assert len(small)==1
                hits.append(dict(source=stem,case=case,marks=marks,small_type='A1',small_nodes=small,
                    D15_nodes=large,large_mark=mark,**role))
    assert count==216 and corner_count==832 and len(hits)==10
    hist=Counter(r['role']for r in hits)
    assert hist=={'fork_center':2,'short_fork_leaf':4,'long_arm':4}
    assert Counter(tuple(r['canonical_mark_images'])for r in hits)=={(12,):2,(13,14):4,(10,):2,(11,):2}
    return dict(records=count,corners=corner_count,keys=keys,D15_cross_corners=hits,
        role_histogram=dict(hist),small_rank3_to5_D15_corners=0,
        scope='Saved Gram connectivity and intrinsic D15 role recognition; inherited sealing/completeness, no fresh inertia or embeddings.')


def later_inventory(doc):
    rows=[];total=0
    for packet in doc['later_finite_output_interfaces']:
        hits=[]
        for corner in packet['corners']:
            total+=1;marks=corner['marks'];parts=corner['whole_components']
            aa,bb=[next(b for b in parts if m in b['nodes'])for m in marks]
            if aa==bb:continue
            for sm,la,m in ((aa,bb,marks[1]),(bb,aa,marks[0])):
                if la['type']=='D15':hits.append(dict(path=corner['path'],small=sm,large=la,large_mark=m))
        assert not hits,packet['source']
        rows.append(dict(source=packet['source'],proof_sections=packet['proof_sections'],
            corners=packet['checked_corners'],D15_cross_corners=hits))
    assert len(rows)==22 and total==564
    return dict(packets=rows,corner_count=total,D15_cross_corners=0,
        scope='Frozen verified whole-component terminal metadata; not a replacement for the bound complete case proofs.')


def costs():
    rows=[]
    for small in C.TYPES:
        r=small[1]
        if r not in (3,4,5):continue
        old=F(C.actual_cross(small,('D',15)),32);new=old-F(r,8)
        path=F(C.actual_cross(small,('A',15)),32);other=F(C.actual_cross(small,('B',15)),32)
        assert old==F(3*r,2) and path==F(r,2) and other==0 and new==F(11*r,8)
        assert new>=max(path,other)
        rows.append(dict(small_type=small,D15_current=old,D15_if_center_three_eighths=new,
            A15_unchanged=path,B15_unchanged=other,shared_rank15_allowance_can_dominate_these_three_types=True))
    assert len(rows)==11
    return dict(rows=rows,changed_rule_proved=False,forest_count_performed=False,
        center=dict(canonical_nodes=[12],old='1/2',candidate='3/8',saving_per_small_mark='1/8'),
        leaves=dict(canonical_nodes=[13,14],old='1/2',candidate='1/2',saving='0'),
        outside_scope=['all A1/D15 roles','rank2/D15 roles','A1/D15 node11 sixteenth','same-component D15','A15 path coupons'])


def run():
    stems=sorted(set(T.STEMS)|{C.TAG,COMPOSE,SENS,LITERAL}|{s for _,ss,_ in ROUTES for s in ss})
    docs={s:json.loads((ROOT/f'results/{s}.json').read_text())for s in stems}
    for s,d in docs.items():
        assert d['status'].startswith('PASS'),s
    composed=docs[COMPOSE]
    assert composed['baseline']=='fresh_order10_low8_ambient58_frontier_20261018'
    assert composed['finite_output_corners']==564 and docs[C.TAG]['finite']['certified_upper']=='1527/4'
    literal=literal_inventory(docs);later=later_inventory(composed)
    expected=docs[LITERAL]['literal_interfaces']
    assert literal['keys']==[[r['source'],r['case']]for r in expected]
    # Each later terminal source is itself a bound input, as are its full
    # arbitrary-support proof report and local certificate.
    stems=sorted(set(stems)|{r['source']for r in later['packets']})
    paths={f'scripts/{TAG}.py'}|{f'{folder}/{s}.{ext}'for s in stems for folder,ext in
        [('scripts','py'),('results','json'),('results','md')]}
    return enc(dict(status='PASS_BOUNDED_SMALL_D15_FORK_COUPON_PROVENANCE',
        sources={p:digest(p)for p in sorted(paths)},baseline=COMPOSE,
        universal_output_inventory=[dict(route=k,sources=ss,conclusion=desc)for k,ss,desc in ROUTES],
        original_geometric_small_D15_type='whole A1',broad_rank3_to5_added_by_counting_enlargement=True,
        literal_inventory=literal,later_finite_inventory=later,candidate_accounting=costs(),
        geometric_replacement_required_by_audited_rank3_to5_D15_outputs=[],
        no_new_payment_theorem=True,no_forest_optimization=True,new_dimension_bound=None,
        scope='Bounded proof-output provenance plus saved-corner and component-pair arithmetic. No universal validity inferred solely from corner absence; no changed rule, full forest certificate or new dimension theorem.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    path=ROOT/f'results/{TAG}.json';assert args.verify or not path.exists()
    out=run()
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['literal_inventory']['records'],out['later_finite_inventory']['corner_count'],flush=True)
