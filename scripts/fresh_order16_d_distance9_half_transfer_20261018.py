#!/usr/bin/env python3
"""Two affine-D15 distance-nine units become halves with an existing coupon.

Complete bounded target-map reconstruction, arbitrary-contact exclusions,
and whole D5/D11 alternatives. No new finite cost or dimension assertion.
"""
import argparse,json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_d16_partial_weight_transfer_20261017 as D
import fresh_finite_rank15_rank16_rank62_budget_20261018 as B
import fresh_finite_rank15_existing_half_transfer_20261018 as R

ROOT=D.ROOT;P=D.P;Q=D.Q;L=D.L
TAG='fresh_order16_d_distance9_half_transfer_20261018'
PAIRS={(0,9),(2,11)}
EXPECTED=[359,363,496,499,548,966,981,1016,1020,1106,1141,1149,
          1221,1232,1267,1271,1286,1295,1315,1625,1655,2021,2030,2530]
STEMS=tuple(dict.fromkeys((D.TAG,D.SMALL,D.OVER,
    'fresh_distance8_adaptive_leaf_overlap_20261017',
    'fresh_distance8_low_core_20261017','fresh_distance8_high_disjoint_20261017',
    'fresh_finite8_ideal5_local_replacement_20261017',
    'fresh_finite8_cycle_replacement_20261017',
    'fresh_finite8_a17_d11_coupon_budget_20261017',
    'fresh_order16_affine_b_singleton_projection_transfer_20261018',
    'fresh_order14_ambient62_ideal_margin_closure_20261018',
    'fresh_order17_e6_unit_coupon_zero_transfer_20261018',
    'fresh_rank15_rank16_dimension61_20261018',
    'fresh_finite_rank15_rank16_rank62_budget_20261018',
    'fresh_depth11_finite8_zero_transfer_20261018',
    'fresh_finite7_ideal5_triangle_payments_20261017',*R.STEMS,R.TAG)))


def encode(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:encode(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [encode(v) for v in x]
    return x


def special_entry(c,core,constraints,target,mp):
    """The two relevant-context maps already force a forbidden first entry."""
    image={mp[j]:i for j,i in enumerate(c['marked_component'])}
    missing=sorted(set(range(16))-set(image));forced=[]
    J=[[F(x) for x in row] for row in core['inverse']]
    for t in missing:
        for j,i in image.items():
            if not target[t][j]:continue
            assert target[t][j]==3 and i in constraints['allowed_contacts']
            if J[i][i]>2:
                forced.append(dict(new_target_node=t,old_target_node=j,old_literal_node=i,
                    native_label=3,physical_contact_square=1,projection_lower_bound=str(J[i][i])))
    assert forced and c['index'] in (359,966)
    return dict(outcome='excluded_forced_entry_projection',forced_entries=forced,
        witness=forced[0],reason='Any required neighbor of the old marked component is outside the entire old core; one-end positivity bounds every arbitrary extra contact. The opposite end gives a retained short8 route.')


def sealed_coupon(c,row,g):
    assert row['outcome']=='sealed_completion'
    assert c['index'] in (1315,2530)
    assert row['signature']==[17,1,0] and len(g)==18
    assert not row['new_roots'] and row['sealing'].startswith('No core-contact root')
    candidates=[]
    for corner in row['corners']:
        if corner['kind']!='finite':continue
        marks=corner['marks'];blocks=corner['components']
        selected=[blocks[i] for i in corner['marked_component_indices']]
        for si,li in ((0,1),(1,0)):
            small,large=selected[si],selected[li]
            if small['type']=='D5' and large['type']=='D11' and large['roles'].get(str(marks[li]))=='fork_center':
                assert small is not large and corner['finite_half']=='1/2'
                candidates.append(dict(marks=marks,small_nodes=small['nodes'],small_type='D5',
                    small_mark=marks[si],small_role='any_node',large_nodes=large['nodes'],
                    large_type='D11',large_mark=marks[li],large_role='fork_center',weight='1/2'))
    assert len(candidates)==1 and candidates[0]['marks']!=c['preferred_corner']
    assert row['finite_payment_total']=='1/2' and row['total_payment']=='1'
    return dict(existing_whole_component_coupon=candidates[0],new_ideal_half='1/2',
        distinct_corner_total='1',no_same_component_payment_used=True)


def run():
    data={s:json.loads((ROOT/f'results/{s}.json').read_text()) for s in STEMS}
    for stem,doc in data.items():
        assert doc['status'].startswith('PASS'),stem
        for path,h in doc['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==h,(stem,path)
    assert data['fresh_rank15_rank16_dimension61_20261018']['dimension_bound']==61
    old=data[D.TAG];assert old['new_cost']=='106'
    rec=next(r for r in data[D.SMALL]['records'] if r['order']==16 and r['kind']=='D')
    target=rec['target'];assert target==old['target']
    autos=P.embeddings(target,target);adj,ds,features=P.distance_data(target)
    weights={tuple(r['pair']):F(r['weight']) for r in old['weighted_pairs']}
    assert all(weights[p]==1 and ds[p[0]][p[1]]==9 for p in PAIRS)
    assert all(sorted([len(adj[i]) for i in p])==[2,3] for p in PAIRS)
    assert all({tuple(sorted(a[i] for i in p)) for p in PAIRS}==PAIRS for a in autos)
    for p in PAIRS:weights[p]=F(1,2)
    assert sum(weights.values())==105
    short8={p for p in weights if ds[p[0]][p[1]]<=8}
    triangles={p for p in weights if ds[p[0]][p[1]]<=5 or
        (ds[p[0]][p[1]]<=8 and any(ds[i][f]<=5 for i in p for f in features))}
    assert all(weights[p]==1 for p in short8|triangles)
    assert all(weights[tuple(p)]==F(1,2) for p in old['half_weight_pairs'])
    assert all(weights[tuple(p)]==0 for p in old['inherited_zero_pairs'])
    over=data[D.OVER]['isolation'];routing=Counter();protected=[]
    for c in over['cases']:
        if c['method']=='finite_leaf9_or_ideal_feature5':
            assert c['protected_marked_leaves']
            routing['protected_marked_leaf']+=1
            protected.append(dict(case=c['index'],protected_marked_leaves=c['protected_marked_leaves'],
                forced_finite=c['forced_finite'],new_pairs_have_no_leaf=True))
        elif c['method']=='nine_context_adaptive_finite_payment':routing['adaptive_finite']+=1
        elif c['component_order_bound']<16:routing['whole_order_below16']+=1
        elif len(c['marked_component'])>16:routing['displayed_order_above16']+=1
        else:routing['eligible']+=1
    assert routing==dict(protected_marked_leaf=74,adaptive_finite=9,whole_order_below16=2045,
                         displayed_order_above16=92,eligible=394)
    eligible=[p for p in data[D.SMALL]['patterns'] if p['order_bound']>=16 and len(p['nodes'])<=16]
    assert len(eligible)==394
    saved={r['case']:r for r in rec['case_embeddings']};reconstructed=[];affected=[]
    for p in eligible:
        maps=P.embeddings(p['labels'],target)
        images=sorted({tuple(sorted(mp[i] for i in p['marks'])) for mp in maps})
        if maps:
            assert maps==saved[p['index']]['embeddings']
            assert list(map(list,images))==saved[p['index']]['marked_images']
        else:assert p['index'] not in saved
        reconstructed.append(dict(case=p['index'],map_count=len(maps),marked_images=list(map(list,images))))
        if set(images)&PAIRS:affected.append((p,maps))
    assert {p['case'] for p in reconstructed if p['map_count']}==set(saved)
    assert [p['index'] for p,maps in affected]==EXPECTED
    records=[];outcomes=Counter();mapcount=0
    for p,maps in affected:
        c=over['cases'][p['index']];core=over['cores'][c['core_key']]
        g,constraints=D.core_constraints(c,core);completions=[]
        for mp in maps:
            pair=tuple(sorted(mp[i] for i in p['marks']))
            if pair not in PAIRS:continue
            assert len(set(mp))==len(mp)
            assert all(p['labels'][i][j]==target[mp[i]][mp[j]] for i in range(len(mp)) for j in range(len(mp)))
            if c['index'] in (359,966):
                q=special_entry(c,core,constraints,target,mp)
            else:
                q,completed=D.complete(c,core,g,constraints,target,mp)
                if q['outcome']=='sealed_completion':q['current_payment']=sealed_coupon(c,q,completed)
            q['embedding']=mp;q['marked_image']=list(pair)
            outcomes[q['outcome']]+=1;mapcount+=1;completions.append(q)
        records.append(dict(case=c['index'],method=c['method'],core_key=c['core_key'],
            preferred_corner=c['preferred_corner'],marked_component=c['marked_component'],
            core_gram=core['gram'],core_inverse=core['inverse'],constraints=constraints,completions=completions))
    assert mapcount==128 and outcomes=={
        'excluded_forced_entry_projection':8,'excluded_multiple_core_contact_roots':48,
        'excluded_projection_above_norm':32,'excluded_lorentz_index':24,'sealed_completion':16}
    # Existing exact total:5 possible D5 marks times3 D11 fork roles at1/2.
    assert B.cross(('D',5),('D',11))==15
    literal=[];literal_affine16=[]
    current={(r['source'],str(r['case'])):r for r in data[R.TAG]['literal_rows']}
    for stem,index,g,white,cs in R.T.literal_inputs(data):
        key=(stem,str(list(index) if isinstance(index,tuple) else index))
        baseline=current[key];assert len(cs)==len(baseline['corners'])
        corners=[]
        for marks,unit,oldweight in cs:
            ids=sorted(set(white)|set(marks));parts=D.components(g,ids)
            shared=next((ns for ns in parts if set(marks)<=set(ns)),None);images=[]
            if shared and len(shared)==16 and Q.inertia(L.sub(g,shared))==[15,0,1]:
                labels=[[int(2-g[i][j]) if i!=j and g[i][j] else 0 for j in shared] for i in shared]
                images=sorted({tuple(sorted(mp[shared.index(i)] for i in marks)) for mp in P.embeddings(labels,target)})
                assert images and not(set(images)&PAIRS)
                literal_affine16.append(dict(source=stem,case=index,marks=marks,images=list(map(list,images))))
            corners.append(dict(marks=marks,affine_D16_images=list(map(list,images)),new_debit='0'))
        assert F(baseline['total'])>=1
        literal.append(dict(source=stem,case=index,corners=corners,unchanged_current_total=baseline['total']))
    assert len(literal)==216
    paths=[Path(__file__),Path(D.__file__),Path(L.__file__),Path(Q.__file__),Path(P.__file__),Path(B.__file__),Path(R.__file__),Path(R.T.__file__)]
    paths += [ROOT/f'{folder}/{s}.{ext}' for s in STEMS for folder,ext in
              (('scripts','py'),('results','json'),('results','md'))]
    return encode(dict(status='PASS_ORDER16_D_DISTANCE9_TWO_HALF_TRANSFER',
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in paths},
        diagram_order=16,affine_type='D15',target=target,automorphisms=autos,
        previous_cost='106',new_cost='105',changed_pairs=[list(p) for p in sorted(PAIRS)],new_pair_weight='1/2',
        new_family_totals=dict(A='80',B='105',C='103',D='105'),new_order16_maximum='105',
        new_normalized_rank15_envelope='105/16',weighted_pairs=[dict(pair=list(p),weight=str(w)) for p,w in sorted(weights.items())],
        complete_routing=dict(routing),protected_leaf_interfaces=protected,
        complete_target_domains=reconstructed,successful_patterns=len(saved),complete_target_maps=sum(len(r['embeddings']) for r in saved.values()),
        affected_contexts=EXPECTED,affected_maps=128,completion_outcomes=dict(outcomes),cases=records,
        existing_coupon=dict(small='whole_D5',small_role='any_node',large='whole_D11',large_role='fork_center',
            marked_weight='1/2',unchanged_total_for_all_three_D11_fork_roles='15/2'),
        complete_current_literal_configurations=216,literal_rows=literal,
        current_literal_affine_D16_interfaces=literal_affine16,
        new_finite_coupons=[],all_finite_weights_unchanged=True,
        short8_pairs_retained=[list(p) for p in sorted(short8)],triangle_pairs_retained=[list(p) for p in sorted(triangles)],
        triangle_minimum=2,quadrilateral_minimum=1,
        scope='Complete bounded affine-D16 target maps and arbitrary-support local transfer under the completed dimension61 rule. Two distance9 pairs, neither marked at a leaf, become halves; only an already paid whole D5/D11 fork-center half is used. Short8, all prior ideal exceptions, finite payments and coupons remain. No dependence on the new depth8 theorem, new forest count, global scalar comparison or ambient dimension claim.'))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:assert not path.exists();path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['affected_maps'],out['completion_outcomes'],'cost',out['new_cost'],flush=True)
