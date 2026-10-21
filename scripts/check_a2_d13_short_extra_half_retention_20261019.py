#!/usr/bin/env python3
"""Fail-closed current A2/D13-short saved-output audit; no enumeration.

Coordinator execution only, bounded parent-lock wrapper --write/--verify.
Imports saved readers and pure graph/matrix helpers; never an old run body.
The direct whole-target seal and actual basic-route priority are written inputs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import check_a15_inner_unit_a1_absence_retention_20261018 as W
from check_d16_d13_physical_output_retention_20261019 import (
    a_order, d_order, inverse, chosen_type, outputs)

ROOT = W.ROOT
CERT = ROOT/'results/a2_d13_short_extra_half_retention_20261019.json'
PLAN = 'results/d16_d13_physical_output_audit_inputs_20261019.json'
OLD = 'fresh_finite_b15_d15_existing_payment_transfer_20261018'
SPARSE = 'fresh_rank15_sparse_payment_exchange_20261018'
PINS = {
    PLAN: '6b6af762943a9279f747fcceca67de75522d6864201e067aba22630b619cab78',
    'scripts/check_a15_inner_unit_a1_absence_retention_20261018.py':
    '328c473c3fa69e547050793485f51cf399dcd667770d642431de25da5f50c908',
    'scripts/check_d16_d13_physical_output_retention_20261019.py':
    '36c3b674e296faf6cea5da48d2c1caaf586ec503763377d27089247bf4f15d70',
    'results/finite_quad_nonsingular_full_seal_rank_obstruction_20261019.md':
    '9d16f3bcdcefe85d3138aae6967a484907f3978badf02eef6f5fde3cc878a2c4',
    'results/a2_d13_short_whole_target_direct_seal_20261019.md':
    '5090ee74aacaca547354fd61c699d224ccb3027e4c48af1ddb1ffe8afb38623a',
    'results/d16_node2_current_weight_correction_and_consumer_scope_20261019.md':
    'f53ef6016f7b64a2aa38b63c476ae275bd8c2d089015cf2a7169c9e7ce512baf',
}


def target(raw, white, marks, ids=None):
    g = W.mat(raw)
    assert isinstance(marks,list) and len(marks) == len(set(marks)) == 2
    assert all(isinstance(i,int) and not isinstance(i,bool) and 0 <= i < len(g) for i in marks)
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    assert len(ids) == len(set(ids)) and set(marks) <= set(ids)
    blocks = W.parts(g,ids)
    first,second = [next(b for b in blocks if m in b) for m in marks]
    if first is second or first == second:
        return None
    for a,d,am,dm in ((first,second,*marks),(second,first,marks[1],marks[0])):
        if len(a) != 2 or len(d) != 13: continue
        ao,do = a_order(g,a),d_order(g,d)
        if ao and do and am in ao and dm in do[1]:
            return dict(A2=a,D13=d,A2_mark=am,D13_leaf=dm,D13_spine=do[0],
                D13_short_leaves=do[1],corner_nodes=ids,finite=W.pd(g,ids))
    return None


def signature(raw):
    """Exact congruence elimination, with a two-by-two offdiagonal pivot."""
    a = [[F(v) for v in row] for row in raw]; positive = negative = 0
    while a:
        n = len(a)
        pivot = next((i for i in range(n) if a[i][i]),None)
        if pivot is not None:
            p = a[pivot][pivot]; positive += p > 0; negative += p < 0
            keep = [i for i in range(n) if i != pivot]
            a = [[a[i][j]-a[i][pivot]*a[pivot][j]/p for j in keep] for i in keep]
        else:
            pair = next(((i,j) for i in range(n) for j in range(i+1,n) if a[i][j]),None)
            if pair is None: return [positive,negative,n]
            i,j = pair; p = a[i][j]; positive += 1; negative += 1
            keep = [k for k in range(n) if k not in pair]
            a = [[a[r][s]-(a[r][i]*a[j][s]+a[r][j]*a[i][s])/p for s in keep] for r in keep]
    return [positive,negative,0]


def models():
    doc = W.load(OLD); selected = {p['index']:p for p in doc['selected_proofs']}
    result = []
    expected = {1162:([2,11],[0,1],[13]),1175:([2,12],[0,1],[13]),
                1910:([1,12],[3],[10,11]),2303:([2,12],[3],[10,11])}
    for index,(marks,left,right) in expected.items():
        p = selected[index]; g = W.mat(p['gram']); white = sorted(p['white'])
        assert len(g) == 18 and all(g[i][i] == 2 for i in range(18))
        assert signature(g) == [17,1,0] and len(white) == 14 and W.pd(g,white)
        assert all(g[i][j] <= 0 for i in range(18) for j in range(18) if i != j)
        assert p['core_key'] == '(25, 25, ((8, 9), (9, 8)))'
        black = sorted(set(range(18))-set(white)); cs = p['corners']
        assert len(black) == 4 and len(cs) == 4
        assert {m for c in cs for m in c['marks']} == set(black)
        assert len({tuple(c['marks']) for c in cs}) == 4
        assert all(W.pd(g,sorted(set(white)|set(c['marks']))) for c in cs)
        assert all(W.pd(g,white+[b]) for b in black)
        row = next(c for c in cs if c['marks'] == marks)
        role = target(p['gram'],white,marks); assert role and role['finite']
        assert F(row['payment']) == F(1,2)
        assert row['reason'] == 'existing_whole_A2_D13_short_fork_leaf_half'
        whole = sorted(role['A2']+role['D13'])
        assert len(whole) == 15 and sorted(set(range(18))-set(whole)) == left+right
        J = inverse(g)
        assert all(sum(g[i][k]*J[k][j] for k in range(18)) == int(i == j)
                   for i in range(18) for j in range(18))
        expected_blocks = ([[30,15],[15,8]],[[4]]) if index < 1900 else ([[4]],[[30,15],[15,8]])
        assert tuple([[J[i][j] for j in group] for i in group] for group in (left,right)) == expected_blocks
        assert min(J[i][j] for group in (left,right) for i in group for j in group) >= 4
        end1,end2 = (0,13) if index < 1900 else (3,10)
        intermediate = 1 if index < 1900 else 11
        near = end1 if index < 1900 else end2
        assert end1 in black and end2 in black and intermediate in white
        assert [end1,end2] in [sorted(c['marks']) for c in cs]
        assert g[near][intermediate] == -1
        # Contacts to both blocks give one actual two/three-edge path through
        # extra white e, never through another black side. Zero contacts seal.
        assert (left,right) == (([0,1],[13]) if index < 1900 else ([3],[10,11]))
        result.append(dict(index=index,gram=g,white=white,marks=sorted(marks),
            corner_nodes=role['corner_nodes'],A2=role['A2'],D13=role['D13'],
            remaining_blocks=[left,right],inverse_blocks=expected_blocks,
            black_shortcut_endpoints=[end1,end2],old_white_intermediate=intermediate,
            independent_short_route_bound=3,core_signature=[17,1,0]))
    del doc, selected
    new = W.load(SPARSE)
    coupon = new['upgraded_coupon']
    assert coupon == dict(small='whole A2',small_marks=[0,1],large='whole D13',
        large_marks=[11,12],old_weight='1/2',new_weight='1',
        component_pair_cost_increase='2',new_total_component_pair_cost='5')
    for p in new['selected_proofs']:
        if p['index'] not in expected: continue
        m = next(r for r in result if r['index'] == p['index'])
        assert W.mat(p['gram']) == m['gram'] and sorted(p['white']) == m['white']
        row = next(c for c in p['corners'] if sorted(c['marks']) == m['marks'])
        assert F(row['payment']) == 1 and row['reason'] == 'upgraded_whole_A2_D13_short_fork_leaf_unit'
    context = W.load('fresh_rank13_variable_context_retention_20261018')
    assert F(context['chosen_mark_weights']['A2_D13_center']) == F(23,64)
    # All subsequent named changes have disjoint roles: this chronology remains
    # an explicit written source premise, not a consequence of the hashes.
    assert 4+2*F(23,64) == F(151,32) and 4*F(1,2)+2*F(23,64) == F(87,32)
    return result


# Aggregate fields are alternative lower-bound statements, not additive terms.
# Historical before/previous/saved values, whole-face totals and intrinsic bases
# are metadata. current_lower/retained_lower and used_coupon are NOT ignored.
AGGREGATES = ('payment','retained_payment','retained_payment_lower_bound',
    'retained_lower','retained_weight','new_lower','new_payment','new_weight',
    'current_weight','current_lower','lower','weight','new','joint_lower')
PARTS = ('used_selected_half','used_coupon','finite_unit','finite_half',
         'unit','half_coupon','old_half','finite_small17_unit','cross_coupon')
SUM_SCHEMAS = (('used_selected_half','used_coupon'),('finite_unit','finite_half'),
    ('unit','half_coupon'),('unit','old_half'),('finite_small17_unit','half_coupon'))
METADATA = ('base','ideal_weight','coupon','coupon_roles','other_same_component_units_used',
    'full_A1_D11_upgrade','parameter_values','coefficient','total','total_lower_bound',
    'payment_total','total_payment','finite_unit_total','finite_coupon_total')


def amounts(row):
    values, malformed = {},[]
    for key in AGGREGATES+PARTS:
        if key not in row: continue
        value = row[key]
        if isinstance(value,bool) or not isinstance(value,(str,int)):
            malformed.append(key); continue
        try: values[key] = F(value)
        except (ValueError,ZeroDivisionError): malformed.append(key)
    if 'parameter_values' in row:
        table = row['parameter_values']
        if not isinstance(table,dict) or '9/64' not in table:
            malformed.append('parameter_values/9/64')
        else:
            try: values['parameter_values/9/64'] = F(table['9/64'])
            except (TypeError,ValueError,ZeroDivisionError): malformed.append('parameter_values/9/64')
    for group in SUM_SCHEMAS:
        if all(key in row for key in group):
            if all(key in values for key in group):
                values['+'.join(group)] = sum(values[k] for k in group)
            else: malformed.append('malformed additive schema '+','.join(group))
    if any(key in row for key in ('used_selected_half','used_coupon')):
        if not all(key in row for key in ('used_selected_half','used_coupon')):
            malformed.append('incomplete additive selected-half/coupon schema')
    for key in row:
        if key in AGGREGATES+PARTS+METADATA: continue
        if key.startswith(('old_','previous','saved_','before','current_before','exact_current_before')): continue
        if any(token in key for token in ('weight','payment','lower','coupon','half','unit')):
            malformed.append('unclassified possible amount '+key)
    return values,malformed


def run():
    for name,expected in PINS.items(): assert W.digest(ROOT/name) == expected,name
    plan = json.loads((ROOT/PLAN).read_text()); paths = W.input_paths()
    assert [str(p.relative_to(ROOT)) for p in paths] == [r['path'] for r in plan['files']]
    for row in plan['files']:
        p = ROOT/row['path']
        assert p.stat().st_size == row['bytes'] and W.digest(p) == row['sha256'],row['path']
    fixed = models(); counts,hits,unknown,nonfinite = Counter(),[],[],0

    def one(stem,where,raw,white,row,index=None):
        nonlocal nonfinite
        counts[stem] += 1
        role = target(raw,white,row['marks'],row.get('corner_nodes'))
        if role is None: return
        if not role['finite']: nonfinite += 1; return
        record = dict(source=stem,path=where,index=index,marks=row['marks'],role=role)
        values,malformed = amounts(row)
        record['operative_amounts'] = {k:str(v) for k,v in values.items()}
        match = [m for m in fixed if W.mat(raw) == m['gram'] and white is not None
            and sorted(white) == m['white'] and sorted(row['marks']) == m['marks']
            and role['corner_nodes'] == m['corner_nodes']
            and role['A2'] == m['A2'] and role['D13'] == m['D13']]
        if not malformed and len(match) == 1 and all(0 <= v <= 1 for v in values.values()):
            record.update(status='complete_whole_target_short3_or_rank_exclusion',native_case=match[0]['index'])
        elif not malformed and values and all(0 <= v <= F(1,2) for v in values.values()):
            record['status'] = 'unchanged_explicit_half_or_less'
        else:
            record.update(status='UNKNOWN_TARGET_DEMAND_OR_ACTUAL_MODEL',malformed_amount_fields=malformed)
            unknown.append(record)
        hits.append(record)

    legacy = 0
    for stem,keys in W.finite_inventory().items():
        doc = W.load(stem)
        for key in keys:
            for where,row,raw,white,index in W.physical_rows(doc[key],key,require_components=True):
                one(stem,where,raw,white,row,index); legacy += 1
        del doc
    assert legacy == 564
    common = W.native(1142)[1]['gram']; literal = 0
    for stem,index,raw,white,cs in W.literal_rows(common):
        for j,row in enumerate(cs): one(stem,f'literal/{index}/{j}',raw,white,row,index)
        literal += 1
    assert literal == 216
    newer = []
    for stem,keys in W.NEWER:
        doc = W.load(stem); section = doc
        for key in keys: section = section[key]
        before = counts[stem]
        for where,row,raw,white,index in W.physical_rows(section,'/'.join(keys)):
            one(stem,where,raw,white,row,index)
        assert counts[stem] > before,('missing physical collection',stem,keys)
        newer.append(dict(source=stem,section=list(keys),corners=counts[stem]-before))
        del doc,section
    chosen = []
    for stem,expected in zip(W.CHOSEN,(280,52)):
        doc = W.load(stem); rs = [r for r in outputs(doc) if F(r['payment']) > 0]
        assert len(rs) == expected
        for row in rs:
            assert 1 <= len(row['roles']) <= 2,'unknown positive chosen role schema'
            types = [chosen_type(r) for r in row['roles']]
            assert not ('A2' in types and 'D13' in types),('unmapped chosen A2/D13',stem,row)
        chosen.append(dict(source=stem,positive_outputs=expected,target_outputs=0))
        del doc,rs
    extensions = []
    for stem in W.EXTENSIONS:
        doc = W.load(stem); before = counts[stem]
        for where,row,raw,white,index in W.physical_rows(doc): one(stem,where,raw,white,row,index)
        extensions.append(dict(source=stem,physical_corners=counts[stem]-before)); del doc
    if unknown:
        raise AssertionError(json.dumps(dict(status='UNRESOLVED_A2_D13_PHYSICAL_OUTPUT',
            unknown_count=len(unknown),unknown=unknown[:12])))
    assert {h['native_case'] for h in hits if h['status'].startswith('complete_')} == {1162,1175,1910,2303}
    assert {h['native_case'] for h in hits if h['source'] == SPARSE
            and h['status'].startswith('complete_')} == {1162,1175,1910,2303}
    sources = {r['path']:r['sha256'] for r in plan['files']} | PINS
    sources[str(Path(__file__).relative_to(ROOT))] = W.digest(Path(__file__))
    return dict(status='PASS_CURRENT_A2_D13_SHORT_SAVED_PHYSICAL_INTERFACE',sources=sources,
        input_files=len(paths),input_bytes=sum(p.stat().st_size for p in paths),
        historical_packets=22,historical_corners=legacy,literal_configurations=literal,
        newer_sections=newer,extensions=extensions,chosen_D12=chosen,physical_corner_counts=dict(counts),
        full_actual_models=[{k:v for k,v in m.items() if k != 'gram'} for m in fixed],
        target_hits=hits,target_summary=dict(Counter(h['status'] for h in hits)),
        nonfinite_target_displays=nonfinite,unknown_targets=0,
        current_short_leaf_weight='1',candidate_short_leaf_weight='1/2',
        current_center_weight='23/64',current_ordinary_pair_cost='151/32',candidate_ordinary_pair_cost='87/32',
        finite_candidate_delta='-2*a2*d13',new_positive_families=0,ideal_weights_changed=False,
        universal_geometry_formalized=False,complete_new_finite_maximum=None,
        scope='Saved-output audit only. Whole-target routing, native contact minimum, actual short2/3 '
        'basic priority, common-white independence and the rank contradiction remain written premises. '
        'Full-model recognition uses the same-corner whole A2+D13 seal, never a transported D15 guard. '
        'No D16 new-orbit refund, census, forest optimization, attainment or dimension improvement.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true'); mode.add_argument('--verify',action='store_true')
    args = parser.parse_args(); signal.alarm(45)
    result = json.loads(json.dumps(run(),default=lambda v:str(v) if isinstance(v,F) else v))
    payload = json.dumps(result,sort_keys=True,indent=2)+'\n'; assert len(payload.encode()) < 1024**2
    if args.verify: assert json.loads(CERT.read_text()) == result,'Saved certificate differs'
    else:
        with CERT.open('x') as stream: stream.write(payload)
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ('status','input_files','input_bytes','target_summary')}))


if __name__ == '__main__':
    main()
