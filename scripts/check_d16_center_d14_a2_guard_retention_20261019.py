#!/usr/bin/env python3
"""Source-prepared saved-output audit for the narrow D16-center/D14-A2 trade.

Coordinator only: bounded parent-lock wrapper, --write or separate --verify.
No producer/census/forest engine is imported or called. No inner flock.
Universal whole-component routing and the direct seal remain written inputs.
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
from check_a2_d13_short_extra_half_retention_20261019 import amounts, signature

ROOT = W.ROOT
CERT = ROOT/'results/d16_center_d14_a2_guard_retention_20261019.json'
PLAN = 'results/d16_d13_physical_output_audit_inputs_20261019.json'
R16 = 'fresh_rank16_existing_payment_exchange_20261018'
R14 = 'fresh_finite_rank14_existing_payment_transfer_20261018'
CURRENT = 'fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018'
PINS = {
    PLAN: '6b6af762943a9279f747fcceca67de75522d6864201e067aba22630b619cab78',
    'scripts/check_a15_inner_unit_a1_absence_retention_20261018.py':
    '328c473c3fa69e547050793485f51cf399dcd667770d642431de25da5f50c908',
    'scripts/check_d16_d13_physical_output_retention_20261019.py':
    '36c3b674e296faf6cea5da48d2c1caaf586ec503763377d27089247bf4f15d70',
    'scripts/check_a2_d13_short_extra_half_retention_20261019.py':
    '59c775bf54e0a21f0e83488813fb2dd15ddb87590ade3585ef6ad3be97d23dbb',
    'scripts/fresh_d12_half_face_composition_20261018.py':
    '0f66fca4eb6084d013ff46a24b491d10dd272332969705c86de2f874a9819155',
    'scripts/run_lean_batches_20261018.py':
    '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
    'results/d16_node2_current_weight_correction_and_consumer_scope_20261019.md':
    'f53ef6016f7b64a2aa38b63c476ae275bd8c2d089015cf2a7169c9e7ce512baf',
    'results/d16_common_a1_countermodel_and_center_consumer_repair_research_20261019.md':
    '9ee57b0c8f2a220dcae63656c9260602238941e371f37f2403149ea273865cda',
    'results/d16_center_direct_seal_and_guarded_surcharge_addendum_20261019.md':
    '91b9d6bacea7cbcd4251410017ba185fdaacf52cc144b285dfdb6d52cfdd4bd6',
    'results/finite_quad_nonsingular_full_seal_rank_obstruction_20261019.md':
    '9d16f3bcdcefe85d3138aae6967a484907f3978badf02eef6f5fde3cc878a2c4',
}


def role(raw, white, marks, rank, ids=None):
    g = W.mat(raw)
    assert isinstance(marks,list) and len(marks) == len(set(marks)) == 2
    assert all(isinstance(i,int) and not isinstance(i,bool) and 0 <= i < len(g) for i in marks)
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    assert len(ids) == len(set(ids)) and set(marks) <= set(ids)
    block = next((b for b in W.parts(g,ids) if set(marks) <= set(b)),None)
    if block is None or len(block) != rank: return None
    order = d_order(g,block)
    if order is None: return None
    spine,leaves = order
    want = [2,13] if rank == 16 else [0,11]
    if not set(marks) <= set(spine) or sorted(spine.index(i) for i in marks) != want:
        return None
    other = [b for b in W.parts(g,ids) if b != block]
    return dict(nodes=block,spine=spine,short_leaves=leaves,canonical_marks=want,
        corner_nodes=ids,A1=[b for b in other if len(b) == 1],
        A2=[b for b in other if len(b) == 2 and a_order(g,b)],finite=W.pd(g,ids))


def basic_model(p, marks, kind):
    g = W.mat(p['gram']); white = sorted(p['white']); n = len(g)
    assert all(g[i][i] == 2 and all(g[i][j] == g[j][i] and g[i][j] <= 0
        for j in range(n) if j != i) for i in range(n))
    black = sorted(set(range(n))-set(white))
    assert len(black) == 4 and W.pd(g,white)
    assert all(W.pd(g,white+[b]) for b in black)
    cs = p['corners']; pairs = [tuple(sorted(c['marks'])) for c in cs]
    assert len(cs) == len(set(pairs)) == 4
    assert all(set(pair) <= set(black) and W.pd(g,sorted(set(white)|set(pair))) for pair in pairs)
    out = role(g,white,marks,16); assert out and out['finite']
    return dict(index=p['index'],kind=kind,gram=g,white=white,marks=sorted(marks),
        corner_nodes=out['corner_nodes'],D16=out['nodes'],displayed_A1=out['A1'])


def initial_models():
    doc = W.load(R16); inv = next(r for r in doc['candidate_inventory'] if r['kind'] == 'D')
    orbit = next(r for r in inv['remaining_long_orbits'] if r['orbit'] == [[2,13]])
    assert (inv['patterns'],inv['maps']) == (394,578)
    assert [r['index'] for r in orbit['selected_contexts']] == [1018,1265,1305,2409]
    assert all(r['component_size'] == 16 for r in orbit['selected_contexts'])
    result = []
    for p in doc['selected_proofs']:
        if p['index'] not in (1018,1265,1305,2409): continue
        row = next(c for c in p['corners'] if c['canonical_images'] == [[2,13]])
        m = basic_model(p,row['marks'],'rank_exclusion' if p['index'] < 1300 else 'retained_whole_A1')
        g,n = m['gram'],len(m['gram']); seal = p['sealing']
        assert F(row['payment']) == F(1,2)
        if m['kind'] == 'rank_exclusion':
            assert n == 18 and signature(g) == [17,1,0] and not m['displayed_A1']
            J = inverse(g); omitted = sorted(set(range(n))-set(m['D16']))
            assert omitted == seal['remaining_contact_nodes'] and len(omitted) == 2
            assert sorted(J[i][i] for i in omitted) == [4,8]
            assert tuple(omitted) in [tuple(sorted(c['marks'])) for c in p['corners']]
            m.update(omitted=omitted,inverse_diagonals=[J[i][i] for i in omitted],short_route_bound=2)
        else:
            assert n == 19 and signature(g) == [17,1,1] and len(m['displayed_A1']) == 1
            assert seal['method'] == 'radical'
            z = list(map(F,seal['radical'])); allowed = seal['allowed_contacts']
            assert all(sum(g[i][j]*z[j] for j in range(n)) == 0 for i in range(n))
            assert allowed and all(z[i] for i in allowed)
            assert seal['mixed_sign_shortcuts'] and all(
                z[r['nodes'][0]]*z[r['nodes'][1]] < 0 and r['length'] <= 8
                for r in seal['mixed_sign_shortcuts'])
            m.update(radical=z,allowed_contacts=allowed,short_route_bound=8)
        result.append(m)
    assert {r['index'] for r in result} == {1018,1265,1305,2409}
    return result


MAPS = {
    (616,None): ([4,3,2,0,5,6,7,8,17,16,15,14,13,9,10,11],1,12),
    (1240,None): ([14,13,12,10,15,16,17,9,8,7,6,5,4,0,1,2],11,3),
    (1284,13): ([17,13,12,10,14,15,16,9,8,7,6,5,4,0,1,2],11,3),
    (1623,3): ([17,3,2,0,4,5,6,7,8,9,16,15,14,10,11,12],1,13),
}


def dot(a,b): return sum(x*y for x,y in zip(a[:-1],b[:-1]))-2*a[-1]*b[-1]


def later_models():
    doc = W.load(R14); selected = {p['index']:p for p in doc['selected_proofs']}
    rule = next(r for r in doc['rules'] if r['kind'] == 'D')
    assert next(F(r['weight']) for r in rule['weighted_pairs'] if r['pair'] == [0,11]) == F(1,2)
    roots = []
    for i in range(15):
        v = [F()]*17; v[i]=1; v[i+1]=-1; roots.append(v)
    roots.append([F(int(i in (14,15))) for i in range(16)]+[F()])
    w3 = [F(int(i < 4)) for i in range(16)]+[F()]
    wl = [F(1,2)]*16+[F()]
    assert [dot(w3,r) for r in roots] == [int(i == 3) for i in range(16)]
    assert [dot(wl,r) for r in roots] == [int(i == 15) for i in range(16)]
    assert (dot(w3,w3),dot(wl,wl),dot(w3,wl)) == (4,4,2)
    c,d = [-x for x in w3],[-x for x in wl]; c[-1]=d[-1]=1
    assert dot(c,c) == dot(d,d) == 2 and dot(c,d) == 0
    result = []
    for (index,extra),(mp,ci,di) in MAPS.items():
        old = selected[index]
        p = old if extra is None else next(b for b in old['extra_white_branches'] if b['extra_contact'] == extra)
        p = dict(p,index=index); marks = [mp[2],mp[13]]
        m = basic_model(p,marks,'D14_A2_repair'); g = m['gram']
        assert len(g) == 18 and signature(g) == [16,1,1] and not m['displayed_A1']
        assert set(mp+[ci,di]) == set(range(18))
        native = [None]*18
        for j,k in enumerate(mp): native[k]=roots[j]
        native[ci],native[di]=c,d
        assert all(dot(native[i],native[j]) == g[i][j] for i in range(18) for j in range(18))
        assert sorted(m['white']+[mp[2],mp[13]]) == sorted(mp)
        target = next(r for r in p['corners'] if sorted(r['marks']) == sorted(marks))
        assert target['reason'] == 'existing_whole_D16_half' and F(target['payment']) == F(1,2)
        donor_marks = [ci,mp[13]]; donor = role(g,m['white'],donor_marks,14)
        assert donor and donor['finite'] and donor['A2'] == [sorted(mp[:2])]
        assert not donor['A1'] and len(donor['corner_nodes']) == 16
        old_donor = next(r for r in p['corners'] if sorted(r['marks']) == sorted(donor_marks))
        assert old_donor['reason'] == 'new_selected_weight' and F(old_donor['payment']) == F(1,2)
        # This vector identity is checked in a full Lorentz realization;
        # its universality is proved by reverse-Cauchy BEFORE sealing.
        assert all(c[j]-d[j]+w3[j]-wl[j] == 0 for j in range(17))
        assert tuple(sorted([ci,di])) in [tuple(sorted(r['marks'])) for r in p['corners']]
        m.update(extra_contact=extra,canonical_D16_nodes=mp,omitted=[ci,di],
            donor_marks=sorted(donor_marks),D14=donor['nodes'],A2=donor['A2'][0],
            current_lower='1/2',retained_lower='0',donor_current_lower='1/2',used_coupon='1/2',
            new_complete_face_lower='1',short_route_bound=2)
        result.append(m)
    current = W.load(CURRENT)['local_proof']['D16_ordinary_weighted_pairs']
    assert next(F(r['weight']) for r in current if r['pair'] == [2,13]) == F(1,2)
    return result


def saved_budget_summary():
    stem = 'fresh_d16_rank17_witness_extension_scout_20261018'
    rows = W.load(stem)['finite']['original_compact_class_maxima']
    assert len(rows) == 38 and all(F(r['maximum']) > 364 for r in rows)
    recipient_classes = []
    for r16,r17 in ((1,0),(0,1),(0,0)):
        cs = [0]*9; cs[5]=1; cs[7]=r16; cs[8]=r17
        for k5 in (0,1):
            assert not any(r['counts'] == cs and r['k5'] == k5 and r['a12'] == 0 for r in rows)
            recipient_classes.append(dict(counts=cs[:],k5=k5,a12=0))
    return dict(source=stem,above_curvature_rows=38,recipient_classes=recipient_classes,
        unlisted_original_class_upper='364',recipient_current_upper='364',new_surcharge='1/2',
        recipient_successor_upper='729/2',unlisted_class_and_chronology_are_written_premises=True)


def run():
    for name,expected in PINS.items(): assert W.digest(ROOT/name) == expected,name
    plan = json.loads((ROOT/PLAN).read_text()); files = W.input_paths()
    assert [str(p.relative_to(ROOT)) for p in files] == [r['path'] for r in plan['files']]
    assert len(files) == 54
    for row in plan['files']:
        p = ROOT/row['path']; assert p.stat().st_size == row['bytes'] and W.digest(p) == row['sha256'],row['path']
    fixed = initial_models()+later_models(); budget = saved_budget_summary()
    counts,hits,recipients,unknown = Counter(),[],[],[]
    nonfinite = Counter()
    def one(stem,where,raw,white,row,index=None):
        counts[stem] += 1
        recipient = role(raw,white,row['marks'],14,row.get('corner_nodes'))
        if recipient and recipient['finite']:
            recipients.append(dict(source=stem,path=where,index=index,marks=row['marks'],role=recipient,
                current_lower='1/2',increment='1/2',operative_guard='a1=0,a2>0,R16/17,D>=18'))
        out = role(raw,white,row['marks'],16,row.get('corner_nodes'))
        if out is None: return
        if not out['finite']: nonfinite[stem] += 1; return
        rec = dict(source=stem,path=where,index=index,marks=row['marks'],role=out)
        values,malformed = amounts(row)
        rec['operative_amounts'] = {k:str(v) for k,v in values.items()}
        match = [m for m in fixed if W.mat(raw) == m['gram'] and white is not None
            and sorted(white) == m['white'] and sorted(row['marks']) == m['marks']
            and out['corner_nodes'] == m['corner_nodes'] and out['nodes'] == m['D16']]
        if not malformed and len(match) == 1 and all(0 <= v <= 1 for v in values.values()):
            m = match[0]; rec.update(status=m['kind'],native_case=m['index'],extra_contact=m.get('extra_contact'),
                current_lower='1/2',retained_lower='1/2' if m['kind'] == 'retained_whole_A1' else '0',
                used_coupon='1/2' if m['kind'] == 'D14_A2_repair' else '0',
                changed_branch_guard='a1=0,ell=1,rho>=2; otherwise retain old D16 half',
                full_face_route='independent short2/8 or the matched complete-current model')
        elif not malformed and values and all(v == 0 for v in values.values()):
            rec.update(status='explicit_zero_no_debit_obligation',current_lower='1/2',retained_lower='0',used_coupon='0')
        else:
            rec.update(status='UNKNOWN_D16_CENTER_AMOUNT_OR_ACTUAL_MODEL',malformed_amount_fields=malformed)
            unknown.append(rec)
        hits.append(rec)
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
        newer.append(dict(source=stem,section=list(keys),corners=counts[stem]-before)); del doc,section
    chosen = []
    for stem,expected in zip(W.CHOSEN,(280,52)):
        doc = W.load(stem); rs = [r for r in outputs(doc) if F(r['payment']) > 0]
        assert len(rs) == expected
        for row in rs:
            assert 1 <= len(row['roles']) <= 2
            types = [chosen_type(r) for r in row['roles']]
            assert 'D16' not in types,('unmapped chosen D16',stem,row)
        chosen.append(dict(source=stem,positive_outputs=expected,D16_target_outputs=0)); del doc,rs
    extensions = []
    for stem in W.EXTENSIONS:
        doc = W.load(stem); before = counts[stem]
        for where,row,raw,white,index in W.physical_rows(doc): one(stem,where,raw,white,row,index)
        extensions.append(dict(source=stem,physical_corners=counts[stem]-before)); del doc
    if unknown:
        raise AssertionError(json.dumps(dict(status='UNRESOLVED_D16_CENTER_OUTPUT',
            unknown_count=len(unknown),unknown=unknown[:12])))
    assert {h['native_case'] for h in hits if h['status'] == 'rank_exclusion'} == {1018,1265}
    assert {h['native_case'] for h in hits if h['status'] == 'retained_whole_A1'} == {1305,2409}
    assert {(h['native_case'],h['extra_contact']) for h in hits if h['status'] == 'D14_A2_repair'} == set(MAPS)
    sources = {r['path']:r['sha256'] for r in plan['files']} | PINS
    sources[str(Path(__file__).relative_to(ROOT))] = W.digest(Path(__file__))
    return dict(status='PASS_CURRENT_D16_CENTER_D14_A2_SAVED_INTERFACE',sources=sources,
        input_files=len(files),input_bytes=sum(r['bytes'] for r in plan['files']),historical_packets=22,
        historical_corners=legacy,literal_configurations=literal,newer_sections=newer,extensions=extensions,
        chosen_D12=chosen,physical_corner_counts=dict(counts),target_hits=hits,
        target_summary=dict(Counter(h['status'] for h in hits)),nonfinite_target_displays=dict(nonfinite),
        recipient_hits=recipients,recipient_canonical_pair=[0,11],target_canonical_pair=[2,13],
        full_models=[{k:v for k,v in m.items() if k != 'gram'} for m in fixed],unknown_targets=0,
        saved_recipient_budget_summary=budget,
        current_center_weight='1/2',recipient_current_weight='1/2',recipient_new_weight='1',
        exact_cost_delta='-(1/2)*d16*[a1=0]*ell*[rho>=2]+(1/2)*d14*[a1=0]*[a2>0]*ell*[rho>=1]',
        rho='a16+d16+a17+d17',ell='Boolean whole D_s>=18',positive_family_maximum_surcharge='1/2',
        safe_new_positive_family_upper='729/2',safe_complete_finite_upper='2965/8',
        budget_is_written_sequential_premise=True,complete_new_finite_maximum=None,ideal_weights_changed=False,
        scope='Saved physical-output audit only. Whole-component routing, native contact minimum, '
        'reverse Lorentz Cauchy direct seal before spectator orthogonality, actual short/basic '
        'priority, the selected-domain completeness/rank argument and the saved-summary sequential '
        'budget remain written premises. No shared-A1 capacity, node2 leaf debit, historical '
        'census, forest optimization, new attained maximum or dimension improvement.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--write',action='store_true'); mode.add_argument('--verify',action='store_true')
    args = parser.parse_args(); signal.alarm(45)
    result = json.loads(json.dumps(run(),default=lambda x:str(x) if isinstance(x,F) else x))
    payload = json.dumps(result,sort_keys=True,indent=2)+'\n'; assert len(payload.encode()) < 1024**2
    if args.verify: assert json.loads(CERT.read_text()) == result,'Saved certificate differs'
    else:
        with CERT.open('x') as stream: stream.write(payload)
    signal.alarm(0)
    print(json.dumps({k:result[k] for k in ('status','input_files','input_bytes','target_summary')}))


if __name__ == '__main__': main()
