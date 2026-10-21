#!/usr/bin/env python3
"""Saved current-output audit for two physical cross-coupon capacities.

Coordinator only: bounded wrapper --parent-lock -- THIS --write/--verify.
No inner flock; 45s child alarm. No old run(), producer, optimizer or census.
Universal routing, actual walls and whole-component seals remain written inputs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import check_a15_inner_unit_a1_absence_retention_20261018 as W
from check_d16_d13_physical_output_retention_20261019 import (
    a_order, d_order, chosen_type, outputs)

ROOT = W.ROOT
TAG = 'a1_rank16_physical_coupon_retention_20261020'
CERT = ROOT/'results'/(TAG+'.json')
PLAN = 'results/d16_d13_physical_output_audit_inputs_20261019.json'
COUPLED = 'fresh_d17_leaf_a16_context_coupled_transfer_20261018'
REPORT = 'results/a1_rank16_physical_coupon_capacity_research_20261020.md'
CORRECTION = 'results/a1_rank16_partial_display_scope_correction_20261020.md'
PINS = {
    PLAN:'6b6af762943a9279f747fcceca67de75522d6864201e067aba22630b619cab78',
    'scripts/run_lean_batches_20261018.py':
        '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    'scripts/check_a15_inner_unit_a1_absence_retention_20261018.py':
        '328c473c3fa69e547050793485f51cf399dcd667770d642431de25da5f50c908',
    'scripts/fresh_d12_half_face_composition_20261018.py':
        '0f66fca4eb6084d013ff46a24b491d10dd272332969705c86de2f874a9819155',
    'scripts/check_d16_d13_physical_output_retention_20261019.py':
        '36c3b674e296faf6cea5da48d2c1caaf586ec503763377d27089247bf4f15d70',
    'results/fresh_ideal7_finite_cross_payments_20261017.md':
        '7c9d08979421838d5620593c74f428977727b55798941b0c57f588f085c41d4a',
    'results/fresh_d17_leaf_a16_context_coupled_transfer_20261018.md':
        '19ed13114b8dddcf5a03febea1180209fa53f52b360ca0ec8c58f2653f4c0738',
    REPORT:'3612f2c2aefc8ff7cc493cadc524805c3ab8144ecf3211c429a1762e4a90636c',
    CORRECTION:'180ec1307f0136ad9568574696eb6978446a71944bf577e20cbeff4e452a340a',
    'results/fresh_small_d16_coupon_and_ideal17_rebalance_20261018.md':
        'c0c10822e322c31689b21a53b2bc64b03522dd7077fbcac1758dfdf19ddacb28',
    'results/review_d16_d13_physical_output_execution_root_20261019.md':
        '978e164f04b6c418c41c5b7e1bd54efb3101bd20826c0f7639ccd91952425262',
    'results/review_d16_center_d14_a2_execution_root_20261019.md':
        '0b4866c836109a52ca816600921ec79bf855600555e15b0406e9c0a9db710841',
    'results/review_a2_d13_short_retention_execution_root_20261019.md':
        '3d76e0c80bc3dbc9e486d5ddf66b58374ddeebf397d51f868bd1830b155b6638',
    'results/d17_node2_a12_d5_single_coupon_capacity_research_20261020.md':
        '60bd1f37cb5afaa30aac5182581b5ffab7b88d94b6f04fd9547be9f96a6c90fa',
    'results/d16_guarded_long_leaf_half_return_with_d17_capacity_20261020.md':
        '665dd7cd3a620bc6eb84455a23c5ea557149e5248159f6c3084022d752443c88',
}
EDGES = [(1,0),(2,0),(3,2),(4,0),(5,4),(6,5),(7,6),(8,7),(9,8),
         (11,10),(12,10),(13,12),(14,10),(15,14),(16,15),(17,16),(18,9),(18,17)]
FORK_CORNERS = {(400,(0,13)),(675,(1,13)),(813,(3,10)),
                (828,(3,11)),(841,(2,13)),(841,(3,12))}
ZERO_PARTIAL = {2367:[3,14],2371:[3,15],2511:[4,13],2539:[5,13]}


def target(raw, white, marks, ids=None):
    g = W.mat(raw); n = len(g)
    assert len(marks) == len(set(marks)) == 2
    assert all(isinstance(i,int) and not isinstance(i,bool) and 0 <= i < n for i in marks)
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    assert len(ids) == len(set(ids)) and set(marks) <= set(ids)
    assert all(isinstance(i,int) and not isinstance(i,bool) and 0 <= i < n for i in ids)
    blocks = W.parts(g,ids)
    aa,bb = [next(b for b in blocks if m in b) for m in marks]
    if aa == bb: return None
    for small,large,am,lm in ((aa,bb,*marks),(bb,aa,marks[1],marks[0])):
        if len(small) != 1 or len(large) != 16: continue
        ao,do = a_order(g,large),d_order(g,large)
        if ao and lm in (ao[0],ao[-1]):
            return dict(kind='A16_endpoint',A1=am,large=large,mark=lm,
                        corner_nodes=ids,finite=W.pd(g,ids))
        if do and lm in [do[0][-1],*do[1]]:
            return dict(kind='D16_fork',A1=am,large=large,mark=lm,
                        corner_nodes=ids,finite=W.pd(g,ids))
    return None


def base_gram():
    g = [[2*int(i == j) for j in range(19)] for i in range(19)]
    for i,j in EDGES: g[i][j] = g[j][i] = -1
    return g


def models(g):
    result = []; old = W.load(W.OLD17)['order17']['replacement_cases']
    for p in old:
        white = sorted(set(range(19))-set(p['black']))
        for row in p['corners']:
            out = target(g,white,row['marks'])
            if not out or out['kind'] != 'D16_fork': continue
            assert out['finite'] and F(row['cross_coupon']) == F(row['weight']) == F(1,2)
            result.append(dict(kind=out['kind'],index=p['case'],marks=sorted(row['marks']),
                               white=white,gram=g))
    assert {(p['index'],tuple(p['marks'])) for p in result} == FORK_CORNERS
    stage = W.load(COUPLED)['local_proof']['old_ideal17_seals']
    for p in stage:
        assert W.mat(p['gram']) == W.mat(g)
        z = list(map(F,p['full_radical']))
        assert all(sum(g[i][j]*z[j] for j in range(19)) == 0 for i in range(19))
        assert z[3]*z[13] < 0
        ns = p['affine_nodes']; kernel = list(map(F,p['positive_affine_kernel']))
        assert ns == sorted(set(range(19))-{3,13}) and all(v > 0 for v in kernel)
        assert all(sum(g[i][j]*kernel[k] for k,j in enumerate(ns)) == 0 for i in ns)
        assert W.pd(W.mat(g),p['white'])
        for row in p['corners']:
            out = target(g,p['white'],row['marks'])
            if not out or out['kind'] != 'A16_endpoint': continue
            assert out['finite'] and F(row['payment']) == F(1,4)
            assert [r['role'] for r in row['roles']] == ['A1_A16_endpoint_quarter']
            result.append(dict(kind=out['kind'],index=p['case'],marks=sorted(row['marks']),
                               white=sorted(p['white']),gram=g))
    assert {(p['index'],tuple(p['marks'])) for p in result if p['kind'] == 'A16_endpoint'} == {
        (675,(3,13)),(828,(3,13))}
    partial = []
    for p in W.load(W.PARTIAL)['cases']:
        if p['index'] not in ZERO_PARTIAL: continue
        marks = ZERO_PARTIAL[p['index']]
        row = next(c for c in p['corners'] if c['marks'] == marks)
        assert target(p['gram'],p['white'],marks) is None
        assert all(F(row[k]) == 0 for k in ('finite_half','finite_unit','ideal_weight'))
        ids = sorted(set(p['white'])|set(marks)); blocks = W.parts(W.mat(p['gram']),ids)
        assert sorted(map(len,blocks)) == [1,16]
        large = next(b for b in blocks if len(b) == 16); do = d_order(W.mat(p['gram']),large)
        marked = next(m for m in marks if m in large)
        assert do and marked not in [do[0][-1],*do[1]]
        partial.append(dict(index=p['index'],marks=marks,status='explicit_zero_nonfork_display'))
    assert len(partial) == 4
    return result,partial


def full_match(raw,white,marks,model):
    """Exact normalized FULL core isomorphism; only orthogonal white extras.

    All edges must be simple after positive root rescaling. Vertex colors
    preserve common whites, both marked blacks and the two omitted blacks.
    Backtracking tries every consistent unused vertex, never just a subcore.
    """
    if white is None: return None
    g = W.mat(raw); nodes = next(b for b in W.parts(g,list(range(len(g)))) if marks[0] in b)
    if len(nodes) != 19 or not set(marks) <= set(nodes): return None
    if not (set(range(len(g)))-set(nodes)) <= set(white): return None
    if any(g[i][i] <= 0 or any(g[i][j] > 0 or (g[i][j] and
        4*g[i][j]**2 != g[i][i]*g[j][j]) for j in nodes if j != i) for i in nodes): return None
    edge = lambda i,j: i != j and bool(g[i][j])
    mg = model['gram']; mw,mm = set(model['white']),set(model['marks'])
    color = {i:(i in white,i in marks,sum(edge(i,j) for j in nodes)) for i in nodes}
    mc = {i:(i in mw,i in mm,sum(bool(mg[i][j]) for j in range(19) if j != i)) for i in range(19)}
    candidates = {i:[j for j in nodes if color[j] == mc[i]] for i in range(19)}
    def visit(mapping,used):
        if len(mapping) == 19: return [mapping[i] for i in range(19)]
        choices = [(i,[j for j in candidates[i] if j not in used and all(
            edge(j,v) == bool(mg[i][k]) for k,v in mapping.items())])
            for i in range(19) if i not in mapping]
        i,options = min(choices,key=lambda r:len(r[1]))
        for j in options:
            found = visit(mapping|{i:j},used|{j})
            if found is not None: return found
        return None
    return visit({},set())


AMOUNTS = ('payment','weight','lower','retained_lower','retained_weight','new_lower',
           'new_payment','new_weight','current_weight','current_lower','new','joint_lower',
           'retained_payment','retained_payment_lower_bound','finite_unit','finite_half',
           'cross_coupon','used_coupon','used_selected_half','unit','half_coupon','old_half')
META = ('base','ideal_weight','coupon','coupon_roles','other_same_component_units_used',
        'full_A1_D11_upgrade','coefficient','total','total_lower_bound','payment_total',
        'total_payment','finite_unit_total','finite_coupon_total')


def amounts(row):
    values,unknown = {},[]
    for key in AMOUNTS:
        if key not in row: continue
        try:
            assert isinstance(row[key],(str,int)) and not isinstance(row[key],bool)
            values[key] = F(row[key])
        except (AssertionError,ValueError,ZeroDivisionError): unknown.append(key)
    for group in (('finite_unit','finite_half'),('used_selected_half','used_coupon'),
                  ('unit','half_coupon'),('unit','old_half')):
        if all(k in values for k in group): values['+'.join(group)] = sum(values[k] for k in group)
    for key in row:
        if key in AMOUNTS+META or key.startswith(('old_','previous','saved_','before','current_before')): continue
        if any(t in key for t in ('weight','payment','lower','coupon','half','unit')): unknown.append(key)
    return values,unknown


def scalar_identities():
    an = lambda i,j: F((min(i,j)+1)*(16-max(i,j)),17)
    assert (an(2,2),an(14,14),an(2,14),an(2,13)) == (F(42,17),F(30,17),F(6,17),F(9,17))
    assert an(2,2)-F(8,17) == an(14,14)+F(1,2)-F(9,34) == 2
    assert an(2,13)-F(8,17) == F(1,17)
    assert an(14,14)-F(9,34) == F(3,2)
    # D16 leaf duals are half-sign vectors; next-to-end dual is e0+e1.
    sigma = [F(1,2)]*16; other = sigma[:-1]+[-F(1,2)]; om = [F(1),F(1)]+[F()]*14
    dot = lambda a,b:sum(x*y for x,y in zip(a,b))
    assert (dot(sigma,sigma),dot(om,om),dot(sigma,om),dot(sigma,other)) == (4,2,1,F(7,2))
    assert dot(sigma,other)-2 == dot(om,om)-F(1,2) == F(3,2)
    assert 3*F(433,4)+36 == F(1443,4)
    assert F(657,2)+28+F(9,2) == 361 and F(657,2)+27+F(1,4) == F(1423,4)
    return dict(A16_endpoint_cap='1/4*[a1*a16>0]',D16_fork_cap='3/2*d16*[a1>0]',
                three_rank16_family_upper='361',global_upper_unchanged='2965/8')


def run():
    own = Path(__file__); own_before = W.digest(own)
    for name,expected in PINS.items(): assert W.digest(ROOT/name) == expected,name
    plan = json.loads((ROOT/PLAN).read_text()); paths = W.input_paths()
    assert len(paths) == 54 and [str(p.relative_to(ROOT)) for p in paths] == [r['path'] for r in plan['files']]
    for r in plan['files']:
        p = ROOT/r['path']; assert p.stat().st_size == r['bytes'] and W.digest(p) == r['sha256'],r['path']
    common = base_gram(); fixed,partial = models(common)
    counts,hits,unknown,nonfinite = Counter(),[],[],Counter()
    def one(stem,where,raw,white,row,index=None):
        counts[stem] += 1
        role = target(raw,white,row['marks'],row.get('corner_nodes'))
        if role is None: return
        if not role['finite']: nonfinite[role['kind']] += 1; return
        vals,bad = amounts(row); limit = F(1,4) if role['kind'] == 'A16_endpoint' else F(1,2)
        rec = dict(source=stem,path=where,index=index,marks=row['marks'],role=role,
                   declared_amounts={k:str(v) for k,v in vals.items()},current_lower=str(limit))
        matches = [(m,full_match(raw,white,row['marks'],m)) for m in fixed if m['kind'] == role['kind']]
        matches = [(m,mp) for m,mp in matches if mp is not None]
        if not bad and vals and all(v == 0 for v in vals.values()):
            rec.update(status='explicit_zero_no_chosen_demand',retained_lower='0',used_coupon='0')
        elif not bad and matches and all(0 <= v <= limit for v in vals.values()):
            rec.update(status='retained_full_physical_coupon_model',retained_lower=str(limit),used_coupon=str(limit),
                       missing_amount_conservatively_activates_known_model=not bool(vals),
                       models=[dict(index=m['index'],marks=m['marks'],physical_map=mp) for m,mp in matches])
        else:
            rec.update(status='UNKNOWN_POSITIVE_TARGET_AMOUNT_OR_FULL_MODEL',malformed_fields=bad); unknown.append(rec)
        hits.append(rec)
    legacy = 0
    for stem,keys in W.finite_inventory().items():
        doc = W.load(stem)
        for key in keys:
            for where,row,raw,white,index in W.physical_rows(doc[key],key,require_components=True):
                one(stem,where,raw,white,row,index); legacy += 1
        del doc
    assert legacy == 564
    literal = 0
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
        assert counts[stem] > before,('missing newer section',stem,keys)
        newer.append(dict(source=stem,section=list(keys),corners=counts[stem]-before)); del doc,section
    chosen = []
    for stem,expected in zip(W.CHOSEN,(280,52)):
        doc = W.load(stem); rows = [r for r in outputs(doc) if F(r['payment']) > 0]
        assert len(rows) == expected
        for row in rows:
            assert 1 <= len(row['roles']) <= 2
            types = [chosen_type(r) for r in row['roles']]
            assert sorted(types) not in [['A1','A16'],['A1','D16']],('foreign unmodeled chosen cross target',stem)
        chosen.append(dict(source=stem,positive_outputs=expected,target_outputs=0)); del doc,rows
    extensions = []
    for stem in W.EXTENSIONS:
        doc = W.load(stem); before = counts[stem]
        for where,row,raw,white,index in W.physical_rows(doc): one(stem,where,raw,white,row,index)
        extensions.append(dict(source=stem,corners=counts[stem]-before)); del doc
    if unknown: raise AssertionError(json.dumps(dict(status='UNKNOWN_CROSS_TARGET',count=len(unknown),first=unknown[:8])))
    positive = [h for h in hits if h['status'] == 'retained_full_physical_coupon_model']
    assert {(h['index'],tuple(h['marks'])) for h in positive if h['source'] == W.OLD17} == FORK_CORNERS
    assert {h['index'] for h in positive if h['source'] == COUPLED and h['role']['kind'] == 'A16_endpoint'} == {675,828}
    sources = {r['path']:r['sha256'] for r in plan['files']}|PINS
    sources[str(own.relative_to(ROOT))] = own_before
    for name,expected in sources.items():
        assert W.digest(ROOT/name) == expected,('source changed during audit',name)
    return json.loads(json.dumps(dict(status='PASS_SAVED_A1_RANK16_PHYSICAL_COUPON_INTERFACE',sources=sources,
        input_files=54,input_bytes=sum(r['bytes'] for r in plan['files']),historical_packets=22,
        historical_corners=legacy,literal_configurations=literal,newer_sections=newer,extensions=extensions,
        chosen_D12=chosen,corner_counts=dict(counts),target_hits=hits,nonfinite_displays=dict(nonfinite),
        target_summary=dict(Counter(h['role']['kind']+'/'+h['status'] for h in hits)),
        explicit_partial_zero_nonfork_displays=partial,unknown_targets=0,scalar_identities=scalar_identities(),
        inherited_written_scope='Universal complete route/output inventory, actual-wall feasibility and full arbitrary-spectator seal; newer monotone/type-disjoint output chronology. Saved-model audit does not prove these premises.',
        normalization='Positive root rescaling to norm2; all19 full-core edges checked simple, zeros preserved, marks/common whites preserved, only orthogonal common-white extras permitted.',
        owns_lock=False,child_deadline_seconds=45)))


def main():
    p = argparse.ArgumentParser(); m = p.add_mutually_exclusive_group(required=True)
    m.add_argument('--write',action='store_true'); m.add_argument('--verify',action='store_true'); args = p.parse_args()
    signal.alarm(45)
    if args.write: assert not CERT.exists(), 'refuse overwrite'
    result = run(); encoded = json.dumps(result,sort_keys=True,separators=(',',':'))+'\n'
    assert len(encoded.encode()) < 1024**2, 'one-MiB retained certificate cap'
    if args.write:
        with CERT.open('x') as stream: stream.write(encoded)
    else:
        with CERT.open() as stream: assert result == json.load(stream), 'saved replay differs'
    print(json.dumps(dict(status=result['status'],mode='write' if args.write else 'verify',
        targets=len(result['target_hits']),summary=result['target_summary'],unknown=0)))


if __name__ == '__main__': main()
