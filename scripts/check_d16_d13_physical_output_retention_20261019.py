#!/usr/bin/env python3
"""Current saved-interface audit; no historical census or forest optimization.

Coordinator only: run_bounded_research_check_20261018.py --parent-lock --
THIS --write/--verify [--mode d16|d13|both]. No inner flock; child alarm45s.
Whole-component/actual-star routing and the two reports remain written inputs.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import json
from pathlib import Path
import signal
import check_a15_inner_unit_a1_absence_retention_20261018 as W

ROOT = W.ROOT
TAG = 'd16_d13_physical_output_retention_20261019'
PLAN = 'results/d16_d13_physical_output_audit_inputs_20261019.json'
D16 = 'fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018'
DISTAL = 'fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018'
PINS = {
    PLAN: '6b6af762943a9279f747fcceca67de75522d6864201e067aba22630b619cab78',
    'scripts/check_a15_inner_unit_a1_absence_retention_20261018.py':
    '328c473c3fa69e547050793485f51cf399dcd667770d642431de25da5f50c908',
    'results/finite_quad_nonsingular_full_seal_rank_obstruction_20261019.md':
    '9d16f3bcdcefe85d3138aae6967a484907f3978badf02eef6f5fde3cc878a2c4',
    'results/order17_distal_d13_physical_coupon_capacity_research_20261019.md':
    'b5e699b35fd434df989b2b0a616f58a7d62753a5879c5858b226ce460f8858a6',
}


def simple_edges(g, nodes):
    return all(g[i][i] > 0 and all(g[i][j] < 0 and
        4*g[i][j]*g[i][j] == g[i][i]*g[j][j] for j in nodes if i != j and g[i][j]) for i in nodes)


def d_order(g, nodes):
    adj = {i: [j for j in nodes if i != j and g[i][j]] for i in nodes}
    centers = [i for i in nodes if len(adj[i]) == 3]
    if len(centers) != 1 or any(len(adj[i]) not in (1, 2, 3) for i in nodes):
        return None
    center = centers[0]
    leaves = sorted(i for i in nodes if adj[i] == [center])
    end = [i for i in nodes if len(adj[i]) == 1 and i not in leaves]
    if len(leaves) != 2 or len(end) != 1:
        return None
    order = [end[0]]
    while order[-1] != center:
        nxt = [j for j in adj[order[-1]] if j not in order]
        if len(nxt) != 1:
            return None
        order.append(nxt[0])
    return (order, leaves) if len(order)+2 == len(nodes) and simple_edges(g, nodes) else None


def a_order(g, nodes):
    adj = {i: [j for j in nodes if i != j and g[i][j]] for i in nodes}
    ends = sorted(i for i in nodes if len(adj[i]) == 1)
    if len(ends) != 2 or any(len(adj[i]) not in (1, 2) for i in nodes) or not simple_edges(g, nodes):
        return None
    order = [ends[0]]
    while len(order) < len(nodes):
        nxt = [j for j in adj[order[-1]] if j not in order]
        if len(nxt) != 1:
            return None
        order.append(nxt[0])
    return order


def roles(raw, white, marks, ids=None):
    g = W.mat(raw)
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    blocks = W.parts(g, ids)
    result = {}
    common = next((b for b in blocks if set(marks) <= set(b)), None)
    if common is not None and len(common) == 16:
        order = d_order(g, common)
        if order and order[0][0] in marks and len(set(marks) & set(order[1])) == 1:
            result['d16'] = dict(nodes=common, spine=order[0], short_leaves=order[1],
                                 displayed_A1=[b[0] for b in blocks if len(b) == 1])
    if common is None:
        first, second = [next(b for b in blocks if m in b) for m in marks]
        for a, d, am, dm in ((first, second, *marks), (second, first, marks[1], marks[0])):
            if len(a) != 4 or len(d) != 13:
                continue
            ao, do = a_order(g, a), d_order(g, d)
            if ao and do and am in (ao[0], ao[-1]) and dm in do[1]:
                result['d13'] = dict(A4=a, D13=d, A4_endpoint=am, D13_leaf=dm,
                                     A4_order=ao, D13_spine=do[0], short_leaves=do[1])
    if result:
        result['corner_nodes'] = ids
        result['finite'] = W.pd(g, ids)
    return result


def inverse(raw):
    n = len(raw)
    a = [[F(x) for x in row]+[F(int(i == j)) for j in range(n)] for i, row in enumerate(raw)]
    for i in range(n):
        k = next(j for j in range(i, n) if a[j][i])
        a[i], a[k] = a[k], a[i]
        p = a[i][i]; a[i] = [x/p for x in a[i]]
        for j in range(n):
            if j != i:
                p = a[j][i]; a[j] = [x-p*y for x, y in zip(a[j], a[i])]
    return [row[n:] for row in a]


def d16_models():
    models = []
    doc = W.load(D16)['local_proof']
    assert (doc['complete_D16_patterns'], doc['validated_D16_maps']) == (394, 578)
    assert doc['selected_D16_contexts'] == [1644, 1664, 1931, 2324]
    for p in doc['whole_proofs']:
        old = p['old_complete_seal']; g = W.mat(old['gram']); white = old['white']
        seal = old['sealing']; z = list(map(F, seal['radical']))
        assert seal['method'] == 'radical' and len(g) == 19
        assert all(sum(g[i][j]*z[j] for j in range(19)) == 0 for i in range(19))
        assert all(z[i] for i in seal['allowed_contacts'])
        assert seal['mixed_sign_shortcuts'] and max(r['length'] for r in seal['mixed_sign_shortcuts']) <= 8
        targets = [(r, roles(old['gram'], white, r['marks'])) for r in old['corners']]
        targets = [(r, out) for r, out in targets if 'd16' in out]
        assert len(targets) == 1
        r, out = targets[0]
        assert out['finite'] and len(out['d16']['displayed_A1']) == 1
        models.append(dict(index=p['index'], gram=g, white=sorted(white), marks=sorted(r['marks']),
                           corner_nodes=out['corner_nodes'], A1=out['d16']['displayed_A1']))
    del doc
    old = W.load(W.OLD16)['exceptional_D15']; g = W.mat(old['gram']); J = inverse(g)
    assert len(g) == 18 and old['signature'] == [17, 1, 0]
    assert all(sum(g[i][k]*J[k][j] for k in range(18)) == int(i == j)
               for i in range(18) for j in range(18))
    assert [[J[i][j] for j in (3, 13)] for i in (3, 13)] == [[4, -4], [-4, 4]]
    ns = old['affine_nodes']; kernel = list(map(F, old['affine_positive_kernel']))
    assert ns == sorted(set(range(18))-{3, 13}) and all(v > 0 for v in kernel)
    assert all(sum(g[i][j]*kernel[k] for k, j in enumerate(ns)) == 0 for i in ns)
    excluded = []
    for case in old['cases']:
        white = sorted(set(range(18))-set(case['black']))
        assert len(white) == 14 and len(set(case['black'])) == 4 and {3, 13} <= set(case['black'])
        assert W.pd(g, white)
        assert any(set(white+c['marks']) == set(ns) and not c['finite'] for c in case['corners'])
        for r in case['corners']:
            out = roles(old['gram'], white, r['marks'])
            if 'd16' in out:
                assert out['finite'] and not out['d16']['displayed_A1']
                excluded.append(dict(index=case['case'], gram=g, white=white, marks=sorted(r['marks']),
                                     corner_nodes=out['corner_nodes']))
    assert [(p['index'], p['marks']) for p in excluded] == [(663, [1, 13]), (663, [3, 11]),
                                                          (676, [3, 12]), (829, [2, 13])]
    return models, excluded


def d13_models():
    # Rational realization: 13 D coordinates +5 A coordinates +time,
    # with metric time contribution -t^2/5. No square-root package is used.
    dot = lambda u, v: sum(a*b for a, b in zip(u[:-1], v[:-1]))-u[-1]*v[-1]/5
    ds, aa = [], []
    for i in range(12):
        v = [F()]*19; v[i]=1; v[i+1]=-1; ds.append(v)
    v = [F()]*19; v[11]=v[12]=1; ds.append(v)
    for i in range(4):
        v = [F()]*19; v[13+i]=1; v[14+i]=-1; aa.append(v)
    c = [-F(1, 2)]*13+[F()]*5+[F(5, 2)]
    d = [-F(1)]+[F()]*12+[-F(2, 5)]*3+[F(3, 5)]*2+[F(1)]
    assert dot(c, c) == dot(d, d) == 2 and dot(c, d) == 0
    mapping = {
        1974: ([15,16,17,18,9,8,7,6,5,4,0,1,2], [13,12,10,11], 3, 14),
        2367: ([15,16,17,18,9,8,7,6,5,4,0,1,2], [13,12,10,11], 3, 14),
        2504: ([5,6,7,8,9,18,17,16,15,14,10,11,12], [3,2,0,1], 13, 4),
        2511: ([5,6,7,8,9,18,17,16,15,14,10,11,12], [3,2,0,1], 13, 4),
    }
    models = []
    doc = W.load(W.PARTIAL)
    for p in doc['cases']:
        if p['index'] not in mapping:
            continue
        g = W.mat(p['gram']); dm, am, ci, di = mapping[p['index']]
        vectors = [None]*19
        for i, j in enumerate(dm): vectors[j] = ds[i]
        for i, j in enumerate(am): vectors[j] = aa[i]
        vectors[ci], vectors[di] = c, d
        assert all(v is not None for v in vectors)
        assert all(dot(vectors[i], vectors[j]) == g[i][j] for i in range(19) for j in range(19))
        z = list(map(F, p['radical']))
        assert all(sum(g[i][j]*z[j] for j in range(19)) == 0 for i in range(19)) and z[3]*z[13] < 0
        ns, kernel = p['affine_nodes'], list(map(F, p['positive_affine_kernel']))
        assert ns == sorted(set(range(19))-{3, 13}) and all(v > 0 for v in kernel)
        assert all(sum(g[i][j]*kernel[k] for k, j in enumerate(ns)) == 0 for i in ns)
        white = p['white']; assert W.pd(g, white) and len(white) == 15
        assert all(W.pd(g, white+[b]) for b in sorted(set(range(19))-set(white)))
        targets = [(r, roles(p['gram'], white, r['marks'])) for r in p['corners']]
        targets = [(r, out) for r, out in targets if 'd13' in out]
        assert len(targets) == 1
        r, out = targets[0]; assert out['finite']
        selected = [r['marks'] for r in p['corners'] if r['kind'] == 'affineD16_order17']
        assert len(selected) == 1
        models.append(dict(index=p['index'], gram=g, white=sorted(white), marks=sorted(r['marks']),
                           corner_nodes=out['corner_nodes'], selected_ideal_marks=sorted(selected[0])))
    assert [p['index'] for p in models] == list(mapping)
    del doc
    stage = W.load(DISTAL)['local']
    assert (stage['complete_patterns'], stage['successful_contexts'], stage['complete_maps']) == (140, 71, 488)
    assert sorted(stage['changed_pairs']) == [[1,15],[1,16],[11,13],[11,14]]
    assert stage['new_weight'] == '7/32' and stage['finite_coupon_new'] == '25/32'
    assert [p['index'] for p in stage['proofs']] == list(mapping)
    for p, model in zip(stage['proofs'], models):
        assert W.mat(p['gram']) == model['gram'] and sorted(p['white']) == model['white']
        assert sorted(p['selected_marks']) == model['selected_ideal_marks'] and p['selected_payment'] == '7/32'
        assert sorted(p['finite_coupon']['marks']) == model['marks'] and p['finite_coupon']['new'] == '25/32'
    # Subsequent ideal changes are order11 only; the unchanged order17 entry
    # is an explicit written chronology premise, not inferred from hashes.
    # Written actual-facet projection supplies the geometric exclusions.
    assert -F(5, 4) == -F(5, 2)*F(5, 2)/5
    assert F(11, 4)-F(5, 4) == F(3, 2)
    assert 1-F(1, 5) == F(4, 5) and 1+F(4, 5)-F(1, 5) == F(8, 5)
    assert F(1, 4) < F(25, 64) < F(1, 2)
    return models


AMOUNT_KEYS = ('payment', 'weight', 'lower', 'finite_half', 'new',
               'new_weight', 'new_payment', 'current_weight', 'new_lower')


def amounts(row):
    out = {}
    for key in AMOUNT_KEYS:
        v = row.get(key)
        if isinstance(v, (str, int)) and not isinstance(v, bool):
            try: out[key] = F(v)
            except (ValueError, ZeroDivisionError): pass
    return out


def outputs(x):
    if isinstance(x, dict):
        if 'payment' in x and isinstance(x.get('roles'), list): yield x
        for v in x.values(): yield from outputs(v)
    elif isinstance(x, list):
        for v in x: yield from outputs(v)


def chosen_type(role):
    assert isinstance(role, dict) and {'type', 'nodes', 'mark'} <= set(role), 'chosen-role schema'
    typ = role['type']
    assert isinstance(typ, str) and typ, 'missing canonical chosen type'
    if typ == 'E6':
        n = 6
    else:
        assert typ[0] in 'ABD' and typ[1:].isdigit(), ('unknown chosen type', typ)
        n = int(typ[1:])
        assert typ == typ[0]+str(n) and n >= {'A':1, 'B':2, 'D':4}[typ[0]], ('invalid chosen rank', typ)
    nodes = role['nodes']
    assert isinstance(nodes, list) and len(nodes) == n and len(set(nodes)) == n
    assert all(isinstance(i, int) and not isinstance(i, bool) and i >= 0 for i in nodes)
    assert isinstance(role['mark'], int) and not isinstance(role['mark'], bool) and role['mark'] in nodes
    return typ


def run(mode):
    for name, expected in PINS.items(): assert W.digest(ROOT/name) == expected, name
    plan = json.loads((ROOT/PLAN).read_text())
    files = W.input_paths()
    assert [str(p.relative_to(ROOT)) for p in files] == [p['path'] for p in plan['files']]
    for row in plan['files']:
        p = ROOT/row['path']
        assert p.stat().st_size == row['bytes'] and W.digest(p) == row['sha256'], row['path']
    m16, excluded = d16_models() if mode != 'd13' else ([], [])
    m13 = d13_models() if mode != 'd16' else []
    # These sources predate the marginal-alpha increment, or only change its
    # disjoint inner-A4 roles. Missing per-angle amounts use this WRITTEN premise.
    half_only = set(W.finite_inventory()) | {W.OLD16, W.OLD17, W.D13, W.D14, W.PARTIAL, W.D16PART, W.NEWD}
    half_only |= {s for s, _ in W.NEWER if s != DISTAL and 'order11_fractional' not in s}
    half_only |= set(W.EXTENSIONS)
    counts, hits, unknown, nonfinite = Counter(), [], [], Counter()
    def matches(raw, white, marks, ids, models):
        return [p['index'] for p in models if W.mat(raw) == p['gram'] and white is not None
                and sorted(white) == p['white'] and sorted(marks) == p['marks'] and ids == p['corner_nodes']]
    def one(stem, where, raw, white, row, index=None):
        counts[stem] += 1
        out = roles(raw, white, row['marks'], row.get('corner_nodes'))
        for kind in ('d16', 'd13'):
            if kind not in out or mode not in (kind, 'both'): continue
            if not out['finite']:
                nonfinite[kind] += 1; continue
            rec = dict(kind=kind, source=stem, path=where, index=index, marks=row['marks'],
                       corner_nodes=out['corner_nodes'], role=out[kind])
            if kind == 'd16':
                good = matches(raw, white, row['marks'], out['corner_nodes'], m16)
                old = matches(raw, white, row['marks'], out['corner_nodes'], excluded)
                if len(good) == 1 and out['d16']['displayed_A1'] == next(p['A1'] for p in m16 if p['index'] == good[0]):
                    rec.update(status='retained_guarded_unit', native_case=good[0])
                elif len(old) == 1 and not out['d16']['displayed_A1']:
                    rec.update(status='full_order16_core_independent_short2', native_case=old[0])
                else: rec['status'] = 'UNKNOWN_D16_WHOLE_A1_OR_ACTUAL_SHORT2_SEAL'; unknown.append(rec)
            else:
                good = matches(raw, white, row['marks'], out['corner_nodes'], m13)
                requested = amounts(row)
                rec['declared_per_corner_amounts'] = {k: str(v) for k, v in requested.items()}
                if len(good) == 1 and all(v <= 1 for v in requested.values()):
                    # Replace this complete face proof by the parametric identity,
                    # not by its historical numeric request at alpha=9/32.
                    rec.update(status='complete_distal_marginal_geometry', native_case=good[0],
                               symbolic_payment='1/2+alpha; 0<=alpha<=1/2')
                elif all(v <= F(1, 2) for v in requested.values()) and (requested or stem in half_only):
                    rec.update(status='unchanged_half_only', missing_amount_uses_written_source_scope=not bool(requested))
                else: rec['status'] = 'UNKNOWN_D13_REQUEST_ABOVE_HALF_OR_UNTYPED_SOURCE'; unknown.append(rec)
            hits.append(rec)
    legacy = 0
    for stem, keys in W.finite_inventory().items():
        doc = W.load(stem)
        for key in keys:
            for where, row, raw, white, index in W.physical_rows(doc[key], key, require_components=True):
                one(stem, where, raw, white, row, index); legacy += 1
        del doc
    assert legacy == 564
    common = W.native(1142)[1]['gram']; literal = 0
    for stem, index, raw, white, cs in W.literal_rows(common):
        for j, row in enumerate(cs): one(stem, f'literal/{index}/{j}', raw, white, row, index)
        literal += 1
    assert literal == 216
    sections = []
    for stem, keys in W.NEWER:
        doc = W.load(stem); section = doc
        for key in keys: section = section[key]
        before = counts[stem]
        for where, row, raw, white, index in W.physical_rows(section, '/'.join(keys)):
            one(stem, where, raw, white, row, index)
        assert counts[stem] > before, ('missing physical collection', stem, keys)
        sections.append(dict(source=stem, section=list(keys), corners=counts[stem]-before))
        del doc, section
    chosen = []
    for stem, expected in zip(W.CHOSEN, (280, 52)):
        doc = W.load(stem); rs = [r for r in outputs(doc) if F(r['payment']) > 0]
        assert len(rs) == expected
        # Their saved roles specify the actual positive terminal components.
        for row in rs:
            assert 1 <= len(row['roles']) <= 2, 'empty/unknown positive chosen role list'
            types = [chosen_type(r) for r in row['roles']]
            if mode != 'd13': assert 'D16' not in types, ('unmapped chosen-D12 D16', stem)
            if mode != 'd16': assert not ('A4' in types and 'D13' in types), ('unmapped chosen-D12 A4/D13', stem)
        chosen.append(dict(source=stem, positive_outputs=expected, target_outputs=0))
        del doc, rs
    extensions = []
    for stem in W.EXTENSIONS:
        doc = W.load(stem); before = counts[stem]
        for where, row, raw, white, index in W.physical_rows(doc): one(stem, where, raw, white, row, index)
        extensions.append(dict(source=stem, physical_corners=counts[stem]-before)); del doc
    if unknown:
        raise AssertionError(json.dumps(dict(status='UNRESOLVED_PHYSICAL_TARGET_OUTPUT', mode=mode,
                                             unknown_count=len(unknown), unknown=unknown[:12])))
    if mode != 'd13':
        assert {h['native_case'] for h in hits if h['kind'] == 'd16' and h['status'] == 'retained_guarded_unit'} == {1644,1664,1931,2324}
        assert sum(h['kind'] == 'd16' and h['source'] == W.OLD16 and h['path'].startswith('literal/') for h in hits) == 4
    if mode != 'd16':
        assert {h['native_case'] for h in hits if h['kind'] == 'd13' and h['status'] == 'complete_distal_marginal_geometry'} == {1974,2367,2504,2511}
    sources = {r['path']: r['sha256'] for r in plan['files']} | PINS
    sources[str(Path(__file__).relative_to(ROOT))] = W.digest(Path(__file__))
    return dict(status='PASS_CURRENT_D16_D13_SAVED_PHYSICAL_INTERFACE', mode=mode, sources=sources,
        input_files=len(files), input_bytes=sum(p.stat().st_size for p in files),
        historical_packets=22, historical_corners=legacy, literal_configurations=literal,
        newer_sections=sections, extensions=extensions, chosen_D12=chosen,
        physical_corner_counts=dict(counts), target_hits=hits, nonfinite_target_displays=dict(nonfinite),
        target_summary=dict(Counter((h['kind']+'/'+h['status']) for h in hits)), unknown_targets=0,
        d16_candidate_delta='-d16*[a1=0]', d13_candidate_extra_upper='2*alpha*[a4*d13>0]',
        d13_parameter_domain='0<=alpha<=1/2', d13_ideal_selected_weight='1/2-alpha',
        d13_ideal_D_total='113-4alpha', d13_rank16_envelope='(113-4alpha)/17',
        whole_half_only_source_registry=sorted(half_only), new_positive_families=0,
        universal_geometry_formalized=False, complete_new_finite_maximum=None,
        scope='Saved physical-output/normalization/model audit only. Whole-factor routing, full contact '
              'domains, rank obstruction, actual facet projection, two-angle capacity and old-half '
              'source semantics are explicit written premises; no census or finite forest optimization.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mode', choices=('d16', 'd13', 'both'), default='both')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true'); group.add_argument('--verify', action='store_true')
    args = parser.parse_args(); signal.alarm(45)
    suffix = '' if args.mode == 'both' else '_'+args.mode
    target = ROOT/'results'/(TAG+suffix+'.json')
    result = json.loads(json.dumps(run(args.mode)))
    encoded = json.dumps(result, indent=2, sort_keys=True)+'\n'
    assert len(encoded.encode()) < 1024**2
    if args.verify: assert json.loads(target.read_text()) == result, 'Saved certificate differs'
    else:
        with target.open('x') as stream: stream.write(encoded)
    signal.alarm(0)
    print(json.dumps({k: result[k] for k in ('status','mode','input_files','input_bytes','target_summary')}))


if __name__ == '__main__':
    main()
