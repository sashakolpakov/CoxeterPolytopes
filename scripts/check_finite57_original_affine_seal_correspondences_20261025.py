#!/usr/bin/env python3
"""SOURCE ONLY: exact old/new D16/D17 affine exception and current receipt join.

A visible whole affine corner gives a positive kernel. Its two omitted-black
contacts then seal the entire old star, or give an independently paid path2.
The code checks this witness and full current receipt at every registered
correspondence. It is not the complete G6/global finite57 check.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
LOADER = 'scripts/index_finite57_funded_selected_whole_providers_20261025.py'
LOADER_SHA = '0a245deec21b283f0d47af0323cf9557fd6ba4465f6e25a2b7c47897f2a46e86'
PINS = 'results/finite57_selected_affine_completion_chain_source_pins_20261024.json'
PINS_SHA = '5f94aff40b8178158e8421c0cbfe1eb377469cd20fdc333aafe6e94e7fe611b5'
SMALL = 'fresh_distance8_ideal5_small_factor_20261017'
MAX_OUTPUT = 500 * 1024


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded original-affine source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded original-affine source read')
    return data


def prepare():
    data = raw(LOADER)
    need(sha256(data).hexdigest() == LOADER_SHA, 'original-affine loader pin')
    obj = types.ModuleType('finite57_original_affine_loader')
    obj.__file__ = str(ROOT / LOADER)
    exec(compile(data, str(ROOT / LOADER), 'exec'), obj.__dict__)
    r, c = obj.prepare()
    data = raw(PINS)
    need(sha256(data).hexdigest() == PINS_SHA, 'selected-affine source register pin')
    rows = json.loads(data, object_pairs_hook=c.unique)['sources']
    need(len(rows) == 68 and len({p['path'] for p in rows}) == 68, 'unique accepted source-register paths')
    for p in rows:
        if p['path'] in c.PINS:
            need(c.PINS[p['path']] == p['sha256'], 'original-affine pin collision')
        c.PINS[p['path']] = p['sha256']
    script = f'scripts/{SMALL}.py'
    need(c.digest(script) == c.PINS[script], 'old/new map producer source pin')
    c.READ.add(script)
    return obj, r, c


def zero_product(c, g, nodes, coefficients):
    for i in nodes:
        value = c.E.Radical.scalar(0)
        for j, q in zip(nodes, coefficients, strict=True):
            value = value + g[i][j] * c.E.Radical.scalar(q)
        need(not value.terms, 'exact positive-kernel/radical annihilator identity')


def seal(c, r, gram, white, corners, core, order, old_record):
    g = c.E.matrix(gram)
    n = len(g)
    need(n == order + 2 and all(g[i][i].rational() == 2 for i in range(n)),
         'original sealed star has norm-two roots and two omitted blacks')
    target = c.W.classical_graph('D', order, affine=True)
    found = None
    for corner in corners:
        marks = corner['marks']
        ids = sorted(set(white) | set(marks))
        parts = c.E.components(g, ids, c.FUEL)
        if parts != [ids] or len(ids) != order:
            continue
        labels = c.E.labels(g, ids)
        witness = next(r.embeddings(labels, target, c.FUEL, c.E), None)
        if witness is not None:
            found = ids, marks, witness
            break
    need(found is not None, 'an actual old corner must contain the full visible affineD factor')
    ids, marks, mapping = found
    kernel = [F(1) if mapping[i] <= order - 5 else F(1, 2) for i in range(order)]
    need(all(q > 0 for q in kernel), 'positive actual affine kernel')
    zero_product(c, g, ids, kernel)
    omitted = sorted(set(range(n)) - set(ids))
    need(len(omitted) == 2 and not set(omitted) & set(white), 'exact two omitted old black walls')
    proof = dict(affine_nodes=ids, actual_affine_marks=marks, full_native_isomorphism_witness=mapping,
                 positive_kernel=list(map(str, kernel)), omitted_black=omitted,
                 common_escape='Contact with both omitted blacks gives the independent current length-two whole-face unit.')
    if order == 16:
        need(core['signature'] == [17, 1, 0] and core['inverse'] is not None,
             'nondegenerate original Lorentz-core signature/inverse premise')
        diagonal = [F(core['inverse'][i][i]) for i in omitted]
        need(all(x > 2 for x in diagonal), 'every sole omitted-black contact projects above norm2')
        proof.update(method='positive_affine_kernel_and_strict_single_contact_projection',
                     omitted_inverse_diagonals=list(map(str, diagonal)),
                     retained_computed_premise='The stored original inverse is the inverse of this physical Gram.')
    else:
        need(order == 17 and core['signature'] == [17, 1, 1], 'original radical Lorentz-core signature premise')
        radical = list(map(F, old_record['radical']))
        need(len(radical) == n and any(radical), 'actual old radical-vector schema')
        zero_product(c, g, list(range(n)), radical)
        need(radical[omitted[0]] * radical[omitted[1]] < 0, 'opposite nonzero omitted-black radical coefficients')
        proof.update(method='positive_affine_kernel_and_actual_zero_radical',
                     actual_zero_radical=list(map(str, radical)),
                     omitted_radical_coefficients=[str(radical[i]) for i in omitted],
                     retained_computed_premise='The stored signature is correct; an actual Lorentz span containing a negative direction has no nonzero null vector in its orthogonal complement.')
    proof['conclusion'] = 'Without the independent path2 escape every additional white is orthogonal to the ENTIRE old star; all its displayed marked factors are whole.'
    return proof


def run(order, c, r):
    high, cores, partition = r.domain(c)
    doc = c.document(SMALL, ('records',), ('patterns',))
    selected = [p for p in doc['records'] if (p['kind'], p['order']) == ('D', order)]
    need(len(selected) == 1, 'unique exact old/new affine target record')
    rec = selected[0]
    target = c.W.classical_graph('D', order, affine=True)
    need(len(rec['target']) == order and all(
        type(row) is list and len(row) == order and all(type(x) is int and x == target[i].get(j, 0)
        for j, x in enumerate(row)) for i, row in enumerate(rec['target'])), 'entire old/new canonical target graph')
    if order == 16:
        stem = c.OLD16
        old = c.field(stem, ('exceptional_D15',))
        original_gram = old['gram']
        originals = {p['case']: p for p in old['cases']}
        expected_old = old['affected_cases']
    else:
        stem = c.OLD17
        old = c.field(stem, ('order17',))
        original_gram = cores['(25, 25, ((9, 9),))']['gram']
        originals = {p['case']: p for p in old['replacement_cases']}
        expected_old = old['dropped_case_indices']
    need(len(originals) == len(expected_old) == 9 and set(originals) == set(expected_old), 'exact nine original physical exception keys')
    correspondence = rec['old_sealed_case_correspondence']
    need(len(correspondence) == 9 and len({p['new_case'] for p in correspondence}) == 9
         and {p['old_case'] for p in correspondence} == set(originals), 'exact bijective nine-case old/new correspondence')
    wanted = {p['new_case'] for p in correspondence}
    need(set(rec['dropped_affected_cases']) == wanted and len(rec['dropped_affected_cases']) == 9,
         'exact current dropped-case support domain')
    saved = {p['case']: p for p in rec['case_embeddings']}
    need(len(saved) == len(rec['case_embeddings']), 'unique current saved map-row keys')
    patterns = {p['index']: p for p in doc['patterns']}
    need(len(patterns) == len(doc['patterns']) == 1153, 'complete pinned1153 pattern domain')
    dropped = {tuple(p) for p in rec['dropped_pairs']}
    need(len(dropped) == 9, 'exact original fork-triple orbit')
    source_affected = set()
    for row in rec['case_embeddings']:
        pattern = patterns[row['case']]
        positions = pattern['marks']
        c.E.indices(positions, len(pattern['nodes']), 2)
        for images in row['embeddings']:
            c.E.indices(images, order, len(pattern['nodes']))
        pairs = sorted({tuple(sorted(images[i] for i in positions)) for images in row['embeddings']})
        need(pairs == sorted(tuple(q) for q in row['marked_images']), 'entire source saved-map/marked-image association')
        if set(pairs) & dropped:
            source_affected.add(row['case'])
    need(source_affected == wanted, 'entire source affected-map domain equals the old/new correspondence keys')
    result = []
    for link in correspondence:
        index, old_index = link['new_case'], link['old_case']
        case = high[index]
        core = cores[case['core_key']]
        base, marks = case['marked_component'], case['preferred_corner']
        need(all(link[name] == case[name] for name in ('left', 'right', 'match')), 'entire marked-core correspondence descriptor')
        need(case['component_order_bound'] >= order and len(base) <= order, 'actual current target eligibility')
        need(c.E.matrix(original_gram) == c.E.matrix(core['gram']), 'old/current full physical core identity')
        p = originals[old_index]
        c.E.indices(p['black'], len(original_gram), 4)
        white = sorted(set(range(len(original_gram))) - set(p['black']))
        current_white = sorted(set(range(len(original_gram))) - set(case['left'][1])
                               - {case['right_copy_map'][i] for i in case['right'][1]})
        need(white == current_white, 'original/current actual white support')
        need({tuple(sorted(q['marks'])) for q in p['corners']} ==
             {tuple(sorted(q['preferred_corner'])) for q in case['all_corners']}, 'complete actual corner-domain join')
        labels = c.E.labels(c.E.matrix(core['gram']), base)
        pattern = patterns[index]
        need(pattern['nodes'] == base and pattern['labels'] == labels
             and pattern['marks'] == [base.index(m) for m in marks], 'current marked-pattern native/base/mark identity')
        maps = list(r.embeddings(labels, target, c.FUEL, c.E))
        need(sorted(maps) == sorted(saved[index]['embeddings']), 'ALL full current native maps equal the saved source list')
        need(len({tuple(mp) for mp in maps}) == len(maps), 'no repeated full map')
        pairs = sorted({tuple(sorted(mp[base.index(m)] for m in marks)) for mp in maps})
        need(pairs == sorted(tuple(q) for q in saved[index]['marked_images']), 'complete saved marked-image association')
        maps = [mp for mp in maps if tuple(sorted(mp[base.index(m)] for m in marks)) in dropped]
        need(maps, 'registered exception must have actual dropped full-map provider keys')
        sealing = seal(c, r, original_gram, white, p['corners'], core, order, p)
        receipt = c.receipt(stem, old_index, original_gram, white, p['corners'], 'old affine positive-kernel full-star seal')
        need(receipt['direct_current_receiver_condition'] and not receipt['legacy_pointwise_bounds_used'], 'direct CURRENT receiver required')
        keys = [dict(high_index=index, core_key=case['core_key'], base_nodes=base, selected_marks=marks,
                     target_kind='D', target_order=order, affine=True, induced_map=mp,
                     marked_image=sorted(mp[base.index(m)] for m in marks)) for mp in maps]
        result.append(dict(new_case=index, old_source=stem, old_case=old_index, complete_dropped_map_keys=keys,
                           sealing_witness=sealing, direct_current_receipt=receipt))
    return dict(order=order, complete_correspondence_keys=correspondence, current_exception_providers=result,
                scope='Exactly these nine original physical stars, modulo independent current path2. Other selected/universal/provider domains and the final G6 union are separate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=('old-d16', 'old-d17'), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    helper, r, c = prepare()
    out = run(16 if args.batch == 'old-d16' else 17, c, r)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'original-affine consumed source/input changed')
    need(raw(relative) == before and sha256(raw(LOADER)).hexdigest() == LOADER_SHA, 'original-affine source/loader changed')
    need(sha256(raw(PINS)).hexdigest() == PINS_SHA, 'selected-affine source register changed')
    need(sha256(helper.raw(helper.COLLECTOR)).hexdigest() == helper.COLLECTOR_SHA, 'direct collector changed')
    need(sha256(helper.raw(helper.ROUTES)).hexdigest() == helper.ROUTES_SHA, 'native route definitions changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'base110 input register changed')
    payload = json.dumps(dict(status='ORIGINAL_AFFINE_CURRENT_SEAL_BINDING_ONLY', batch=args.batch,
        source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
        consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0], result=out), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'original-affine output cap; refusal not truncated maps')
    print(payload)


if __name__ == '__main__':
    main()
