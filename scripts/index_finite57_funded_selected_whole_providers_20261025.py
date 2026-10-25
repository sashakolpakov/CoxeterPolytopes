#!/usr/bin/env python3
"""SOURCE ONLY: current-funded whole-model provider keys for the22 finite packets.

Physical model keys and complete native maps are reconstructed. The accepted
source-specific full-variant/exhaustiveness theorem remains an explicit premise.
This is a provider index, not the final selected-completion coverage check.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
ROUTES = 'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTES_SHA = '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
COLLECTOR = 'scripts/check_finite57_terminal_collectors_direct_receivers_20261025.py'
COLLECTOR_SHA = 'a2a9d61a7637e45d5a4d7af1c3283e120ca7685d07e9f686ae7b6ee585520eca'
SHARDS = 8
MAX_OUTPUT = 500 * 1024


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded whole-provider source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded whole-provider source read')
    return data


def module(relative, digest, name):
    data = raw(relative)
    need(sha256(data).hexdigest() == digest, 'whole-provider exact source identity')
    obj = types.ModuleType(name)
    obj.__file__ = str(ROOT / relative)
    exec(compile(data, str(ROOT / relative), 'exec'), obj.__dict__)
    return obj


def prepare():
    r = module(ROUTES, ROUTES_SHA, 'finite57_whole_provider_routes')
    c = module(COLLECTOR, COLLECTOR_SHA, 'finite57_whole_provider_direct_collector')
    c.parent_lock()
    manifest = c.raw(c.MANIFEST)
    need(sha256(manifest).hexdigest() == c.MANIFEST_SHA, 'whole-provider manifest identity')
    doc = json.loads(manifest, object_pairs_hook=c.unique)
    c.PINS = {row['path']: row['sha256'] for row in doc['inputs']}
    need(len(c.PINS) == len(doc['inputs']), 'unique whole-provider input paths')
    c.W = c.module(c.LIB, 'finite57_whole_provider_literal')
    c.E = c.module(c.ENGINE, 'finite57_whole_provider_engine')
    c.J = c.module(c.PROJECTOR, 'finite57_whole_provider_projection')
    c.J.json = types.SimpleNamespace(loads=lambda value: json.loads(
        value, object_pairs_hook=c.unique,
        parse_constant=lambda value: (_ for _ in ()).throw(ValueError('nonfinite JSON number'))))
    primary = 'manuscripts/finite_volume57.tex'
    need(c.digest(primary) == c.PINS[primary] == c.W.TEX_SHA256, 'literal/primary source binding')
    c.READ.add(primary)
    return r, c


def original_white(case, order):
    black = case['left'][1] + [case['right_copy_map'][i] for i in case['right'][1]]
    return sorted(set(range(order)) - set(black))


def completion_key(case, typ, images):
    base = case['marked_component']
    return dict(high_index=case['index'], core_key=case['core_key'], base_nodes=base,
                selected_marks=case['preferred_corner'], target_kind=typ[0],
                target_order=typ[1], affine=False, induced_map=images,
                marked_image=sorted(images[base.index(m)] for m in case['preferred_corner']))


def record(c, r, stem, section, ordinal, proof, case, core):
    base = case['marked_component']
    if 'core_key' in proof:
        need(proof['core_key'] == case['core_key'], 'whole proof authoritative core key')
    models = list(c.complete_models(proof, f'{section}/{ordinal}'))
    if not models:
        return None
    original = c.E.matrix(core['gram'])
    n = len(original)
    white = original_white(case, n)
    funded = {}
    providers = {}
    for path, encoded, current_white, corners, gp, wp in models:
        physical = c.E.matrix(encoded)
        need(len(physical) >= n and all(physical[i][j] == original[i][j]
                                       for i in range(n) for j in range(n)),
             'whole model retains original physical core with original numbering')
        need(sorted(i for i in current_white if i < n) == white,
             'whole model original white-set intersection')
        full_receipt = c.receipt(stem, path, encoded, current_white, corners, dict(gram=gp, white=wp))
        need(full_receipt['direct_current_receiver_condition']
             and not full_receipt['legacy_pointwise_bounds_used'], 'DIRECT current receiver conclusion required')
        model_key = sha256(c.canonical([encoded, current_white, [q['marks'] for q in corners]]).encode()).hexdigest()
        if model_key not in funded:
            funded[model_key] = dict(model_key=model_key, source_paths=[],
                context_sha256=full_receipt['context_sha256'], current_total=full_receipt['total'],
                current_corners=[dict(marks=q['marks'], kind=q['kind'], lower=q['weight'])
                                 for q in full_receipt['corners']])
        funded[model_key]['source_paths'].append(path)
        if not any(set(q['marks']) == set(case['preferred_corner']) for q in corners):
            continue
        ids = sorted(set(current_white) | set(case['preferred_corner']))
        parts = c.E.components(physical, ids, c.FUEL)
        selected = next((p for p in parts if set(case['preferred_corner']) <= set(p)), None)
        need(selected is not None and set(base) <= set(selected), 'whole model preferred factor contains original selected base')
        if not 10 <= len(selected) <= 17:
            continue
        typ, marked_images, affine = c.E.classify(physical, selected, case['preferred_corner'], c.W, c.FUEL)
        if affine:
            continue  # This registry makes no implicit selected-affine theorem claim.
        need(typ[0] in ('A', 'B', 'D'), 'only finite classical factors can contain this long selected base')
        target = c.W.classical_graph(*typ)
        labels = c.E.labels(physical, selected)
        for full_map in r.embeddings(labels, target, c.FUEL, c.E):
            restricted = [full_map[selected.index(node)] for node in base]
            key = completion_key(case, typ, restricted)
            name = c.canonical(key)
            entry = providers.setdefault(name, dict(key=key, funded_model_keys=[]))
            if model_key not in entry['funded_model_keys']:
                entry['funded_model_keys'].append(model_key)
    return dict(source=stem, source_section=section, source_row=ordinal,
                high_index=case['index'], core_key=case['core_key'],
                all_current_funded_source_models=[funded[k] for k in sorted(funded)],
                full_completion_keys=[providers[k] for k in sorted(providers)],
                geometric_applicability_premise='For this authoritative high row and matching whole finite completion, the accepted selected theorem of this exact source proof exhausts its listed variants and proves its marked/witness seals, modulo independent current S8.',
                automatic_exhaustiveness_claim=False)


def run(shard, c, r):
    high, cores, partition = r.domain(c)
    packets = []
    for packet, stem in enumerate(c.LATER):
        if packet % SHARDS != shard:
            continue
        rows = []
        source_domain = []
        no_models = []
        for section in c.LATER[stem]:
            section_rows = c.field(stem, (section,))
            need(type(section_rows) is list, 'exact later selected source section')
            for ordinal, proof in enumerate(section_rows):
                index = proof['index']
                need(type(index) is int and 0 <= index < len(high), 'selected source authoritative high index')
                case = high[index]
                need(case.get('marked_component') is not None, 'selected proof uses a bounded component')
                core = cores[case['core_key']]
                source_domain.append([section, ordinal, index])
                got = record(c, r, stem, section, ordinal, proof, case, core)
                if got is None:
                    no_models.append([section, ordinal, index])
                else:
                    rows.append(got)
        packets.append(dict(source=stem, complete_source_row_keys=source_domain,
                            current_funded_whole_model_providers=rows,
                            rows_without_whole_model=no_models))
    return dict(shard=shard, shards=SHARDS, complete_packet_keys=list(c.LATER), packets=packets,
                scope='The22 later finite packets only. Missing whole keys, original/newer/affine providers, S8/triangle/protected routes and final G6 union remain explicit.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(f'whole-{i}' for i in range(SHARDS)), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    r, c = prepare()
    out = run(int(args.batch[-1]), c, r)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'whole-provider consumed input changed')
    need(raw(relative) == before and sha256(raw(ROUTES)).hexdigest() == ROUTES_SHA,
         'whole-provider/index source changed')
    need(sha256(raw(COLLECTOR)).hexdigest() == COLLECTOR_SHA, 'direct collector source changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'whole-provider manifest changed')
    payload = json.dumps(dict(status='CURRENT_FUNDED_WHOLE_PROVIDER_INDEX_ONLY', batch=args.batch,
                              source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
                              consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)},
                              fuel_used=c.FUEL[0], result=out), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'whole-provider output cap; refusal not truncated coverage')
    print(payload)


if __name__ == '__main__':
    main()
