#!/usr/bin/env python3
"""SOURCE ONLY: exact proper finite-map exclusion provider index for G6.

The five accepted source families supply geometric exclusion implications at
explicit computed inverse/residual/contact premises. This adapter binds their
FULL map keys and original native contacts; it does not prove final coverage.
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
ROUTES = 'scripts/index_finite57_current_terminal_routes_syntax_safe_20261025.py'
ROUTES_SHA = '43a0c4e986cbdb92d112771bf0edc1d128ce2d0feced30fffd9d9825dda7c29d'
MAX_OUTPUT = 500 * 1024
FAMILIES = (
    ('fresh_finite_a16_d16_d17_sparse_coupon_transfer_20261018', 'proper_completion_cases', 'D', 16),
    ('fresh_finite_rank14_existing_payment_transfer_20261018', 'selected_proofs', 'row-kind', 14),
    ('fresh_finite_b15_d15_existing_payment_transfer_20261018', 'selected_proofs', 'D', 15),
    ('fresh_finite_d16_fork_exchange_20261018', 'proper_completion_proofs', 'D', 16),
    ('fresh_finite_d16_second_pair_zero_transfer_20261018', 'proper_completion_proofs', 'D', 16),
)
OUTCOMES = {
    'projection_exceeds_two', 'two_orthogonal_fork_leaves_cauchy',
    'core_orthogonal_two_leaf_fork', 'arbitrary_contact_projection_exceeds_norm2',
    'two_distinct_short_leaves_would_coincide', 'first_entry_projection_above_norm2',
    'zero_residual_fork_center', 'forced_projection_above_two',
    'two_orthogonal_leaves_cauchy',
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded source read')
    return data


def prepare():
    data = raw(ROUTES)
    need(sha256(data).hexdigest() == ROUTES_SHA, 'route source identity')
    obj = types.ModuleType('finite57_proper_finite_routes')
    obj.__file__ = str(ROOT / ROUTES)
    exec(compile(data, str(ROOT / ROUTES), 'exec'), obj.__dict__)
    return obj, obj.prepare()


def completion_key(high, kind, order, images):
    base = high['marked_component']
    return dict(high_index=high['index'], core_key=high['core_key'], base_nodes=base,
                selected_marks=high['preferred_corner'], target_kind=kind,
                target_order=order, affine=False, induced_map=images,
                marked_image=sorted(images[base.index(m)] for m in high['preferred_corner']))


def actual_entries(target, images, base, missing):
    return sorted((new, old, base[i], target[new][old])
                  for new in missing for i, old in enumerate(images) if old in target[new])


def bind_contacts(c, row, target, images, base, core, method):
    missing = sorted(set(range(len(target))) - set(images))
    name = next((name for name in ('missing_target_nodes', 'missing_nodes', 'missing') if name in row), None)
    need(name is not None and row[name] == missing, 'exact full-map missing target nodes')
    required = actual_entries(target, images, base, missing)
    need(required, 'proper completion must have an entry from its displayed base')
    inverse = core['inverse']
    need(inverse is not None, 'exact original inverse input required')
    if method == 'impossible_proper_completion':
        need(len(required) == 1, 'unique required contact in D15 proper branch')
        new, old, node, label = required[0]
        need((row['forced_old_target_node'], row['core_contact'], row['label']) == (old, node, label),
             'D15 proper forced contact must match its actual induced-map edge')
        need(F(row['projection_lower_bound']) == F(inverse[node][node]) > 2,
             'D15 proper projection bound at actual original contact')
        return required
    if 'forced_entries' in row:
        entries = []
        for e in row['forced_entries']:
            new = e.get('new_target_node', e.get('new_node'))
            old = e.get('old_target_node', e.get('old_node'))
            node = e['core_node']
            need(type(new) is int and type(old) is int and type(node) is int,
                 'typed forced-entry schema')
            need(F(e['inverse_diagonal']) == F(inverse[node][node]),
                 'forced-entry inverse diagonal at actual physical node')
            entries.append((new, old, node, target[new].get(old, 0)))
        need(sorted(entries) == required, 'complete exact forced-entry contact list')
    elif 'required_contacts' in row:
        entries = []
        for new, triples in row['required_contacts'].items():
            need(str(int(new)) == new, 'canonical target node key')
            entries.extend((int(new), old, node, label) for old, node, label in triples)
        need(sorted(entries) == required, 'complete exact rank14 required-contact list')
        node = row['forced_node']
        need(node in [q[2] for q in required], 'rank14 forced node is an actual required contact')
        if row['outcome'] == 'arbitrary_contact_projection_exceeds_norm2':
            need(F(row['projection_lower_bound']) == F(inverse[node][node]) > 2,
                 'rank14 strict projection bound at its actual forced contact')
        else:
            need(F(row['projection_square']) == F(inverse[node][node]) == 2,
                 'rank14 saturated fork contact at its actual inverse entry')
    else:
        raise ValueError('unsupported proper contact provider schema')
    need(all(q[3] == 3 for q in required), 'accepted finite proper providers require simple native entry bonds')
    return required


def run(shard, c, r):
    high, cores, partition = r.domain(c)
    packets = []
    for family, (stem, section, kind, order) in enumerate(FAMILIES):
        if family % 4 != shard:
            continue
        rows = c.field(stem, (section,))
        need(type(rows) is list and rows, 'complete pinned provider section list')
        out = []
        exact_rows = []
        for ordinal, p in enumerate(rows):
            method = p.get('method')
            if section == 'selected_proofs' and method not in ('proper_completion_exclusion', 'impossible_proper_completion'):
                continue
            index = p['index']
            need(type(index) is int and 0 <= index < len(high), 'authoritative high row index')
            case = high[index]
            core = cores[case['core_key']]
            need(p['core_key'] == case['core_key'], 'provider original core key')
            a, b = c.E.matrix(p['gram']), c.E.matrix(core['gram'])
            need(a == b, 'provider original physical Gram equality')
            base = case['marked_component']
            target_kind = p['kind'] if kind == 'row-kind' else kind
            target = c.W.classical_graph(target_kind, order)
            need(len(base) < order, 'proper provider strictly enlarges displayed component')
            labels = c.E.labels(b, base)
            entries = []
            keys = set()
            for j, row in enumerate(p['maps']):
                images = row.get('map', row.get('images'))
                c.E.indices(images, order, len(base))
                need(all(labels[i][k] == target[images[i]].get(images[k], 0)
                         for i in range(len(base)) for k in range(i)),
                     'FULL native induced map, including all edges/nonedges')
                key = completion_key(case, target_kind, order, images)
                if 'marked_images' in row:
                    need(row['marked_images'] == key['marked_image'], 'ordered base/mark image association')
                label = 'required_norm2_root_projects_above_two' if method == 'impossible_proper_completion' else row.get('outcome', row.get('method'))
                need(label in OUTCOMES or label == 'required_norm2_root_projects_above_two',
                     'accepted exact proper source outcome')
                contacts = bind_contacts(c, row, target, images, base, core, method)
                need(all(node in case['possible_contacts'] for new, old, node, bond in contacts),
                     'all forced contacts belong to actual high contact domain')
                encoded = c.canonical(key)
                need(encoded not in keys, 'duplicate full proper map key')
                keys.add(encoded)
                entries.append(dict(key=key, source_path=f'{section}/{ordinal}/maps/{j}',
                                    outcome=label, forced_entries=contacts,
                                    computed_record_sha256=sha256(c.canonical(row).encode()).hexdigest(),
                                    geometric_premise='Accepted source theorem: saved inverse/residual/contact identities imply this exact proper completion is impossible, modulo its explicitly independent current short8 escape.'))
            need(entries, 'nonempty proper map domain')
            out.append(dict(source=stem, high_index=index, maps=entries))
            exact_rows.append([section, ordinal, index])
        need(out, 'named proper family must contain its exact source proofs')
        packets.append(dict(source=stem, source_section=section, complete_source_row_keys=exact_rows,
                            proper_map_providers=out))
    return dict(shard=shard, shards=4, packets=packets,
                scope='Explicit five-family proper finite provider index. Final G6 union/exhaustiveness and all whole/affine/current-short providers remain separate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(f'proper-{i}' for i in range(4)), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    r, c = prepare()
    out = run(int(args.batch[-1]), c, r)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'consumed provider source/input changed')
    need(raw(relative) == before and sha256(raw(ROUTES)).hexdigest() == ROUTES_SHA, 'provider/index source changed')
    need(sha256(c.raw(r.COLLECTOR)).hexdigest() == r.COLLECTOR_SHA, 'collector source changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'input manifest changed')
    payload = json.dumps(dict(status='EXACT_PROPER_FINITE_PROVIDER_INDEX_ONLY', batch=args.batch,
                              source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
                              consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)},
                              fuel_used=c.FUEL[0], result=out), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'provider output cap; refuse rather than truncate keys')
    print(payload)


if __name__ == '__main__':
    main()
