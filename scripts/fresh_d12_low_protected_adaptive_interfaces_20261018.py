#!/usr/bin/env python3
"""D12 long/short half: low66 and protected74/adaptive9 interfaces only.

Inherit arbitrary-support closure, freshly identify every changed selected
finite output and its replacement. No universal rule or dimension claim.
"""
import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import sympy as s

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_d12_low_protected_adaptive_interfaces_20261018'
LOW = 'fresh_distance8_low_core_20261017'
HIGH = 'fresh_distance8_high_overlap_20261017'
ADAPT = 'fresh_distance8_adaptive_leaf_overlap_20261017'
HALF = 'fresh_d12_remaining_selected_half_replacements_20261018'


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def key(row):
    return json.dumps([row[k] for k in ('low', 'high', 'match')])


def path(g, ids, marks):
    parents = {marks[0]: None}
    todo = [marks[0]]
    for a in todo:
        for b in ids:
            if b != a and g[a, b] != 0 and b not in parents:
                parents[b] = a
                todo.append(b)
    assert marks[1] in parents
    out = [marks[1]]
    while out[-1] != marks[0]:
        out.append(parents[out[-1]])
    return list(reversed(out))


def shape(g, ids, marks):
    adj = {i: [j for j in ids if i != j and g[i, j] != 0] for i in ids}
    assert set(marks) <= set(ids)
    pp = path(g, ids, marks)
    simple = all(g[i, i] == 2 for i in ids) and all(
        g[i, j] == -1 for i in ids for j in adj[i])
    degrees = sorted(map(len, adj.values()))
    d12 = False
    arms = []
    if simple and len(ids) == 12 and degrees == [1]*3 + [2]*8 + [3]:
        branch = next(i for i in ids if len(adj[i]) == 3)
        for start in adj[branch]:
            previous, current, arm = branch, start, [start]
            while len(adj[current]) == 2:
                nxt = next(j for j in adj[current] if j != previous)
                previous, current = current, nxt
                assert current not in arm and current != branch
                arm.append(current)
            assert len(adj[current]) == 1
            arms.append(arm)
        assert len(set().union(*map(set, arms))) == 11
        d12 = sorted(map(len, arms)) == [1, 1, 9]
    changed = d12 and all(len(adj[i]) == 1 for i in marks) and len(pp) == 11
    return dict(order=len(ids), simple=simple, degrees=degrees,
                marked_degrees=[len(adj[i]) for i in marks], path=pp,
                fork_arms=arms, D12=d12, changed_long_short_pair=changed)


def a_endpoint(g, ids, mark):
    adj = {i: [j for j in ids if i != j and g[i, j] != 0] for i in ids}
    assert all(g[i, i] == 2 for i in ids)
    assert all(g[i, j] == -1 for i in ids for j in adj[i])
    assert sorted(map(len, adj.values())) == [1, 1] + [2]*(len(ids)-2)
    assert len(path(g, ids, [mark, next(i for i in ids if len(adj[i]) == 1 and i != mark)])) == len(ids)
    assert len(adj[mark]) == 1


def run():
    docs = {stem: json.loads((ROOT / f'results/{stem}.json').read_text())
            for stem in (LOW, HIGH, ADAPT, HALF)}
    for doc in docs.values():
        for p, h in doc['sources'].items():
            assert digest(p) == h, p
    low = docs[LOW]
    geom = {key(row): row for row in low['extension_geometry']['cases']}
    assert len(geom) == len(low['extension_closure']['cases']) == 66
    rows, affected = [], []
    for index, row in enumerate(low['extension_closure']['cases']):
        geo = geom[key(row)]
        corner = next(c for c in geo['corners'] if c['marks'] == row['corner'])
        core = low['extension_geometry']['cores'][row['core_key']]
        g = s.Matrix(core['gram'])
        ids = corner['component']
        assert row['component_order_bound'] <= 12
        variants = [shape(g, ids, row['corner'])]
        if row['method'] == 'unique_sealed_label3_root':
            n, q = g.rows, row['contact_node']
            h = g.row_join(s.Matrix([-int(i == q) for i in range(n)]))
            h = h.col_join(s.Matrix([[-int(i == q) for i in range(n)] + [2]]))
            variants.append(shape(h, ids+[n], row['corner']))
        changed = any(v['changed_long_short_pair'] for v in variants)
        rows.append(dict(index=index, inherited=row, variants=variants, changed=changed))
        if changed:
            affected.append(index)
    assert affected == [23, 25, 59, 61]

    # These four inherited closure mechanisms seal the entire displayed core,
    # not merely the originally chosen finite factor. Extra orthogonal factors
    # cannot connect to a displayed factor without a first nonzero core contact.
    repairs = []
    for index in affected:
        row = low['extension_closure']['cases'][index]
        geo = geom[key(row)]
        corner = next(c for c in geo['corners'] if c['marks'] == row['corner'])
        core = low['extension_geometry']['cores'][row['core_key']]
        g = s.Matrix(core['gram'])
        if core['radical']:
            z = s.Matrix(core['radical'][0])
            assert g*z == s.zeros(g.rows, 1)
            assert corner['isolated'] and corner['zero_nodes'] == []
            assert corner['radical_shortcut'] <= 8
            assert all(z[i] != 0 for i in corner['possible_contacts'])
            sealing = 'radical annihilation or independent short8'
        else:
            J = s.Matrix(core['inverse'])
            assert corner['candidates'] == [] and corner['end_partitions']
            assert all(J[i, i] > 2 for i in corner['possible_contacts'])
            sealing = 'all allowed single-contact projections exceed2 in positive end blocks'
        other = next(c for c in geo['corners'] if c['marks'] == [row['corner'][0], 8])
        assert other['marks'] != row['corner']
        if index in (23, 59):
            ids = other['component']
            assert len(ids) == 12 and len(other['blocks']) == 1
            for mark in other['marks']:
                a_endpoint(g, ids, mark)
            payment = dict(kind='retained whole A12 endpoint pair', value='1', nodes=ids)
        else:
            one = next(b['nodes'] for b in other['blocks'] if 8 in b['nodes'])
            eleven = next(b['nodes'] for b in other['blocks'] if row['corner'][0] in b['nodes'])
            assert one == [8] and len(eleven) == 11 and len(other['blocks']) == 2
            a_endpoint(g, eleven, row['corner'][0])
            assert all(g[i, 8] == 0 for i in eleven)
            payment = dict(kind='new whole A1/A11 endpoint half plus selected D12 half',
                           value='1/2+1/2', A1=one, A11=eleven)
        repairs.append(dict(index=index, sealing=sealing,
                            replacement_corner=other['marks'], payment=payment))

    high = docs[HIGH]['isolation']
    protected, adaptive = [], []
    for row in high['cases']:
        if row['method'] not in ('finite_leaf9_or_ideal_feature5', 'nine_context_adaptive_finite_payment'):
            continue
        g = s.Matrix(high['cores'][row['core_key']]['gram'])
        marks = row['preferred_corner'] if row['method'] == 'finite_leaf9_or_ideal_feature5' else row['payment']['corner']
        corner = next(c for c in row['all_corners'] if c['preferred_corner'] == marks)
        pp = path(g, corner['marked_component'], marks)
        assert len(pp) == 10
        out = dict(index=row['index'], marks=marks, retained_path=pp, distance=9)
        (protected if row['method'] == 'finite_leaf9_or_ideal_feature5' else adaptive).append(out)
    assert len(protected) == 74 and len(adaptive) == 9
    paths = {f'scripts/{TAG}.py'} | {f'{folder}/{stem}.{ext}'
        for stem in docs for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))}
    return dict(status='PASS_D12_LOW66_PROTECTED74_ADAPTIVE9_INTERFACES',
        sources={p: digest(p) for p in sorted(paths)}, low=rows, repairs=repairs,
        protected=protected, adaptive=adaptive,
        summary=dict(low_contexts=66, unchanged_low=62, repaired_low=4,
            repair_kinds=dict(Counter(r['payment']['kind'] for r in repairs)),
            protected_paths=74, adaptive_paths=9),
        scope='No-short8 inherited low/protected/adaptive routes for changing whole D12 long/short unit to half. The universal ideal/short-path interfaces and global composition/count remain separate. No new dimension bound.')


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    path_out = ROOT / f'results/{TAG}.json'
    assert args.verify or not path_out.exists()
    out = run()
    if args.verify:
        assert out == json.loads(path_out.read_text())
    else:
        path_out.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['summary'], flush=True)
