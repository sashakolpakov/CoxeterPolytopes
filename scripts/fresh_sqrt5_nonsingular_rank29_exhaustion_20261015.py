#!/usr/bin/env python3
"""Atomic, dimension-specific necessary attachment closure of311 finite roots.

All singular stopping choices are exported, irrespective of prior weights.
The31-node limit is the full quadrilateral-star size for ambient dimension29.
No all-rank conclusion, actual realization, or weight payment is asserted.
"""
import argparse
import json
import sqlite3
import time
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
from pathlib import Path
import fresh_sqrt5_singular_extension_frontier_20261015 as F
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_nonsingular_rank29_exhaustion_20261015'
Q, T = F.Q, F.T
DB = ROOT / f'data/enumerations/{TAG}_running.sqlite'
OUT = ROOT / f'results/{TAG}.json'
MAX_NODES = 31


def load(tag):
    return json.loads((ROOT / f'results/{tag}.json').read_text())


def manifest():
    names = [F.TAG, J.TAG, 'fresh_sqrt5_sealed_subcore_closure_20261015',
             'fresh_sqrt5_quad_exhaustion_20261015',
             'fresh_sqrt5_quad_contacts_20261015_review',
             'fresh_sqrt5_quad_contacts_20261015',
             'fresh_sqrt5_sealed_weights_20261015']
    paths = [f'scripts/{TAG}.py']
    for name in names:
        paths.extend([f'scripts/{name}.py', f'results/{name}.json'])
    paths.append('results/arithmetic_dimension_bounds_and_effective_high_dimension_search.md')
    Q.manifest_check(paths)
    return {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def initial_rows():
    base = load('fresh_sqrt5_quad_exhaustion_20261015')['labelled_union_supports']
    rows = [r for r in load('fresh_sqrt5_sealed_subcore_closure_20261015')['records']
            if r['status'] == 'OPEN_NONSINGULAR']
    assert len(rows) == 311
    return [(r, F.extension(base[r['source_union_index']], r['contacts'])) for r in rows]


def corner_patterns(row):
    g = Q.edges(row)
    answer = []
    for a, b in product((0, 2), (1, 3)):
        H = T.graph(row['white'] + [a, b], g)
        remaining = set(H)
        components = []
        while remaining:
            ids = T.component(H, min(remaining))
            remaining.difference_update(ids)
            components.append({i: {j: m for j, m in H[i].items() if j in ids} for i in ids})
        marked = next(h for h in components if a in h)
        connected = b in marked
        answer.append(dict(marked=[a, b], forest=[T.typename(h) for h in components],
                           a1_count=sum(len(h) == 1 for h in components),
                           marked_component=T.typename(marked) if connected else None,
                           pair_orbit=T.keystr(T.tree_key(marked, (a, b))) if connected else None))
    return answer


def evaluate(row):
    n = row['nodes']
    assert n <= MAX_NODES and row['white'] == list(range(4, n))
    assert F.corners(row) and not F.short(row) and not F.OLD.h4_payments(row)
    M = Q.tree_scaled(n, Q.edges(row))
    assert M is not None
    physical = Q.inertia(M)
    conjugate = Q.inertia([[x.conj() for x in r] for r in M])
    assert physical[1] == 1 and physical[2] in (0, 1)
    assert conjugate[1] == 0 and conjugate[2] == physical[2]
    singular = bool(physical[2])
    out = dict(singular=singular, physical=physical, conjugate=conjugate)
    if singular:
        out['zero_touch_corner_patterns'] = corner_patterns(row)
        out['orthogonal_spectator_rank_budget'] = {str(d): d + 2 - n for d in (28, 29) if d + 2 >= n}
    if n == MAX_NODES:
        return dict(**out, status='DIMENSION29_SINGULAR_WEIGHT_OBLIGATION' if singular
                    else 'REJECT_RANK_EXCEEDS_DIMENSION29', domain_counts={}, trials=[], accepted=[])
    domain, dc = F.profiles(row)
    trials, accepted = [], []
    for contact in domain:
        ext = F.extension(row, contact)
        r = dict(contacts=contact)
        if (w := F.OLD.h4_payments(ext)):
            r.update(status='PAID_ENTIRE_H4', witness=w)
        elif (G := Q.tree_scaled(n + 1, Q.edges(ext))) is None:
            r['status'] = 'REJECT_FINITE_CYCLE_FIELD'
        else:
            # Use the gauge of the enlarged matrix for the inherited relation.
            ker = Q.kernel([z[:-1] for z in G[:-1]]) if singular else []
            products = [sum((x * G[i][-1] for i, x in enumerate(v)), Q.K()) for v in ker]
            if any(products):
                r.update(status='REJECT_INHERITED_RELATION', products=[x.enc() for x in products])
            else:
                p = Q.inertia(G)
                c = Q.inertia([[x.conj() for x in z] for z in G])
                r.update(physical=p, conjugate=c)
                if p[1] > 1:
                    r['status'] = 'REJECT_TWO_NEGATIVE'
                elif c[1]:
                    r['status'] = 'REJECT_CONJUGATE'
                elif p[2] > 1:
                    r['status'] = 'REJECT_MULTIPLE_LORENTZ_RELATIONS'
                else:
                    assert p[1] == 1 and p[2] == c[2]
                    if singular:
                        assert p[2] == 1
                    k = Q.kernel(G) if p[2] else []
                    if k and (all(x.sign() >= 0 for x in k[0]) or all(x.sign() <= 0 for x in k[0])):
                        r.update(status='REJECT_NONNEGATIVE_ROOT_DEPENDENCE', witness=[x.enc() for x in k[0]])
                    else:
                        r['status'] = 'ACCEPT_NECESSARY_UNPAID'
                        accepted.append(contact)
        trials.append(r)
    if accepted:
        status = 'EXPANDED_SINGULAR' if singular else 'EXPANDED_NONSINGULAR'
    else:
        status = 'SEALED_SINGULAR_WEIGHT_OBLIGATION' if singular else 'REJECT_SEALED_NONSINGULAR_RANK'
    return dict(**out, status=status, domain_counts=dc, trials=trials, accepted=accepted)


class Registry:
    def __init__(self, con):
        self.con = con
        self.buckets = defaultdict(list)
        self.next_id = 0
        for i, encoded in con.execute('SELECT id,row FROM cores ORDER BY id'):
            assert i == self.next_id
            M = J.matrix(json.loads(encoded))
            self.buckets[J.invariant(M)].append((i, M))
            self.next_id += 1

    def add(self, row):
        M = J.matrix(row)
        key = J.invariant(M)
        for i, target in self.buckets[key]:
            mp = J.isomorphism(M, target)
            if mp is not None:
                J.check_map(M, target, mp)
                return i, mp
        i = self.next_id
        self.next_id += 1
        self.con.execute('INSERT INTO cores VALUES(?,?,NULL)', (i, json.dumps(row, separators=(',', ':'))))
        self.buckets[key].append((i, M))
        return i, list(range(row['nodes']))


def export(con, state):
    assert state['sources'] == manifest()
    records, counts, trial_counts, obligations = [], Counter(), Counter(), []
    for i, row, payload in con.execute('SELECT id,row,result FROM cores ORDER BY id'):
        rec = dict(index=i, row=json.loads(row), result=json.loads(payload) if payload else None)
        records.append(rec)
        if rec['result']:
            z = rec['result']
            counts[z['status']] += 1
            trial_counts.update(t['status'] for t in z['trials'])
            if z['singular']:
                obligations.append(i)
        else:
            counts['PENDING'] += 1
    assert [r['index'] for r in records] == list(range(len(records)))
    assert all((r['result'] is not None) == (r['index'] < state['next_core']) for r in records)
    complete = state['next_core'] == len(records)
    out = dict(status='COMPLETE_FINITE_SEED_NECESSARY_ATTACHMENT_QUEUE_THROUGH_DIMENSION29' if complete
               else 'PARTIAL_RESUMABLE_FINITE_SEED_ATTACHMENT_QUEUE_THROUGH_DIMENSION29',
               sources=state['sources'], state=state, dimension_bound=29, maximum_star_nodes=MAX_NODES,
               database=str(DB.relative_to(ROOT)), database_sha256=sha256(DB.read_bytes()).hexdigest(),
               marked_states=len(records), complete_queue=complete, counts=dict(counts), trial_counts=dict(trial_counts),
               root_maps=[json.loads(x[0]) for x in con.execute('SELECT payload FROM roots ORDER BY literal_index')],
               weight_obligations=obligations, records=records,
               scope='Necessary finite-label attachment closure from311 literal nonsingular first extensions, for ambient dimensions28/29 only. All proper cores are tested with complete one-white finite-corner/no-short domains and exact Qsqrt5 field, Lorentz and conjugate conditions. No prior paid-subcore pruning is used. Every singular state, even an expanded one or a state previously paid by old weights, is an explicit zero-further-contact weight obligation, with its full corner forests and remaining orthogonal-spectator rank budget. At31nodes no further white fits ambient dimension29, and nonsingular rank31 exceeds ambient rank30. Whole-core transactions and all literal child maps are retained; batch limits only pause. Short/H4 payments are inherited universally and must be respected by any final rule. We do not claim all affine subdiagrams excluded, actual quadrilateral realization, all-rank sealing, payment by a new rule, coverage of dotted or other seed families, or a new global dimension bound.')
    OUT.write_text(json.dumps(out, indent=2) + '\n')
    return out


def verify(sources):
    saved = json.loads(OUT.read_text())
    assert saved['sources'] == sources
    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    assert con.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    assert sha256(DB.read_bytes()).hexdigest() == saved['database_sha256']
    state = json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])
    assert state == saved['state'] and state['sources'] == sources
    roots = [json.loads(z[0]) for z in con.execute('SELECT payload FROM roots ORDER BY literal_index')]
    assert roots == saved['root_maps'] and len(roots) == 311
    reached = set()
    for i, ((raw, row), mp) in enumerate(zip(initial_rows(), roots)):
        assert mp['literal_index'] == i
        assert all(mp[k] == raw[k] for k in ('source_union_index', 'source_class_index', 'trial_index', 'contacts'))
        J.check_map(J.matrix(row), J.matrix(saved['records'][mp['target']]['row']), mp['map'])
        reached.add(mp['target'])
    assert reached == set(range(state['initial_marked_roots']))
    counts, trials, obligations = Counter(), Counter(), []
    for rec in saved['records']:
        i = rec['index']
        assert i in reached
        stored = con.execute('SELECT row,result FROM cores WHERE id=?', (i,)).fetchone()
        assert json.loads(stored[0]) == rec['row']
        assert (rec['result'] is not None) == (i < state['next_core'])
        if rec['result']:
            z = evaluate(rec['row'])
            old = rec['result']
            children = old['children']
            assert z == {k: v for k, v in old.items() if k != 'children'}
            assert len(children) == len(z['accepted'])
            for c, child in zip(z['accepted'], children):
                assert child['contacts'] == c
                J.check_map(J.matrix(F.extension(rec['row'], c)), J.matrix(saved['records'][child['index']]['row']), child['map'])
                reached.add(child['index'])
            assert json.loads(stored[1]) == old
            counts[z['status']] += 1
            trials.update(t['status'] for t in z['trials'])
            if z['singular']:
                obligations.append(i)
        else:
            assert stored[1] is None
            counts['PENDING'] += 1
        if i and i % 100 == 0:
            print('verified', i, flush=True)
    assert reached == set(range(len(saved['records'])))
    assert dict(counts) == saved['counts'] and dict(trials) == saved['trial_counts']
    assert obligations == saved['weight_obligations']
    assert saved['complete_queue'] == (state['next_core'] == len(saved['records']))
    con.close()
    print('PASS_SAVED_RANK29_ATTACHMENT_BATCH', state['next_core'], saved['marked_states'], len(obligations))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--max-cores', type=int, default=100)
    ap.add_argument('--max-seconds', type=float, default=600)
    args = ap.parse_args()
    sources = manifest()
    if args.verify:
        verify(sources)
        return
    DB.parent.mkdir(parents=True, exist_ok=True)
    if not args.resume:
        assert not DB.exists(), 'Use --resume; never overwrite an existing frontier.'
    con = sqlite3.connect(DB)
    con.execute('PRAGMA journal_mode=DELETE')
    con.execute('PRAGMA synchronous=FULL')
    if not args.resume:
        con.execute('BEGIN IMMEDIATE')
        con.execute('CREATE TABLE cores(id INTEGER PRIMARY KEY,row TEXT,result TEXT)')
        con.execute('CREATE TABLE roots(literal_index INTEGER PRIMARY KEY,payload TEXT)')
        con.execute('CREATE TABLE meta(id INTEGER PRIMARY KEY,value TEXT)')
        registry = Registry(con)
        for i, (raw, row) in enumerate(initial_rows()):
            j, mp = registry.add(row)
            root = dict(literal_index=i, target=j, map=mp, **{k: raw[k] for k in ('source_union_index', 'source_class_index', 'trial_index', 'contacts')})
            con.execute('INSERT INTO roots VALUES(?,?)', (i, json.dumps(root)))
        state = dict(sources=sources, next_core=0, initial_literal_roots=311,
                     initial_marked_roots=registry.next_id, batches=[])
        assert state['initial_marked_roots'] == 199
        con.execute('INSERT INTO meta VALUES(1,?)', (json.dumps(state),))
        con.commit()
    else:
        state = json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])
        assert state['sources'] == sources
        registry = Registry(con)
    start, began = state['next_core'], time.monotonic()
    stop = 'MAX_CORES'
    for k in range(args.max_cores):
        if k and time.monotonic() - began >= args.max_seconds:
            stop = 'TIME_GUARD_AFTER_COMPLETE_CORE'
            break
        i = state['next_core']
        record = con.execute('SELECT row FROM cores WHERE id=?', (i,)).fetchone()
        if record is None:
            stop = 'QUEUE_EXHAUSTED'
            break
        row = json.loads(record[0])
        z = evaluate(row)
        con.execute('BEGIN IMMEDIATE')
        children = []
        for contact in z['accepted']:
            j, mp = registry.add(F.extension(row, contact))
            children.append(dict(index=j, contacts=contact, map=mp))
        z['children'] = children
        con.execute('UPDATE cores SET result=? WHERE id=?', (json.dumps(z, separators=(',', ':')), i))
        state['next_core'] = i + 1
        con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state),))
        con.commit()
        print('core', i, 'nodes', row['nodes'], z['status'], 'children', len(children), 'states', registry.next_id, flush=True)
    state['batches'].append(dict(start=start, end=state['next_core'], stop_reason=stop))
    con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state),))
    con.commit()
    out = export(con, state)
    con.close()
    print(out['status'], state['next_core'], out['marked_states'], out['counts'], 'seconds', round(time.monotonic() - began, 3))


if __name__ == '__main__':
    main()
