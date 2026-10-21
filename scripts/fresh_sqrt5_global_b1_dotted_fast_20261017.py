#!/usr/bin/env python3
"""Exact B1 dotted closure using the conjugate radical lemma.

Safe payment pruning is conditional on each adjacent discount being <=1/2.
Singular states are both stopping obligations and possible growth states.
"""
import argparse
import json
import sqlite3
import time
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path

import fresh_sqrt5_global_b1_dotted_parameters_20261016 as P
import fresh_sqrt5_all_quad_contacts_20261015 as C
import fresh_sqrt5_adjacent_discount_frontier_20261016 as W
import fresh_sqrt5_nonsingular_rank29_fast_20261015 as A
import fresh_sqrt5_sealed_subcore_closure_20261015_review as J

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_global_b1_dotted_fast_20261017'
INPUT = f'results/{P.TAG}.json'
DB = ROOT / f'data/enumerations/{TAG}.sqlite'
OUT = ROOT / f'results/{TAG}.json'
D, Q, K = P.D, P.Q, P.K


def edge_key(i, j):
    return ','.join(map(str, sorted((i, j))))


def field_matrix(row, parameters):
    """Tree-product gauge for arbitrary finite and fixed dotted edges."""
    n = row['nodes']
    edges = Q.edges(row)
    squares = {3: K(1), 4: K(2), 5: Q.PHI * Q.PHI, 6: K(3), 10: 2 + Q.PHI}
    adj = {i: [] for i in range(n)}
    values = {}
    assert set(parameters) == {edge_key(i, j) for (i, j), m in edges.items() if m == 0}
    for (i, j), m in edges.items():
        if m == 2:
            continue
        z = K(*parameters[edge_key(i, j)]) if m == 0 else squares[m]
        assert z.sign() > 0 and z.conj().sign() > 0
        values[i, j] = z
        adj[i].append(j)
        adj[j].append(i)
    norms = {}
    for start in reversed(range(n)):
        if start in norms:
            continue
        norms[start] = K(2)
        stack = [start]
        while stack:
            i = stack.pop()
            for j in sorted(adj[i], reverse=True):
                if j not in norms:
                    norms[j] = norms[i] * values[tuple(sorted((i, j)))]
                    stack.append(j)
    M = [[K() for _ in range(n)] for _ in range(n)]
    for i in range(n):
        M[i][i] = norms[i]
    for (i, j), z in values.items():
        root = D.square_root(z * norms[i] * norms[j])
        if root is None:
            return None
        M[i][j] = M[j][i] = -root / 2
        assert 4 * M[i][j] ** 2 / (norms[i] * norms[j]) == z
    return M


def classify(row, parameters):
    M = field_matrix(row, parameters)
    if M is None:
        return dict(status='REJECT_CYCLE_FIELD')
    tuples = [[(x.a, x.b) for x in r] for r in M]
    p, c, vectors = A.analyze(tuples, row['white'])
    out = dict(physical=p, conjugate=c)
    if p[1] > 1:
        return dict(**out, status='REJECT_TWO_NEGATIVE')
    if c[1]:
        return dict(**out, status='REJECT_CONJUGATE')
    if p[2] > 1:
        return dict(**out, status='REJECT_MULTIPLE_LORENTZ_RELATIONS')
    assert p[1] == 1 and c == [row['nodes'] - p[2], 0, p[2]]
    if row['nodes'] == 23 and not p[2]:
        return dict(**out, status='REJECT_AMBIENT_RANK')
    if vectors:
        assert len(vectors) == 1
        signs = [A.R.sign(x) for x in vectors[0]]
        out.update(kernel=[A.R.enc(x) for x in vectors[0]], kernel_signs=signs)
        if all(x >= 0 for x in signs) or all(x <= 0 for x in signs):
            return dict(**out, status='REJECT_NONNEGATIVE_ROOT_RELATION')
        out['black_kernel_alternates'] = (signs[0] == signs[2] and
                                          signs[1] == signs[3] and
                                          signs[0] * signs[1] == -1)
    cycles = D.cycles(M)
    bad = next((r for r in cycles if not r['integral']), None)
    if bad:
        return dict(**out, status='REJECT_CYCLIC_INTEGRALITY', witness=bad)
    out['integral_cycles'] = len(cycles)
    return dict(**out, status='ACCEPT_SINGULAR' if p[2] else 'ACCEPT_NONSINGULAR')


def evaluate(row, parameters):
    cl = classify(row, parameters)
    assert cl['status'] in ('ACCEPT_SINGULAR', 'ACCEPT_NONSINGULAR'), cl
    analysis = W.analysis(row)
    assert analysis['baseline'] == 1 and analysis['potentially_affected']
    singular = cl['status'] == 'ACCEPT_SINGULAR'
    out = dict(classification=cl, analysis=analysis, singular=singular,
               orthogonal_spectator_rank_budget=23 - row['nodes'])
    if row['nodes'] == 23:
        assert singular
        return dict(**out, status='RANK23_SINGULAR_STOP', trials=[], accepted=[],
                    domain_counts=dict(rank_boundary=True, profiles=0))
    domain, counts = C.domains(row)
    trials, accepted = [], []
    for contacts in domain:
        child = C.extension(row, contacts)
        ca = W.analysis(child)
        W.transition_check(analysis, ca)
        if ca['baseline'] >= 2:
            z = dict(status='HEREDITARY_SAFE_BASELINE_AT_LEAST_TWO', baseline=ca['baseline'])
        elif not ca['potentially_affected']:
            z = dict(status='HEREDITARY_SAFE_NO_POSSIBLE_ADJACENT_DISCOUNT')
        else:
            z = classify(child, parameters)
            if z['status'] in ('ACCEPT_SINGULAR', 'ACCEPT_NONSINGULAR'):
                if singular:
                    assert z['status'] == 'ACCEPT_SINGULAR'
                accepted.append(contacts)
        trials.append(dict(contacts=contacts, **z))
    status = (('EXPANDED_SINGULAR' if singular else 'EXPANDED_NONSINGULAR') if accepted
              else ('SEALED_SINGULAR_STOP' if singular else 'REJECT_SEALED_NONSINGULAR_RANK'))
    return dict(**out, status=status, trials=trials, accepted=accepted,
                domain_counts=dict(**counts, profiles=len(domain)))


def marked_matrix(row, parameters):
    n = row['nodes']
    M = [[(2, '', '') for _ in range(n)] for _ in range(n)]
    for (i, j), m in Q.edges(row).items():
        z = (0, *parameters[edge_key(i, j)]) if m == 0 else (m, '', '')
        M[i][j] = M[j][i] = z
    return M


class Registry:
    def __init__(self, con):
        self.con, self.buckets, self.next_id = con, defaultdict(list), 0
        for i, row, params in con.execute('SELECT id,row,params FROM cores ORDER BY id'):
            assert i == self.next_id
            M = marked_matrix(json.loads(row), json.loads(params))
            self.buckets[J.invariant(M)].append((i, M))
            self.next_id += 1

    def add(self, row, parameters):
        M = marked_matrix(row, parameters)
        key = J.invariant(M)
        for i, T in self.buckets[key]:
            mp = J.isomorphism(M, T)
            if mp is not None:
                J.check_map(M, T, mp)
                return i, mp
        i = self.next_id
        self.next_id += 1
        self.con.execute('INSERT INTO cores VALUES(?,?,?,NULL)',
                         (i, json.dumps(row), json.dumps(parameters, sort_keys=True)))
        self.buckets[key].append((i, M))
        return i, list(range(row['nodes']))


def manifest():
    paths = [f'scripts/{TAG}.py', 'results/fresh_conjugate_radical_lemma_20261017.md',
             'scripts/fresh_sqrt5_global_b1_dotted_migration_20261017.py',
             'results/fresh_sqrt5_global_b1_dotted_closure_20261016.json',
             'scripts/fresh_sqrt5_global_b1_dotted_closure_20261016.py', INPUT, f'scripts/{P.TAG}.py', f'scripts/{C.TAG}.py',
             f'scripts/{W.TAG}.py', f'scripts/{A.TAG}.py', f'scripts/{J.TAG}.py',
             f'scripts/{D.TAG}.py']
    return {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def seeds():
    doc = json.loads((ROOT / INPUT).read_text())
    assert doc['accepted_parameter_seeds'] == 377
    for r in doc['records']:
        for i, a in enumerate(r['result']['accepted']):
            ratios = a.get('dotted_ratios', {'0,2': a.get('ratio')})
            params = {k: (4 * K(*v)).enc() for k, v in ratios.items()}
            yield r, i, params


def export(con, state):
    rows, counts, trials, obligations = [], Counter(), Counter(), []
    for i, row, params, payload in con.execute('SELECT id,row,params,result FROM cores ORDER BY id'):
        z = json.loads(payload) if payload else None
        assert (z is not None) == (i < state['next_core'])
        rows.append(dict(index=i, row=json.loads(row), parameters=json.loads(params), result=z))
        if z:
            counts[z['status']] += 1
            trials.update(t['status'] for t in z['trials'])
            if z['singular']:
                obligations.append(i)
        else:
            counts['PENDING'] += 1
    complete = state['next_core'] == len(rows)
    doc = dict(status='COMPLETE_DOTTED_B1_ATTACHMENT_CLOSURE' if complete
               else 'PARTIAL_ATOMIC_DOTTED_B1_ATTACHMENT_CLOSURE',
               sources=state['sources'], state=state, dimension=21, maximum_support_nodes=23,
               affine_exclusion='Conjugate-PSD radical transport plus rejection of nonnegative full relations; see fresh_conjugate_radical_lemma_20261017.md',
               conditional_weight_hypotheses=dict(distance2_discount='0',
                   each_adjacent_discount='0 <= delta <= 1/2', bonuses='nonnegative'),
               database=str(DB.relative_to(ROOT)), database_sha256=sha256(DB.read_bytes()).hexdigest(),
               complete_queue=complete, marked_states=len(rows), counts=dict(counts), trial_counts=dict(trials),
               root_maps=[json.loads(r[0]) for r in con.execute('SELECT payload FROM roots ORDER BY id')],
               singular_stopping_obligations=obligations, records=rows,
               scope='All377 exact necessary minimal B1 dotted parameter seeds, including both dots and nonsingular starting cores. Dotted values are included in every checked isomorphism. All finite-corner nonorthogonal white profiles are expanded through23nodes, subject only to exact arithmetic/compactness/rank gates and displayed hereditary weight pruning. Every singular state remains a stopping obligation even when expanded. Pending states are explicit. Orthogonal spectators and exact payments are separate; no dimension theorem or global compact realization.')
    OUT.write_text(json.dumps(doc, indent=2) + '\n')
    return doc


def verify():
    doc = json.loads(OUT.read_text())
    assert doc['sources'] == manifest()
    assert sha256(DB.read_bytes()).hexdigest() == doc['database_sha256']
    con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
    assert con.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]) == doc['state']
    stored = list(con.execute('SELECT id,row,params,result FROM cores ORDER BY id'))
    con.close()
    rows = doc['records']
    assert len(rows) == len(stored)
    reached = set()
    roots = list(seeds())
    assert len(roots) == len(doc['root_maps']) == 377
    for (r, a, params), m in zip(roots, doc['root_maps']):
        assert m['source_record'] == r['index'] and m['accepted_index'] == a
        target = rows[m['target']]
        J.check_map(marked_matrix(r['row'], params),
                    marked_matrix(target['row'], target['parameters']), m['map'])
        reached.add(m['target'])
    counts, trials, obligations = Counter(), Counter(), []
    for rec, saved in zip(rows, stored):
        i = rec['index']
        assert i == saved[0] and i in reached
        assert json.loads(saved[1]) == rec['row'] and json.loads(saved[2]) == rec['parameters']
        assert (json.loads(saved[3]) if saved[3] else None) == rec['result']
        if rec['result'] is None:
            counts['PENDING'] += 1
            continue
        old = rec['result']
        new = json.loads(json.dumps(evaluate(rec['row'], rec['parameters'])))
        assert new == {k: v for k, v in old.items() if k != 'children'}
        counts[new['status']] += 1
        trials.update(t['status'] for t in new['trials'])
        if new['singular']:
            obligations.append(i)
        children = {tuple(c['contacts']): c for c in old['children']}
        assert set(children) == set(map(tuple, new['accepted'])) and len(children) == len(old['children'])
        for contact in new['accepted']:
            ch = children[tuple(contact)]
            child, target = C.extension(rec['row'], contact), rows[ch['index']]
            J.check_map(marked_matrix(child, rec['parameters']),
                        marked_matrix(target['row'], target['parameters']), ch['map'])
            reached.add(ch['index'])
        if i % 25 == 0:
            print('FAST_DOTTED_B1_REPLAY', i, flush=True)
    assert reached == set(range(len(rows))) and dict(counts) == doc['counts']
    assert dict(trials) == doc['trial_counts'] and obligations == doc['singular_stopping_obligations']
    print('PASS_SAVED_DOTTED_B1_ATTACHMENT_CLOSURE', doc['state']['next_core'],
          doc['marked_states'], doc['complete_queue'], flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--max-cores', type=int, default=10000)
    ap.add_argument('--max-seconds', type=float, default=300)
    args = ap.parse_args()
    if args.verify:
        verify()
        return
    sources = manifest()
    if not args.resume:
        assert not DB.exists(), 'Use --resume; never overwrite a queue.'
    con = sqlite3.connect(DB)
    con.execute('PRAGMA journal_mode=DELETE')
    con.execute('PRAGMA synchronous=FULL')
    if not args.resume:
        con.execute('BEGIN IMMEDIATE')
        con.execute('CREATE TABLE cores(id INTEGER PRIMARY KEY,row TEXT,params TEXT,result TEXT)')
        con.execute('CREATE TABLE roots(id INTEGER PRIMARY KEY,payload TEXT)')
        con.execute('CREATE TABLE meta(id INTEGER PRIMARY KEY,value TEXT)')
        registry = Registry(con)
        for literal, (r, a, params) in enumerate(seeds()):
            j, mp = registry.add(r['row'], params)
            m = dict(source_record=r['index'], accepted_index=a, target=j, map=mp)
            con.execute('INSERT INTO roots VALUES(?,?)', (literal, json.dumps(m)))
        state = dict(sources=sources, next_core=0, initial_literal_roots=377,
                     initial_marked_roots=registry.next_id, batches=[])
        con.execute('INSERT INTO meta VALUES(1,?)', (json.dumps(state),))
        con.commit()
    else:
        state = json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])
        assert state['sources'] == sources
        registry = Registry(con)
    start, began, reason = state['next_core'], time.monotonic(), 'MAX_CORES'
    for k in range(args.max_cores):
        if k and time.monotonic() - began >= args.max_seconds:
            reason = 'TIME_GUARD_AFTER_COMPLETE_CORE'
            break
        i = state['next_core']
        rec = con.execute('SELECT row,params FROM cores WHERE id=?', (i,)).fetchone()
        if rec is None:
            reason = 'QUEUE_EXHAUSTED'
            break
        row, params = map(json.loads, rec)
        z = evaluate(row, params)
        con.execute('BEGIN IMMEDIATE')
        children = []
        for contacts in z['accepted']:
            j, mp = registry.add(C.extension(row, contacts), params)
            children.append(dict(index=j, contacts=contacts, map=mp))
        z['children'] = children
        con.execute('UPDATE cores SET result=? WHERE id=?', (json.dumps(z), i))
        state['next_core'] = i + 1
        con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state),))
        con.commit()
        if i % 10 == 0:
            print('FAST_DOTTED_B1_CORE', i, row['nodes'], registry.next_id, z['status'], flush=True)
    state['batches'].append(dict(start=start, end=state['next_core'], stop_reason=reason,
                                elapsed_seconds=round(time.monotonic() - began, 3)))
    con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state),))
    con.commit()
    doc = export(con, state)
    con.close()
    print(doc['status'], state['next_core'], doc['marked_states'],
          len(doc['singular_stopping_obligations']), flush=True)


if __name__ == '__main__':
    main()
