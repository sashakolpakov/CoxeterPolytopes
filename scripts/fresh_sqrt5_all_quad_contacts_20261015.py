#!/usr/bin/env python3
"""Resumable complete one-white necessary domains for ALL retained quad cores.

No no-short test and no weight/payment pruning are used. Completed per-core
records and the next class index are committed atomically; a partial core is
never claimed complete. Finite extensions use a cached spherical-corner Schur
complement over Q(phi). Dotted distances remain symbolic.
"""
import argparse
import json
import sqlite3
import time
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import fresh_sqrt5_all_quad_unions_20261015 as U
import fresh_sqrt5_quad_contacts_20261015 as OLD
import fresh_sqrt5_quad_filters_20261015 as A

ROOT = Path(__file__).resolve().parents[1]
BASE = 'fresh_sqrt5_all_quad_unions_20261015'
TAG = 'fresh_sqrt5_all_quad_contacts_20261015'
DB = ROOT / 'data/enumerations' / (TAG + '_running.sqlite')
REPORT = ROOT / 'results' / (TAG + '.json')
ALPHABET = (2, 3, 4, 5, 6, 10)


@lru_cache(maxsize=50000)
def finite_word(n, word):
    return U.E.finite_word.__wrapped__(n, word)


def edges(row):
    return dict(zip(combinations(range(row['nodes']), 2), row['labels']))


def corner_ok(n, g):
    return all(finite_word(len(ids), tuple(g.get(tuple(sorted((i, j))), 2)
                                           for i, j in combinations(ids, 2)))
               for a, b in product((0, 2), (1, 3))
               for ids in [tuple(sorted(list(range(4, n)) + [a, b]))])


def domains(row):
    n = row['nodes']
    g = edges(row)
    contacts = [2] * n
    out = []
    counts = Counter()

    def visit(k):
        counts['nodes'] += 1
        if k == n:
            if any(m != 2 for m in contacts):
                out.append(contacts[:])
            return
        for m in ALPHABET:
            contacts[k] = m
            g[k, n] = m
            counts['label_trials'] += 1
            if m == 2 or corner_ok(n + 1, g):
                visit(k + 1)
            else:
                counts['corner_prunes'] += 1
        del g[k, n]
        contacts[k] = 2

    assert corner_ok(n, g)
    visit(0)
    assert len(out) == len(set(map(tuple, out)))
    return out, dict(counts)


def extension(row, contacts):
    n = row['nodes']
    g = edges(row)
    g.update({(i, n): m for i, m in enumerate(contacts)})
    return dict(nodes=n + 1, white=row['white'] + [n],
                labels=[g[p] for p in combinations(range(n + 1), 2)])


def inverse_pd(M):
    n = len(M)
    T = [r[:] + [A.K(i == j) for j in range(n)] for i, r in enumerate(M)]
    for k in range(n):
        d = T[k][k]
        assert d.sign() > 0
        T[k] = [x / d for x in T[k]]
        for i in range(n):
            if i != k:
                q = T[i][k]
                T[i] = [x - q * y for x, y in zip(T[i], T[k])]
    assert all(T[i][j] == (i == j) for i in range(n) for j in range(n))
    return [r[n:] for r in T]


def setup(row):
    colors = A.coloring(row)
    assert colors['passes']
    M = A.scaled_matrix(row, colors['colors'])
    I = row['white'] + [0, 1]
    Cinv = inverse_pd(A.sub(M, I))
    B = [[M[i][j] for i in I] for j in (2, 3)]
    T = [[A.dot(r, b) for r in Cinv] for b in B]
    S = [[M[i][j] - A.dot(B[ii], T[jj])
          for jj, j in enumerate((2, 3))] for ii, i in enumerate((2, 3))]
    assert S == A.schur(M, I)[0]
    determinant = S[0][0] * S[1][1] - S[0][1] * S[0][1]
    assert determinant.sign() <= 0
    return dict(colors=colors['colors'], I=I, Cinv=Cinv, B=B, T=T, S=S,
                singular=not determinant, scaled_gram=A.encode(M))


def finite_test(ext, cached, check_full=False):
    colors = A.coloring(ext)
    if not colors['passes']:
        return dict(status='REJECT_CYCLE_FIELD', conflict=colors['conflict'])
    assert colors['colors'][:-1] == cached['colors']
    G = A.scaled_matrix(ext, colors['colors'])
    v = [r[-1] for r in G[:-1]]
    w = [v[i] for i in cached['I']]
    z = [A.dot(r, w) for r in cached['Cinv']]
    rho = G[-1][-1] - A.dot(w, z)
    assert rho.sign() > 0 and rho.conj().sign() > 0
    q = [v[j] - A.dot(cached['T'][i], w) for i, j in enumerate((2, 3))]
    S = [[cached['S'][i][j] - q[i] * q[j] / rho for j in range(2)] for i in range(2)]
    det = S[0][0] * S[1][1] - S[0][1] * S[1][0]
    signs = [S[0][0].sign(), S[1][1].sign(), det.sign()]
    cs = [S[0][0].conj().sign(), S[1][1].conj().sign(), det.conj().sign()]
    assert signs[:2] == [-1, -1]
    status = ('REJECT_TWO_NEGATIVE_DIRECTIONS' if signs[-1] > 0 else
              'REJECT_CONJUGATE_NOT_PSD' if min(cs) < 0 else
              'PASS_EXTRA_WHITE_NECESSARY_TESTS')
    out = dict(status=status, physical_signs=signs, conjugate_signs=cs,
               determinant=det.enc())
    if check_full:
        old = A.conditions(ext, G)
        assert old['physical']['signs'] == signs and old['conjugate']['signs'] == cs
        assert old['physical']['schur'] == A.encode(S)
    return out


def core_record(index, row, union_index, check_full=False):
    out = dict(class_index=index, union_index=union_index, nodes=row['nodes'],
               dotted=0 in row['labels'])
    witness = OLD.affine_witness(row)
    if witness:
        out.update(status='REJECT_CORE_AFFINE_PRINCIPAL_SUBDIAGRAM', affine_witness=witness)
        return out
    cached = None if out['dotted'] else setup(row)
    out['singular'] = None if cached is None else cached['singular']
    profiles, counts = domains(row)
    out['domain_counts'] = counts
    out['structural_profile_count'] = len(profiles)
    trials, accepted = [], []
    for contact in profiles:
        ext = extension(row, contact)
        if not U.cycle_field(ext['nodes'], edges(ext)):
            trial = dict(status='REJECT_FIXED_FINITE_CYCLE_FIELD')
        elif cached is not None:
            trial = finite_test(ext, cached, check_full)
        else:
            trial = dict(status='DOTTED_PARAMETER_UNFILTERED')
        if trial['status'] in ('PASS_EXTRA_WHITE_NECESSARY_TESTS', 'DOTTED_PARAMETER_UNFILTERED'):
            witness = OLD.affine_witness(ext)
            if witness:
                trial = dict(status='REJECT_AFFINE_PRINCIPAL_SUBDIAGRAM', affine_witness=witness)
            else:
                accepted.append(contact)
        trial['contacts'] = contact
        trials.append(trial)
    out.update(trials=trials, accepted=accepted)
    if out['dotted']:
        out['status'] = 'DOTTED_STRUCTURAL_ONE_WHITE_DOMAIN'
    elif not accepted:
        out['status'] = ('SEALED_SINGULAR_FINITE_CORE' if out['singular'] else
                         'REJECT_SEALED_NONSINGULAR_FULL_STAR_RANK')
    else:
        out['status'] = 'FINITE_CORE_WITH_ONE_WHITE_EXTENSIONS'
    return out


def manifest():
    paths = [f'scripts/{TAG}.py', f'scripts/{BASE}.py', f'results/{BASE}.json',
             f'results/{BASE}_classes.json', f'results/{BASE}_review_complete.json',
             'scripts/fresh_sqrt5_quad_contacts_20261015.py',
             'results/fresh_sqrt5_quad_contacts_20261015_review.json',
             'scripts/fresh_sqrt5_quad_filters_20261015.py',
             'results/fresh_sqrt5_quad_filters_20261015_review.json']
    return {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in paths}


def load_inputs():
    base = json.loads((ROOT / 'results' / (BASE + '.json')).read_text())
    classes = json.loads((ROOT / 'results' / (BASE + '_classes.json')).read_text())
    assert base['status'] == 'COMPLETE_ALL_QUAD_UNIONS' and classes['marked_classes'] == 8298
    return base['labelled_union_supports'], classes['marked_isomorphism_classes']


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--resume', action='store_true')
    ap.add_argument('--verify', action='store_true')
    ap.add_argument('--max-cores', type=int, default=25)
    ap.add_argument('--check-full', action='store_true')
    args = ap.parse_args()
    rows, classes = load_inputs()
    if args.verify:
        saved = json.loads(REPORT.read_text())
        assert saved['sources'] == manifest()
        con = sqlite3.connect(f'file:{DB}?mode=ro', uri=True)
        assert con.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
        assert json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]) == saved['state']
        stored = [json.loads(r[0]) for r in con.execute('SELECT value FROM records ORDER BY class_index')]
        assert stored == saved['records']
        assert [r['class_index'] for r in stored] == list(range(saved['state']['next_class']))
        assert sha256(DB.read_bytes()).hexdigest() == saved['database_sha256']
        con.close()
        for rec in saved['records']:
            j = rec['class_index']
            ui = classes[j]['representative']
            assert rec == core_record(j, rows[ui], ui, args.check_full)
        print('PASS_SAVED_COMPLETE_PER_CORE_CONTACT_DOMAINS', len(saved['records']))
        return
    DB.parent.mkdir(parents=True, exist_ok=True)
    if not args.resume:
        assert not DB.exists(), 'Use --resume; an existing per-core frontier is never overwritten.'
    con = sqlite3.connect(DB)
    con.execute('PRAGMA journal_mode=DELETE')
    con.execute('PRAGMA synchronous=FULL')
    if not args.resume:
        con.execute('CREATE TABLE meta (id INTEGER PRIMARY KEY CHECK(id=1), value TEXT NOT NULL)')
        con.execute('CREATE TABLE records (class_index INTEGER PRIMARY KEY, value TEXT NOT NULL)')
        state = dict(sources=manifest(), next_class=0, batches=[])
        con.execute('INSERT INTO meta VALUES(1,?)', (json.dumps(state),))
        con.commit()
    else:
        state = json.loads(con.execute('SELECT value FROM meta WHERE id=1').fetchone()[0])
        assert state['sources'] == manifest()
    start = state['next_class']
    began = time.monotonic()
    stop = min(len(classes), start + args.max_cores)
    for j in range(start, stop):
        ui = classes[j]['representative']
        rec = core_record(j, rows[ui], ui, args.check_full)
        con.execute('INSERT INTO records VALUES(?,?)', (j, json.dumps(rec, separators=(',', ':'))))
        state['next_class'] = j + 1
        con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state, sort_keys=True),))
        con.commit()
        print('core', j, rec['status'], 'profiles', rec.get('structural_profile_count', 0),
              'accepted', len(rec.get('accepted', [])), flush=True)
    state['batches'].append(dict(first_class=start, next_class=state['next_class']))
    con.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state, sort_keys=True),))
    con.commit()
    records = [json.loads(r[0]) for r in con.execute('SELECT value FROM records ORDER BY class_index')]
    assert [r['class_index'] for r in records] == list(range(state['next_class']))
    con.close()
    counts = Counter(r['status'] for r in records)
    profile_counts = Counter(t['status'] for r in records for t in r.get('trials', []))
    out = dict(status='COMPLETE_ALL_QUAD_ONE_WHITE_DOMAINS' if len(records) == len(classes) else
               'PARTIAL_RESUMABLE_ALL_QUAD_ONE_WHITE_DOMAINS', sources=manifest(), state=state,
               database=str(DB.relative_to(ROOT)), database_sha256=sha256(DB.read_bytes()).hexdigest(),
               counts=dict(counts), profile_counts=dict(profile_counts), records=records,
               scope='Complete per-core domains for the processed prefix. No short-path or H4-payment pruning. Finite exact corner/field/Lorentz/conjugate tests; dotted distances remain symbolic. Nonzero contacts only; arbitrary mutually compatible collections and geometry remain separate.')
    REPORT.write_text(json.dumps(out, indent=2) + '\n')
    print(out['status'], state['next_class'], 'elapsed', round(time.monotonic() - began, 3),
          out['counts'], out['profile_counts'])


if __name__ == '__main__':
    main()
