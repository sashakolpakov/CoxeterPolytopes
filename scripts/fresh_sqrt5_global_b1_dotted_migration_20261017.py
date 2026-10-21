#!/usr/bin/env python3
"""Recompute every stopped prefix decision before migrating to the fast queue."""
import argparse
import json
import sqlite3
from collections import Counter
from hashlib import sha256

import fresh_sqrt5_global_b1_dotted_closure_20261016 as OLD
import fresh_sqrt5_global_b1_dotted_fast_20261017 as NEW

ROOT = NEW.ROOT
TAG = 'fresh_sqrt5_global_b1_dotted_migration_20261017'
OUT = ROOT / f'results/{TAG}.json'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def normalize(value):
    return json.loads(json.dumps(value))


def run(verify=False):
    old = json.loads(OLD.OUT.read_text())
    assert old['sources'] == OLD.manifest()
    assert digest(OLD.DB) == old['database_sha256']
    previous = sqlite3.connect(f'file:{OLD.DB}?mode=ro', uri=True)
    assert previous.execute('PRAGMA integrity_check').fetchone()[0] == 'ok'
    assert json.loads(previous.execute('SELECT value FROM meta WHERE id=1').fetchone()[0]) == old['state']
    raw = list(previous.execute('SELECT id,row,params,result FROM cores ORDER BY id'))
    assert len(raw) == len(old['records'])
    replacements, comparisons, changes = [], [], Counter()
    accepted = 0
    for rec, db in zip(old['records'], raw):
        i = rec['index']
        assert db[0] == i and json.loads(db[1]) == rec['row']
        assert json.loads(db[2]) == rec['parameters']
        assert (json.loads(db[3]) if db[3] else None) == rec['result']
        prior = rec['result']
        if prior is None:
            assert i >= old['state']['next_core']
            continue
        assert i < old['state']['next_core']
        fresh = normalize(NEW.evaluate(rec['row'], rec['parameters']))
        unchanged = {k: v for k, v in prior.items() if k not in ('children', 'trials')}
        assert unchanged == {k: v for k, v in fresh.items() if k != 'trials'}, i
        before = {tuple(t['contacts']): t for t in prior['trials']}
        after = {tuple(t['contacts']): t for t in fresh['trials']}
        assert before.keys() == after.keys()
        for contacts, t in before.items():
            n = after[contacts]
            if t['status'] != 'REJECT_AFFINE_PRINCIPAL':
                assert t == n, (i, contacts)
            else:
                assert n['status'].startswith('REJECT_'), (i, contacts, n)
                witness = t['witness']
                child = NEW.C.extension(rec['row'], list(contacts))
                ids, edges = witness['nodes'], NEW.Q.edges(child)
                labels = [edges[a, b] for a, b in OLD.combinations(ids, 2)]
                assert labels == witness['labels']
                block = dict(nodes=len(ids), labels=labels)
                coloring = OLD.C.OLD.A.coloring(block)
                assert coloring['passes']
                M = OLD.C.OLD.A.scaled_matrix(block, coloring['colors'])
                assert OLD.C.OLD.A.inertia(M) == [len(ids) - 1, 0, 1]
                changes[n['status']] += 1
        fresh['children'] = prior['children']
        payload = json.dumps(fresh)
        replacements.append((payload, i))
        accepted += len(fresh['accepted'])
        comparisons.append(dict(index=i, fresh_result_sha256=sha256(payload.encode()).hexdigest(),
                                trials=len(before), accepted=len(fresh['accepted'])))
        if i % 100 == 0:
            print('EXACT_DOTTED_PREFIX_COMPARISON', i, flush=True)
    assert len(replacements) == old['state']['next_core']
    sources = NEW.manifest()
    sources.update({str(OLD.DB.relative_to(ROOT)): digest(OLD.DB)})
    report = dict(status='PASS_ALL_STOPPED_PREFIX_DECISIONS_RECOMPUTED_WITH_RADICAL_LEMMA',
                  sources=sources, predecessor_processed=old['state']['next_core'],
                  predecessor_marked_states=old['marked_states'],
                  accepted_attachments_compared=accepted,
                  affine_rejections_replaced=dict(changes), comparisons=comparisons,
                  scope='Every processed predecessor core is freshly evaluated; complete contact domains, retained attachments, unchanged tests and child maps agree. Affine rejection witnesses are checked separately, and the new exact arithmetic tests also reject each such child. This verifies a partial prefix and the migration, not exhaustion of the remaining queue or face payments.')
    if verify:
        assert report == json.loads(OUT.read_text())
    else:
        assert not NEW.DB.exists(), 'Never overwrite or migrate into an existing queue.'
        successor = sqlite3.connect(NEW.DB)
        previous.backup(successor)
        state = dict(sources=NEW.manifest(), next_core=old['state']['next_core'],
                     initial_literal_roots=old['state']['initial_literal_roots'],
                     initial_marked_roots=old['state']['initial_marked_roots'], batches=[],
                     predecessor=dict(path=str(OLD.OUT.relative_to(ROOT)), sha256=digest(OLD.OUT),
                                      processed=old['state']['next_core'],
                                      marked_states=old['marked_states']))
        successor.execute('BEGIN IMMEDIATE')
        successor.executemany('UPDATE cores SET result=? WHERE id=?', replacements)
        successor.execute('UPDATE meta SET value=? WHERE id=1', (json.dumps(state),))
        successor.commit()
        NEW.export(successor, state)
        successor.close()
        OUT.write_text(json.dumps(report, indent=2) + '\n')
    previous.close()
    print(report['status'], len(replacements), dict(changes), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    run(parser.parse_args().verify)
