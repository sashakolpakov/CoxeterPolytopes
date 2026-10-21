#!/usr/bin/env python3
"""Full profile reconstruction and terminal SQLite/export consistency replay."""
import argparse,json,sqlite3
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_triangle_extension_complete_20261016 as P
ROOT=P.ROOT;TAG='fresh_sqrt5_triangle_extension_complete_20261016_replay'
RESULT='results/fresh_sqrt5_triangle_extension_complete_20261016_rank22.json'
DATABASE='data/enumerations/fresh_sqrt5_triangle_extension_complete_20261016_rank22.sqlite'
def run():
    P.verify(ROOT/RESULT);d=json.loads((ROOT/RESULT).read_text())
    assert d['status']=='PASS_COMPLETE_CONNECTED_TRIANGLE_GRAPH_DOMAIN_DIMENSION21'
    assert d['complete_dimension21_connected_domain']and d['complete_through_bound']and not d['pending']
    db=sqlite3.connect(f'file:{ROOT/DATABASE}?mode=ro',uri=True)
    assert db.execute('PRAGMA integrity_check').fetchone()==('ok',)
    assert db.execute("SELECT value FROM meta WHERE key='source_hash'").fetchone()[0]==P.digest(f'scripts/{P.TAG}.py')
    assert db.execute("SELECT value FROM meta WHERE key='max_nodes'").fetchone()[0]=='22'
    stored=list(db.execute('SELECT id,row,status,payload FROM states ORDER BY id'));assert len(stored)==len(d['records'])
    for (idx,row,status,payload),rec in zip(stored,d['records']):
        expanded=dict(index=idx,row=json.loads(row),status=status,**json.loads(payload));assert expanded==rec
    assert [json.loads(r[0])for r in db.execute('SELECT payload FROM seeds ORDER BY id')]==d['seed_maps']
    assert db.execute("SELECT count(*) FROM states WHERE status='PENDING'").fetchone()[0]==0
    db.close()
    files=[f'scripts/{TAG}.py',f'scripts/{P.TAG}.py',RESULT,DATABASE]
    return dict(status='PASS_FRESH_COMPLETE_TRIANGLE_PROFILE_MAP_AND_DATABASE_REPLAY',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in files},connected_states=d['connected_states'],complete_dimension21_connected_domain=True,pending=0,counts=d['counts'],trial_counts=d['trial_counts'],literal_maps=76+sum(t['status']=='ACCEPT'for r in d['records']for t in r.get('trials',[])),scope='Same-producer exact reconstruction of every profile, arithmetic decision and literal map, plus terminal database equality/integrity. Independent mathematical review and all spectator/facet-mark payments remain separate obligations.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['connected_states'],d['literal_maps'],flush=True)
