#!/usr/bin/env python3
"""Immutable batch610 extraction of20 complete unpaid sealed domains.

The archived frontier is provenance, not a claim of whole-queue closure.
The source replays every extracted domain and every saved reachability map.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_nonsingular_attachment_search_20261015 as S

ROOT=Path(__file__).resolve().parents[1];TAG=S.TAG+'_batch610'
ARCHIVE='data/enumerations/'+TAG+'.json'

def run():
    D=json.loads((ROOT/ARCHIVE).read_text());assert D['sources']==S.manifest()
    assert D['state']['next_core']==610 and D['marked_states']==762
    assert D['counts']['PENDING']==152
    rows=D['records'];assert len(rows)==762
    base,lookup,weight=S.context();reach=set();maps=0
    for r in D['root_maps']:
        raw=S.F.extension(base[r['source_union_index']],r['contacts'])
        S.C.check_map(raw,rows[r['target']]['row'],r['map']);reach.add(r['target']);maps+=1
    assert len(D['root_maps'])==311 and reach==set(range(199))
    for i,rec in enumerate(rows):
        assert rec['index']==i
        if i>=610:assert rec['result']is None;continue
        assert i in reach
        result=rec['result'];children=result.get('children',[])
        assert len(children)==len(result.get('accepted',[]))
        for child,contact in zip(children,result.get('accepted',[])):
            assert child['contacts']==contact
            raw=S.F.extension(rec['row'],contact)
            S.C.check_map(raw,rows[child['index']]['row'],child['map']);reach.add(child['index']);maps+=1
    assert reach==set(range(762)) and maps==1477
    selected=[r for r in rows if r['result']and r['result']['status']=='UNPAID_SEALED_FRONTIER']
    assert [r['index']for r in selected]==D['unpaid_stopping_cases']and len(selected)==20
    out=[];trials=Counter()
    for r in selected:
        result=S.evaluate(r['row'],lookup,weight)
        assert result==r['result']and not result['accepted']and result['singular']
        assert result['payment_if_no_further_contact']['sums_without_with_external_a1']==['0','0']
        trials.update(t['status']for t in result['trials'])
        out.append(dict(index=r['index'],row=r['row'],result=result))
    paths=[f'scripts/{TAG}.py',f'scripts/{S.TAG}.py',ARCHIVE]
    S.Q.manifest_check(paths)
    return dict(status='PASS_IMMUTABLE_BATCH610_TWENTY_UNPAID_SEALED_DOMAINS',sources={p:sha256((ROOT/p).read_bytes()).hexdigest()for p in paths},
      archived_processed=610,archived_states=762,archived_pending=152,checked_literal_maps=maps,selected_indices=[r['index']for r in out],trial_counts=dict(trials),records=out,
      scope='Twenty literal sealed unpaid singular principal stars from the immutable610-core snapshot. All extracted nonzero one-white domains are replayed exactly and fail a necessary test or create persistent baseline payment. Both no-further-contact weight sums are0 under the bound current rule. This is primary shared-engine replay, not an independent arithmetic review or geometric realization. The whole search remains partial with152pending states; the archive database path/hash is historical storage metadata and is not reopened by this extraction.')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');args=p.parse_args();out=run();path=ROOT/'results'/(TAG+'.json')
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],len(out['records']),out['trial_counts'])
