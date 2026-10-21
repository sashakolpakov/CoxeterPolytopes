#!/usr/bin/env python3
"""Immutable31-star input for the enlarged forest-sensitive weight model."""
import argparse,json
from hashlib import sha256
from pathlib import Path
import fresh_sqrt5_all_quad_unions_20261015_classes as C

ROOT=Path(__file__).resolve().parents[1]
BASE='fresh_sqrt5_two_chain_star_closure_20261015';TAG=BASE+'_export'


def run():
    data=json.loads((ROOT/f'results/{BASE}.json').read_text())
    assert data['open_attachment_states']==0
    for p,h in data['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    geometry={r['index']:r for r in data['unpaid_geometry']}
    obligations={r['index']:r for r in data['new_weight_obligations']}
    rows=[];keys=set()
    for rec in data['records']+data['child_closures']:
        if rec['status'] not in ('NEW_SEALED_PAID','NEW_SEALED_UNPAID','SEALED_PAID_CHILD'):continue
        assert not rec['accepted']
        key,order=C.canonical(rec['row']);assert key not in keys;keys.add(key)
        child='parent' in rec
        out=dict(index=len(rows),source_kind='child' if child else 'initial',
                 source_index=rec['parent'] if child else rec['index'],
                 source_status=rec['status'],row=rec['row'],canonical_order=list(order),
                 current_payment=rec['payment'])
        if rec['status']=='NEW_SEALED_UNPAID':
            out['exact_geometry']=geometry[rec['index']]
            out['necessary_weight_obligation']=obligations[rec['index']]
        rows.append(out)
    assert len(rows)==31 and sum(not r['current_payment']['paid'] for r in rows)==3
    paths=[f'scripts/{TAG}.py',f'scripts/{BASE}.py',f'results/{BASE}.json',
           'scripts/fresh_sqrt5_all_quad_unions_20261015_classes.py']
    return dict(status='PASS_IMMUTABLE_31_TWO_CHAIN_SEALED_STAR_EXPORT',
                sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                marked_classes=31,currently_paid=28,currently_unpaid=3,records=rows,
                scope='Literal marked diagrams for31 newly checked sealed classes from the two-chain closure. All sources are bound. This is an input export, not a new independent domain or geometry replay; arbitrary further unpaid white attachments are orthogonal by the source sealing proof, with externalA1 contexts retained.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],out['marked_classes'])
