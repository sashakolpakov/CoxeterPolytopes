#!/usr/bin/env python3
"""Exact preferred-core images in order16 affine factors.

This is a local census; universal replacement of the antipodal A15 cycle
corner and any further finite coupons are separate obligations.
"""
import argparse,json
from collections import Counter
from hashlib import sha256
from pathlib import Path
import fresh_ideal7_finite_cross_payments_20261017 as P

ROOT=P.ROOT
TAG='fresh_order16_preferred_embeddings_20261017'

def run():
    old=json.loads((ROOT/f'results/{P.OLD}.json').read_text())
    diagrams={r['index']:r for r in json.loads((ROOT/f'results/{P.CENSUS}.json').read_text())['diagrams']}
    low=json.loads((ROOT/'results/fresh_distance9_low_core_20261017.json').read_text())
    assert max(r['selected_corner']['component_order_bound']for r in low['extensions']['cases'])<=12
    patterns=[]
    for index,r in enumerate(old['isolation']['cases']):
        if r['component_order_bound']<16:continue
        di,u=r['left'];dj,v=r['right']
        g,mp=P.R.unite(P.R.bond_graph(diagrams[di]),P.R.bond_graph(diagrams[dj]),tuple(map(tuple,r['match'])))
        ids=r['marked_component'];labels=[[int(2-g[i][j])if i!=j and g[i][j]else 0 for j in ids]for i in ids]
        patterns.append(dict(index=index,nodes=ids,labels=labels,marks=[ids.index(x)for x in r['preferred_corner']],
            method=r['method'],order_bound=r['component_order_bound']))
    assert len(patterns)==260
    records={}
    for kind in'ABCD':
        target=P.labelled(16,P.A.P.graph(kind,16)[0]);basic=P.basic_pairs(target);extra=set();rows=[]
        for pattern in patterns:
            maps=P.embeddings(pattern['labels'],target)if len(pattern['nodes'])<=16 else[]
            pairs={tuple(sorted(m[i]for i in pattern['marks']))for m in maps};missing=pairs-basic
            extra|=missing
            if maps:rows.append(dict(case=pattern['index'],embeddings=maps,
                marked_images=[list(p)for p in sorted(pairs)],extra_images=[list(p)for p in sorted(missing)]))
        paid=basic|extra
        records[kind]=dict(target=target,basic_pairs=[list(p)for p in sorted(basic)],extra_pairs=[list(p)for p in sorted(extra)],
            paid_pairs=[list(p)for p in sorted(paid)],cost=len(paid),case_embeddings=rows,embedding_count=sum(len(r['embeddings'])for r in rows))
        print('ORDER16',kind,'basic',len(basic),'extra',len(extra),'cost',len(paid),'cases',len(rows),flush=True)
    stems=[P.OLD,P.CENSUS,'fresh_distance9_low_core_20261017',P.TAG]
    paths=[Path(__file__),Path(P.__file__)]+[ROOT/f'{folder}/{stem}.{ext}'for stem in stems for folder,ext in(('scripts','py'),('results','json'),('results','md'))]
    return dict(status='PASS_COMPLETE_PREFERRED_ORDER16_EMBEDDINGS',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        patterns=patterns,targets=records,component_bound_counts=dict(Counter(r['order_bound']for r in patterns)),
        low_core_order_bound=12,
        scope='Complete induced labelled embeddings of all260 inherited preferred components with extension bound>=16 into each classical order16 affine diagram. Pattern images plus the retained basic pairs are complete under the inherited no-short9 hypothesis. Universal treatment of the exceptional ideal antipodal cycle corner is a separate prerequisite; no new dimension bound asserted.')

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();d=json.loads(json.dumps(run()))
    p=ROOT/f'results/{TAG}.json'
    if args.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],flush=True)
