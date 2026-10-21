#!/usr/bin/env python3
"""All possible additional cross-core edges in the finite7 overlap domain.

Monotone actual-edge/corner witnesses exclude any larger native contacts.
Remaining positions are a necessary domain, not asserted realizations.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_distance7_high_overlap_domain_20261017 as B
import fresh_distance8_high_disjoint_20261017 as D

H=B.H;F=H.F;ROOT=H.ROOT
TAG='fresh_distance7_high_cross_domain_20261017'


def normalized(a):
    try:
        g,scale=H.native_to_rational(a)
        return g,scale,'positive_diagonal_congruence'
    except AssertionError:
        g=[[F(2 if i==j else -int(bool(a[i][j])))for j in range(len(a))]for i in range(len(a))]
        return g,[1]*len(a),'unit_bond_relaxation'


def corner_bad(native,u,v):
    white=[i for i in range(len(native))if i not in u+v]
    for p,q in product(u,v):
        ids=sorted(white+[p,q]);g,scale,mode=normalized(H.extract(native,ids))
        o=D.strict_obstruction(g)
        if o:
            x=list(map(F,o['vector']))
            assert all(z>=0 for z in x)and D.R.square(g,x)<0
            return dict(corner=[p,q],nodes=ids,gauge_squared=scale,comparison=mode,**o)
    return None


def run():
    base=json.loads((ROOT/f'results/{B.TAG}.json').read_text())
    for rel,h in base['sources'].items():assert sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
    data=json.loads((ROOT/f'results/{H.CENSUS}.json').read_text())
    diagrams={d['index']:d for d in data['diagrams']}
    graphs={i:H.bond_graph(d)for i,d in diagrams.items()}
    @lru_cache(None)
    def edge_test(native,black):
        return H.edge_bad_native(native,list(black))
    @lru_cache(None)
    def corner_test(native,u,v):
        return corner_bad(native,list(u),list(v))
    proofs=[];proof_ids={}
    def record(proof):
        text=json.dumps(proof,sort_keys=True,separators=(',',':'))
        if text not in proof_ids:proof_ids[text]=len(proofs);proofs.append(proof)
        return proof_ids[text]
    counts=Counter();rows=[]
    for index,b in enumerate(base['survivors']):
        di,u=b['left'];dj,v0=b['right'];match=tuple(map(tuple,b['match']))
        native,mp=H.unite(graphs[di],graphs[dj],match);v=[mp[i]for i in v0]
        assert min(H.graph_dist(native,u)[i]for i in v)==b['minimum_marked_distance']>=8
        key=tuple(map(tuple,native));black=tuple(sorted(u+v))
        ob=corner_test(key,tuple(u),tuple(v))
        row=dict(index=index,base_context=b,right_copy_map=mp)
        if ob:
            row.update(status='base_actual_corner_not_PSD',base_corner_witness=record(ob))
            counts['base_corner_rejected']+=1
        else:
            left=[i for i in range(len(graphs[di]))if i not in{a for a,z in match}]
            right=[mp[i]for i in range(len(graphs[dj]))if i not in{z for a,z in match}]
            shortcuts=[];edge_bad=[];corner_rejected=[];remaining=[]
            for i,j in product(left,right):
                assert i!=j and native[i][j]==0
                g=[r.copy()for r in native];g[i][j]=g[j][i]=F(-1)
                ds=min(H.graph_dist(g,u)[p]for p in v)
                if ds<=7:
                    shortcuts.append([i,j,ds]);continue
                literal=tuple(map(tuple,g));ob=edge_test(literal,black)
                if ob:
                    edge_bad.append([i,j,record(ob)]);continue
                ob=corner_test(literal,tuple(u),tuple(v))
                if ob:corner_rejected.append([i,j,record(ob)])
                else:remaining.append([i,j,ds])
            row.update(status='remaining_cross_positions'if remaining else'base_gram_forced',
                short7_positions=shortcuts,edge_PD_rejections=edge_bad,
                corner_PSD_rejections=corner_rejected,remaining_positions=remaining)
            counts.update(contexts_with_remaining_positions=int(bool(remaining)),forced_base_grams=int(not remaining),
                short7_positions=len(shortcuts),edge_PD_rejections=len(edge_bad),
                corner_PSD_rejections=len(corner_rejected),remaining_positions=len(remaining))
        rows.append(row)
        if (index+1)%250==0:print('CROSS',index+1,'/6691',dict(counts),flush=True)
    assert len(rows)==6691
    assert counts['base_corner_rejected']+counts['contexts_with_remaining_positions']+counts['forced_base_grams']==6691
    paths=[Path(__file__),Path(B.__file__),Path(H.__file__),Path(D.__file__),Path(D.R.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths.extend(ROOT/f'{folder}/{tag}.{ext}'for tag in(B.TAG,H.CENSUS)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths=list(dict.fromkeys(paths))
    out=dict(status='PASS_COMPLETE_DISTANCE7_HIGH_CROSS_POSITION_DOMAIN',sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        base_context_count=len(rows),counts=dict(counts),witnesses=proofs,cases=rows,
        remaining_position_multiplicities=dict(Counter(len(r.get('remaining_positions',[]))for r in rows if r['status']=='remaining_cross_positions')),
        scope='All additional individual cross-core positions in the complete6691 necessary base-overlap contexts. Shortcuts and nonnegative strict-edge/strict-negative-corner witnesses exclude any simultaneous larger contacts containing the tested position. Remaining positions require their complete subset/label and supporting-extension analysis; they are not asserted realizations. No payment or dimension theorem.')
    return json.loads(json.dumps(out))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    out=run();path=ROOT/f'results/{TAG}.json'
    if a.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['counts'],out['remaining_position_multiplicities'],flush=True)
