#!/usr/bin/env python3
"""Complete base-overlap domain with no marked path of length at most seven.

This is a necessary-domain diagnostic. Extra cross-core edges and arbitrary
supporting whites have not been exhausted, and no payment theorem is claimed.
"""
import argparse
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import fresh_distance8_high_overlap_20261017 as H

ROOT=H.ROOT
TAG='fresh_distance7_high_overlap_domain_20261017'


def matchings(ga,gb,left,right):
    """Every nonempty compatible partial injection, each exactly once."""
    choices={i:[j for j,e in enumerate(right['distances'])if d+e>=8]
             for i,d in enumerate(left['distances'])}
    nodes=[i for i in choices if choices[i]]
    def visit(pos,used,pairs):
        if pos==len(nodes):
            if pairs:yield tuple(pairs)
            return
        i=nodes[pos]
        yield from visit(pos+1,used,pairs)
        for j in choices[i]:
            if j in used or any(ga[i][a]!=gb[j][b]for a,b in pairs):continue
            yield from visit(pos+1,used|{j},pairs+[(i,j)])
    yield from visit(0,set(),[])


def run():
    data=json.loads((ROOT/f'results/{H.CENSUS}.json').read_text())
    for rel,h in data['sources'].items():assert sha256((ROOT/rel).read_bytes()).hexdigest()==h,rel
    diagrams={d['index']:d for d in data['diagrams']}
    graphs={i:H.bond_graph(d)for i,d in diagrams.items()}
    contexts=H.all_contexts(data['diagrams'])
    assert max(max(c['distances'])for c in contexts)==7
    @lru_cache(None)
    def diagonals(di,black,node):
        out=[]
        for b in black:
            ids=[i for i in range(len(graphs[di]))if i!=b]
            g,scale=H.native_to_rational(H.extract(graphs[di],ids));H.psd(g,True)
            iv=H.inverse(g);i=ids.index(node);out.append(scale[i]*iv[i][i])
        return tuple(out)
    @lru_cache(None)
    def union_test(di,dj,match,bl,br):
        native,mp=H.unite(graphs[di],graphs[dj],match)
        black=list(bl)+[mp[i]for i in br]
        # Identifications can create shorter routes than the one-vertex test.
        distance=min(H.graph_dist(native,bl)[mp[i]]for i in br)
        if distance<=7:return 'created_shortcut',distance,None
        bad=H.edge_bad_native(native,black)
        if bad:return 'actual_edge_not_PD',None,bad
        return 'survives_base_checks',distance,None
    counts=Counter();sizes=Counter();survivors=[];stream=sha256();number=0;next_progress=10000
    for left,right in product(contexts,repeat=2):
        ga,gb=graphs[left['diagram']],graphs[right['diagram']]
        for match in matchings(ga,gb,left,right):
            number+=1;sizes[len(match)]+=1
            key=[left['diagram'],list(left['black']),right['diagram'],list(right['black']),match]
            reject=None
            for i,j in match:
                for side,node,othernode,c,other,shared in (
                    ('left',i,j,left,gb,{b for a,b in match}),
                    ('right',j,i,right,ga,{a for a,b in match})):
                    ext=[k for k in range(len(other))if k!=othernode and other[othernode][k]and k not in shared]
                    if ext:
                        value=max(diagonals(c['diagram'],tuple(c['black']),node))
                        other_black=right['black']if side=='left'else left['black']
                        # A white extension sits in a strictly PD edge. A
                        # black extension sits only in a PSD corner: equality
                        # at2 is allowed and must not be discarded.
                        bad=[k for k in ext if value>2 or(value==2 and k not in other_black)]
                        if bad:
                            reject=[side,node,bad,str(value),
                                    ['corner_PSD'if k in other_black else'edge_PD'for k in bad]];break
                if reject:break
            if reject:kind='inverse_rejection';distance=proof=None
            else:
                kind,distance,proof=union_test(left['diagram'],right['diagram'],match,
                                              tuple(left['black']),tuple(right['black']))
            counts[kind]+=1
            stream.update((json.dumps([key,kind,reject,proof],separators=(',',':'))+'\n').encode())
            if kind=='survives_base_checks':
                survivors.append(dict(left=[left['diagram'],list(left['black'])],
                    right=[right['diagram'],list(right['black'])],match=[list(p)for p in match],
                    minimum_marked_distance=distance))
        if number>=next_progress:
            print('DOMAIN',number,dict(counts),flush=True);next_progress=number+10000
    assert sum(counts.values())==sum(sizes.values())==number
    old=json.loads((ROOT/f'results/{H.TAG}.json').read_text())
    def key(c):return(c['left'][0],tuple(c['left'][1]),c['right'][0],tuple(c['right'][1]),tuple(map(tuple,c['match'])))
    old_keys={key(c)for c in old['isolation']['cases']};new_keys={key(c)for c in survivors}
    assert len(new_keys)==len(survivors) and old_keys<=new_keys
    paths=[Path(__file__),Path(H.__file__),ROOT/'scripts/fresh_selective_distance10_dimension116_20261017_review.py']
    paths.extend(ROOT/f'{folder}/{s}.{ext}'for s in(H.CENSUS,H.TAG)for folder,ext in(('scripts','py'),('results','json'),('results','md')))
    paths=list(dict.fromkeys(paths))
    return dict(status='PASS_COMPLETE_DISTANCE7_BASE_OVERLAP_DOMAIN',
        marked_contexts=len(contexts),consistent_matchings=number,classification=dict(counts),
        matching_sizes={str(k):v for k,v in sorted(sizes.items())},
        deterministic_classification_digest=stream.hexdigest(),survivors=survivors,
        old_distance8_survivors=len(old_keys),new_survivors=len(new_keys-old_keys),
        sources={str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest()for p in paths},
        scope='Every label-consistent partial white identification under the necessary distance-sum>=8 condition; inverse/actual-edge exclusions and created-shortcut tests. Only base unions are classified. Extra cross edges, arbitrary-support extensions, low/disjoint cores, ideal compression and quadrilateral payments remain open; no global bound.')


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    out=json.loads(json.dumps(run()));path=ROOT/f'results/{TAG}.json'
    if args.verify:assert out==json.loads(path.read_text())
    else:path.write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(out['status'],out['consistent_matchings'],out['classification'],
          'new',out['new_survivors'],flush=True)
