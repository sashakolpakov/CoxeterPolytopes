#!/usr/bin/env python3
"""Read-only independent literal transport check; no label search replay."""
import argparse,json
from hashlib import sha256
from itertools import combinations
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_star_pi4_resume_20261015_gates_review'
P='results/fresh_h9_star_pi4_resume_20261015_gates.json'
N='results/fresh_h9_star_third_four_20261015.json'
M='results/fresh_h9_star_core_transfers_20261014.json'
G='results/fresh_h9_tree_cube_wedges_20261014.json'
OLD='results/fresh_h9_star_single_third_opposite_20261015_gates.json'
BATCH='results/fresh_h9_star_single_third_opposite_20261015_search_batch_014.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def read(p):return json.loads((ROOT/p).read_text())
def key(r):return tuple(tuple(sorted((min(a,b),max(a,b),m)for a,b,m in r[k]))for k in('fixed','bounds'))
def run():
    p,n,m,g,old,batch=map(read,[P,N,M,G,OLD,BATCH])
    for d in[p,n,m,g,old,batch]:
        for f,h in d.get('sources',{}).items():assert sha(f)==h,f
    review=read('results/fresh_h9_star_third_four_20261015_review.json')
    mapsreview=read('results/fresh_h9_star_core_transfers_20261014_review.json')
    assert review['status'].startswith('PASS_')and review['sources'][N]==sha(N)
    assert mapsreview['status'].startswith('PASS_')and mapsreview['sources'][M]==sha(M)
    geometry=g['representatives'][0]['geometry'];vertices=list(map(lambda x:set(x['incident']),geometry['vertices']))
    oldkeys={key(f)for f in old['families']};assert len(oldkeys)==397
    newkeys={key(f)for f in p['families']}-oldkeys
    assert len(p['families'])==493 and len(newkeys)==96
    supportlist=[r['support']for r in n['positive_kernel_guards']+n['projection_budgets']]
    assert len(supportlist)==114 and len(n['families'])==8
    core=set(n['rank_core']);assert core==set(m['source_core'])and len(core)==11
    assert n['complete_maximal_core_faces']==m['source_maximal_faces']
    transported=set();counts=[]
    for placement in m['placements']:
        f=dict(placement['core_map']);assert set(f)==core and len(set(f.values()))==11
        image=set(f.values());traces={frozenset(v&image)for v in vertices}
        maxima={v for v in traces if not any(v<w for w in traces)}
        assert maxima=={frozenset(f[x]for x in v)for v in m['source_maximal_faces']}
        assert all(any(set(t)<=v for v in vertices)for t in combinations(image,2))
        count=0
        for S in supportlist:
            target={f[x]for x in S};assert any(target<=v for v in vertices);count+=1
        for family in n['families']:
            mapped={k:[[f[a],f[b],v]for a,b,v in family[k]]for k in('fixed','bounds')}
            kk=key(mapped);assert kk not in oldkeys and kk not in transported;transported.add(kk)
        counts.append(count)
    assert transported==newkeys and counts==[114]*12
    labels={tuple(x[:2]):x[2]for x in batch['solution']['pair_labels']}
    def hit(k):return all(labels[a,b]==v for a,b,v in k[0])and all(labels[a,b]>=v for a,b,v in k[1])
    assert not any(hit(k)for k in oldkeys)
    hits=[f['id']for f in p['families']if key(f)in newkeys and hit(key(f))]
    assert sorted(hits)==sorted(p['saved_survivor_new_gate_hits'])and len(hits)==3
    assert p['saved_cursor']==batch['cumulative_additional_nodes']==60015
    paths=[P,N,M,G,OLD,BATCH,'results/fresh_h9_star_third_four_20261015_review.json',
        'results/fresh_h9_star_core_transfers_20261014_review.json',
        'scripts/fresh_h9_star_pi4_resume_20261015_gates.py',
        'scripts/fresh_h9_star_pi4_resume_20261015_search.py',f'scripts/{TAG}.py']
    return dict(status='PASS_INDEPENDENT_PI4_BINDER_LITERAL_TRANSPORT_SCOPE',sources={f:sha(f)for f in paths},
        placements=12,support_occurrences=sum(counts),new_predicates=96,old_predicates=397,total_predicates=493,
        saved_cursor=60015,saved_survivor_hits=hits,
        scope='Independent literal maps, face premises, exact predicate transport and saved-survivor check. Third4 mathematical exclusions and prior397 predicates remain reviewed dependencies. No resumed DFS replay, whole-incidence exclusion, or metric construction.')
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run();out=ROOT/'results'/f'{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text());print('PASS saved independent pi4 binder review')
    else:out.write_text(json.dumps(d,indent=2)+'\n');print(d['status'])
