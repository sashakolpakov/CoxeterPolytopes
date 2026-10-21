#!/usr/bin/env python3
"""Small complete finite-vertex incidence gate for the three new walls."""
import argparse,json
from fractions import Fraction as F
from itertools import combinations
from pathlib import Path
from hashlib import sha256
from fresh_h20_bridge_edge_closure_20261014_review import Field,parse
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h20_bridge_edge_closure_20261014_review_census'
BASE='results/fresh_h20_e7_edge_cusp_repair_20261014_exact.json'
PRE='results/fresh_h20_e7_bridge_face_20261014_exact.json'
NEW='results/fresh_h20_bridge_edge_closure_20261014.json'
REVIEW='results/fresh_h20_bridge_edge_closure_20261014_review.json'
def sha(p):return sha256((ROOT/p).read_bytes()).hexdigest()
def enc(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,(list,tuple)):return[enc(y)for y in x]
    if isinstance(x,dict):return{str(k):enc(v)for k,v in x.items()}
    return x
def run():
    old,pre,new,review=[json.loads((ROOT/p).read_text())for p in[BASE,PRE,NEW,REVIEW]]
    for d in[old,pre,new,review]:
        for p,h in d['sources'].items():assert sha(p)==h,p
    f=Field([parse(old['field']['discriminant'],0)])
    roots={i:tuple(parse(x,1)for x in old['normals'][i])for i in[61,28,87,0]}
    Q=[list(map(F,row))for row in old['form']]
    def dot(v,w):return f.total(f.scale(f.mul(v[i],w[j]),q)for i,row in enumerate(Q)for j,q in enumerate(row)if q)
    records=[]
    for i,j in[(61,28),(87,0)]:
        g=dot(roots[i],roots[j]);p=F(old['norms'][i])*F(old['norms'][j])
        margin=f.add(f.mul(g,g),f.const(-p))
        assert f.sign(g)<0 and f.sign(margin)>0
        if i==61:assert g==f.const(F(-11,2))
        records.append(dict(pair=[i,j],pairing=g,squared_dotted_margin=margin))
    S=pre['recipe']['right_support']
    finite=sorted(row['wall']for row in pre['new_pairs']if row['relation']in['right','pi3','pi4','pi5','pi6'])
    assert finite==sorted(S+[0,28])and len(finite)==20
    bad=[{61,28},{87,0}];assert bad[0].isdisjoint(bad[1])
    checks=0
    for selected in combinations(finite,19):
        assert any(pair<=set(selected)for pair in bad);checks+=1
    assert checks==20
    assert review['new_pairs_checked']==193
    assert review['new_pair_counts']=={'right':34,'pi3':4,'dotted':155}
    assert review['actual_distinct_finite_endpoints']==2
    # The bound full replay verifies precisely17 right neighbors per new cap,
    # pi3 at61/95, and strictly dotted relations to every remaining wall.
    assert new['finite_edge']['shared_support']==S+[95]
    assert [row['incident']for row in new['finite_endpoints']]==[S+[95,96],S+[95,97]]
    paths=[BASE,PRE,NEW,REVIEW,f'scripts/{TAG}.py',
           'scripts/fresh_h20_bridge_edge_closure_20261014_review.py',
           'scripts/fresh_h20_bridge_faces_20261014_review.py']
    return enc(dict(status='PASS_COMPLETE_FINITE_VERTEX_CENSUS_ON_WALLS95_96_97',sources={p:sha(p)for p in paths},
        independently_reconstructed_disjoint_dotted_pairs=records,wall95_old_finite_neighbors=finite,
        old19_subsets_checked=checks,wall95_finite_vertices_in96=0,
        each_new_wall_finite_neighbor_set=S+[95],each_new_wall_unique_finite_vertex=True,
        finite_vertices_on_union_of_three_walls_in98=2,
        scope='Two old pairings independently reconstructed. Uses the bound complete95-pair and193-pair replays for finite-neighbor coverage. No global old-chamber vertex census.'))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');a=ap.parse_args();d=run()
    p=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(p.read_text())
    else:p.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'])
