#!/usr/bin/env python3
"""Two exact general endpoint caps, retaining the non-quasi-arithmetic seed.

The caps lie in Q(sqrt5), although the retained source requires the quartic
field. They are mutually right; neither is an old extreme-ray polar.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
from fresh_h9_mixed_d4_coupling_20261014 import (
    L,K,Z,O,sm,pair,mv,inertia,cone,relation,escape,enc,rank)
ROOT=Path(__file__).resolve().parents[1]
TAG='fresh_h9_mixed_d4_endpoint_pair_20261014'
BASE='results/fresh_h9_mixed_d4_coupling_20261014.json'
QUEUE='results/fresh_h9_mixed_d4_coupling_20261014_frontier.json'

def dec(x):return L(K(*x[0]),K(*x[1]))
def dv(v):return [dec(x) for x in v]

def build():
    old=json.loads((ROOT/BASE).read_text());queue=json.loads((ROOT/QUEUE).read_text())
    for doc in (old,queue):
        for p,h in doc['sources'].items():assert sha256((ROOT/p).read_bytes()).hexdigest()==h
    Q=[dv(r) for r in old['form']];N=[dv(r) for r in old['final_normals']]
    coeffs=[
        [('65/4','29/4'),('20','9'),('57/4','27/4'),('19/2','9/2'),('19/4','9/4'),
         ('2','1'),('-9/2','-3/2'),('-4','-3/2'),('-11/4','-3/4'),('-9/4','-3/4')],
        [('33/4','15/4'),('21/2','9/2'),('15/2','3'),('5','2'),('5/2','1'),
         ('1/2','1/2'),('-3','-1'),('-2','-1'),('-3/2','-1/2'),('-3/2','-1/2')]]
    caps=[[L(K(a,b)) for a,b in row] for row in coeffs]
    recipes=[dict(edge_id=6,contacts=[[2,3],[9,3]],anchor=1,anchor_label=5),
             dict(edge_id=78,contacts=[[5,3],[7,3]],anchor=2,anchor_label=5)]
    assert pair(caps[0],caps[1],Q)==Z
    oldE=dv(old['strict_escape']['vector']);cut_values=[pair(r,oldE,Q) for r in caps]
    assert all(x.sign()>0 for x in cut_values)
    intended=[]
    for recipe,r in zip(recipes,caps):
        assert pair(r,r,Q)==L(2)
        entry=queue['queue'][recipe['edge_id']];S=entry['support'];start=dv(entry['start_point'])
        labels=dict(recipe['contacts']);labels[recipe['anchor']]=recipe['anchor_label']
        for i in S+[recipe['anchor']]:
            expected=Z if labels.get(i,2)==2 else -O if labels[i]==3 else -L(K(F(1,2),F(1,2)))
            assert pair(r,N[i],Q)==expected
        assert pair(r,start,Q).sign()<0
        intended.append(dict(**recipe,support=S,start=start,start_pairing=pair(r,start,Q)))
    # All old finite vertices remain strict at both caps, with no heuristic test.
    oldfinite=[dv(r['vector']) for r in old['final_complete_ray_census']['rays'] if r['kind']=='finite']
    assert all(pair(v,r,Q).sign()<0 for v in oldfinite for r in caps)
    N+=caps
    h=[sm(v[j] for v in oldfinite)+dv(old['center'])[j]/1000 for j in range(10)]
    assert pair(h,h,Q).sign()<0 and all(pair(h,n,Q).sign()<0 for n in N)
    pairs=[];counts=Counter()
    for i,j in combinations(range(len(N)),2):
        rel=relation(pair(N[i],N[j],Q),pair(N[i],N[i],Q),pair(N[j],N[j],Q))
        assert rel['type'] in ('right','finite','parallel','dotted')
        if j>=12:pairs.append(dict(pair=[i,j],**rel));counts[rel['type']+str(rel.get('label',''))]+=1
    census=cone(Q,N);print('complete rays',census['counts'],flush=True)
    finite=[row for row in census['rays'] if row['kind']=='finite'];newfinite=[]
    for row in finite:
        ids=row['incident'];assert len(ids)==9
        assert inertia([[pair(N[i],N[j],Q) for j in ids] for i in ids])==[9,0]
        if any(i>=12 for i in ids):newfinite.append(row)
    for item in intended:
        found=[row for row in finite if set(item['support'])<=set(row['incident']) and any(i>=12 for i in row['incident'])]
        assert len(found)==1;item['finite_endpoint']=found[0]
    facets=[]
    for i,n in enumerate(N):
        v=[a-pair(h,n,Q)*b/pair(n,n,Q) for a,b in zip(h,n)]
        assert pair(v,v,Q).sign()<0 and all(pair(v,w,Q).sign()==(0 if i==j else -1) for j,w in enumerate(N))
        facets.append(v)
    # Recompute the complete finite-corner queue from the complete cone census.
    rows=[mv(Q,n) for n in N];directions=[];cached={};edge_counts=Counter()
    for vi,v in enumerate(census['rays']):
        if v['kind']!='finite':continue
        for released in v['incident']:
            support=sorted(set(v['incident'])-{released});key=tuple(support)
            if key not in cached:cached[key]=rank([rows[i] for i in support])
            assert cached[key]==8
            other=[j for j,r in enumerate(census['rays']) if j!=vi and set(support)<=set(r['incident'])]
            assert len(other)==1;j=other[0];kind=census['rays'][j]['kind'];assert kind!='ideal'
            edge_counts[kind]+=1
            directions.append(dict(start_ray=vi,released_wall=released,support=support,
                other_cone_ray=j,kind='finite_endpoint' if kind=='finite' else 'open_to_light_cone'))
    files=[BASE,QUEUE,f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014.py']
    return enc(dict(status='PASS_EXACT_TWO_GENERAL_ENDPOINT_CAPS_MIXED_D4_H9',
        sources={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        form=Q,normals=N,added_normals=caps,added_norms=[L(2),L(2)],mutual_product=Z,
        old_escape_cut_pairings=cut_values,center=h,center_norm=pair(h,h,Q),actual_facet_points=facets,
        new_pairs=pairs,new_pair_counts=dict(counts),targeted_edge_endpoints=intended,
        preserved_old_finite_vertex_count=len(oldfinite),new_finite_vertices=newfinite,
        complete_ray_census=census,finite_corner_directions=directions,finite_corner_counts=dict(edge_counts),
        strict_escape=escape(Q,N,h,census),
        scope='Coherent fourteen-wall H9 chamber with two mutually right general endpoint caps. Complete exact cone census, all finite vertices and existing finite-corner directions. Both caps cut the prior escape; the new strict rational escape proves infinite volume. No arbitrary-new-mirror enumeration or new finite-volume dimension.'))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=build();out=ROOT/f'results/{TAG}.json'
    if a.verify:assert d==json.loads(out.read_text())
    else:out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['complete_ray_census']['counts'],len(d['new_finite_vertices']),d['finite_corner_counts'],flush=True)
