#!/usr/bin/env python3
"""Independent H8 facet support census, current witnesses and map checks."""
import argparse
import json
from itertools import combinations
import fresh_h9_beta131_pilot_20261016_review as I

ROOT=I.ROOT
TAG='fresh_h9_beta131_h100_witnesses_20261016_review'
PREFIX='results/fresh_h9_beta131_h100_20261016.json'
GEOMETRY='results/fresh_h9_beta131_h100_20261016_geometry.json'
FACET='results/fresh_h9_beta131_h100_20261016_facet.json'
MAPS='results/fresh_h9_beta131_h100_20261016_maps.json'
KNOWN='data/known_examples/bugaenko_compact_h8.json'


def matrix(rows): return [I.decode(row) for row in rows]
def transpose(rows): return [list(row) for row in zip(*rows)]
def multiply(a,b): return [[I.dot(row,col) for col in zip(*b)] for row in a]


def rref(rows):
    a=[list(row) for row in rows]
    n=len(a[0]);pivots=[]
    for j in range(n):
        i=len(pivots)
        k=next((k for k in range(i,len(a)) if a[k][j]),None)
        if k is None: continue
        a[i],a[k]=a[k],a[i]
        pivot=a[i][j]
        a[i]=[v/pivot for v in a[i]]
        for k in range(len(a)):
            if k!=i:
                factor=a[k][j]
                if factor: a[k]=[x-factor*y for x,y in zip(a[k],a[i])]
        pivots.append(j)
    return a,pivots


def run():
    paths=[PREFIX,GEOMETRY,FACET,MAPS,KNOWN]
    data,geometry,facet,maps,known=[json.loads((ROOT/p).read_text()) for p in paths]
    for document in (data,geometry,facet,maps):
        for p,h in document['sources'].items(): assert I.sha(p)==h,p
    form=matrix(data['form']);roots=matrix(data['roots'])
    center=I.decode(data['strict_center']);base=I.decode(data['basepoint'])
    assert len(roots)==41
    assert I.pair(center,center,form)<0 and all(I.pair(center,r,form)<0 for r in roots)
    escape=I.decode(geometry['strict_escape']['point'])
    assert I.pair(escape,escape,form)>0 and I.pair(escape,center,form)<0
    assert all(I.pair(escape,r,form)<0 for r in roots)
    seen=set()
    for row in geometry['finite_vertices']:
        v=I.decode(row['point']);active=tuple(row['incident'])
        assert len(active)==9 and active not in seen
        seen.add(active)
        assert I.pair(v,v,form)<0 and I.pair(v,center,form)<0
        values=[I.pair(v,r,form) for r in roots]
        assert all(z<=0 for z in values) and tuple(i for i,z in enumerate(values) if not z)==active
        assert len(rref([I.mv(form,roots[i]) for i in active])[1])==9
    assert len(seen)==45
    # Facet t=0 is literally the archived compact H8 example.
    ids=facet['source_to_current_wall_map'];r=roots[8]
    delta=I.K(1,12)
    assert I.mv(form,r)==[I.K()]*9+[-delta]
    q=I.pair(r,r,form)
    assert q==2*I.K(12,13)
    restricted=[row[:9] for row in form[:9]]
    source_roots=[roots[i][:9] for i in ids]
    assert restricted==matrix(known['lattice_gram']) and source_roots==matrix(known['simple_roots'])
    assert all(roots[i][9]==0 and I.pair(roots[i],r,form)==0 for i in ids)
    assert {i for i,u in enumerate(roots) if I.pair(u,r,form)==0}==set(ids)
    interior=I.decode(known['geometry']['interior_point']);interior10=interior+[I.K()]
    assert I.pair(interior10,interior10,form)<0 and I.pair(interior10,base,form)<0
    redundant=[]
    for i,u in enumerate(roots):
        if i==8 or i in ids: continue
        pairing=I.pair(u,r,form)
        projection=[a-pairing*b/q for a,b in zip(u,r)]
        assert projection[9]==0 and I.pair(projection,r,form)==0
        assert I.pair(projection,projection,form)<0 and I.pair(projection,interior10,form)<0
        redundant.append(i)
    assert len(redundant)==29
    covectors=[I.mv(restricted,u) for u in source_roots]
    assert len(rref(covectors)[1])==9
    assert all(I.dot(interior,c)<0 for c in covectors)
    finite={};supports=0
    for support in combinations(range(11),8):
        supports+=1
        rows,pivots=rref([covectors[i] for i in support])
        if len(pivots)!=8: continue
        free=next(i for i in range(9) if i not in pivots)
        v=[I.K()]*9;v[free]=I.K(1)
        for i,j in enumerate(pivots): v[j]=-rows[i][free]
        values=[I.dot(v,c) for c in covectors]
        if any(z>0 for z in values):
            if any(z<0 for z in values): continue
            v=[-z for z in v];values=[-z for z in values]
        assert I.pair(v,v,restricted)<0 and I.pair(v,interior,restricted)<0
        active=tuple(i for i,z in enumerate(values) if not z)
        assert len(active)==8
        finite[active]=v
    assert supports==165 and len(finite)==41
    assert set(finite)=={tuple(v) for v in known['geometry']['vertices']}
    assert {tuple(sorted([8]+[ids[i] for i in active])) for active in finite}=={a for a in seen if 8 in a}
    # Check each saved integral isometry and their common timelike fixed point;
    # this does not independently reproduce the corner-map search domain.
    fixed=I.decode(maps['common_fixed_timelike']['point'])
    assert I.pair(fixed,fixed,form)==I.K(I.F(13,29),I.F(-9,29))<0
    assert I.pair(fixed,base,form)<0
    identity=[[I.K(i==j) for j in range(10)] for i in range(10)]
    for row in maps['integral_maps']:
        m=matrix(row['matrix']);inverse=matrix(row['inverse'])
        assert all(x.is_integral for r0 in m+inverse for x in r0)
        assert multiply(m,inverse)==identity
        assert multiply(transpose(m),multiply(form,m))==form
        assert I.mv(m,fixed)==fixed
        assert all(I.mv(m,roots[i])==roots[j] for i,j in zip(row['source_ids'],row['target_ids']))
    assert len(maps['integral_maps'])==66
    paths += [f'scripts/{TAG}.py','scripts/fresh_h9_beta131_pilot_20261016_review.py',
              'src/coxeter_polytopes/exact/golden.py']
    return I.enc(dict(status='PASS_INDEPENDENT_BETA131_H100_WITNESSES_AND_COMPACT_H8_FACET',
                      sources={p:I.sha(p) for p in paths},finite_vertex_witnesses=45,
                      source_to_current_wall_map=ids,facet_supports=supports,
                      complete_facet_vertices=[dict(active=list(a),point=v) for a,v in sorted(finite.items())],
                      redundant_walls=redundant,integral_maps_with_common_timelike_fixed_point=66,
                      fixed_point=fixed,strict_escape=escape,
                      scope='Complete independent H8 facet census, all29 remaining wall projections,45 finite-vertex witnesses, strict ambient escape and66 saved isometries. Does not independently establish completeness of the41-wall ambient cone census or corner-map search.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true')
    args=parser.parse_args();result=run();path=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['finite_vertex_witnesses'],len(result['complete_facet_vertices']))
