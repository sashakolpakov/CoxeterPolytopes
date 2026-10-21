#!/usr/bin/env python3
"""Independent exhaustive support enumeration, rather than double description."""
import argparse
import json
from itertools import combinations
import fresh_h9_beta131_pilot_20261016_review as I

ROOT=I.ROOT
TAG='fresh_h9_beta131_geometry_20261016_review'
PRIMARY='results/fresh_h9_beta131_pilot_20261016.json'
GEOMETRY='results/fresh_h9_beta131_pilot_20261016_geometry.json'


def reduced(matrix):
    rows=[list(row) for row in matrix]
    pivots=[]
    for j in range(10):
        k=next((i for i in range(len(pivots),len(rows)) if rows[i][j]),None)
        if k is None: continue
        i=len(pivots)
        rows[i],rows[k]=rows[k],rows[i]
        pivot=rows[i][j]
        rows[i]=[v/pivot for v in rows[i]]
        for k in range(len(rows)):
            if k!=i:
                factor=rows[k][j]
                if factor: rows[k]=[a-factor*b for a,b in zip(rows[k],rows[i])]
        pivots.append(j)
    return rows,pivots


def run():
    data=json.loads((ROOT/PRIMARY).read_text())
    geometry=json.loads((ROOT/GEOMETRY).read_text())
    for item in (data,geometry):
        for path,digest in item['sources'].items(): assert I.sha(path)==digest,path
    form=[I.decode(row) for row in data['form']]
    roots=[I.decode(row) for row in data['roots']]
    center=I.decode(data['strict_center'])
    covectors=[I.mv(form,r) for r in roots]
    assert len(roots)==12 and len(reduced(covectors)[1])==10
    assert I.pair(center,center,form)<0 and all(I.dot(center,c)<0 for c in covectors)
    census={name:{} for name in ('finite_vertices','cusps','hyperideal_rays','opposite_sheet_rays')}
    supports=rank_nine=feasible=0
    for support in combinations(range(12),9):
        supports+=1
        rows,pivots=reduced([covectors[i] for i in support])
        if len(pivots)!=9: continue
        rank_nine+=1
        free=next(i for i in range(10) if i not in pivots)
        vector=[I.K() for _ in range(10)]
        vector[free]=I.K(1)
        for i,j in enumerate(pivots): vector[j]=-rows[i][free]
        values=[I.dot(vector,c) for c in covectors]
        if any(v>0 for v in values):
            if any(v<0 for v in values): continue
            vector=[-v for v in vector]
            values=[-v for v in values]
        assert all(v<=0 for v in values) and any(v<0 for v in values)
        feasible+=1
        active=tuple(i for i,v in enumerate(values) if not v)
        square=I.pair(vector,vector,form)
        sheet=I.pair(vector,center,form)
        kind='hyperideal_rays' if square>0 else 'opposite_sheet_rays' if sheet>=0 else 'finite_vertices' if square<0 else 'cusps'
        assert active not in census[kind] or all(
            a*next(v for v in census[kind][active] if v)==b*next(v for v in vector if v)
            for a,b in zip(vector,census[kind][active]))
        census[kind][active]=vector
    assert supports==220
    assert {k:set(v) for k,v in census.items()}=={
        k:{tuple(s) for s in v} for k,v in geometry['complete_census'].items()}
    escape=I.decode(geometry['strict_escape']['point'])
    assert I.pair(escape,escape,form)>0 and I.pair(escape,center,form)<0
    assert all(I.dot(escape,c)<0 for c in covectors)
    paths=[PRIMARY,GEOMETRY,f'scripts/{TAG}.py',
           'scripts/fresh_h9_beta131_pilot_20261016_review.py',
           'src/coxeter_polytopes/exact/golden.py']
    return I.enc(dict(status='PASS_INDEPENDENT_EXHAUSTIVE_BETA131_CONE',sources={p:I.sha(p) for p in paths},
                      supports=supports,rank_nine_supports=rank_nine,feasible_supports=feasible,
                      counts={k:len(v) for k,v in census.items()},
                      rays={k:[dict(active=list(s),point=v) for s,v in sorted(rows.items())]
                            for k,rows in census.items()},strict_escape=escape,
                      completeness='The full-dimensional cone is pointed: all wall covectors have rank10. Every extreme ray has a rank9 active support. All220 nine-wall supports are enumerated independently of incremental double description.',
                      scope='Complete cone geometry and infinite volume of the height20 prefix only. No conclusion on the full reflection chamber.'))


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    result=run()
    path=ROOT/'results'/f'{TAG}.json'
    if args.verify: assert result==json.loads(path.read_text())
    else: path.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['counts'])
