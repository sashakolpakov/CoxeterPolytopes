#!/usr/bin/env python3
"""Standalone necessary-star replay; no inference of a geometric RA5 realization."""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path
import json

R=Path(__file__).resolve().parents[1]
D=json.loads((R/'results/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.json').read_text())
assert D['source_sha256']==sha256((R/'scripts/fresh_ideal_ra5_known_carrier_cusp_compatibility_20261018.py').read_bytes()).hexdigest()
for path,h in D['sources'].items():assert sha256((R/path).read_bytes()).hexdigest()==h


def bil(x,y):
    return sum((-1 if i==0 else 1)*x[i]*y[i] for i in range(5))


base=json.loads((R/'results/fresh_ideal_ra4_double24_volume_20261018.json').read_text())
double=json.loads((R/'results/fresh_ideal_ra4_positive_kappa_doubling_20261018.json').read_text())
inputs={1:base['original'],2:double['stages'][0],16:double['stages'][-1]}
assert set(D['facet_types'])=={'1','2','16'}
records={}
for typ,inp in inputs.items():
    ns=[tuple(map(Q,x)) for x in inp['walls']];vs=[tuple(map(Q,x)) for x in inp['vertices']]
    assert all(bil(n,n)==1 for n in ns) and all(v[0]==1 and bil(v,v)==0 for v in vs)
    sizes=[sum(bil(n,v)==0 for v in vs) for n in ns]
    got=D['facet_types'][str(typ)]
    assert got['vertices']==len(vs) and got['facets']==len(ns) and len(got['records'])==len(vs)
    assert got['ridge_size_histogram']=={str(k):v for k,v in sorted(Counter(sizes).items())}
    shapes=Counter(); rr=[]
    for vi,v in enumerate(vs):
        ids=[i for i,n in enumerate(ns) if bil(n,v)==0]
        pairs=sorted((i,j) for i in ids for j in ids if i<j and bil(ns[i],ns[j])==-1)
        assert len(ids)==6 and len(pairs)==3 and sorted(x for p in pairs for x in p)==ids
        axes=[]
        for i,j in pairs:
            width=ns[i][0]+ns[j][0];assert width>0
            assert all(ns[i][k]+ns[j][k]==width*v[k] for k in range(5))
            axes.append({'width':str(width),'ridge_vertices':[sizes[i],sizes[j]],'facet_indices':[i,j]})
        record={'vertex':vi,'axes':axes};assert got['records'][vi]==record;rr.append(record)
        ww=sorted(Q(a['width']) for a in axes);shapes[tuple(str(w/ww[0]) for w in ww)]+=1
    assert got['cusp_shape_histogram']==[{'shape':list(k),'count':v} for k,v in sorted(shapes.items())]
    records[typ]=rr

nodes=tuple(product(range(4),range(2)))
adj={n:tuple(m for m in nodes if m[0]!=n[0]) for n in nodes}
indices={n:{m:j for j,m in enumerate(adj[n])} for n in nodes}
special=[r for r in records[16] if sorted(x for a in r['axes'] for x in a['ridge_vertices'])==[14,14,24,24,24,24]]
assert len(special)==1 and special[0]['vertex']==D['special_vertex']==85
assert sorted((Q(a['width']),sorted(a['ridge_vertices'])) for a in special[0]['axes'])==[(4,[24,24]),(4,[24,24]),(8,[14,14])]
for r in records[16]:
    ww=sorted(Q(a['width']) for a in r['axes'])
    if tuple(w/ww[0] for w in ww)==(1,2,2):
        assert all(a['ridge_vertices']==[9,9] for a in r['axes'] if Q(a['width'])==ww[-1])


def solve(dom):
    # Arc-consistency elimination followed by complete reverse-order branching.
    while True:
        changed=False
        for n in reversed(nodes):
            keep={p for p in dom[n] if all(any(p[indices[n][m]]==q[indices[m][n]] for q in dom[m]) for m in adj[n])}
            if not keep:return False
            if keep!=dom[n]:changed=True;dom[n]=keep
        if not changed:break
    if all(len(dom[n])==1 for n in nodes):return True
    n=min((n for n in reversed(nodes) if len(dom[n])>1),key=lambda n:len(dom[n]))
    for p in sorted(dom[n],reverse=True):
        if solve({m:({p} if m==n else set(dom[m])) for m in nodes}):return True
    return False


assert len(D['local_cases'])==4
assert {(r['new_width'],r['opposite_type']) for r in D['local_cases']}==set(product((4,8),(2,16)))
for row in D['local_cases']:
    widths=[row['new_width'],4,4,8];dom={}
    for node in nodes:
        typ=row['opposite_type'] if node==(0,1) else 16
        remain=[a for a in range(4) if a!=node[0]];vals=set()
        for record in records[typ]:
            for mapping in permutations(remain):
                axes=record['axes'];ratio={Q(q['width'])/widths[a] for a,q in zip(mapping,axes)}
                if len(ratio)!=1:continue
                for signs in product((0,1),repeat=3):
                    by_neighbor={}
                    for a,q,sgn in zip(mapping,axes,signs):
                        by_neighbor[(a,0)]=q['ridge_vertices'][sgn]
                        by_neighbor[(a,1)]=q['ridge_vertices'][1-sgn]
                    vals.add(tuple(by_neighbor[m] for m in adj[node]))
        dom[node]=vals
    dom[(0,0)]={(24,24,24,24,14,14)}
    assert len(row['domains'])==8 and {tuple(x['node']) for x in row['domains']}==set(nodes)
    for saved in row['domains']:
        assert len(saved['labels'])==len({tuple(x) for x in saved['labels']})
        assert dom[tuple(saved['node'])]=={tuple(x) for x in saved['labels']}
    compatible=solve({n:set(x) for n,x in dom.items()})
    assert compatible==row['compatible']==(row['new_width']==4 and row['opposite_type']==16)
    if compatible:
        assert len(row['witness'])==8
        w={tuple(x['node']):tuple(x['labels']) for x in row['witness']};assert set(w)==set(nodes)
        assert all(w[n] in dom[n] for n in nodes)
        assert all(w[n][indices[n][m]]==w[m][indices[m][n]] for n in nodes for m in adj[n])
    else:assert row['witness'] is None

assert D['necessary_counts']=={'Q16_multiple':6,'Q16_minimum':30,'double24_minimum':66}
# Universal nonnegative counts: below30, positive multiples6 violate the 22-neighbor/rank bound.
assert all(11*z>Q(2,5)*z*z for z in (6,12,18,24))
w=D['abstract_incidence_witness'];x,y,z=w['facet_counts'];a,b,c,d=w['ambient_cusp_types']
assert (x,y,z,a,b,c,d)==(361,66,30,2051,241,105,0)
assert 6*z-2*y==48 and 8*a+2*b+2*d==24*x+36*y+195*z
assert 6*b+4*c==6*y+49*z and 4*c+6*d==14*z
assert x+y+z==w['facets']==457 and a+b+c+d==w['vertices']==2397
assert 24*x+38*y+198*z==2*(3*w['vertices']+3*w['facets']-6)
assert w['is_geometric_realization'] is False
assert D['global_claims']=={'ideal_RA5_example':False,'spectrum_excluded':False,'unrestricted_RA5_excluded':False}
print(json.dumps({'status':'PASS','stars':324,'independent_domain_cases':4,
                  'compatible_cases':1,'Q16_minimum':30,'double24_minimum':66,
                  'scope':'saved necessary-star and arithmetic replay, not geometry realization'}))
