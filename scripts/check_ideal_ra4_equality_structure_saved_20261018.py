#!/usr/bin/env python3
"""Independent finite-support checks; the universal equality proof is prose."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
data=json.loads((ROOT/'results/fresh_ideal_ra4_equality_structure_20261018.json').read_text())
assert data['source_sha256']==sha256((ROOT/'scripts/fresh_ideal_ra4_equality_structure_20261018.py').read_bytes()).hexdigest()
for path,expected in data['sources'].items():assert sha256((ROOT/path).read_bytes()).hexdigest()==expected
assert data['status']=='EXACT_FINITE_SUPPORTS_FOR_WRITTEN_EQUALITY_PROOF'
assert data['scope']=={'historical_classification_replayed':False,'four_polytope_census':False,
                      'universal_geometry_proved_by_this_computation':False}

def elimination(rows,rhs=None):
    # Reverse columns and bottom-up pivot search, independent of producer.
    width=len(rows[0]);a=[[F(x) for x in row] for row in rows]
    if rhs is not None:
        for row,b in zip(a,rhs):row.append(F(b))
    pivots=[];r=0
    for c in reversed(range(width)):
        p=next((i for i in reversed(range(r,len(a))) if a[i][c]),None)
        if p is None:continue
        a[r],a[p]=a[p],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        pivots.append(c);r+=1
        if r==len(a):break
    if rhs is None:return len(pivots)
    if len(pivots)!=width:return None
    out=[F(0)]*width
    for i,c in enumerate(pivots):out[c]=a[i][-1]
    return tuple(out)

def inner(a,b,metric):return sum(F(m)*F(x)*F(y) for m,x,y in zip(metric,a,b))

def boundary(c):return {tuple(sorted((c[i],c[(i+1)%len(c)]))) for i in range(len(c))}

t=data['ten_vertex'];xy={v:tuple(p) for v,p in t['vertices'].items()}
edges={tuple(e) for e in t['edges']};faces=t['faces'];quads=t['quadrilaterals']
assert len(xy)==10 and len(edges)==20 and len(faces)==12 and len({tuple(f) for f in faces})==12
assert Counter(x for e in edges for x in e)==Counter({v:4 for v in xy})
assert Counter(map(len,faces))=={3:8,4:4}
assert Counter(e for f in faces for e in boundary(f))==Counter({e:2 for e in edges})
def cross(a,b,c):return (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
for (a,b),(c,d) in combinations(edges,2):
    if {a,b}&{c,d}:continue
    assert not(cross(xy[a],xy[b],xy[c])*cross(xy[a],xy[b],xy[d])<0 and
               cross(xy[c],xy[d],xy[a])*cross(xy[c],xy[d],xy[b])<0)
def area(c):return abs(sum(xy[c[i]][0]*xy[c[(i+1)%len(c)]][1]-xy[c[i]][1]*xy[c[(i+1)%len(c)]][0] for i in range(len(c))))
outer=list('ABCD');assert outer in faces
assert sum(area(f) for f in faces if f!=outer)==area(outer)
assert {tuple(q) for q in quads}=={tuple(f) for f in faces if len(f)==4}
matching=[list(p) for p in combinations(range(4),2) if boundary(quads[p[0]])&boundary(quads[p[1]])]
assert matching==t['quadrilateral_edge_adjacency'] and sorted(x for p in matching for x in p)==list(range(4))
assert t['historical_classification_replayed'] is False

old=json.loads((ROOT/'results/fresh_ideal_ra4_double24_volume_20261018.json').read_text())['double']
walls=[[F(x) for x in a] for a in old['walls']];vertices=[[F(x) for x in a] for a in old['vertices']]
metric=[-1,1,1,1,1];facet=[{i for i,v in enumerate(vertices) if inner(v,n,metric)==0} for n in walls]
ridges=[set(r) for r in old['faces']['2']]
local=data['double_oct_local'];assert len(local)==8 and len({r['facet'] for r in local})==8
assert {r['facet'] for r in local}=={i for i,f in enumerate(facet) if len(f)==9}
for row in local:
    i=row['facet'];ns=set(row['nonsquare_vertices']);f=facet[i]
    qs=[r for r in ridges if r<=f and len(r)==4];ts=[r for r in ridges if r<=f and len(r)==3]
    assert {frozenset(q) for q in row['quads']}==set(map(frozenset,qs))
    assert len(ns)==3 and len(qs)==3 and len(ts)==8
    assert all(len(q&ns)==2 for q in qs) and all(len(t&ns)<=1 for t in ts)
    assert all(len(q&r)==1 and q&r<=ns for q,r in combinations(qs,2))
    assert row['triangle_nonsquare_intersections']==sorted(len(t&ns) for t in ts)
    assert {w['vertex'] for w in row['widths']}==ns and len(row['widths'])==3
    for vi in f:
        qn=[j for j,g in enumerate(facet) if j!=i and f&g in qs and vi in g]
        tn=[j for j,g in enumerate(facet) if j!=i and f&g in ts and vi in g]
        assert len(qn)==(2 if vi in ns else 1)
        if vi in ns:
            assert len(tn)==2 and inner(walls[qn[0]],walls[qn[1]],metric)==-1
            assert inner(walls[tn[0]],walls[tn[1]],metric)==-1
            w=next(w for w in row['widths'] if w['vertex']==vi)
            assert F(w['quad_pair_width'])==sum(walls[j][0] for j in qn)
            assert F(w['triangle_pair_width'])==sum(walls[j][0] for j in tn)==2*F(w['quad_pair_width'])

c=data['cuboctahedron'];m=list(map(F,c['scaled_metric']))
assert m==[-1,F(1,2),F(1,2),F(1,2)]
ns=[tuple(map(F,n)) for n in c['walls']];vs=[tuple(map(F,v)) for v in c['vertices']]
expected_ns={(F(1),*p) for p in product((-1,1),repeat=3)}
for axis in range(3):
    for sign in (-1,1):
        p=[F(0)]*3;p[axis]=2*sign;expected_ns.add((F(1),*p))
assert set(ns)==expected_ns
expected_vs={(F(1),*p) for p in product((-1,0,1),repeat=3) if sum(x*x for x in p)==2}
assert set(vs)==expected_vs and len(vs)==12 and len(set(ns))==len(ns)==14
assert all(inner(v,v,m)==0 for v in vs)
assert all(inner(v,n,m)<=0 for v in vs for n in ns)
assert all(inner(n,k,m)<=0 for n,k in combinations(ns,2))
found=set()
for chosen in combinations(ns,3):
    x=elimination([[n[i]*m[i] for i in range(1,4)] for n in chosen],[1]*3)
    if x is not None and all(inner((1,*x),n,m)<=0 for n in ns):found.add((F(1),*x))
assert found==expected_vs and c['all_active_triples']==364
assert c['facet_sizes']==[sum(inner(v,n,m)==0 for v in vs) for n in ns]==[4]*6+[3]*8
assert len(c['cusps'])==12 and {r['vertex'] for r in c['cusps']}==set(range(12))
for row in c['cusps']:
    v=vs[row['vertex']];ids=[i for i,n in enumerate(ns) if inner(v,n,m)==0]
    pairs=[list(p) for p in combinations(ids,2) if inner(ns[p[0]],ns[p[1]],m)<0]
    assert pairs==row['opposite_pairs'] and len(pairs)==2 and sorted(i for p in pairs for i in p)==ids
    widths=[]
    for a,b in pairs:
        square=inner(ns[a],ns[a],m);assert square==inner(ns[b],ns[b],m)==-inner(ns[a],ns[b],m)
        total=tuple(x+y for x,y in zip(ns[a],ns[b]));assert total==tuple(total[0]*x for x in v)
        widths.append(total[0]**2/square)
    assert widths==list(map(F,row['width_squares'])) and max(widths)==2*min(widths)
assert c['all_cusp_aspect_squared']=='2'

assert data['remaining_defect_partitions']==[[3]*6+[6],[3]*5+[4,5],[3]*4+[4,4,4],[3]*8]
graphs=data['six_cusp_graphs'];assert len(graphs)==15
domain={tuple(p) for p in combinations(list(combinations(range(6),2)),3) if sorted(x for e in p for x in e)==list(range(6))}
assert {tuple(map(tuple,g['missing_matching'])) for g in graphs}==domain
for g in graphs:
    miss=set(map(tuple,g['missing_matching']));edges=set(combinations(range(6),2))-miss
    tris=[t for t in combinations(range(6),3) if set(combinations(t,2))<=edges]
    assert list(map(list,tris))==g['triangles'] and len(tris)==8
    dual=[list(p) for p in combinations(range(8),2) if len(set(tris[p[0]])&set(tris[p[1]]))==2]
    assert dual==g['dual_edges'] and len(dual)==12
    signatures=g['cube_signatures'];assert {tuple(s) for s in signatures}==set(product((0,1),repeat=3))
    for t,s in zip(tris,signatures):assert all(g['missing_matching'][i][s[i]] in t for i in range(3))
    assert all(sum(a!=b for a,b in zip(signatures[i],signatures[j]))==1 for i,j in dual)

a=data['cube_algebra'];ss=list(product((-1,1),repeat=3));assert a['signs']==list(map(list,ss))
equations=[]
for i in range(3):
    j,k=[v for v in range(3) if v!=i]
    for sign in (-1,1):equations.append([s[j]*s[k] if s[i]==sign else 0 for s in ss])
assert equations==a['square_equations'] and elimination(equations)==a['equation_rank']==4
modes=[[1]*8]+[[s[i] for s in ss] for i in range(3)]
assert a['affine_modes']==modes and elimination(modes)==4
assert all(sum(x*y for x,y in zip(e,f))==0 for e in equations for f in modes)
gram=[[F(-1+sum(x*y for x,y in zip(s,t)),2) for t in ss] for s in ss]
assert gram==[list(map(F,r)) for r in a['normal_gram']] and elimination(gram)==4
assert a['orthogonal_affine_mode_squares']==['-1/2','1/2','1/2','1/2']
print(json.dumps({'status':'PASS','scope':'independent finite supports; no universal geometry or classification replay',
  'ten_vertex_faces':12,'double_oct_facets':8,'cuboctahedron_active_triples':364,
  'six_cusp_graphs':15,'cube_square_rank':4}))
