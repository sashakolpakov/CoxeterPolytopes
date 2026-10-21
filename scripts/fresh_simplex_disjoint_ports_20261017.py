#!/usr/bin/env python3
"""Export one-port simplex charges and test the preceding exact sparse dual.

The general geometric proof is in the companion report. This is a single
exact arithmetic implementation, with no optimizer or new face coverage.
"""
import argparse
import hashlib
import json
from collections import Counter, defaultdict
from fractions import Fraction as F
from pathlib import Path
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_simplex_disjoint_ports_20261017'
LP = 'fresh_sqrt5_rank13_simplex_packing_repair_20261017'
STARS = 'fresh_sqrt5_rank13_full_triangle_simplex_separation_20261017'


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def read(tag):
    doc = json.loads((ROOT / f'results/{tag}.json').read_text())
    for path, h in doc.get('sources',{}).items():
        assert sha(path) == h, (tag,path)
    return doc


def code(pieces, marks=()):
    return tuple(sorted(C.code(p,marks) for p in pieces))


def forest_graph(forest, models):
    g = {}
    for name in forest:
        offset = len(g)
        for i,row in models[name].items():
            g[offset+i] = {offset+j:m for j,m in row.items()}
    assert len(g) == 13
    return g


def run():
    old = read(LP)
    stars = read(STARS)['stars']
    dual = next(c for c in old['certificates'] if c['kind']=='EXACT_RESTRICTED_MODEL_DUAL')
    alpha = {i:F(q) for i,q in dual['face_multipliers']}
    mu = {i:F(q) for i,q in dual['vertex_multipliers']}
    models = dict(C.models())
    finite = {C.code(g):name for name,g in models.items() if len(g)<=13}
    vertices = {tuple(v['forest']):v['index'] for v in old['vertices']}
    columns = []
    lookup = {}
    for v in old['vertices']:
        g = forest_graph(v['forest'],models)
        pieces = C.pieces(g)
        orbits = defaultdict(list)
        for port in g:
            orbits[code(pieces,(port,))].append(port)
        assert sum(map(len,orbits.values())) == 13
        for marked, ports in sorted(orbits.items()):
            j = len(columns)
            columns.append(dict(index=j,vertex=v['index'],forest=v['forest'],
                                marked_forest_code=marked,literal_port_orbit=ports,
                                vertex_cost=len(ports),face_coefficients=[]))
            lookup[v['index'],marked] = j
    loads = defaultdict(F)
    face_types = Counter()
    for face in old['faces']:
        fi = face['index'] if 'index' in face else None
        # The authoritative domain stores rows in index order.
        if fi is None:
            raise AssertionError('explicit row index required')
        assert old['faces'][fi] == face and len(face['provenance']) == 1
        si = face['provenance'][0]['star']
        star = stars[si]
        S = {r['omitted'] for r in star['spherical_deletions']}
        d = star['simplex_dimension']
        assert len(S) == d+1
        g = C.graph(range(14),{(i,j):m for i,j,m in star['edges']})
        coefficients = Counter()
        for row in star['spherical_deletions']:
            o = row['omitted']
            h = {i:{j:m for j,m in rr.items() if j!=o} for i,rr in g.items() if i!=o}
            pieces = C.pieces(h)
            forest = tuple(sorted(finite[C.code(p)] for p in pieces))
            assert forest == tuple(row['forest'])
            v = vertices[forest]
            for port in sorted(S-{o}):
                coefficients[lookup[v,code(pieces,(port,))]] += 1
        assert sum(coefficients.values()) == d*(d+1)
        for j,n in sorted(coefficients.items()):
            columns[j]['face_coefficients'].append([fi,n])
            loads[j] += alpha.get(fi,F()) * n
        face_types[d] += 1
    numerator = sum((a*old['faces'][i]['required'] for i,a in alpha.items()),F())
    normalizer = sum(mu.values(),F())
    assert numerator == F(dual['objective_numerator'])
    assert normalizer == F(dual['normalizer'])
    assert numerator/normalizer == F(dual['bound'])
    summaries = []
    best = {}
    for col in columns:
        j,v,cost = col['index'],col['vertex'],col['vertex_cost']
        load = loads[j]
        floor = load/cost
        if floor > best.get(v,(F(),None))[0]:
            best[v] = floor,j
        if floor > mu.get(v,F()):
            summaries.append(dict(column=j,vertex=v,forest=col['forest'],
                                  literal_port_orbit=col['literal_port_orbit'],
                                  vertex_mass=str(mu.get(v,F())),
                                  total_port_load=str(load),orbit_size=cost,
                                  required_vertex_mass=str(floor),
                                  ratio=str(floor/mu[v]) if mu.get(v) else None,
                                  positive_dual_rows=[[i,n,str(alpha[i])] for i,n in col['face_coefficients'] if i in alpha]))
    extra = sum((max(F(),x-mu.get(v,F())) for v,(x,j) in best.items()),F())
    witness = min(summaries,key=lambda r:(len(r['positive_dual_rows']),r['vertex'],r['column'])) if summaries else None
    paths = [f'scripts/{TAG}.py',f'scripts/{C.TAG}.py',f'results/{LP}.json',f'results/{STARS}.json']
    return dict(status='PASS_EXACT_DISJOINT_SIMPLEX_PORT_EXPORT_AND_DUAL_SEPARATION',
        sources={p:sha(p) for p in paths},
        counts=dict(vertices=len(old['vertices']),faces=len(old['faces']),columns=len(columns),
                    nonzero_columns=sum(bool(c['face_coefficients']) for c in columns),
                    nonzero_coefficients=sum(len(c['face_coefficients']) for c in columns),
                    violating_columns=len(summaries),violating_forests=len({r['vertex'] for r in summaries})),
        simplex_dimension_counts=dict(sorted(face_types.items())),
        old_dual_bound=str(numerator/normalizer),old_dual_numerator=str(numerator),old_vertex_mass=str(normalizer),
        additional_vertex_mass_with_fixed_face_multipliers=str(extra),
        old_demand_over_necessary_new_vertex_mass=str(numerator/(normalizer+extra)),
        minimal_violation=witness,violations=summaries,successor_columns=columns,
        theorem='At a vertex of a simple convex polytope, maximal simplex faces of dimension at least2 use pairwise disjoint port sets. Each port has total maximal-simplex incidence at most1. The exported orbit column pays the number of its ports used at each simplex corner and has vertex cost equal to the orbit size.',
        scope='All10185 old flexible simplex rows and1870 old full forests, with exact one-port automorphism orbits. No optimizer, no full face coverage and no dimension exclusion. New charges pay actual maximal-simplex flags only, so quadrilateral rows receive zero. The fixed-alpha mass correction rejects a saved dual distribution when positive; it is not an augmented feasibility certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--verify',action='store_true');a=p.parse_args()
    d=json.loads(json.dumps(run()));out=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert d==json.loads(out.read_text())
    else:
        out.write_text(json.dumps(d,indent=2)+'\n')
    print(d['status'],d['counts'],d['old_demand_over_necessary_new_vertex_mass'],flush=True)
