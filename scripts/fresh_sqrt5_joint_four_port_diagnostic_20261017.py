#!/usr/bin/env python3
"""Exact four-port facet scan after the frozen three-port repair; no LP."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations
from math import comb
from pathlib import Path
import json

import fresh_sqrt5_joint_triangle_square_ports_20261017 as J

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_joint_four_port_diagnostic_20261017'
THEOREM = 'fresh_sqrt5_joint_four_port_domain_20261017'
C = J.C


def dump(x):
    return json.dumps(x, sort_keys=True, separators=(',', ':'))


def quartet_orbits(g):
    pieces = C.pieces(g)
    groups = defaultdict(list)
    for U in combinations(sorted(g), 4):
        groups[tuple(sorted(C.code(p, U) for p in pieces))].append(U)
    assert sum(map(len, groups.values())) == 715
    return groups


def run():
    joint = J.read(J.TAG)
    theorem = J.read(THEOREM)
    assert theorem['status'] == 'PASS_EXACT_SELECTED_DOTTED_PRISM_FOUR_PORT_FACET'
    local = J.read(J.LOCAL)
    t = next(x for x in local['tests'] if x['name'] == 'legacy_export_seed')
    assert t['complete_positive_simplex_forest_coverage'] and not t['open_forests']
    assert joint['face_multipliers'] == t['face_multipliers']
    alpha = {i:F(a) for i,a in joint['face_multipliers']}
    mu = {v:F(a) for v,a in joint['vertex_multipliers']}
    oldmu = {v:F(a) for v,a in t['vertex_multipliers']}
    assert len(alpha) == 2219 and all(mu[v] >= a for v,a in oldmu.items())
    model = json.loads((ROOT / f'results/{J.MODEL}.json').read_text())
    for p,h in model['prepared']['sources'].items():
        assert J.B.P.sha(p) == h, p
    pmap = dict(model['prepared']['provenance_to_face'])
    inverse = {f:j for j,f in pmap.items()}
    selected_q = sorted(set(alpha) & set(inverse))
    assert len(selected_q) == 355
    actual = J.read(J.ACTUAL)
    quad = J.read(J.B.Q.TAG)
    raw = J.read(J.B.Q.CORES)
    spectators = J.read(J.B.Q.SPECTATORS)
    forests = J.read(J.FORESTS)['vertices']
    forest_lookup = {tuple(v['forest']):v['index'] for v in forests}
    models = dict(C.models())
    finite = {C.code(g):name for name,g in models.items()}
    assert len(forest_lookup) == 6558
    graphs, orbit_data, code_to_oi, representatives = {}, {}, {}, {}

    def ensure(v):
        if v in graphs:
            return
        g = {}
        for name in forests[v]['forest']:
            offset = len(g)
            for i,row in models[name].items():
                g[offset+i] = {offset+j:m for j,m in row.items()}
        assert set(g) == set(range(13))
        graphs[v] = g
        rows = []
        for oi,(code,sets) in enumerate(sorted(quartet_orbits(g).items())):
            rows.append(dict(orbit=oi,code=json.loads(dump(code)),representative=list(sets[0]),multiplicity=len(sets)))
            code_to_oi[v,code] = oi
            representatives[v,oi] = sum(1<<i for i in sets[0])
        orbit_data[v] = rows

    blockcols = {c['index']:c for c in t['literal_block_orbits']}
    blockweights = defaultdict(F)
    sfeatures_by_vertex = defaultdict(list)
    for f in t['positive_simplex_features']:
        touched = set()
        for j,m in f['block_orbits']:
            blockweights[j] += m * alpha[f['face']]
            touched.add(blockcols[j]['vertex'])
        for v in touched:
            sfeatures_by_vertex[v].append(f)
    assert set(blockweights) == set(blockcols)
    assert len(blockcols) == 3485 and len(t['positive_simplex_features']) == 1864
    block_coeff, sload = {}, defaultdict(F)
    simplex_identities = []
    for serial,(j,c) in enumerate(sorted(blockcols.items())):
        v = c['vertex']
        ensure(v)
        masks = c['literal_masks']
        assert masks and len(masks) == len(set(masks))
        assert all(0 < m < 1<<13 and m.bit_count() == c['block_size'] for m in masks)
        for row in orbit_data[v]:
            oi = row['orbit']
            u = representatives[v,oi]
            value = F(sum(max(0,(m&u).bit_count()-1) for m in masks),len(masks))
            if value:
                block_coeff[j,oi] = value
                sload[v,oi] += blockweights[j] * value
        k = c['block_size']
        total = sum((block_coeff.get((j,r['orbit']),F()) * r['multiplicity'] for r in orbit_data[v]),F())
        expected = sum((a-1)*comb(k,a)*comb(13-k,4-a) for a in range(2,min(k,4)+1))
        assert total == expected
        simplex_identities.append([j,k,str(total)])
        if serial % 500 == 0:
            print('FOUR_PORT_SIMPLEX_ORBITS',serial+1,'/',len(blockcols),'forests',len(graphs),flush=True)

    qload, qcoeff = defaultdict(F), defaultdict(Counter)
    qrecords, quad_identities = [], []
    oldq = {r['face']:r for r in joint['actual_quad_features']}
    assert set(oldq) == set(selected_q)
    for serial,fi in enumerate(selected_q):
        j = inverse[fi]
        rec = actual['provenance_masses'][j]
        ci = rec['canonical_prism_core']
        core = quad['canonical_cores'][ci]
        g0,_ = J.B.Q.graph(raw['records'][core['source_core']])
        g = {i:{k:m[0] for k,m in row.items()} for i,row in g0.items()}
        sp = spectators['contexts'][rec['source_context']]['spectator_forest']
        for name in sp:
            offset = len(g)
            for i,row in models[name].items():
                g[offset+i] = {offset+k:m for k,m in row.items()}
        assert len(g) == 15
        corners = J.X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
        assert len(corners) == 2*rec['prism_dimension']
        previous = oldq[fi]
        assert (previous['provenance'],previous['canonical_core'],previous['spectator_forest']) == (j,ci,sp)
        corner_records = []
        for c,oldc in zip(corners,previous['corners']):
            omitted = set(c['vertex_omitted'])
            h = {i:{k:m for k,m in row.items() if k not in omitted} for i,row in g.items() if i not in omitted}
            assert len(h) == 13 and all(m>0 for row in h.values() for m in row.values())
            pieces = C.pieces(h)
            forest = tuple(sorted(finite[C.code(p)] for p in pieces))
            v = forest_lookup[forest]
            ensure(v)
            edges = set(map(tuple,c['occupied_pairs']))
            e = len(edges)
            assert e == rec['prism_dimension']-1
            assert oldc['omitted'] == sorted(omitted) and oldc['vertex'] == v
            assert oldc['occupied_pairs'] == c['occupied_pairs']
            hubs = set.intersection(*(set(edge) for edge in edges))
            assert hubs
            hits = Counter()
            hitsets = []
            # All715 subsets are included: those with fewer than two occupied
            # pairs have coefficient zero, and every other subset is encoded.
            for U in combinations(sorted(h),4):
                if sum(set(edge) <= set(U) for edge in edges) < 2:
                    continue
                hitsets.append(U)
                code = tuple(sorted(C.code(p,U) for p in pieces))
                oi = code_to_oi[v,code]
                hits[oi] += 1
            expected = comb(e,2)*(12-e)+comb(e,3)
            assert len(hitsets) == expected
            coeff = {}
            for oi,hit in sorted(hits.items()):
                assert 0 < hit <= orbit_data[v][oi]['multiplicity']
                value = F(hit,orbit_data[v][oi]['multiplicity'])
                coeff[oi] = value
                qcoeff[v,oi][fi] += value
                qload[v,oi] += alpha[fi] * value
            total = sum((value*orbit_data[v][oi]['multiplicity'] for oi,value in coeff.items()),F())
            assert total == expected
            quad_identities.append([fi,sorted(omitted),e,str(total)])
            corner_records.append(dict(omitted=sorted(omitted),vertex=v,occupied_pairs=c['occupied_pairs'],
                four_port_orbit_coefficients=[[oi,str(value)] for oi,value in coeff.items()]))
        qrecords.append(dict(face=fi,provenance=j,canonical_core=ci,prism_dimension=rec['prism_dimension'],
            spectator_forest=sp,corners=corner_records))
        if serial % 50 == 0:
            print('FOUR_PORT_ACTUAL_QUADS',serial+1,'/',len(selected_q),flush=True)

    updates, checks, critical = dict(mu), [], []
    for v,rows in sorted(orbit_data.items()):
        for row in rows:
            oi = row['orbit']
            s,q = sload[v,oi],qload[v,oi]
            load = s+q
            need = load/3
            checks.append([v,oi,str(s),str(q),str(load-3*mu.get(v,F()))])
            updates[v] = max(updates.get(v,F()),need)
            if need > mu.get(v,F()):
                terms = Counter(qcoeff[v,oi])
                for f in sfeatures_by_vertex[v]:
                    value = sum((m*block_coeff.get((j,oi),F()) for j,m in f['block_orbits']
                                 if blockcols[j]['vertex']==v),F())
                    if value:
                        terms[f['face']] += value
                assert sum((alpha[fi]*a for fi,a in terms.items()),F()) == load
                critical.append(dict(index=len(critical),vertex=v,four_port_orbit=oi,
                    representative=row['representative'],orbit_multiplicity=row['multiplicity'],vertex_cost=3,
                    positive_face_coefficients=[[fi,str(a)] for fi,a in sorted(terms.items())],
                    simplex_load=str(s),quad_load=str(q),load=str(load),old_vertex_mass=str(mu.get(v,F()))))
    numerator = F(t['objective_numerator'])
    assert numerator == F(joint['objective_numerator'])
    den = sum(updates.values(),F())
    assert all(updates[v] >= a for v,a in mu.items())
    assert set(graphs) == set(mu) and len(graphs)==593
    for v,oi,s,q,_ in checks:
        assert F(s)+F(q) <= 3*updates[v]
    paths = [f'scripts/{TAG}.py',f'scripts/{THEOREM}.py',f'results/{THEOREM}.json',f'results/{THEOREM}.md',
        f'scripts/{J.TAG}.py',f'results/{J.TAG}.json',f'results/{J.LOCAL}.json',
        f'results/{J.MODEL}.json',f'results/{J.ACTUAL}.json',f'results/{J.FORESTS}.json',
        f'scripts/{J.X.TAG}.py',f'results/{J.B.Q.TAG}.json',f'results/{J.B.Q.CORES}.json',f'results/{J.B.Q.SPECTATORS}.json']
    return dict(status='PASS_EXACT_FOUR_PORT_FACET_AFTER_STRENGTHENED_THREE_PORT_REPAIR',
        sources={p:J.B.P.sha(p) for p in paths},dimension=13,target=12,optimizer_invocations=0,
        initial_bound=joint['bound'],bound=str(numerator/den),objective_numerator=str(numerator),normalizer=str(den),
        counts=dict(simplex_faces=1864,simplex_block_orbits=3485,actual_quad_provenances=len(qrecords),
            quad_vertex_occurrences=sum(len(r['corners']) for r in qrecords),tested_forests=len(graphs),
            literal_four_port_subsets=715*len(graphs),four_port_orbits=len(checks),
            violated_columns=len(critical),increased_vertex_masses=sum(updates[v]>mu.get(v,F()) for v in updates)),
        finite_forest_four_port_orbits=[dict(vertex=v,forest=forests[v]['forest'],orbits=rows)
            for v,rows in sorted(orbit_data.items())],actual_quad_features=qrecords,
        simplex_coefficient_sum_identities=simplex_identities,quad_coefficient_sum_identities=quad_identities,
        all_orbit_loads_and_slacks=checks,critical_columns=critical,
        face_multipliers=t['face_multipliers'],initial_vertex_multipliers=joint['vertex_multipliers'],
        vertex_multipliers=[[v,str(a)] for v,a in sorted(updates.items())],
        scope='Same1864 simplex rows and355 actual selected dotted-Q provenances as the frozen restricted whole-local-cone and strengthened three-port diagnostics. All715 four-subsets are covered for each contributing forest, using exact whole-forest automorphism averages. The face/provenance masses are unchanged and mu is raised by a per-vertex maximum, preserving every inherited positive-cost/signed flag condition, full local-simplex cover and strengthened triple-port condition. Unexported face coefficients can be set to zero as valid weaker nonnegative payments. No optimizer, general non-dotted-Q application, complete joint packing cone, ambient realization or dimension theorem.')


if __name__ == '__main__':
    p=argparse.ArgumentParser()
    p.add_argument('--verify',action='store_true')
    a=p.parse_args()
    result=json.loads(dump(run()))
    out=ROOT/f'results/{TAG}.json'
    if a.verify:
        assert result==json.loads(out.read_text())
    else:
        out.write_text(dump(result)+'\n')
    print(result['status'],result['counts'],result['bound'],float(F(result['bound'])),flush=True)
