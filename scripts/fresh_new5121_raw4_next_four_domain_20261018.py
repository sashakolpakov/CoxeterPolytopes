#!/usr/bin/env python3
"""Reconstruct the next four positive raw4 classes and bound their local benefit.

This is an exact necessary corner domain, not a new packing bound or chamber.
The current324-class compatibility graphs and multipliers are not changed.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_raw125_raw205_actual_corner_frontier_20261017 as A
import fresh_new5121_raw4_two_face_local_covers_20261018 as Z

ROOT = A.ROOT
R, M, V, K, C, X, B, P = A.R, A.M, A.V, A.K, A.C, A.X, A.B, A.P
TAG = 'fresh_new5121_raw4_next_four_domain_20261018'
BASE = Z.TAG
MIGRATION = 'fresh_cross_vertex_support_migration_20261018'
CURRENT = 'fresh_cross_vertex_historical32_local_covers_20261018'
FACES = (14935, 14368, 14684, 14670)


def normalize(x):
    return json.loads(json.dumps(x, sort_keys=True))


def run():
    inputs = {}
    def read(stem):
        if stem not in inputs:
            inputs[stem] = json.loads((ROOT/f'results/{stem}.json').read_text())
        return inputs[stem]

    migration, current, baseline = read(MIGRATION), read(CURRENT), read(BASE)
    assert baseline["represented_classes"] == 324
    represented = set(migration['current290_classes']) | {
        r['selected_face']['face'] for r in current['added_faces']} | {r['face'] for r in baseline['selected_faces']}
    assert len(represented) == 324
    missing = [r for r in migration['quad_rows'] if r['face'] not in represented]
    assert len(missing) == 1954
    missing_types = Counter(r['provenance_kind'] for r in missing)
    assert missing_types == dict(explicit=286, unique_saved_context=1210, merged_unexpanded=458)
    family = sorted((r for r in missing if r.get('raw_family') == 4),
                    key=lambda r: (-F(r['demand_mass']), r['face']))
    assert len(family) == 89 and tuple(r['face'] for r in family[:4]) == FACES
    selected = family[:4]
    alpha = {i: F(a) for i, a in baseline['face_multipliers']}

    raw = read(M.TAGS['raw'])
    quads = read(M.TAGS['quad'])
    spectators = read(B.Q.SPECTATORS)
    forests = read(P.FORESTS)
    models = dict(C.models())
    finite = {C.code(g): name for name, g in models.items()}
    lookup = {tuple(v['forest']): v['index'] for v in forests['vertices']}
    names = {v['index']: v['forest'] for v in forests['vertices']}
    core = quads['canonical_cores'][4]
    assert (core['index'], core['source_core'], core['nodes'], core['required']) == (4, 4, 5, 3)
    graph = V.raw_graph(raw['records'][4])
    assert len(graph) == 5 and all(e[0] in (0, 2, 3, 5) for row in graph for e in row)
    corners = X.quad_features(dict(quadrilateral_black_partitions=core['quad_partitions']))
    assert len(corners) == len({tuple(c['vertex_omitted']) for c in corners}) == 6
    targets, inverses = {}, {}
    occurrences, stars = defaultdict(list), defaultdict(dict)
    face_records = []
    maps_total = orientations_total = 0
    for q in selected:
        fi = q['face']
        assert q['provenance_kind'] == 'unique_saved_context' and q['provenance'] is None
        assert q['canonical_core'] == 4 and q['required'] == 3
        assert alpha[fi] == F(q['mass'])
        face = forests['faces'][fi]
        assert face['kind'] == 'canonical_dotted_quad_aggregate'
        assert face['source_contexts'] == [[q['source_context'], 4]]
        ctx = spectators['contexts'][q['source_context']]
        assert quads['source_core_maps'][ctx['core']]['canonical_core'] == 4
        spectator = ctx['spectator_forest']
        assert sum(len(models[n]) for n in spectator) == 10
        records = []
        for corner in corners:
            omitted = tuple(corner['vertex_omitted'])
            ids = [i for i in range(5) if i not in omitted]
            h = {i: {j: graph[i][j][0] for j in ids
                     if i != j and graph[i][j] != A.RIGHT} for i in ids}
            assert all(m in (3, 5) for row in h.values() for m in row.values())
            vf = tuple(sorted(spectator + [finite[C.code(a)] for a in C.pieces(h)]))
            vertex = lookup[vf]
            if vertex not in targets:
                target = P.forest_graph(vf, models)
                assert len(target) == 13
                targets[vertex] = target
                inverses[vertex] = R.projection_data(target)
            maps = R.active_embeddings(h, targets[vertex])
            orientations = {}
            for mp in maps:
                reverse = {j: i for i, j in mp.items()}
                full = tuple(sorted(tuple(graph[x][reverse[j]] if j in reverse else A.RIGHT
                                          for j in range(13)) for x in omitted))
                assert full[0] != full[1]
                occupied = tuple(sorted(tuple(sorted(mp[z] for z in pair))
                                        for pair in corner['occupied_pairs']))
                assert len(occupied) == len(set(occupied)) == 2
                if full in orientations:
                    assert orientations[full][0] == occupied
                else:
                    orientations[full] = (occupied, sorted(mp.items()))
            assert orientations
            rec = dict(face=fi, provenance=None, raw_family=4, omitted=omitted,
                       vertex=vertex, active_embeddings=len(maps), literal_orientations=len(orientations))
            maps_total += len(maps)
            orientations_total += len(orientations)
            J, active, components = inverses[vertex]
            indices = []
            for full, (occupied, mp) in sorted(orientations.items()):
                covectors = tuple(R.sparse_covector(es) for es in full)
                assert all(i in active for p in covectors for i, e in p)
                native = graph[omitted[0]][omitted[1]]
                field = A.pair_field(covectors, (native,))
                heights = []
                for p in covectors:
                    norm = R.projection(p, p, J, field)
                    assert not any(norm[1:])
                    height = norm[0] - 2
                    assert height.sign() > 0
                    heights.append(height)
                projection = R.projection(covectors[0], covectors[1], J, field)
                prescribed = field.negative_contact(native)
                assert field.signs_difference(tuple(a-b for a, b in zip(projection, prescribed)),
                                              heights[0]*heights[1]) == 0
                key = A.H.frozen(full)
                data = normalize(dict(contacts=[[[i, list(e)] for i, e in p] for p in covectors],
                    occupied_pairs=occupied, squared_timelike_heights=[x.out() for x in heights],
                    native_omitted_product=list(native)))
                if key in stars[vertex]:
                    assert stars[vertex][key]['native'] == data
                else:
                    stars[vertex][key] = dict(index=len(stars[vertex]), native=data,
                        active_embedding=mp, pair_field_radicals=[field.tau.out(), field.sigma.out()],
                        projection=R.qout(projection))
                indices.append(stars[vertex][key]['index'])
            assert len(indices) == len(set(indices)) == rec['literal_orientations']
            records.append(rec)
            occurrences[vertex].append(dict(origin=rec, literal_star_indices=sorted(indices)))
        face_records.append(dict(face=fi, source_context=q['source_context'], canonical_prism_core=4,
            provenance=None, raw_family=4, mass=q['mass'], demand_mass=q['demand_mass'], required=3,
            spectator_forest=spectator, corners=records))

    rows = []
    for vertex in sorted(targets):
        rows.append(dict(vertex=vertex, forest=names[vertex], target=targets[vertex],
            active_projection_blocks=inverses[vertex][2], occurrences=occurrences[vertex],
            literal_stars=list(stars[vertex].values())))
    inventory = cover_inventory(rows, baseline)
    paths = {str(Path(mod.__file__).relative_to(ROOT)) for mod in (A, R, M, V, C, X, B, P, R.T, Z, Z.U, Z.T, Z.H, Z.W)}
    paths.add(str(Path(__file__).relative_to(ROOT)))
    for stem, data in inputs.items():
        assert data['status'].startswith('PASS'), (stem, data['status'])
        for path, expected in data.get('sources', data.get('prepared', {}).get('sources', {})).items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, (stem, path)
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')):
            path = f'{folder}/{stem}.{ext}'
            if (ROOT/path).exists():
                paths.add(path)
    result = dict(status='PASS_NEW5121_RAW4_NEXT_FOUR_EXACT_DOMAIN_AND_BENEFIT_CEILING',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)},
        baseline_represented_classes=324, baseline_represented_forests=600,
        missing_positive_Q_classes=len(missing), missing_Q_provenance_counts=dict(missing_types),
        missing_raw4_classes=len(family), missing_raw4_demand_mass=str(sum((F(r['demand_mass']) for r in family), F())),
        selection='Four largest exact demand masses among currently unrepresented positive unambiguous raw4 classes.',
        selected_faces=face_records, forest_records=rows,
        cover_inventory=inventory,
        counts=dict(faces=len(selected), actual_corner_occurrences=sum(len(r['corners']) for r in face_records), vertex_types=len(rows),
            active_embeddings=maps_total, literal_orientation_occurrences=orientations_total,
            distinct_native_stars=sum(len(s) for s in stars.values())),
        current_compatibility_graphs_changed=False, pair_compatibility_tests=0,
        new_local_cover=False, new_global_optimization=False, new_dimension_bound=None,
        scope='Exact necessary corner reconstruction for four positive previously unrepresented raw4 classes, with full whole-component embeddings and physical norm/projection/native-product checks. Inactive spectator permutations do not alter contact data. Current324 native domains and frozen covers are transported solely for inventory and a singleton-cover upper bound on possible multiplier increases. No new physical pair test, exclusion, new complete cone, revised actual witness, optimizer, realization or dimension theorem.')
    return normalize(result)


def cover_inventory(raw_rows, baseline):
    """Extend frozen covers by literal singletons; no physical pair is tested."""
    p = Z.prepare()
    latest = {r['vertex']: r for r in baseline['forest_records']}
    known = set(p['previous']['latest']) | set(p['rows32'])
    assert len(known) == 600
    mu = {v: F(x) for v, x in baseline['vertex_multipliers']}
    alpha = {i: F(x) for i, x in baseline['face_multipliers']}
    records, increments = [], F()
    for raw in raw_rows:
        v = raw['vertex']
        if v in known:
            nodes, groups, graph, owner, forest, prior = Z.current(v, p)
            if v in latest:
                row = latest[v]; domain = row['domain']; oldn = len(nodes)
                assert domain['literal_stars'][:oldn] == nodes
                assert domain['corner_orbits'][:len(groups)] == groups
                nodes = Z.normal(domain['literal_stars'])
                groups = Z.normal(domain['corner_orbits'])
                graph.add_nodes_from(range(oldn, len(nodes)))
                from itertools import combinations
                pairs = [(i,j) for i,j in combinations(range(len(nodes)),2) if j >= oldn]
                assert [tuple(z[:2]) for z in row['newly_owned_pair_domain']] == pairs
                for i,j,code in row['newly_owned_pair_domain']:
                    if code is None: graph.add_edge(i,j)
                assert Z.L.jhash(sorted(tuple(sorted(e)) for e in graph.edges())) == row['graph_edges_sha256']
                demand = [sum((alpha.get(f,F())*F(m) for f,m in z['positive_face_coefficients']),F()) for z in groups]
                value = Z.W.check_certificate(row['maximal_count_vectors'], demand, row['certificate'])
                assert value == F(row['exact_local_optimum'])
                prior = dict(source=BASE, record_sha256=Z.L.jhash(row), exact_local_optimum=str(value),
                    graph_edges_sha256=row['graph_edges_sha256'])
            assert forest == raw['forest']
            oldvalue = F(prior['exact_local_optimum'])
            assert oldvalue <= mu.get(v,F())
        else:
            nodes, groups, oldvalue, prior = [], [], F(), None
        assert normalize(raw['target']) == normalize(R.P.forest_graph(raw['forest'],dict(R.C.models())))
        oldn, oldgroups = len(nodes), len(groups)
        keys = {Z.H.frozen(s['contacts']): i for i,s in enumerate(nodes)}
        assert len(keys) == len(nodes)
        mapping = []
        for s in raw['literal_stars']:
            native = s['native']; key = Z.H.frozen(native['contacts'])
            if key in keys:
                i = keys[key]
                assert all(nodes[i][k] == native[k] for k in Z.H.NATIVE_FIELDS)
            else:
                i = len(nodes); keys[key] = i; nodes.append(dict(index=i, **normalize(native)))
            mapping.append(i)
        assert len(mapping) == len(set(mapping))
        orbitkeys = {frozenset(z['literal_stars']): z['index'] for z in groups}
        occupied = {i for z in groups for i in z['literal_stars']}
        origins = {Z.H.origin_key(o) for z in groups for o in z['origins']}
        extras = Counter(); occurrence_rows = []
        for o in raw['occurrences']:
            origin = o['origin']; key = Z.H.origin_key(origin)
            assert key not in origins; origins.add(key)
            ids = sorted(mapping[i] for i in o['literal_star_indices']); sk = frozenset(ids)
            assert len(ids) == origin['literal_orientations']
            if sk in orbitkeys: gi = orbitkeys[sk]
            else:
                assert not occupied.intersection(ids), 'partial native orbit overlap'
                gi = len(groups); orbitkeys[sk] = gi; occupied.update(ids)
                groups.append(dict(index=gi, literal_stars=ids))
            extras[gi] += alpha[origin['face']]
            occurrence_rows.append(dict(origin=origin, orbit=gi, literal_star_indices=ids))
        assert occupied == set(range(len(nodes)))
        additional = sum(extras.values(), F())
        upper = max(mu.get(v,F()), oldvalue + additional)
        increments += upper - mu.get(v,F())
        nnew = len(nodes)-oldn
        row = dict(vertex=v, forest=raw['forest'], baseline=prior,
            old_literal_stars=oldn, new_literal_stars=nnew,
            old_corner_orbits=oldgroups,new_corner_orbits=len(groups)-oldgroups,
            raw4_to_merged_star_indices=mapping,actual_occurrences=occurrence_rows,
            extra_orbit_demands=[[i,str(x)] for i,x in sorted(extras.items())],
            candidate_new_pair_positions=oldn*nnew+nnew*(nnew-1)//2,
            singleton_extension_mass=str(additional), baseline_local_cover=str(oldvalue),
            baseline_vertex_mass=str(mu.get(v,F())), possible_vertex_mass_upper=str(upper),
            maximum_possible_mu_increase=str(upper-mu.get(v,F())))
        print('NEXT_RAW4_DOMAIN',v,'old',oldn,'new',nnew,'pairs',row['candidate_new_pair_positions'],flush=True)
        records.append(row)
    total = sum((F(r['singleton_extension_mass']) for r in records),F())
    selected = {o['origin']['face'] for r in raw_rows for o in r['occurrences']}
    assert total == 6*sum((alpha[f] for f in selected),F())
    num, den = F(baseline['objective_numerator']), F(baseline['normalizer'])
    assert den == sum(mu.values(),F())
    return dict(forests=records,affected_existing_forests=sum(r['vertex'] in known for r in records),
        potential_new_forests=sum(r['vertex'] not in known for r in records),
        new_literal_stars=sum(r['new_literal_stars'] for r in records),
        new_corner_orbits=sum(r['new_corner_orbits'] for r in records),
        candidate_new_pair_positions=sum(r['candidate_new_pair_positions'] for r in records),
        singleton_extension_mass=str(total),maximum_possible_mu_increase=str(increments),
        current_objective_numerator=str(num),current_normalizer=str(den),current_ratio=str(num/den),
        coarse_ratio_lower_bound=str(num/(den+total)),
        slack_aware_ratio_lower_bound=str(num/(den+increments)),
        claim='Every weaker graph obtained by testing only pairs involving these added stars admits the transported old covers plus singleton covers of all added actual-corner demands. Consequently its minimal monotone mu repair is bounded above as saved; ratios are bounded below. No physical exclusion or actual repaired witness is asserted.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
    else:
        assert not path.exists()
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'], flush=True)
