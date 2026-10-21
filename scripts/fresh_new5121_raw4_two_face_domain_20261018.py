#!/usr/bin/env python3
"""Reconstruct two newly positive raw4 classes in full labelled vertex frames.

This is an exact necessary corner domain, not a new packing bound or chamber.
The existing 322-class compatibility graphs are not changed by this packet.
"""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import fresh_raw125_raw205_actual_corner_frontier_20261017 as A

ROOT = A.ROOT
R, M, V, K, C, X, B, P = A.R, A.M, A.V, A.K, A.C, A.X, A.B, A.P
TAG = 'fresh_new5121_raw4_two_face_domain_20261018'
MIGRATION = 'fresh_cross_vertex_support_migration_20261018'
CURRENT = 'fresh_cross_vertex_historical32_local_covers_20261018'
FACES = (14792, 14927)


def normalize(x):
    return json.loads(json.dumps(x, sort_keys=True))


def run():
    inputs = {}
    def read(stem):
        if stem not in inputs:
            inputs[stem] = json.loads((ROOT/f'results/{stem}.json').read_text())
        return inputs[stem]

    migration, current = read(MIGRATION), read(CURRENT)
    represented = set(migration['current290_classes']) | {
        r['selected_face']['face'] for r in current['added_faces']}
    assert len(represented) == 322
    missing = [r for r in migration['quad_rows'] if r['face'] not in represented]
    assert len(missing) == 1956
    missing_types = Counter(r['provenance_kind'] for r in missing)
    assert missing_types == dict(explicit=286, unique_saved_context=1212, merged_unexpanded=458)
    family = sorted((r for r in missing if r.get('raw_family') == 4),
                    key=lambda r: (-F(r['demand_mass']), r['face']))
    assert len(family) == 91 and tuple(r['face'] for r in family[:2]) == FACES
    selected = family[:2]
    alpha = {i: F(a) for i, a in current['face_multipliers']}

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
    paths = {str(Path(mod.__file__).relative_to(ROOT)) for mod in (A, R, M, V, C, X, B, P, R.T)}
    paths.add(str(Path(__file__).relative_to(ROOT)))
    for stem, data in inputs.items():
        assert data['status'].startswith('PASS'), (stem, data['status'])
        for path, expected in data.get('sources', data.get('prepared', {}).get('sources', {})).items():
            assert sha256((ROOT/path).read_bytes()).hexdigest() == expected, (stem, path)
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md')):
            path = f'{folder}/{stem}.{ext}'
            if (ROOT/path).exists():
                paths.add(path)
    result = dict(status='PASS_NEW5121_RAW4_TWO_FACE_EXACT_CORNER_DOMAIN',
        sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sorted(paths)},
        baseline_represented_classes=322, baseline_represented_forests=600,
        missing_positive_Q_classes=1956, missing_Q_provenance_counts=dict(missing_types),
        missing_raw4_classes=91, missing_raw4_demand_mass=str(sum((F(r['demand_mass']) for r in family), F())),
        selection='Two largest exact demand masses among currently unrepresented positive unambiguous raw4 classes.',
        selected_faces=face_records, forest_records=rows,
        counts=dict(faces=2, actual_corner_occurrences=12, vertex_types=len(rows),
            active_embeddings=maps_total, literal_orientation_occurrences=orientations_total,
            distinct_native_stars=sum(len(s) for s in stars.values())),
        current_compatibility_graphs_changed=False, pair_compatibility_tests=0,
        new_local_cover=False, new_global_optimization=False, new_dimension_bound=None,
        scope='Exact necessary corner reconstruction for two positive previously unrepresented raw4 classes, with full whole-component embeddings and physical norm/projection/native-product checks. Inactive spectator permutations do not alter contact data. No pair-compatibility graph, local cover, revised dual ratio, ambient polytope realization or dimension theorem.')
    return normalize(result)


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
