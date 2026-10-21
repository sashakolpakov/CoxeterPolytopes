#!/usr/bin/env python3
"""All orthogonal finite rank13 spectator contexts of the1399 B1 dot cores."""
import argparse
import json
from collections import Counter
from hashlib import sha256
from pathlib import Path

import fresh_sqrt5_rank13_dotted_b1_quad_cores_20261017 as P

C = P.C
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank13_dotted_b1_spectators_20261017'


def dump(x):
    return json.dumps(x, separators=(',', ':'))


def sha(path):
    return sha256((ROOT/path).read_bytes()).hexdigest()


def run():
    core_path = f'results/{P.TAG}.json'
    prior = json.loads((ROOT/core_path).read_text())
    for path, digest in prior['sources'].items():
        assert sha(path) == digest, path
    assert prior['counts']['cores'] == 1399
    models = [(name, g) for name, g in C.models() if len(g) <= 11]
    codes = {name: json.loads(dump(C.code(g))) for name, g in models}
    ranks = {name: len(g) for name, g in models}
    certified = {r['name'] for r in prior['finite_component_certificates']}
    assert set(codes) <= certified

    def multisets(rank, first=0):
        if rank == 0:
            yield ()
            return
        for j in range(first, len(models)):
            name, graph = models[j]
            if len(graph) <= rank:
                for tail in multisets(rank-len(graph), j):
                    yield tuple(sorted((name,)+tail))

    spectator_sets = {n: sorted(multisets(n)) for n in range(12)}
    assert [len(spectator_sets[n]) for n in range(12)] == [1,1,6,9,29,47,117,196,419,704,1366,2276]
    for n, lists in spectator_sets.items():
        assert len(lists) == len(set(lists))
        assert all(sum(ranks[name] for name in names) == n for names in lists)
    corner_types, corner_lookup = [], {}
    contexts, faces, face_lookup = [], [], {}
    forests = set()
    for core in prior['records']:
        ci = core['index']
        for names in spectator_sets[core['spectator_rank']]:
            corner_indices = []
            for corner in core['corners']:
                forest = tuple(sorted(corner['forest']+list(names)))
                assert sum(int(next(r['rank'] for r in prior['finite_component_certificates'] if r['name'] == name)) for name in forest) == 13
                code = sorted(corner['marked_forest_code']+[codes[name] for name in names])
                key = (forest, dump(code))
                if key not in corner_lookup:
                    corner_lookup[key] = len(corner_types)
                    corner_types.append(dict(index=len(corner_types), forest=list(forest), marked_forest_code=code))
                    forests.add(forest)
                corner_indices.append(corner_lookup[key])
            index = len(contexts)
            coefficients = tuple(sorted(Counter(corner_indices).items()))
            if coefficients not in face_lookup:
                face_lookup[coefficients] = len(faces)
                faces.append(dict(index=len(faces), required=1, corner_coefficients=[list(x) for x in coefficients], contexts=[]))
            face = face_lookup[coefficients]
            faces[face]['contexts'].append(index)
            contexts.append(dict(index=index, core=ci, spectator_forest=list(names), corner_type_indices=corner_indices, face=face))
        if ci % 100 == 0:
            print('DOTTED_SPECTATOR_CONTEXTS', ci, len(contexts), len(corner_types), len(faces), flush=True)
    assert len(contexts) == 119698
    assert sum(len(r['contexts']) for r in faces) == len(contexts)
    paths = [f'scripts/{TAG}.py', core_path, f'scripts/{P.TAG}.py', f'scripts/{C.TAG}.py', P.LEMMA]
    return dict(status='PASS_COMPLETE_ORTHOGONAL_SPECTATOR_CONTEXTS_OF_RANK13_DOTTED_B1_CORES',
                sources={p: sha(p) for p in paths},
                counts=dict(cores=1399, contexts=len(contexts), corner_types=len(corner_types),
                            corner_forests=len(forests), distinct_quad_rows=len(faces)),
                spectator_counts={str(n): len(v) for n, v in spectator_sets.items()},
                corner_types=corner_types, contexts=contexts, face_rows=faces,
                proof='Every core has a certified local arithmetic compact quadrilateral. Every component multiset of the exact remaining rank15-n is enumerated, from all finite Qsqrt5 Coxeter types. Its orthogonal sum with the core yields an admissible rank14 form,15normals and11white roots; the full common-white face remains precisely the original four finite rays. Each exported corner is the original marked corner forest plus the unmarked spectator forest. Equal corner codes and duplicate four-corner inequalities are identified literally.',
                scope='Complete for orthogonal finite spectator extensions of these1399 dotted B1 cores only. This supplies119698 actual local quadrilateral contexts, not a complete rank13 face domain, a global finite-volume example, or a dimension bound. No weight payment is assumed.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(path.read_text())
        print('PASS_SAVED_DOTTED_B1_SPECTATOR_REPLAY')
    else:
        path.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], result['counts'], flush=True)
