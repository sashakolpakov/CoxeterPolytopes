#!/usr/bin/env python3
"""Stream and identify nine old sealed E8 quad consumers; no queue rebuild."""
import argparse
import fcntl
from hashlib import file_digest
from itertools import combinations
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact29_small_faces_20261015 as P

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_compact21_sealed_e8_quad_triangle_activation_20261018'
UNION = 'results/fresh_sqrt5_quad_exhaustion_20261015.json'


def selected_union_rows(indices):
    """Read a pretty-printed top-level array with one bounded item in memory."""
    active, index, item = False, -1, []
    with (ROOT/UNION).open() as stream:
        for line in stream:
            if not active:
                if line.startswith('  "labelled_union_supports": ['):
                    active = True
                continue
            if line.startswith('  ],'):
                break
            if line.startswith('    {'):
                index += 1
                item = []
            if index in indices:
                item.append(line)
                assert sum(map(len, item)) < 65536
                if line.startswith('    }'):
                    yield index, json.loads(''.join(item).rstrip().rstrip(','))


def canonical_e8_map(vertices, edges):
    graph = {i: [j for j in vertices if j != i and edges[tuple(sorted((i, j)))] != 2]
             for i in vertices}
    assert len(vertices) == 8
    assert all(edges[tuple(sorted((i,j)))] == 3 for i in vertices for j in graph[i])
    branches = [i for i in vertices if len(graph[i]) == 3]
    assert len(branches) == 1
    branch = branches[0]
    arms = {}
    for start in graph[branch]:
        path, previous, current = [], branch, start
        while True:
            path.append(current)
            following = [j for j in graph[current] if j != previous]
            if not following:
                break
            assert len(following) == 1
            previous, current = current, following[0]
        assert len(path) not in arms
        arms[len(path)] = path
    assert set(arms) == {1, 2, 4}
    answer = {branch: 2}
    for length, names in ((1, [7]), (2, [1, 0]), (4, [3, 4, 5, 6])):
        answer.update(zip(arms[length], names))
    assert set(answer.values()) == set(range(8))
    return answer


def run():
    lp_path = 'results/fresh_sqrt5_sealed_weights_20261015.json'
    lp = json.loads((ROOT/lp_path).read_text())  # 162 KiB, not the union queue.
    classes = {940,958,982,984,993,995,1107,1143,1145}
    wanted = {r['union_index']: r for r in lp['face_rows'] if r['class_index'] in classes}
    assert len(wanted) == 9
    records = []
    for index, row in selected_union_rows(wanted):
        saved = wanted[index]
        connected = [c for c in saved['corners'] if c.get('connected')]
        assert len(connected) == 1 and connected[0]['type'] == 'E8'
        marked = connected[0]['marked']
        vertices = row['white']+marked
        edges = dict(zip(combinations(range(row['nodes']), 2), row['labels']))
        renaming = canonical_e8_map(vertices, edges)
        outside = sorted(set(range(row['nodes']))-set(vertices))
        assert len(outside) == 2 and edges[tuple(outside)] == 2
        contacts = []
        for root in outside:
            contact = [(j, edges[tuple(sorted((root,j)))]) for j in vertices
                       if edges[tuple(sorted((root,j)))] != 2]
            assert len(contact) == 1 and contact[0][1] == 5
            contacts.append(renaming[contact[0][0]])
        assert sorted(contacts) == [0,6]
        records.append({'class_index': saved['class_index'], 'union_index': index,
                        'canonical_marked_pair': sorted(renaming[i] for i in marked),
                        'outside_roots': outside, 'canonical_phi_contacts': contacts,
                        'new_wall_mutual_label': 2, 'canonical_E8_renaming': renaming})
    assert len(records) == 9
    pairs = {tuple(r['canonical_marked_pair']) for r in records}
    expected = {(i,j) for i in (0,1,2) for j in (3,4,5,6)}-{(2,3),(1,3),(2,4)}
    assert pairs == expected
    s = P.K((0,0,1,0))
    phi, phi_prime = (1+s)/2, (1-s)/2
    f0_conjugate_room = 2-4*phi_prime*phi_prime
    assert f0_conjugate_room.sign() > 0 and (1-f0_conjugate_room).sign() > 0
    h0, h6, projection = 4+2*s, 1+s, 3+s
    assert h0*h6 == projection*projection
    assert h0.sign() > 0
    # Every extra nonzero positive D5 projection increases only H6 and
    # breaks this equality, since the F0 wall has no D5 projection.
    sources = [f'scripts/{TAG}.py', 'scripts/fresh_compact29_small_faces_20261015.py', UNION, lp_path,
               'results/fresh_sqrt5_sealed_weights_20261015_geometry.json',
               'results/fresh_compact21_d5_two_e8_triangle_coupling_20261018.json']
    pins = {}
    for path in sources:
        with (ROOT/path).open('rb') as stream:
            pins[path] = file_digest(stream, 'sha256').hexdigest()
    return {'status': 'PASS_NINE_SEALED_E8_QUADS_FORCE_BOTH_TRIANGLE_BUNDLES',
            'sources': pins, 'records': records,
            'F0_conjugate_D5_projection_room': P.enc(f0_conjugate_room),
            'minimum_nonzero_spherical_D5_extension_projection': '1',
            'F0_F6_heights': P.enc([h0,h6]), 'mutual_old_E8_projection': P.enc(projection),
            'additional_D5_contacts_excluded_in_proof': True,
            'inactive_E8_triangle_state_excludes_these_nine_consumers': True,
            'all_E8_long_pair_quad_consumers_classified': False,
            'new_dimension_bound': None}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    rss = lambda: resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    assert rss() < 1024**3
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
        assert rss() < 1024**3
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == json.loads(json.dumps(result))
        else:
            assert not target.exists()
            target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'sealed_consumers': 9,
                      'peak_rss_bytes': rss(), 'new_dimension_bound': None}))
