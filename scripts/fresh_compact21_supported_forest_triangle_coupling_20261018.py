#!/usr/bin/env python3
"""Complete E7/E8 nonadjacent triangle domains beside whole D5 and A1.

No D5/A1 spectator contact is assumed zero during enumeration. A whole
untouched E7/E8 is excluded by the written finite-extension argument.
First finite-corner trees give a small exact complete contact domain.
"""
import argparse
from collections import Counter
import fcntl
from fractions import Fraction as Q
from hashlib import file_digest, sha256
from itertools import combinations, product
import json
from pathlib import Path
import resource
import signal
import sys
import fresh_compact21_d5_a1_triangle_coupling_20261018 as A

P = A.P
ROOT = A.ROOT
TAG = 'fresh_compact21_supported_forest_triangle_coupling_20261018'


def guard():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if rss*(1 if sys.platform == 'darwin' else 1024) > 1024**3:
        raise MemoryError('One-GiB RSS guard exceeded; no certificate written')


def initial_extensions(old, edges, omitted):
    nodes = [i for i in range(old) if i != omitted]
    blocks, _ = A.parts(nodes, edges)
    where = {node: i for i, block in enumerate(blocks) for node in block}
    options = []
    ledger = sha256()
    trials = 0
    for size in range(4):
        for support in combinations(nodes, size):
            if len({where[i] for i in support}) != size:
                continue
            labels = [(3, 3, 3)] if size == 3 else product((3, 4, 5), repeat=size)
            for values in labels:
                contacts = dict(zip(support, values))
                full = edges | {(i, old): m for i, m in contacts.items()}
                ok = A.finite(nodes+[old], full)
                ledger.update(json.dumps([sorted(contacts.items()), ok], separators=(',', ':')).encode()+b'\n')
                trials += 1
                if ok:
                    options.append(contacts)
    return options, {'trials': trials, 'accepted': len(options), 'digest': ledger.hexdigest()}


def domain(kind, n):
    edges = {(i, i+1): 3 for i in range(n-2)} | {(2, n-1): 3}
    edges |= {(i+n, j+n): m for (i, j), m in A.EDGES.items()}
    old = n+6  # Full marked exceptional component, whole D5 and whole A1.
    gram = P.base_gram(old, edges, old)
    inverse = P.inverse(gram)
    assert P.psd_status(gram) == 1
    blocks, _ = A.parts(list(range(old)), edges)
    where = {node: i for i, block in enumerate(blocks) for node in block}
    rows, roots, principal_cache = [], {}, {}
    counts = Counter()
    for a, b in combinations(range(n), 2):
        if (a, b) in edges:
            continue
        guard()
        options, first = initial_extensions(old, edges, a)
        ledger = sha256()
        eligible = []
        for initial in options:
            for label in range(2, 6):
                contacts = initial | ({a: label} if label > 2 else {})
                full = edges | {(i, old): m for i, m in contacts.items()}
                ok = A.finite([i for i in range(old+1) if i != b], full)
                ledger.update(json.dumps([sorted(contacts.items()), ok], separators=(',', ':')).encode()+b'\n')
                if not ok:
                    continue
                h = sum((P.COS[mi]*inverse[i][j]*P.COS[mj]
                         for i, mi in contacts.items() for j, mj in contacts.items()), P.K(-2))
                row = {'contacts': sorted(contacts.items()), 'height_squared': h}
                if h.sign() <= 0:
                    row['status'] = 'nonpositive_height'
                else:
                    cycles = [P.COS[contacts[i]]*P.COS[contacts[j]]
                              for i, j in combinations(sorted(contacts), 2) if where[i] == where[j]]
                    if any(z.v[1] or z.v[3] for z in cycles):
                        row['status'] = 'cycle_outside_Qsqrt5'
                    else:
                        assert not (h.v[1] or h.v[3])
                        conjugate = P.K((h.v[0], 0, -h.v[2], 0))
                        row['conjugate_height_squared'] = conjugate
                        if conjugate.sign() > 0:
                            row['status'] = 'conjugate_not_PSD'
                        else:
                            key = tuple(sorted(contacts.items()))
                            if key not in principal_cache:
                                groups, _ = A.parts(list(range(old+1)), full)
                                component = next(group for group in groups if old in group)
                                others = sorted(set(component)-{old})
                                g = P.base_gram(old, full, old+1)
                                witness = None
                                for size in range(len(others)+1):
                                    for subset in combinations(others, size):
                                        indices = list(subset)+[old]
                                        cc, _ = A.parts(indices, full)
                                        if len(cc) != 1 or any(m == 5 and i in indices and j in indices
                                                               for (i, j), m in full.items()):
                                            continue
                                        if A.finite(indices, full):
                                            continue
                                        if P.psd_status(P.sub(g, indices)) == 0:
                                            witness = indices
                                            break
                                    if witness is not None:
                                        break
                                principal_cache[key] = witness
                            row['parabolic_subset'] = principal_cache[key]
                            row['status'] = ('parabolic_obstruction' if principal_cache[key] is not None
                                             else 'necessary_arithmetic_triangle_root')
                            if principal_cache[key] is None:
                                assert all(i < n or i == n+5 for i in contacts)
                                record = roots.setdefault(key, {'kind': kind, 'contacts': sorted(contacts.items()),
                                    'A1_label': contacts.get(n+5, 2), 'height_squared': h, 'marked_pairs': []})
                                record['marked_pairs'].append([a, b])
                counts[row['status']] += 1
                eligible.append(row)
        rows.append({'marked': [a, b], 'first_corner': first,
                     'second_corner_trials': 4*len(options),
                     'second_corner_digest': ledger.hexdigest(), 'eligible': eligible})
    expected = ({((0, 5),), ((5, 5),), ((5, 3), (12, 5))} if n == 7 else
                {((0, 5),), ((6, 5),)})
    assert set(roots) == expected
    return {'kind': kind, 'old_nodes': old, 'D5_nodes': list(range(n, n+5)), 'A1_node': n+5,
            'domains': rows, 'status_counts': dict(counts), 'roots': list(roots.values())}, inverse


def integral5(x):
    a, b, c, d = x.v
    assert not b and not d
    return (2*c).denominator == 1 and (a-c).denominator == 1


def run():
    exceptional, inverse_by_kind = [], {}
    for kind, n in (('E7', 7), ('E8', 8)):
        result, inverse = domain(kind, n)
        exceptional.append(result)
        inverse_by_kind[kind] = inverse
    previous_path = ROOT/'results/fresh_compact21_d5_a1_triangle_coupling_20261018.json'
    previous = json.loads(previous_path.read_text())
    assert previous['status'] == 'PASS_COMPACT21_D5_A1_COMPLETE_NONADJACENT_TRIANGLE_DOMAIN'
    inverse_by_kind['D5'] = P.inverse(P.base_gram(6, A.EDGES, 6))
    roots = []
    for row in previous['surviving_roots']:
        contacts = dict(row['contacts'])
        labels = [contacts.get(i, 2) for i in range(6)]
        marked = [r['marked'] for r in previous['domains'] if any(
            q['status'] == 'necessary_arithmetic_triangle_root' and q['labels'] == labels for q in r['eligible'])]
        roots.append({'kind': 'D5', 'contacts': row['contacts'], 'A1_label': contacts.get(5, 2),
                      'height_squared': P.K(tuple(map(Q, row['height_squared']))), 'marked_pairs': marked})
    for data in exceptional:
        roots.extend(data['roots'])
    assert len(roots) == 11
    pairs = []
    compatible = {}
    for i, j in combinations(range(len(roots)), 2):
        a, b = roots[i], roots[j]
        p = P.K()
        if a['kind'] == b['kind']:
            n = int(a['kind'][1:])
            inv = inverse_by_kind[a['kind']]
            p += sum((P.COS[mi]*inv[x][y]*P.COS[mj]
                      for x, mi in a['contacts'] if x < n
                      for y, mj in b['contacts'] if y < n), P.K())
        p += P.COS[a['A1_label']]*P.COS[b['A1_label']]/2
        hprod = a['height_squared']*b['height_squared']
        square_margin = p*p-hprod
        if p and square_margin.sign() > 0:
            assert p.sign() > 0
            ok, status = False, 'positive_mutual_product'
        elif p:
            assert not square_margin
            ok, status = True, 'orthogonal_necessary_pair'
        elif not integral5(hprod):
            ok, status = False, 'nonintegral_length_two_cycle'
        else:
            if (hprod-4).sign() < 0:
                assert hprod in [P.K(1), P.K(2), P.COS[5]*P.COS[5], P.K(3), P.K((Q(5,2),0,Q(1,2),0))]
            ok, status = True, 'native_or_dotted_necessary_pair'
        compatible[i, j] = ok
        pairs.append({'roots': [i, j], 'projection_product': p, 'height_product': hprod,
                      'square_margin': square_margin, 'status': status, 'compatible_necessary': ok})
    maxima, best = [], -1
    for bits in range(1 << len(roots)):
        selected = [i for i in range(len(roots)) if bits >> i & 1]
        if any(not compatible[i, j] for i, j in combinations(selected, 2)):
            continue
        marked = {(roots[i]['kind'], tuple(pair)) for i in selected for pair in roots[i]['marked_pairs']}
        score = len(marked)
        if score > best:
            best, maxima = score, []
        if score == best:
            maxima.append(selected)
    assert best == 7 and len(maxima) == 1
    maximizing = [roots[i] for i in maxima[0]]
    assert sorted((r['kind'], r['contacts']) for r in maximizing) == [
        ('E7', [(5, 5)]), ('E8', [(0, 5)]), ('E8', [(6, 5)])]
    source_paths = [f'scripts/{TAG}.py', 'scripts/fresh_compact21_d5_a1_triangle_coupling_20261018.py',
                    'scripts/fresh_compact29_small_faces_20261015.py',
                    'results/fresh_compact21_d5_a1_triangle_coupling_20261018.json',
                    'results/fresh_compact21_single_field_20261016.json',
                    'results/fresh_compact21_single_field_20261016.md']
    sources = {}
    for path in source_paths:
        with (ROOT/path).open('rb') as stream:
            sources[path] = file_digest(stream, 'sha256').hexdigest()
    return P.enc({'status': 'PASS_COMPACT21_A1_D5_E7_E8_COMPLETE_NONADJACENT_TRIANGLE_COUPLING',
                  'sources': sources, 'forest': ['A1', 'D5', 'E7', 'E8'], 'ground_field': 'Q(sqrt5)',
                  'exceptional_domains': exceptional, 'roots': roots, 'pairs': pairs,
                  'root_subsets_checked': 1 << len(roots), 'nonadjacent_faces': 42,
                  'triangle_maximum': best, 'maximizing_necessary_subsets': maxima,
                  'D5_spectator_contacts_were_included': True,
                  'unchanged_old_face_payments': True, 'new_charge_saving': False,
                  'new_dimension_bound': None, 'Lean_verified': False})


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.alarm(59)
    lock_path = ROOT/'formal/.proof-batches/exclusive.lock'
    lock_path.parent.mkdir(exist_ok=True)
    with lock_path.open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        result = run()
        target = ROOT/f'results/{TAG}.json'
        if args.verify:
            assert json.loads(target.read_text()) == result
        else:
            assert not target.exists()
            target.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'roots': len(result['roots']),
                      'pair_checks': len(result['pairs']), 'triangle_maximum': result['triangle_maximum'],
                      'new_dimension_bound': None, 'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                          (1 if sys.platform == 'darwin' else 1024)}))
