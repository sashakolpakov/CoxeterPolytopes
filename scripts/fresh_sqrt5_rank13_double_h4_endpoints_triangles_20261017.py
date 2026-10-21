#!/usr/bin/env python3
"""Complete local triangles at three dual-active double-H4 rank13 vertex types."""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path

import fresh_sqrt5_adjacent_discount_witness_20261016 as A
import fresh_sqrt5_forest_sensitive_weights_20261015_review as C

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_rank13_double_h4_endpoints_triangles_20261017'
K, Z = A.K, A.Z
PHI = K(F(1, 2), F(1, 2))
TARGETS = [('D5', 'H4', 'H4'), ('A1', 'D4', 'H4', 'H4'),
           ('A1', 'A4', 'H4', 'H4')]


def baseline(g):
    return sum((C.base(piece, a, b) for piece in C.pieces(g)
                for a, b in combinations(piece, 2)), F(0))


def run():
    models = [(name, g) for name, g in C.models() if len(g) <= 13]
    finite = {C.code(g): name for name, g in models}
    assert len(finite) == len(models)
    selected_models = [(name, g) for name, g in models
                       if all(m in (3, 5) for row in g.values() for m in row.values())]
    graphs = dict(selected_models)
    entries = [(name, len(g), baseline(g)) for name, g in selected_models]
    def forests(rank, first=0):
        if not rank:
            yield (), F(0)
            return
        for j in range(first, len(entries)):
            name, n, b = entries[j]
            if n <= rank:
                for tail, cost in forests(rank-n, j):
                    yield tuple(sorted((name,)+tail)), b+cost
    census = sorted((names, b) for names, b in forests(13))
    selected = {names: b for names, b in census if b >= 12}
    assert all(names in selected for names in TARGETS)

    def classify(g):
        names = []
        for piece in C.pieces(g):
            if sum(map(len, piece.values())) != 2*(len(piece)-1):
                return None
            if max(map(len, piece.values())) > 3:
                return None
            name = finite.get(C.code(piece))
            if name is None:
                return None
            names.append(name)
        return tuple(sorted(names))

    targets, counts, third_types = [], Counter(), {}
    for names in TARGETS:
        g = {}
        for name in names:
            offset = len(g)
            for i, row in graphs[name].items():
                g[offset+i] = {offset+j: m for j, m in row.items()}
        assert len(g) == 13 and classify(g) == names
        gram = [[K(2) if i == j else -({3: K(1), 5: PHI}[g[i][j]])
                 if j in g[i] else Z for j in g] for i in g]
        inv = A.inverse(gram)
        ports, occurrences, pairs = [], [], set()
        for port in range(13):
            attachments = []
            remaining = [i for i in g if i != port]
            for size in range(4):
                labelings = {0: [()], 1: [(3,), (5,)],
                             2: [(3, 3), (3, 5), (5, 3)], 3: [(3, 3, 3)]}[size]
                for nodes in combinations(remaining, size):
                    for labels in labelings:
                        counts['endpoint_contact_tests'] += 1
                        contacts = dict(zip(nodes, labels))
                        end = {i: {j: m for j, m in row.items() if j != port}
                               for i, row in g.items() if i != port}
                        end[13] = dict(contacts)
                        for i, m in contacts.items():
                            end[i][13] = m
                        destination = classify(end)
                        if destination not in selected:
                            continue
                        index = len(attachments)
                        record = dict(contacts=[[i, m] for i, m in contacts.items()],
                                      destination=list(destination), destination_baseline=str(selected[destination]))
                        attachments.append(record)
                        counts['accepted_finite_endpoints'] += 1
                        magnitudes = {i: K(1) if m == 3 else PHI for i, m in contacts.items()}
                        a = inv[port][port]
                        b = -sum((inv[port][i]*t for i, t in magnitudes.items()), Z)
                        c = 2-sum((ti*inv[i][j]*tj for i, ti in magnitudes.items()
                                   for j, tj in magnitudes.items()), Z)
                        assert a.sign() > 0 and b.sign() <= 0
                        for label in (2, 3, 4, 5, 6, 10):
                            counts['framing_tests'] += 1
                            if label in (4, 6, 10):
                                if b:
                                    counts['nonfield_framing_cycle'] += 1
                                    continue
                                sq = {4: K(2), 6: K(3), 10: 2+PHI}[label]
                                residual = c-a*sq
                            else:
                                t = {2: Z, 3: K(1), 5: PHI}[label]
                                residual = c+2*b*t-a*t*t
                            if residual.sign() != -1 or residual.conj().sign() != 1:
                                counts['edge_signature_rejection'] += 1
                                continue
                            counts['admitted_edge_labels'] += 1
                            full = {i: dict(row) for i, row in g.items()}
                            full[13] = dict(contacts)
                            for i, m in contacts.items():
                                full[i][13] = m
                            if label != 2:
                                full[port][13] = full[13][port] = label
                            for other in remaining:
                                counts['third_corner_tests'] += 1
                                third = {i: {j: m for j, m in row.items() if j != other}
                                         for i, row in full.items() if i != other}
                                forest = classify(third)
                                if forest is None:
                                    continue
                                counts['accepted_triangle_occurrences'] += 1
                                counts['label_' + str(label)] += 1
                                counts['all_selected_triangles'] += forest in selected
                                pair = tuple(sorted((port, other)))
                                pairs.add(pair)
                                pieces = C.pieces(third)
                                marked = next((p for p in pieces if port in p and 13 in p), None)
                                corner = F(0) if marked is None else C.base(marked, port, 13)
                                third_base = baseline(third)
                                assert sum((baseline(p) for p in pieces), F(0)) == third_base
                                key = '|'.join(forest)
                                rec = third_types.setdefault(key, dict(forest=list(forest), baseline=str(third_base),
                                                                      selected=forest in selected, occurrences=0))
                                assert rec['baseline'] == str(third_base)
                                rec['occurrences'] += 1
                                occurrences.append(dict(port=port, other=other, attachment=index, label=label,
                                    destination=list(destination), third_forest=list(forest),
                                    third_is_selected=forest in selected, third_baseline=str(third_base),
                                    third_corner_baseline=str(corner),
                                    third_marked_forest_code=sorted(C.code(p, (port, 13)) for p in pieces),
                                    schur_coefficients=[x.enc() for x in (a, b, c)],
                                    physical_residual=residual.enc(), conjugate_residual=residual.conj().enc()))
            ports.append(dict(port=port, attachments=attachments))
        targets.append(dict(forest=list(names), baseline=str(selected[names]),
                            edges=[[i, j, m] for i, row in g.items() for j, m in row.items() if i < j],
                            ports=ports, triangle_occurrences=occurrences,
                            possible_triangle_corners=[list(p) for p in sorted(pairs)]))
        print('COMPENSATING_TARGET', names, sum(len(p['attachments']) for p in ports),
              len(occurrences), len(pairs), flush=True)
    paths = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.py', 'results/fresh_sqrt5_rank13_compensating_endpoints_triangles_20261017.json', 'scripts/fresh_sqrt5_adjacent_discount_witness_20261016.py',
             'scripts/fresh_sqrt5_forest_sensitive_weights_20261015_review.py',
             'results/fresh_sqrt5_rank13_small_d_third_types_20261017.json']
    return dict(status='PASS_COMPLETE_THREE_DOUBLE_H4_ENDPOINT_AND_TRIANGLE_DOMAINS',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths},
                dimension=13, allowed_selected_components=[name for name, _ in selected_models],
                all_35_forest_count=len(census), selected_forest_count=len(selected),
                selected_forests=[dict(forest=list(names), baseline=str(b)) for names, b in sorted(selected.items())],
                counts=dict(counts), targets=targets,
                third_forests=sorted(third_types.values(), key=lambda x: x['forest']),
                completeness='The selected set comprises every rank13 finite forest with ADE,H3,H4,I2(5) components and baseline>=12. At a destination in this set the new root has degree<=3; a degree3 root has only3-bonds, and a degree2 root has at most one5-bond. Every such contact subset and assignment is tested. The framing finite-label list2,3,4,5,6,10 is complete because its squared product lies in Q(sqrt5). Every released other root is tested against the full finite Coxeter classification.',
                field_gate='All fixed products are0,-1,-phi. A nonzero b is equivalent to a fixed nonzero-contact path from the old port to the new root. Its cycle forces the framing product intoK, excluding labels4,6,10. When b=0 their squared products suffice. No rational-product exclusion is made: the old H3/H4-containing Gram has nonrational inverse.',
                signatures='The original finite Gram is positive definite in both embeddings. The edge Schur residual c-2*b*g-a*g^2 must be negative physically and positive at the conjugate embedding. All arithmetic uses exact Fractions inQ(sqrt5).',
                scope='Complete necessary finite endpoint and triangular-face domains at three specified double-H4 vertices, with the other endpoint in the explicitly enlarged selected set. This does not enumerate dotted edge lengths or quadrilaterals, and does not claim simultaneous realization, a valid payment rule or a dimension bound. In particular all-selected triangles may now occur; the previous ADE-only all-selected-triangle exclusion is not reused.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    data = json.loads(json.dumps(run()))
    path = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data, indent=2)+'\n')
    print(data['status'], data['selected_forest_count'], data['counts'], flush=True)
