#!/usr/bin/env python3
"""All-length rank completion for a single attached simply-laced white path.

Only its first node touches a nonsingular finite unpaid core. Other external
components, branching, and paths with two attached ends are outside scope.
"""
import argparse
import json
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import fresh_sqrt5_quad_contacts_20261015 as C
import fresh_sqrt5_quad_filters_20261015 as A

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_sqrt5_single_chain_completion_20261015'
K = A.K


def inverse(M):
    n = len(M)
    G = [r[:] + [K(int(i==j)) for j in range(n)] for i, r in enumerate(M)]
    for j in range(n):
        p = next(i for i in range(j, n) if G[i][j])
        G[j], G[p] = G[p], G[j]
        pivot = G[j][j]
        G[j] = [v/pivot for v in G[j]]
        for i in range(n):
            if i != j:
                factor = G[i][j]
                G[i] = [x-factor*y for x, y in zip(G[i], G[j])]
    out = [r[n:] for r in G]
    assert all(A.dot(M[i], [out[k][j] for k in range(n)]) == int(i==j)
               for i in range(n) for j in range(n))
    return out


def chain(row, contacts, length):
    n = row['nodes']
    old = C.edges(row)
    g = {p: old.get(p, 2) for p in combinations(range(n+length), 2)}
    for i, m in enumerate(contacts):
        g[i, n] = m
    for j in range(n, n+length-1):
        g[j, j+1] = 3
    return dict(nodes=n+length, white=list(range(4, n+length)),
                labels=[g[p] for p in combinations(range(n+length), 2)])


def run():
    contacts = json.loads((ROOT/'results/fresh_sqrt5_quad_contacts_20261015.json').read_text())
    base = json.loads((ROOT/'results/fresh_sqrt5_quad_exhaustion_20261015.json').read_text())
    records = []
    counts = Counter()
    lengths = Counter()
    for rec in contacts['records']:
        if rec['status']!='UNPAID_FINITE_CONTACT_DOMAIN' or rec['singular'] or not rec['accepted']:
            continue
        row = base['labelled_union_supports'][rec['union_index']]
        colors = A.coloring(row)
        H = A.scaled_matrix(row, colors['colors'])
        Hi = inverse(H)
        for ti, trial in enumerate(rec['trials']):
            if trial['status']!='PASS_EXTRA_WHITE_NECESSARY_TESTS':
                continue
            one = C.extension(row, trial['contacts'])
            ec = A.coloring(one)
            assert ec['passes'] and ec['colors'][:-1] == colors['colors']
            G = A.scaled_matrix(one, ec['colors'])
            v = [r[-1] for r in G[:-1]]
            d = G[-1][-1]
            c = A.dot(v, [A.dot(r, v) for r in Hi])
            assert d.sign()>0 and d.conj().sign()>0
            denominator = 2*c-d
            out = dict(source_class_index=rec['class_index'], trial_index=ti, contacts=trial['contacts'],
                       first_norm=d.enc(), projection_cost=c.enc(), denominator=denominator.enc())
            if not denominator:
                out['status'] = 'NO_LENGTH_ZERO_DENOMINATOR'
            else:
                length = d/denominator
                out['length_in_field'] = length.enc()
                if length.b or length.a.denominator != 1 or length.a < 1:
                    out['status'] = 'NO_POSITIVE_INTEGER_LENGTH'
                else:
                    length = int(length.a)
                    assert (length+1)*d == 2*length*c
                    out['length'] = length
                    lengths[length] += 1
                    whole = chain(row, trial['contacts'], length)
                    g = C.edges(whole)
                    if not C.E.corners_finite(whole['white'], g):
                        out['status'] = 'UNIQUE_LENGTH_REJECT_NONFINITE_CORNER'
                    else:
                        cc = A.coloring(whole)
                        assert cc['passes']
                        W = A.scaled_matrix(whole, cc['colors'])
                        assert all(W[row['nodes']+i][row['nodes']+i] == d for i in range(length))
                        assert all(W[row['nodes']+i][row['nodes']+i+1] == -d/2 for i in range(length-1))
                        inertia = A.inertia(W)
                        conjugate = A.inertia([[x.conj() for x in r] for r in W])
                        assert inertia == [whole['nodes']-2, 1, 1]
                        assert conjugate == [whole['nodes']-1, 0, 1]
                        out['physical_inertia'] = inertia
                        out['conjugate_inertia'] = conjugate
                        witness = C.affine_witness(whole)
                        if witness:
                            out.update(status='UNIQUE_LENGTH_REJECT_AFFINE', affine_witness=witness)
                        elif not C.E.no_short(whole['white'], g):
                            out['status'] = 'UNIQUE_LENGTH_PAID_SHORT_PATH'
                        elif C.h4_payments(whole):
                            out['status'] = 'UNIQUE_LENGTH_PAID_H4'
                        else:
                            out.update(status='UNIQUE_LENGTH_UNPAID_PRINCIPAL_STAR', row=whole)
            counts[out['status']] += 1
            records.append(out)
    assert len(records) == 458
    files = [f'scripts/{TAG}.py', 'scripts/fresh_sqrt5_quad_contacts_20261015.py',
             'scripts/fresh_sqrt5_quad_filters_20261015.py',
             'scripts/fresh_sqrt5_quad_exhaustion_20261015.py',
             'results/fresh_sqrt5_quad_contacts_20261015.json',
             'results/fresh_sqrt5_quad_exhaustion_20261015.json']
    return dict(status='PASS_EXACT_ALL_LENGTH_SINGLE_WHITE_CHAIN_COMPLETION',
                sources={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
                profiles=len(records), counts=dict(sorted(counts.items())),
                unique_positive_integer_lengths={str(k): v for k, v in sorted(lengths.items())}, records=records,
                determinant_identity='For H nonsingular, v its first-chain-node contact column, d the common scaled chain norm and c=v^T H^-1 v, the external A_l block is(d/2)Cartan(A_l), whose first inverse entry is2l/(d(l+1)). The determinant of the full block is det(H)det(A_l)(1−2lc/(d(l+1))). Corank1 forces(l+1)d−2lc=0, hence l=d/(2c−d). If2c=d there is no solution because d>0. An actual integer length must be the same rational integer in both embeddings. All458 values are checked exactly; no length cutoff is used.',
                scope='Complete rank-completion necessity only for one external simply-laced path with just its first node touching a nonsingular unpaid finite core. Additional white components entirely orthogonal to both core and path do not change this condition. Multiple touching components, branching, paths joining two attachments, dotted cores and arbitrary further completions remain open. Surviving principal stars are not claimed to be new polytopes or a complete face classification.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    out = run()
    path = ROOT/'results'/f'{TAG}.json'
    if args.verify:
        assert out == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], out['counts'], out['unique_positive_integer_lengths'])
