#!/usr/bin/env python3
"""Independent quadratic-field and binary-status check of side-lift rigidity."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
from fresh_h9_mixed_d4_coupling_20261014_review import Q
ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_h9_h8_side_tilts_20261014_review'
PROOF = 'results/fresh_h9_h8_side_tilts_20261014.json'
SOURCE = 'data/known_examples/bugaenko_compact_h8.json'
def sha(p): return sha256((ROOT/p).read_bytes()).hexdigest()
def run():
    d = json.loads((ROOT/PROOF).read_text())
    for p, h in d['sources'].items(): assert sha(p) == h
    src = json.loads((ROOT/SOURCE).read_text())
    def sf(x): return Q(F(x['a'])+F(x['b'])/2, F(x['b'])/2)
    G = [[sf(x) for x in row] for row in src['wall_gram']]
    roots = [[sf(x) for x in row] for row in src['simple_roots']]
    form = [[sf(x) for x in row] for row in src['lattice_gram']]
    assert all(sum((roots[i][a]*form[a][b]*roots[j][b] for a in range(9) for b in range(9)), Q()) == G[i][j]
               for i in range(11) for j in range(11))
    assert all(G[i][i] == 2 for i in range(11))
    edges = [(i, j) for i in range(11) for j in range(i+1, 11) if G[i][j] == -1]
    rights = [(i, j) for i in range(11) for j in range(i+1, 11) if not G[i][j]]
    assert [list(p) for p in edges] == d['pi3_edges']
    phi = Q(F(1, 2), F(1, 2))
    assert G[0][1] == -phi and G[8][7] == -phi
    cosine_square = phi*phi/4
    gaps = [3*cosine_square/4-Q(F(1, 4)), Q(F(1, 2))-3*cosine_square/4,
            cosine_square/2-Q(F(1, 4)), Q(F(1, 2))-cosine_square/2,
            Q(F(1, 4))-cosine_square*(1-cosine_square)]
    assert all(x.sign() > 0 for x in gaps)
    count = dict(pi3_mismatch=0,positive_right_pair=0,pi5_gap=0,all_right=0)
    for mask in product((0, 1), repeat=11):
        if any(mask[i] != mask[j] for i, j in edges): count['pi3_mismatch'] += 1; continue
        if any(mask[i] and mask[j] for i, j in rights): count['positive_right_pair'] += 1; continue
        if any(mask[i] != mask[j] for i, j in ((0, 1), (8, 7))): count['pi5_gap'] += 1; continue
        assert not any(mask); count['all_right'] += 1
    assert sum(count.values()) == 2048 and count['all_right'] == 1
    paths = [PROOF, SOURCE, f'scripts/{TAG}.py','scripts/fresh_h9_mixed_d4_coupling_20261014_review.py']
    return dict(status='PASS_INDEPENDENT_ALL_LABEL_SOURCE_SIDE_LIFT_REVIEW',
                sources={p:sha(p) for p in paths},binary_status_census=count,
                exact_positive_angle_gap_margins=[x.enc() for x in gaps],
                scope='Separate Fraction quadratic arithmetic reconstructs all121source Gram entries and the exact five angle gaps. Exhaustive2048 right/nonright masks verify the logical propagation. The normalized side-lift formula is reviewed algebraically; uniform sine monotonicity covers all finite labels>=5. Additional walls disjoint from the source facet remain outside scope.')
if __name__ == '__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--verify',action='store_true');args=parser.parse_args()
    result=run();out=ROOT/f'results/{TAG}.json'
    if args.verify:assert result==json.loads(out.read_text())
    else:out.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],result['binary_status_census'])
