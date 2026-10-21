#!/usr/bin/env python3
"""Exact facet screen of Vinberg's known nonarithmetic H18 mixture.

The two finite-volume parent chambers are I_(18,1) (37 facets) and
II_(17,1)+<2> (24 facets).  This is an exploratory search for Coxeter
faces, not a new-polytope claim.  All decisions use Q(sqrt(2)).
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from pathlib import Path

import networkx as nx
import numpy as np
from sympy import Matrix

from coxeter_polytopes.exact.silver import Silver, silver_inertia
from scripts.certify_rational_i1n_inverse_design import enumerate_simple_roots, lorentz_pair


def graph(gram, indices):
    result = nx.Graph()
    for i in indices:
        result.add_node(i, norm=gram[i][i])
    for j, i in enumerate(indices):
        for k in indices[:j]:
            if gram[i][k]:
                result.add_edge(i, k, pairing=gram[i][k])
    return result


def isolated(gram):
    return [i for i in range(len(gram)) if all(
        not gram[i][j] or gram[i][j]**2 >= gram[i][i]*gram[j][j]
        for j in range(len(gram)) if j != i)]


def parents():
    roots = enumerate_simple_roots(18)
    left = [[lorentz_pair(x, y) for y in roots] for x in roots]
    case = next(c for c in json.loads(Path('data/known_examples/k3_high_rank_nef_chambers.json').read_text())['cases'] if c['rho'] == 19)
    right = [[2 if i == j else 0 for j in range(24)] for i in range(24)]
    for pairing, pairs in case['intersections'].items():
        for i, j in pairs:
            right[i-1][j-1] = right[j-1][i-1] = -int(pairing)
    return left, right


def mixture():
    left, right = parents()
    li, ri = isolated(left), isolated(right)
    a, b = li[0], ri[0]
    assert left[a][a] == 1 and right[b][b] == 2
    common_left = [i for i in range(len(left)) if left[a][i] == 0]
    common_right = [i for i in range(len(right)) if right[b][i] == 0]
    matcher = nx.algorithms.isomorphism.GraphMatcher(
        graph(left, common_left), graph(right, common_right),
        node_match=lambda x, y: x['norm'] == y['norm'],
        edge_match=lambda x, y: x['pairing'] == y['pairing'])
    mapping = next(matcher.isomorphisms_iter())
    common = Matrix([[left[i][j] for j in common_left] for i in common_left])
    basis_positions = list(common.rref()[1])
    basis = [common_left[i] for i in basis_positions]
    inverse = Matrix([[left[i][j] for j in basis] for i in basis]).inv()
    assert len(basis) == 18
    retained_left = [i for i in range(len(left)) if i != a]
    retained_right = [j for j in range(len(right)) if j != b and j not in common_right]
    labels = [f'L{i+1}' for i in retained_left] + [f'R{j+1}' for j in retained_right]
    gram = [[Silver(0) for _ in labels] for _ in labels]
    for x, i in enumerate(retained_left):
        for y, j in enumerate(retained_left):
            gram[x][y] = Silver(left[i][j])
    shift = len(retained_left)
    for x, i in enumerate(retained_right):
        for y, j in enumerate(retained_right):
            gram[shift+x][shift+y] = Silver(right[i][j])
    for x, i in enumerate(retained_left):
        lpair = Matrix([[left[i][k] for k in basis]])
        for y, j in enumerate(retained_right):
            rpair = Matrix([right[j][mapping[k]] for k in basis])
            rational = Fraction(lpair.dot(inverse*rpair))
            # The gluing normals have norms 1 and 2 and opposite directions.
            radical = -Fraction(left[i][a]*right[j][b], 2)
            gram[x][shift+y] = gram[shift+y][x] = Silver(rational, radical)
    assert len(gram) == 40
    assert all(gram[i][j] <= 0 for i in range(40) for j in range(i))
    allowed = {Silver(0), Silver(Fraction(1,4)), Silver(Fraction(1,2)), Silver(Fraction(3,4)), Silver(Fraction(1,2), Fraction(1,4))}
    assert all(gram[i][j]**2 >= gram[i][i]*gram[j][j] or gram[i][j]**2/(gram[i][i]*gram[j][j]) in allowed for i in range(40) for j in range(i))
    return gram, labels, {'left_isolated_zero_based': li, 'right_isolated_zero_based': ri, 'gluing_facets_zero_based': [a,b], 'common_facet_count': len(common_left), 'common_basis_rank': len(basis), 'matching': {str(k):v for k,v in mapping.items()}}


def facet_screen(gram, labels):
    allowed = {Silver(0), Silver(Fraction(1,4)), Silver(Fraction(1,2)), Silver(Fraction(3,4)), Silver(Fraction(1,2), Fraction(1,4))}
    output = []
    for k in range(len(gram)):
        indices = [i for i in range(len(gram)) if i != k and gram[i][i]*gram[k][k] > gram[i][k]**2]
        restricted = [[gram[i][j]-gram[i][k]*gram[j][k]/gram[k][k] for j in indices] for i in indices]
        failures = []
        for ii, i in enumerate(indices):
            for jj, j in enumerate(indices[:ii]):
                pairing = restricted[ii][jj]
                square = pairing**2 / (restricted[ii][ii]*restricted[jj][jj])
                if pairing > 0 or (square < 1 and square not in allowed):
                    failures.append({'walls': [labels[i], labels[j]], 'pairing': pairing.to_json(), 'normalized_pairing_square': square.to_json(), 'obtuse': pairing > 0})
        output.append({'facet': labels[k], 'dimension': 17, 'wall_count': len(indices), 'forbidden_pair_count': len(failures), 'coxeter_angle_gate': not failures, 'first_forbidden_pair': failures[0] if failures else None})
    return output


def normalized_graph(gram):
    result = nx.Graph()
    result.add_nodes_from(range(len(gram)))
    for i in range(len(gram)):
        for j in range(i):
            if gram[i][j]:
                value = gram[i][j]**2/(gram[i][i]*gram[j][j])
                result.add_edge(i, j, value=value)
    return result


def identify_facets(gram, labels, facets):
    roots = enumerate_simple_roots(17)
    i17 = [[Silver(lorentz_pair(x, y)) for y in roots] for x in roots]
    case = next(c for c in json.loads(Path('data/known_examples/k3_high_rank_nef_chambers.json').read_text())['cases'] if c['rho'] == 18)
    ii17 = [[Silver(2 if i == j else 0) for j in range(19)] for i in range(19)]
    for pairing, pairs in case['intersections'].items():
        for i,j in pairs:
            ii17[i-1][j-1] = ii17[j-1][i-1] = Silver(-int(pairing))
    models = [('known I_(17,1) 22-wall chamber', normalized_graph(i17)), ('known II_(17,1) 19-wall chamber', normalized_graph(ii17))]
    for facet in facets:
        if not facet['coxeter_angle_gate']:
            continue
        k = labels.index(facet['facet'])
        indices = [i for i in range(len(gram)) if i != k and gram[i][i]*gram[k][k] > gram[i][k]**2]
        restricted = [[gram[i][j]-gram[i][k]*gram[j][k]/gram[k][k] for j in indices] for i in indices]
        candidate = normalized_graph(restricted)
        matches = [name for name, graph in models if nx.is_isomorphic(candidate, graph, edge_match=lambda x,y:x['value']==y['value'])]
        assert len(matches) == 1
        facet['exact_weighted_gram_identification'] = matches[0]


def seam_face_screen(gram, labels, construction, codimensions):
    """Screen every face specified by codimension r walls crossing the seam.

    Numeric arithmetic only suggests a forbidden pair. Each exclusion is
    proved by exact Sylvester/Schur computations; failures to suggest a pair
    trigger complete exact pair checking. The list of facets of a face is
    precisely its elliptic one-node extensions.
    """
    seam_labels = [f'L{int(i)+1}' for i in construction['matching']]
    seam_indices = sorted(labels.index(x) for x in seam_labels)
    approximate = np.array([[float(x.a)+2**.5*float(x.b) for x in row] for row in gram])
    allowed = {Silver(0), Silver(Fraction(1,4)), Silver(Fraction(1,2)), Silver(Fraction(3,4)), Silver(Fraction(1,2), Fraction(1,4))}
    rows=[]
    for codim in codimensions:
        total=0; elliptic=0; rejected=0; survivors=[]; witnesses=[]
        for selection in combinations(seam_indices, codim):
            total += 1
            # Seam walls have norm two and rational pairings.
            block = Matrix([[gram[i][j].a for j in selection] for i in selection])
            if any(block[:k,:k].det() <= 0 for k in range(1,codim+1)):
                continue
            elliptic+=1
            inverse=block.inv()
            inv=[[Silver(Fraction(inverse[i,j])) for j in range(codim)] for i in range(codim)]
            cache={}
            def entry(i,j):
                if (i,j) not in cache:
                    cache[i,j]=gram[i][j]-sum((gram[i][selection[a]]*inv[a][b]*gram[selection[b]][j] for a in range(codim) for b in range(codim)), Silver(0))
                return cache[i,j]
            indices=[i for i in range(len(gram)) if i not in selection and entry(i,i)>0]
            # Prefer obvious numerical violations, then certify exactly.
            proj=approximate-np.matmul(np.matmul(approximate[:,selection],np.asarray(inverse).astype(float)),approximate[selection,:])
            proposals=[]
            for i,j in combinations(indices,2):
                val=proj[i,j]
                sq=val*val/(proj[i,i]*proj[j,j])
                if val > 1e-9 or (sq < 1-1e-9 and min(abs(sq-v) for v in (0,.25,.5,.75,(2+2**.5)/4))>1e-8):
                    proposals.append((i,j))
            def failure(i,j):
                val=entry(i,j)
                sq=val**2/(entry(i,i)*entry(j,j))
                return (val>0 or (sq<1 and sq not in allowed)),val,sq
            bad=None
            for i,j in proposals:
                fails,val,sq=failure(i,j)
                if fails:
                    bad=(i,j,val,sq);break
            if bad is None:
                for i,j in combinations(indices,2):
                    fails,val,sq=failure(i,j)
                    if fails:
                        bad=(i,j,val,sq);break
            if bad is None:
                survivors.append({'walls':[labels[i] for i in selection],'facet_count':len(indices)})
            else:
                rejected+=1
                i,j,val,sq=bad
                witnesses.append({'face':[labels[x] for x in selection], 'forbidden_pair':[labels[i],labels[j]], 'pairing':val.to_json(),'normalized_pairing_square':sq.to_json()})
        row={'codimension':codim,'dimension':18-codim,'subsets':total,'elliptic_subsets':elliptic,'exact_rejections':rejected,'survivors':survivors,'witnesses':witnesses}
        print(json.dumps({k:v for k,v in row.items() if k!='witnesses'}),flush=True)
        rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, default=Path('data/enumerations/vinberg_h18_mixture_facets_2026_10_13.json'))
    parser.add_argument('--seam-codimensions', type=int, nargs='*', default=[])
    args = parser.parse_args()
    gram, labels, construction = mixture()
    facets = facet_screen(gram, labels)
    identify_facets(gram, labels, facets)
    inertia = silver_inertia(gram)
    result = {
        'status': 'EXPLORATORY_EXACT_FACET_SCREEN',
        'source': 'https://sfb701.math.uni-bielefeld.de/preprints/sfb14047.pdf',
        'parent_dimensions': [18,18], 'parent_facets': [37,24],
        'mixture_dimension': 18, 'mixture_facets': 40,
        'construction': construction, 'labels': labels,
        'gram': [[x.to_json() for x in row] for row in gram],
        'facets': facets,
        'exact_inertia': [inertia.positive, inertia.negative, inertia.zero],
        'seam_faces': seam_face_screen(gram, labels, construction, args.seam_codimensions),
        'coxeter_facet_count': sum(x['coxeter_angle_gate'] for x in facets),
        'forbidden_pair_histogram': dict(sorted(Counter(x['forbidden_pair_count'] for x in facets).items())),
        'scope': 'No worldwide novelty or new-dimensional example is asserted. The finite-volume H18 mixture is Vinberg 2015. Facet rejection uses exact projected Gram angles.'}
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in {'gram','facets','construction','scope','source','labels','seam_faces'}}, indent=2))


if __name__ == '__main__':
    main()
