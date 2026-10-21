#!/usr/bin/env python3
"""Exact composition of the mixed-simple 3n weight proof; no new face census.

The geometric lemmas are explicit bound inputs, not facts proved by a status
string.  This script checks their interfaces, source integrity, the complete
case partition, and the rational h-vector identities and all-rank thresholds.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
from math import comb
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_all_simple_dimension29_20261017'
WEIGHT = 'fresh_simple_cusp_vinberg_weight_frontier_20261017'
AFFINE = 'fresh_simple_cusp_affine_exclusion_20261017'
ACYCLE = 'fresh_simple_cusp_affine_a_quadrilateral_20261017'
AREVIEW = 'review_simple_cusp_affine_a_quadrilateral_20261017'
TREE = 'fresh_simple_cusp_affine_tree_completion_20261017'
TREVIEW = 'fresh_simple_cusp_affine_tree_terminal_audit_20261017'
FINITE = 'fresh_compact_dimension_bound_20261015'
FREVIEW = 'fresh_compact_dimension_bound_20261015_review'
TRANSFER = 'fresh_simple_finite_quad_vinberg_transfer_20261017'


def sha(path):
    return hashlib.sha256((ROOT/path).read_bytes()).hexdigest()


def read(stem, status):
    data = json.loads((ROOT/f'results/{stem}.json').read_text())
    assert data['status'] == status, stem
    return data


def check_manifests(initial):
    """Recheck every recursively named sources/source_hashes manifest.

    Hashing an inherited computation is integrity checking, not rerunning it.
    No historical open-scope field is overwritten or treated as closed by itself.
    """
    todo = list(initial)
    visited = set()
    manifests = []
    checked = {}
    while todo:
        path = todo.pop()
        if path in visited:
            continue
        visited.add(path)
        doc = json.loads((ROOT/path).read_text())
        docs = [doc]
        if isinstance(doc.get('prepared'), dict):
            docs.append(doc['prepared'])
        refs = 0
        for d in docs:
            for name in ('sources', 'source_hashes'):
                entries = d.get(name, {})
                assert isinstance(entries, dict), (path, name)
                for p, h in entries.items():
                    assert isinstance(h, str) and len(h) == 64, (path, p)
                    assert sha(p) == h, (path, p)
                    assert p not in checked or checked[p] == h
                    checked[p] = h
                    refs += 1
                    if p.endswith('.json'):
                        todo.append(p)
        manifests.append(dict(path=path, references=refs))
    return dict(manifests=sorted(manifests, key=lambda r:r['path']),
                checked_files=len(checked),
                checked_hash_references=sum(r['references'] for r in manifests))


def beta(n):
    return F(n*(n-6), 8) if n % 2 == 0 else F((n-1)*(n-5), 8)


def paired_slack(n, i):
    """Coefficient after averaging the two equal h_i,h_(n-i) entries."""
    return F(5, 2)*(comb(i, 2)+comb(n-i, 2))-comb(n, 2)-beta(n)


def run():
    weight = read(WEIGHT, 'PASS_AFFINE_3N_WEIGHT_INTERFACE_QUADRILATERAL_PAYMENT_OPEN')
    affine = read(AFFINE, 'PASS_EXACT_ALL_SIMPLE_AFFINE_EXCLUSION_DIMENSION61')
    cycle = read(ACYCLE, 'PASS_ALL_RANK_AFFINE_A_QUADRILATERAL_PAYMENT')
    areview = read(AREVIEW, 'PASS_INDEPENDENT_AFFINE_A_RECTANGLE_AND_KERNEL_REVIEW')
    tree = read(TREE, 'PASS_ALL_RANK_CLASSICAL_AFFINE_TREE_QUADRILATERAL_PAYMENT_N18')
    treview = read(TREVIEW, 'PASS_INDEPENDENT_TERMINAL_AND_TWO_E8_LITERAL_AUDIT')
    finite = read(FINITE, 'PASS_EXACT_DIMENSION29_EXTREMAL_VERTEX_REDUCTION_NOT_EXCLUSION')
    freview = read(FREVIEW, 'PASS_INDEPENDENT_COMPACT29_REDUCTION_ARITHMETIC_AND_SCOPE')
    assert weight['minimum_dimension'] == cycle['minimum_dimension'] == 14
    assert tree['minimum_dimension'] == 18
    assert affine['global_affine_order'] == 'Every connected affine subdiagram has exactly n nodes.'
    assert affine['minimal_indefinite_core_order_bound'] == 5
    assert areview['rectangle_checks'] == 256 and treview['fork_E8_exception']
    assert freview['typed_forests'] == finite['extremal_typed_count'] == 2126

    affine_costs = []
    for row in weight['affine_graph_checks']:
        assert row['order'] >= 14 and row['weight'] == row['budget'] == 3*row['order']
        affine_costs.append([row['type'], row['order'], row['weight']])
    assert len(affine_costs) == 348
    assert weight['all_rank_extra_end_costs'] == dict(cycle=0, path=6, single_fork=5, double_fork=4)
    finite_costs = []
    for row in finite['component_table']:
        r, w = row['rank'], row['weight']
        assert w <= 3*r and row['deficit'] == 3*r-w
        finite_costs.append([r, w, 3*r-w])
    assert len(finite_costs) == 29
    # Published finite classical count at rank r>=10 is <=(r-1)+20.
    # The all-rank deficit is 3r-(r+19)=2r-19>=1, not a cutoff inference.
    tail = dict(minimum_rank=10, upper_weight='r+19', deficit='2r-19',
                shifted_deficit_coefficients=[1, 2])  # r=10+t.

    h_checks = []
    for n in range(3, 101):
        m = n//2
        coefficients = []
        for i in range(n+1):
            c = paired_slack(n, i)
            expected = F(5, 2)*((i-m)**2 if n % 2 == 0 else (i-m)*(i-m-1))
            assert c == expected and c >= 0
            coefficients.append(str(c))
        assert F(coefficients[0]) == F(coefficients[-1]) > 0
        h_checks.append(dict(dimension=n, coefficients=coefficients,
                             endpoint_slack=str(2*F(coefficients[0]))))
    # Polynomial identities are also derived coefficientwise, for every t>=0.
    # 8(beta(n)-3n) = n^2-30n (even), n^2-30n+5 (odd).
    polynomials = dict(even=[0, 60, 4], odd=[36, 64, 4])  # ascending powers of t.
    for parity, n0, constant in [('even', 30, 0), ('odd', 31, 5)]:
        derived = [n0*n0-30*n0+constant, 4*n0-60, 4]
        assert derived == polynomials[parity] and all(x >= 0 for x in derived)
    threshold_checks = []
    for t in range(100):
        for parity, n in [('even', 30+2*t), ('odd', 31+2*t)]:
            p = polynomials[parity]
            assert 8*(beta(n)-3*n) == p[0]+p[1]*t+p[2]*t*t >= 0
            threshold_checks.append([n, str(beta(n)-3*n)])
    assert beta(30) == 90 and beta(29) < 87

    stems = [WEIGHT, AFFINE, ACYCLE, AREVIEW, TREE, TREVIEW, FINITE, FREVIEW]
    paths = [f'scripts/{TAG}.py', f'results/{TRANSFER}.md',
             f'results/{TRANSFER}_review.md',
             'results/fresh_simple_dimension29_literature_scope_20261017.md',
             'results/fresh_simple_cusp_affine_exclusion_20261017_review.md']
    for s in stems:
        paths.extend([f'scripts/{s}.py', f'results/{s}.json', f'results/{s}.md'])
    manifests = check_manifests([f'results/{s}.json' for s in stems])
    return dict(
        status='PASS_ALL_SIMPLE_FINITE_VOLUME_DIMENSION29_COMPOSITION',
        sources={p:sha(p) for p in paths}, recursive_integrity=manifests,
        dimension_bound=29, contradiction_from_dimension=30,
        minimum_dimension_for_combined_face_rule=18,
        arithmetic_hypothesis=False,
        scope='Finite-sided finite-volume hyperbolic Coxeter polytopes simple at every vertex, including ideal vertices. This theorem does not apply to polytopes having a nonsimple ideal vertex; unrestricted163 and compact-facet127 unchanged. No novelty/current-best claim. Internally reviewed proof composition, not a formalized proof or external refereeing.',
        geometric_inputs=[
            dict(case='finite vertex', result='Total original Vinberg corner weight <=3n',
                 source=f'results/{FINITE}.md'),
            dict(case='simple ideal vertex', result='Classical connected affine diagram; mixed rule total3n for n>=14',
                 source=f'results/{WEIGHT}.md'),
            dict(case='triangle', result='Total corner weight >=2, finite and ideal corners',
                 source=f'results/{WEIGHT}.md'),
            dict(case='quadrilateral, four finite corners', result='Total corner weight >=1 for n>=18',
                 source=f'results/{TRANSFER}.md'),
            dict(case='quadrilateral with affine-cycle corner', result='Total corner weight >=1 for n>=14',
                 source=f'results/{ACYCLE}.md'),
            dict(case='quadrilateral with classical affine-tree corner', result='Total corner weight >=1 for n>=18',
                 source=f'results/{TREE}.md'),
            dict(case='k-gon with k>=5', result='Nonnegative weights give total >=0>=5-k'),
        ],
        case_partition='At n>=18 every ideal vertex has a connected classical affine n-node diagram, of cycle or tree type; hence the three quadrilateral branches are exhaustive.',
        historical_open_obligations='The frozen weight interface and individual affine reports correctly record their then-open branches. This successor closes those branches by the separately bound transfer and completion proofs; it does not mutate the historical reports.',
        finite_transfer_review_scope=dict(
            root='Independent mathematical/source read of Vinberg sections6.5-11, geometric replacements, carrier sizes, high-bond connectedness, and determinant coefficient; no duplicate local-table census.',
            frontier='Separate primary-source review of the same transfer, including section6.11 core locality versus finite shared-branch classification; no duplicate local-table census.',
            transfer_sha256=sha(f'results/{TRANSFER}.md'),
            review_note=f'results/{TRANSFER}_review.md'),
        all_rank_finite_tail=tail, finite_component_budget_checks=finite_costs,
        inherited_affine_budget_checks=affine_costs,
        h_vector_argument=dict(f0='sum h_i', f2='sum binom(i,2) h_i',
            symmetry='h_i=h_(n-i)', nonnegative=True, endpoints='h_0=h_n=1',
            corner_incidence='sum_F k_F=binom(n,2) f_0',
            curvature='D5=5 f_2-binom(n,2)f_0',
            strict_lower_even='D5>n(n-6)f0/8',
            strict_lower_odd='D5>(n-1)(n-5)f0/8',
            strictness='The paired slack coefficients are nonnegative, and both endpoint h entries are1 with strictly positive coefficients.',
            all_rank_identity_even='D5-beta*f0=(5/2)sum(i-n/2)^2*h_i',
            all_rank_identity_odd='D5-beta*f0=(5/2)sum((i-n/2)^2-1/4)*h_i'),
        h_vector_identity_checks=h_checks,
        threshold_polynomials_times8_ascending=polynomials,
        threshold_checks=threshold_checks,
        contradiction='3n*f0 >= total corner payment >= D5 > beta_n*f0 >=3n*f0 for every n>=30.',
        verification_scope='Exact interface, source-integrity and rational-composition verification. Published finite classifications and separately reviewed geometric proofs remain explicit mathematical inputs; no duplicate full affine-tree or finite-face census is asserted.'
    )


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = json.loads(json.dumps(run()))
    out = ROOT/f'results/{TAG}.json'
    if args.verify:
        assert result == json.loads(out.read_text())
    else:
        out.write_text(json.dumps(result, sort_keys=True, indent=2)+'\n')
    print(result['status'], 'bound', result['dimension_bound'],
          'manifests', len(result['recursive_integrity']['manifests']),
          'files', result['recursive_integrity']['checked_files'])
