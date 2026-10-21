#!/usr/bin/env python3
"""Small independent scalar/quantifier audit, not a new geometric queue replay."""
import argparse
import json
from fractions import Fraction as F
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TAG = 'fresh_arithmetic_dimension23_scope_audit_20261016'


def digest(path):
    return sha256((ROOT / path).read_bytes()).hexdigest()


def contexts(rank):
    # Residual components exclude A1 and I2(5). Their possible total ranks
    # are exactly zero and every integer >=2 (A2 and A3 suffice).
    return {(a, b) for a in range(rank + 1) for b in range(rank // 2 + 1)
            if rank - a - 2*b == 0 or rank - a - 2*b >= 2}


def threshold(n):
    paired = [F((n-i)*(n-i-1) + i*(i-1), 4) for i in range(n+1)]
    minimum = min(paired)
    assert paired[0] > minimum
    answer = 5*minimum - F(n*(n-1), 2)
    closed = F(n*(n-6), 8) if n % 2 == 0 else F((n-1)*(n-5), 8)
    assert answer == closed
    return minimum, answer


def run():
    paths = {f'scripts/{TAG}.py'}

    def read(name):
        path = f'results/{name}.json'
        paths.add(path)
        return json.loads((ROOT/path).read_text())

    terminal = read('fresh_sqrt5_nonsingular_rank29_fast_20261015_review')
    assert terminal['status'] == 'PASS_INDEPENDENT_COMPLETE_FINITE_ATTACHMENT_QUEUE_DIMENSION29'
    assert (terminal['processed'], terminal['marked_states'], terminal['pending'], terminal['complete_queue']) == (5172, 5172, 0, True)
    assert (terminal['exact_profiles'], terminal['literal_maps'], terminal['corner_patterns'], terminal['weight_obligations']) == (51599, 10085, 8280, 2070)
    for suffix in ('.py', '.md'):
        paths.add(('scripts/' if suffix == '.py' else 'results/') + 'fresh_sqrt5_nonsingular_rank29_fast_20261015_review' + suffix)

    # These exact reviews are named computational premises, not rerun here.
    inherited = [
        'fresh_sqrt5_quad_exhaustion_20261015_incidence_review',
        'fresh_sqrt5_quad_contacts_20261015_review',
        'fresh_sqrt5_triangle_weights_20261015_review',
        'fresh_sqrt5_dotted_parameters_20261015_review',
        'fresh_sqrt5_dotted_rank29_exhaustion_20261015_review',
        'fresh_sqrt2_quad_exhaustion_20261015_review',
        'fresh_sqrt2_weight_bound_20261015_review',
        'fresh_compact_global_bound_20261015_review',
    ]
    dependencies = []
    for name in inherited:
        data = read(name)
        assert data['status'].startswith('PASS')
        paths.add(f'results/{name}.md')
        dependencies.append(dict(name=name, status=data['status'], arithmetic_repeated_here=False))

    # Type count includes A2=I2(3), B2=I2(4) only once.
    type_counts = {}
    for n in (23, 24, 29):
        types = ([f'A{r}' for r in range(1, n+1)] +
                 [f'B{r}' for r in range(2, n+1)] +
                 [f'D{r}' for r in range(4, n+1)] +
                 ['E6', 'E7', 'E8', 'F4', 'H3', 'H4', 'I2(5)', 'I2(6)', 'I2(10)'])
        assert len(types) == len(set(types)) == 3*n+5
        type_counts[str(n)] = len(types)

    audits = []
    cases_checked = 0
    for dimension, weight_name, binder_name, expected in [
        (24, 'fresh_sqrt5_dimension24_i25_cross_20261015_review',
         'fresh_sqrt5_dimension23_coverage_20261015', (1897, 60735, F(267527,5000))),
        (23, 'fresh_sqrt5_hub_count_weights_20261016_review_n23_cap2',
         'fresh_sqrt5_hub_count_coverage_20261016_bound22', (1764, 52029, F(485211,10000))),
    ]:
        weight, binder = read(weight_name), read(binder_name)
        assert weight['status'].startswith('PASS_INDEPENDENT_EXACT_')
        assert binder['complete_coverage']
        assert binder['compact_dimension_upper_bound'] == dimension-1
        assert binder['excluded_dimensions'] == list(range(dimension, 30))
        assert binder['arithmetic_ground_field'] == 'Q(sqrt5)'
        assert weight['dimension'] == dimension
        expected_contexts = contexts(dimension)
        fronts = weight['forest_maxima']
        assert len(fronts) == len(expected_contexts)
        assert {tuple(row['context']) for row in fronts} == expected_contexts
        maximum = max(F(row['maximum']) for row in fronts)
        minimum, target = threshold(dimension)
        assert maximum == F(weight['vertex_maximum']) == expected[2] <= target
        assert F(weight['strict_face_threshold']) == target
        assert weight['finite_component_types'] == type_counts[str(dimension)]
        payments = binder['payments']
        assert len(payments) == binder['zero_touch_rows'] == 2661
        actual_cases = feasible = 0
        min_payment = None
        for row in payments:
            rank = dimension + 2 - row['nodes']
            required = contexts(rank) if rank >= 0 else set()
            seen = [(c['external_A1'], c['external_I25']) for c in row['cases']]
            assert len(seen) == len(set(seen)) and set(seen) == required
            feasible += bool(seen)
            for case in row['cases']:
                assert case['dimension'] == dimension
                assert len(case['corners']) == 4
                corners = list(map(F, case['corners']))
                assert min(corners) >= 0
                value = sum(corners, F(0))
                assert value == F(case['total']) >= 1
                min_payment = value if min_payment is None else min(min_payment, value)
                actual_cases += 1
        assert (feasible, actual_cases) == expected[:2]
        assert binder['dimension_feasible_rows'] == feasible
        assert binder['dimension_feasible_cases'] == actual_cases
        assert binder['over_rank_rows'] == len(payments)-feasible
        assert binder['literal_map_checks'] == 471
        vertex = binder['baseline_and_vertex_bounds']
        assert F(vertex['maximum_vertex_weight']) == maximum
        assert F(vertex['strict_face_weight_threshold']) == target
        assert len(vertex['baseline_totals']) == 76
        assert min(map(F, vertex['baseline_totals'])) >= 2
        cases_checked += actual_cases
        audits.append(dict(dimension=dimension, component_types=type_counts[str(dimension)],
                           exact_vertex_contexts=len(fronts), zero_touch_rows=len(payments),
                           dimension_feasible_rows=feasible, spectator_cases=actual_cases,
                           minimum_saved_payment=str(min_payment), vertex_maximum=str(maximum),
                           strict_threshold=str(target), gap=str(target-maximum)))

    # Independently enumerate the logical field/dimension cartesian product.
    compositions = []
    for bound, name in [(23, 'fresh_arithmetic_dimension_composition_20261015_bound23'),
                        (22, 'fresh_arithmetic_hub_dimension_composition_20261016_bound22')]:
        data = read(name)
        assert data['status'] == 'PASS_ARITHMETIC_REFLECTION_DIMENSION_COMPOSITION'
        assert data['finite_volume_arithmetic_reflection_dimension_upper_bound'] == bound
        field_bounds = {'Q(sqrt2)': 20, 'Q(sqrt5)': bound, 'Q(cos(2pi/7))': 20}
        rows = data['compact_high_dimension_field_branches']
        assert len(rows) == 3*(29-bound)
        assert {(r['dimension'], r['field']) for r in rows} == {
            (n, field) for n in range(bound+1,30) for field in field_bounds}
        for row in rows:
            assert row['excluded_by_upper_bound'] == field_bounds[row['field']] < row['dimension']
        assert data['noncompact_arithmetic']['upper_bound'] == 21
        assert data['noncompact_arithmetic']['excluded_dimension'] == 20
        compositions.append(dict(bound=bound, branches=len(rows)))

    scalar = []
    for n in range(21, 30):
        minimum, target = threshold(n)
        scalar.append(dict(dimension=n, paired_f2_minimum=str(minimum), strict_weight_threshold=str(target)))
    # Pure counting identity: capped context never changes these actual denominators.
    for a in range(30):
        for b in range(15):
            if a:
                assert sum((F(1,a) for _ in range(a)), F(0)) == 1
            if b:
                assert sum((F(1,2*b) for _ in range(2*b)), F(0)) == 1

    paths.update({
        'scripts/fresh_sqrt5_dimension27_coverage_20261015.py',
        'scripts/fresh_sqrt5_hub_count_weights_20261016_review.py',
        'scripts/fresh_sqrt5_hub_count_coverage_20261016.py',
        'scripts/fresh_arithmetic_hub_dimension_composition_20261016.py',
        'results/arithmetic_dimension23_and_complete_face_coverage_2026_10_15.md',
        'results/fresh_arithmetic_dimension_composition_20261015.md',
        'results/fresh_sqrt5_dimension23_coverage_20261015.md',
    })
    return dict(
        status='PASS_SCOPED_HOSTILE_GEOMETRY_QUANTIFIER_AND_EXACT_SCALAR_AUDIT',
        sources={p:digest(p) for p in sorted(paths)},
        independently_checked_saved_payment_sums=cases_checked,
        coverage_scalar_audits=audits, dimension_thresholds=scalar,
        finite_component_counts=type_counts, field_compositions=compositions,
        inherited_computational_premises=dependencies,
        terminal_replay_confirmation=dict(processed=5172, pending=0, singular_obligations=2070,
            prior_independent_generation_and_fresh_saved_replay='PASS; exit code 0; original session23308 closed'),
        provenance_correction='The fast producer had generation and terminal structural DB/map checks. Its independent reviewer, not the fast producer, had a fresh complete saved replay. The frozen dimension27 MD wording conflates these; the final checkpoint is accurate.',
        scope='Hostile mathematical/source audit plus independent standard-library Fraction scalar, spectator-quantifier, saved-payment-sum and composition checks. No new Lanner/union catalogue, queue traversal, graph-weight reconstruction, or independent forest optimizer is claimed. Those remain the explicitly bound reviewed inputs. No novelty, external validation or nonarithmetic dimension claim.')


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    data=run();path=ROOT/f'results/{TAG}.json'
    if args.verify:
        assert data == json.loads(path.read_text())
    else:
        path.write_text(json.dumps(data,indent=2)+'\n')
    print(data['status'],data['independently_checked_saved_payment_sums'])
