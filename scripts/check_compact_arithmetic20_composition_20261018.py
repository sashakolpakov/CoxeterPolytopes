#!/usr/bin/env python3
"""Stream-safe integrity and exact arithmetic for the new compact20 composition.

This binds the already reviewed field/coverage proofs. It does not reconstruct
the historical Gram queues or formalize the universal geometric arguments.
"""
import argparse
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
import json
from pathlib import Path
import re
import signal
import resource
import sys
import time
import run_lean_batches_signal_safe_20261018 as safe

ROOT = Path(__file__).resolve().parents[1]
TAG = 'compact_arithmetic20_composition_20261018'
HIGH = 'results/fresh_arithmetic_dimension21_composition_20261016_bound21.json'
FIELDS = 'results/fresh_compact21_single_field_20261016.json'
PRIMAL = 'results/fresh_compact21_self_hub_exact40_reconstruction_20261018.json'
INDEPENDENT = 'results/fresh_compact21_self_hub_exact40_independent_check_20261018.json'
REVIEW = 'results/review_compact_arithmetic20_composition_root_20261018.md'
SMALL_LIMIT = 2*1024**2


def small_json(relative):
    path = ROOT/relative
    assert path.stat().st_size < 8*1024**2, relative
    with path.open() as stream:
        return json.load(stream)


def source_maps(path):
    if path.stat().st_size < SMALL_LIMIT:
        with path.open() as stream:
            doc = json.load(stream)
        return [doc.get(key, {}) for key in ('sources', 'source_hashes')]
    # Some old files contain a single enormous JSON line. Scan fixed-size
    # byte chunks, then decode only bounded source-map objects. Any embedded
    # named hash maps are checked too; no scientific array is materialized.
    pattern = re.compile(rb'"(?:sources|source_hashes)"\s*:\s*\{')
    result, seen, offset, tail = [], set(), 0, b''
    with path.open('rb') as stream:
        while True:
            block = stream.read(65536)
            if not block:
                break
            data, origin = tail+block, offset-len(tail)
            for match in pattern.finditer(data):
                start = origin+match.end()-1
                if start in seen:
                    continue
                seen.add(start)
                with path.open('rb') as object_stream:
                    object_stream.seek(start)
                    payload = b''
                    while True:
                        chunk = object_stream.read(65536)
                        assert chunk and len(payload)+len(chunk) < SMALL_LIMIT
                        payload += chunk
                        try:
                            mapping, _ = json.JSONDecoder().raw_decode(payload.decode(errors='ignore'))
                            break
                        except json.JSONDecodeError:
                            continue
                assert isinstance(mapping, dict)
                result.append(mapping)
            offset += len(block)
            tail = data[-256:]
    return result


def run():
    high, fields, primal, independent = map(small_json, (HIGH, FIELDS, PRIMAL, INDEPENDENT))
    assert high['status'] == 'PASS_ARITHMETIC_REFLECTION_DIMENSION_COMPOSITION'
    assert high['compact_arithmetic_dimension_upper_bound'] == 21
    expected_high = {(n, k) for n in range(22, 30)
                     for k in ('Q(sqrt2)', 'Q(sqrt5)', 'Q(cos(2pi/7))')}
    assert {(r['dimension'], r['field']) for r in high['compact_high_dimension_field_branches']} == expected_high
    assert all(r['excluded_by_upper_bound'] < r['dimension'] for r in high['compact_high_dimension_field_branches'])
    assert fields['status'] == 'PASS_COMPACT_DIMENSION21_REDUCTION_TO_SINGLE_OPEN_FIELD'
    assert fields['remaining_fields'] == ['Q(sqrt5)']
    expected_fields = {'Q', 'Q(sqrt2)', 'Q(sqrt3)', 'Q(sqrt5)', 'Q(sqrt6)',
                       'Q(sqrt2,sqrt3)', 'Q(sqrt2,sqrt5)',
                       'K7', 'K9', 'K11', 'K15', 'K16', 'K20'}
    assert len(fields['branches']) == 13 and {r['field'] for r in fields['branches']} == expected_fields
    for row in fields['branches']:
        assert row['dimension'] == 21
        if row['field'] == 'Q':
            assert row['status'] == 'EXCLUDED_PUBLISHED'
        elif row['field'] == 'Q(sqrt5)':
            assert row['status'] == 'OPEN'
        else:
            assert row['status'] == 'EXCLUDED_RESEARCH_CANDIDATE' and row['upper_bound'] == 20
    assert primal['status'] == 'PASS_EXACT_SELF_HUB_COMPLETE_MODEL_AT_MOST_40'
    assert F(primal['exact_primal_upper']) == 40 and primal['primal_at_most_40']
    assert primal['all_enlarged_native_face_rows_exactly_verified'] and primal['all_finite_forests_exactly_verified']
    assert primal['finite_component_types'] == 68 and primal['exact_hub_contexts'] == 729
    assert independent['status'] == 'PASS_INDEPENDENT_SELF_HUB_FULL_PRIMAL'
    assert independent['native_faces'] == 1979 and independent['native_spectator_cases'] == 102917
    assert independent['exact_hub_tuples'] == 729 and independent['forest_offers'] == 850
    assert F(independent['exact_primal_upper']) == 40 and independent['primal_at_most_40']
    assert F(independent['minimum_payment']) >= 1 and not independent['new_model_helpers_imported']
    pair_coefficients = []
    for i in range(11):
        choose = lambda k: k*(k-1)//2
        coefficient = 5*(choose(i)+choose(21-i))-420
        assert coefficient == 80+5*(i-10)*(i-11) and coefficient >= 80
        pair_coefficients.append(coefficient)
    assert pair_coefficients[0]-80 == 550
    sources = [HIGH, FIELDS, PRIMAL, INDEPENDENT, REVIEW,
               'results/review_compact21_self_hub_pool_model_root_20261018.md',
               'results/review_compact21_self_hub_pool_model_batch_guard_20261018.md',
               'scripts/check_compact21_self_hub_primal_independent_20261018.py',
               str(Path(__file__).relative_to(ROOT))]
    hashes, visited, references, hashed_bytes = {}, set(), 0, 0

    def digest(relative):
        nonlocal hashed_bytes
        assert isinstance(relative, str)
        path = (ROOT/relative).resolve()
        assert path.is_relative_to(ROOT) and path.is_file(), relative
        if relative not in hashes:
            with path.open('rb') as stream:
                hashes[relative] = file_digest(stream, 'sha256').hexdigest()
            hashed_bytes += path.stat().st_size
        return hashes[relative]

    def check(relative):
        nonlocal references
        if not relative.endswith('.json') or relative in visited:
            return
        visited.add(relative)
        for mapping in source_maps(ROOT/relative):
            assert isinstance(mapping, dict)
            for source, expected in mapping.items():
                assert isinstance(expected, str) and digest(source) == expected, (relative, source)
                references += 1
                check(source)
        if len(visited) % 8 == 0:
            safe.resources(start=False)
            peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
            assert peak < 1024**3

    for relative in sources:
        digest(relative); check(relative)
    assert independent['candidate_sha256'] == digest(PRIMAL)
    return {
        'status': 'PASS_COMPACT_ARITHMETIC_DIMENSION20_COMPOSITION',
        'compact_arithmetic_dimension_upper_bound': 20,
        'finite_volume_arithmetic_dimension_upper_bound': 21,
        'noncompact_arithmetic_dimension_upper_bound': 21,
        'noncompact_dimension20_excluded_published': True,
        'dimension21_field_branches': 13,
        'higher_dimension_field_branches': len(expected_high),
        'golden_field_exact_vertex_upper': '40',
        'universal_compact21_total_weight_lower': '40*v+550',
        'symmetric_h_pair_coefficients': pair_coefficients,
        'recursive_manifest_count': len(visited),
        'source_reference_count': references,
        'distinct_hashed_files': len(hashes), 'streamed_hashed_bytes': hashed_bytes,
        'sources': {p: digest(p) for p in sources},
        'scope': 'Exact final composition with separately reviewed internal field, '
                 'whole-contact coverage, angle and forest proofs. Large manifests '
                 'are streamed. Historical Gram queues are not rerun, no external '
                 'human review or full geometric Lean verification is claimed.',
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = safe.resources(); signal.alarm(59); started = time.monotonic()
        data = run(); safe.resources(start=False)
        target = ROOT/'results'/(TAG+'.json')
        if args.verify:
            assert small_json(str(target.relative_to(ROOT))) == data
        else:
            assert not target.exists()
            target.write_text(json.dumps(data, indent=2)+'\n')
        signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'compact_arithmetic_dimension_upper_bound',
        'dimension21_field_branches', 'higher_dimension_field_branches',
        'recursive_manifest_count', 'distinct_hashed_files', 'streamed_hashed_bytes')}))
