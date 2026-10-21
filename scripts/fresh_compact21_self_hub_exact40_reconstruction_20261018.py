#!/usr/bin/env python3
"""Solver-free reconstruction of five nearly integral values, with exact replay."""
import argparse
from fractions import Fraction as F
import fcntl
from hashlib import file_digest
import json
import os
from pathlib import Path
import resource
import signal
import sys
import threading
import time
import fresh_compact21_self_hub_epigraph_minimum_20261018 as S

M, E, C, ROOT, Z = S.M, S.E, S.C, S.ROOT, S.Z
TAG = 'fresh_compact21_self_hub_exact40_reconstruction_20261018'
INPUT = 'results/fresh_compact21_self_hub_epigraph_minimum_20261018.json'


def sha(path):
    with (ROOT/path).open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def reconstruct():
    with (ROOT/INPUT).open() as stream:
        original = json.load(stream)
    for path, expected in original['sources'].items():
        assert sha(path) == expected
    raw = {j: float.fromhex(v) for j, v in original['retained_numeric_bonus_primal']}
    dropped = sorted(j for j, v in raw.items() if 0 < v < 1e-8)
    assert len(dropped) == 3
    near_integer = {j: round(v) for j, v in raw.items() if abs(v-round(v)) < 1e-8}
    assert len(near_integer) == 5 and sum(v != 0 for v in near_integer.values()) == 2
    values = {j: F(v) for field in ('positive_bonus_values', 'positive_pool_values')
              for j, v in original[field] if j not in dropped}
    values.update({j: F(v) for j, v in near_integer.items() if v})
    base, coverage, _, models, profiles, grouped, _, pools, _ = M.build()
    rule = E.rule_from_data(original)
    rule.values = {j: v for j, v in values.items() if j < Z}
    pool_columns = set(pools.values())
    pool_values = {j: v for j, v in values.items() if j in pool_columns}
    violations = []
    for i, (row, rhs, kind) in enumerate(base):
        if kind != 'face':
            continue
        lhs = sum((a*values.get(j, F()) for j, a in row.items()), F())
        if lhs > rhs:
            violations.append((i, str(lhs-rhs)))
    assert not violations, violations[:12]
    maximum, attainer = M.exact_forest_maximum(rule, models, profiles, grouped, pools, pool_values)
    assert maximum <= 40, (maximum, attainer)
    return {
        'status': 'PASS_EXACT_SELF_HUB_COMPLETE_MODEL_AT_MOST_40',
        'exact_primal_upper': str(maximum), 'attaining_forest': attainer,
        'primal_at_most_40': True, 'primal_strictly_below_40': maximum < 40,
        'transfer_value': original['transfer_value'],
        'positive_bonus_values': [[j, str(v)] for j, v in sorted(rule.values.items())],
        'positive_pool_values': [[j, str(v)] for j, v in sorted(pool_values.items())],
        'pool_variable_order': original['pool_variable_order'],
        'native_coverage': coverage,
        'finite_component_types': len(models),
        'exact_hub_contexts': sum(map(len, grouped.values())),
        'capped_context_groups': len(grouped),
        'dropped_numerical_zero_entries': [
            {'global_column': j, 'raw_hex': raw[j].hex(),
             'original_rounded_value': dict(original['positive_bonus_values']+original['positive_pool_values'])[j]}
            for j in dropped],
        'integral_reconstructions': [
            {'global_column': j, 'raw_hex': raw[j].hex(), 'exact_value': str(v)}
            for j, v in sorted(near_integer.items())],
        'reconstruction': 'Reconstruct only the five raw values within 1e-8 of an integer: '
                          'three zeros and pool values 36 and 39. '
                          'All other exact rounded bonuses and transfer are unchanged. '
                          'The threshold selects a candidate only; every face inequality '
                          'and finite-forest maximum is checked over rational numbers.',
        'all_enlarged_native_face_rows_exactly_verified': True,
        'all_finite_forests_exactly_verified': True,
        'universal_coverage_composition_reviewed': False,
        'new_dimension_bound': None,
        'sources': original['sources'] | {INPUT: sha(INPUT), 'scripts/'+TAG+'.py': sha('scripts/'+TAG+'.py'),
                    'scripts/compact21_self_hub_pool_model_20261018.py':
                    sha('scripts/compact21_self_hub_pool_model_20261018.py'),
                    'results/compact21_self_hub_full_model_formulation_20261018.md':
                    sha('results/compact21_self_hub_full_model_formulation_20261018.md')},
    }


def worker():
    owner = int(os.environ['COXETER_SUPPORT_SCOUT_OWNER'])
    assert os.getppid() == owner
    def watch():
        while True:
            if os.getppid() != owner:
                os.killpg(os.getpgrp(), signal.SIGKILL)
            time.sleep(.2)
    threading.Thread(target=watch, daemon=True).start()
    data = reconstruct()
    data['worker_peak_rss_bytes'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*(1 if sys.platform == 'darwin' else 1024)
    path = Path(os.environ['COXETER_FULL_FACE_TEMP'])
    size = C.bounded_write(path, data)
    with path.open('rb') as stream:
        digest = file_digest(stream, 'sha256').hexdigest()
    S.event('certificate_ready', status=data['status'], certificate_bytes=size,
            certificate_sha256=digest, exact_primal_upper=data['exact_primal_upper'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker', action='store_true', help=argparse.SUPPRESS)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    if args.worker:
        worker(); return
    for sig in (signal.SIGTERM, signal.SIGINT, signal.SIGALRM):
        signal.signal(sig, C.H.S.G.interrupt)
    signal.alarm(59)
    target = ROOT/f'results/{TAG}.json'
    temporary = target.with_suffix('.candidate.tmp')
    output_tmp = target.with_suffix('.json.tmp')
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            if args.verify:
                C.H.S.G.resources(); C.H.S.process_tree_rss()
                with target.open() as stream:
                    data = json.load(stream)
                for path, digest in data['sources'].items():
                    assert sha(path) == digest
                assert reconstruct() == {k: v for k, v in data.items()
                                         if k not in ('guard', 'worker_peak_rss_bytes')}
                S.verify(data)
                C.H.S.G.resources(start=False)
            else:
                assert not target.exists() and not temporary.exists() and not output_tmp.exists()
                # Reuse the reviewed signal-safe process-tree/RSS guard, but
                # direct its child to this solver-free worker, never an LP.
                S.__file__, S.TAG = __file__, TAG
                data = S.guarded(lock.fileno(), temporary)
                C.bounded_write(output_tmp, data); os.replace(output_tmp, target)
        finally:
            temporary.unlink(missing_ok=True); output_tmp.unlink(missing_ok=True)
    signal.alarm(0)
    print(json.dumps({k: data[k] for k in ('status', 'exact_primal_upper', 'attaining_forest', 'new_dimension_bound')}))


if __name__ == '__main__':
    main()
