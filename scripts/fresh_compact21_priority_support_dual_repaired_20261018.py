#!/usr/bin/env python3
"""Repair retained numerical multipliers by monotone exact local clipping.

No numerical package or optimizer is imported or run.  The saved 2,805-row
support is small; literal weak duality is replayed on all 6,121 columns.
"""
import argparse
from collections import defaultdict
import fcntl
from fractions import Fraction as F
from hashlib import file_digest
import json
import os
from pathlib import Path
import resource
import signal
import sys
import time
import fresh_compact21_priority_support_exact_dual_20261018 as D

ROOT = D.ROOT
TAG = 'fresh_compact21_priority_support_dual_repaired_20261018'
NUMERIC = 'results/fresh_compact21_priority_support_exact_dual_20261018.json'


def repair():
    started = time.monotonic()
    source = json.loads((ROOT / NUMERIC).read_text())
    rows = list(D.S.rows())
    grid = 10**12
    face, forest = {}, {}
    for i, encoded in source['retained_numeric_dual']:
        mass = max(F(), -F.from_float(float.fromhex(encoded)))
        scaled = mass * grid
        numerator = scaled.numerator // scaled.denominator
        if rows[i][2] == 'forest' and numerator * scaled.denominator < scaled.numerator:
            numerator += 1
        if numerator:
            (face if rows[i][2] == 'face' else forest)[i] = F(numerator, grid)
    capacity, demand = defaultdict(F), defaultdict(F)
    incident = defaultdict(list)
    for i, mass in forest.items():
        row = rows[i][0]
        assert row[93045] == -1
        for j, a in row.items():
            if j != 93045:
                assert a >= 0
                capacity[j] += mass * a
    for i, mass in face.items():
        for j, a in rows[i][0].items():
            assert j != 93045 and a <= 0
            if a:
                demand[j] -= mass * a
                incident[j].append(i)
    initial = face.copy()
    violations = sorted((j for j in demand if demand[j] > capacity[j]),
                        key=lambda j: (-(demand[j] - capacity[j]), j))
    for j in violations:
        # Decreasing any face mass only decreases every column demand.
        # Prefer the smallest mass so isolated numerical noise disappears first.
        for i in sorted(incident[j], key=lambda i: (face[i], i)):
            excess = demand[j] - capacity[j]
            if excess <= 0:
                break
            take = min(face[i], excess / -rows[i][0][j])
            face[i] -= take
            for k, a in rows[i][0].items():
                demand[k] += take * a
        assert demand[j] <= capacity[j]
    assert all(demand[j] <= capacity[j] for j in demand)
    normalizer = max(F(1), sum(forest.values(), F()))
    weights = {i: -mass / normalizer for i, mass in face.items() if mass}
    weights.update({i: -mass / normalizer for i, mass in forest.items()})
    bound = sum((y * rows[i][1] for i, y in weights.items()), F())
    metadata = D.forest_metadata()
    masses = sorted(({'support_row': i, 'positive_mass': str(-y), **metadata[i]}
                     for i, y in weights.items() if i in metadata),
                    key=lambda r: (-F(r['positive_mass']), r['support_row']))
    removed = {i: initial[i] - face[i] for i in initial if initial[i] != face[i]}
    data = {
        'status': 'PASS_EXACT_REPAIRED_PRIORITY_SUPPORT_DUAL',
        'exact_lower_bound': str(bound), 'excess_above_40': str(bound - 40),
        'proves_support_infeasible_at_40': bound > 40,
        'exact_dual_rows': [[i, str(y)] for i, y in sorted(weights.items())],
        'normalizer': str(normalizer), 'rounding_grid': grid,
        'initial_violated_columns': len(violations),
        'clipped_face_rows': len(removed),
        'objective_loss_before_normalization': str(sum(
            (mass * -rows[i][1] for i, mass in removed.items()), F())),
        'largest_forest_masses': masses[:16],
        'A1_D5_E7_E8_forest_mass': str(sum(
            (F(r['positive_mass']) for r in masses
             if sorted(r['types']) == ['A1', 'D5', 'E7', 'E8']), F())),
        'all_nonzero_forest_masses': masses,
        'support_rows': len(rows), 'support_columns': 6121,
        'support_nonzeros': sum(len(row) for row, _, _ in rows),
        'full_model_primal_verified': False, 'new_dimension_bound': None,
        'scope': 'Exact lower certificate for the stated support successor only; '
                 'an obstruction to this restricted weight family, not to polytopes.',
        'algorithm': 'Round face masses down and forest masses up, monotonically '
                     'clip overloaded column demands, normalize forest mass, '
                     'then check every literal weak-duality column.',
    }
    D.verify_dual(data)
    data['resources'] = {
        'elapsed_seconds': round(time.monotonic() - started, 3),
        'peak_rss_bytes': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss *
                          (1 if sys.platform == 'darwin' else 1024),
        'numeric_optimizer_run': False,
    }
    assert data['resources']['peak_rss_bytes'] < 1024**3
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    signal.signal(signal.SIGTERM, D.S.G.interrupt)
    signal.signal(signal.SIGINT, D.S.G.interrupt)
    signal.signal(signal.SIGALRM, D.S.G.interrupt)
    signal.alarm(59)
    target = ROOT / f'results/{TAG}.json'
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        before = D.S.G.resources()
        if args.verify:
            data = json.loads(target.read_text())
            D.verify_dual(data)
            for path, digest in data['sources'].items():
                with (ROOT / path).open('rb') as stream:
                    assert file_digest(stream, 'sha256').hexdigest() == digest
        else:
            assert not target.exists()
            data = repair()
            data['initial_resources'] = before
            data['sources'] = {}
            for path in (f'scripts/{TAG}.py', NUMERIC, f'scripts/{D.TAG}.py',
                         f'scripts/{D.S.TAG}.py', D.S.M.MODEL, D.S.SEED):
                with (ROOT / path).open('rb') as stream:
                    data['sources'][path] = file_digest(stream, 'sha256').hexdigest()
            encoded = json.dumps(data, indent=2) + '\n'
            assert len(encoded.encode()) < 512 * 1024
            temporary = target.with_suffix('.json.tmp')
            try:
                with temporary.open('x') as stream:
                    stream.write(encoded)
                    stream.flush()
                    os.fsync(stream.fileno())
                os.replace(temporary, target)
            finally:
                temporary.unlink(missing_ok=True)
        D.S.G.resources(start=False)
    signal.alarm(0)
    print(json.dumps({
        'status': data['status'], 'exact_lower_bound': data['exact_lower_bound'],
        'decimal_lower_bound': float(F(data['exact_lower_bound'])),
        'above_40': data['proves_support_infeasible_at_40'],
        'clipped_face_rows': data['clipped_face_rows'],
        'dual_terms': len(data['exact_dual_rows']),
        'new_dimension_bound': None,
    }))


if __name__ == '__main__':
    main()
