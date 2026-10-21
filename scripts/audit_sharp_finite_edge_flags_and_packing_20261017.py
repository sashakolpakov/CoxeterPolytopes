#!/usr/bin/env python3
"""Bind the sharp-coefficient, flag-capacity and local-packing checkpoint."""
import argparse
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/sharp_finite_edge_flags_and_packing_20261017.md'
PREVIOUS = 'results/finite_edge_coefficient10_and_refined_arithmetic_20261017_audit.json'
OUT = ROOT / 'results/sharp_finite_edge_flags_and_packing_20261017_audit.json'
STEMS = (
    'fresh_six_simplex_ideal_frontier_20261017',
    'fresh_finite_edge_cusp_tradeoff_28over3_20261017',
    'review_finite_edge_28over3_20261017',
    'fresh_finite_vertex_t2_port_counterexample_20261017',
    'fresh_sqrt5_all_subset_port_packing_20261017',
    'fresh_sqrt5_intersecting_block_cliques_20261017',
    'fresh_sqrt5_pointed_simplex_base_capacity_20261017',
    'fresh_sqrt5_local_simplex_packing_cones_20261017',
    'fresh_sqrt5_all_local_simplex_packing_cones_20261017',
    'review_sqrt5_all_local_packing_cover_20261017',
    'fresh_sqrt5_flag_refined_base_capacity_20261017',
    'fresh_sqrt5_flag_incremental_local_cover_20261017',
    'fresh_sqrt5_flag_local_packing_full_model_20261017',
    'fresh_sqrt5_flag_local_packing_dual_support_20261017',
    'fresh_nonarithmetic_b_disk_bridge_20261017',
)


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run():
    prior = json.loads((ROOT / PREVIOUS).read_text())
    preserved = {}
    for section in ('source_hashes', 'top_artifact_hashes', 'report_hashes'):
        for name, expected in prior[section].items():
            assert digest(ROOT / name) == expected, ('frozen predecessor', name)
        preserved[section] = len(prior[section])
    top = [ROOT / f'results/{stem}.json' for stem in STEMS]
    pending = top[:] + [ROOT / A.RECOVERED_PROVENANCE]
    sources, manifests, seen, bibliography = {}, [], set(), []
    references = records_count = 0
    while pending:
        path = pending.pop()
        name = str(path.relative_to(ROOT))
        if name in seen:
            continue
        seen.add(name)
        doc = json.loads(path.read_text())
        if not isinstance(doc, dict):
            continue
        bound = []
        if 'database_sha256' in doc and name not in A.HISTORICAL_DATABASE_METADATA:
            assert isinstance(doc.get('database'), str), name
            bound.append((doc['database'], doc['database_sha256']))
        # Prepared model manifests have the same path-to-hash semantics as
        # top-level manifests; numerical solver arrays are not manifests.
        manifest_docs = [doc]
        if isinstance(doc.get('prepared'), dict):
            manifest_docs.append(doc['prepared'])
        for md in manifest_docs:
            for field in ('sources', 'source_hashes'):
                for source, value in md.get(field, {}).items():
                    if isinstance(value, dict):
                        assert field == 'sources' and 'title' in value, (name, source)
                        bibliography.append(dict(document=name, citation=source))
                        for key, expected in value.items():
                            if key.endswith('_sha256'):
                                key_path = key[:-7]
                                if isinstance(value.get(key_path), str) and (ROOT / value[key_path]).is_file():
                                    bound.append((value[key_path], expected))
                    else:
                        bound.append((source, value))
        records = list(A.record_refs(doc))
        records_count += len(records)
        for source, expected in bound + records:
            assert re.fullmatch('[0-9a-f]{64}', expected), (name, source)
            source = A.EXTERNAL_SOURCE_ALIASES.get(source, source)
            target = ROOT / source
            assert target.is_file() and digest(target) == expected, (name, source)
            assert source not in sources or sources[source] == expected
            sources[source] = expected
            references += 1
            if target.suffix == '.json':
                pending.append(target)
        if bound or records:
            manifests.append(dict(path=name, status=doc.get('status'),
                                  source_references=len(bound), record_references=len(records)))
    helper = 'scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
    sources[helper] = digest(ROOT / helper)
    lp = json.loads((ROOT / 'results/fresh_sqrt5_flag_local_packing_full_model_20261017.json').read_text())
    assert lp['solver_invocations'] == 1 and len(lp['certificates']) == 2
    assert lp['prepared']['dimension'] == 13 and lp['prepared']['target'] == 12
    primal = next(c for c in lp['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL')
    dual = next(c for c in lp['certificates'] if c['kind'] == 'EXACT_PARTIAL_MODEL_DUAL')
    counts = lp['prepared']['counts']
    for rows, limit in ((dual['face_multipliers'], counts['faces']),
                        (dual['vertex_multipliers'], counts['vertex_forests']),
                        (primal['rounded_variable_numerators'], counts['variables'])):
        assert len(dict(rows)) == len(rows)
        assert all(0 <= i < limit and F(a) >= 0 for i, a in rows)
    assert F(dual['normalizer']) > 0 and primal['integer_denominator'] > 0
    assert F(dual['bound']) == F(dual['objective_numerator']) / F(dual['normalizer'])
    assert F(primal['maximum']) == max(F(a) for _, a in primal['vertex_budgets'])
    assert F(primal['repair_delta']) >= 0
    assert 12 < F(dual['bound']) <= F(primal['maximum'])
    assert dual['excludes_target12'] and not primal['meets_target12']
    small_checks = dict(scope='Independent sign, index, exact-ratio and target checks; the full matrix certificate is separately replayed by its producer.',
                        exact_dual_bound=dual['bound'], exact_primal_maximum=primal['maximum'],
                        nonnegative_face_multipliers=len(dual['face_multipliers']),
                        nonnegative_vertex_multipliers=len(dual['vertex_multipliers']))
    reports = [ROOT / CHECKPOINT,
               ROOT / 'results/fresh_sqrt5_all_subset_port_packing_20261017_review.md',
               ROOT / 'results/fresh_sqrt5_flag_local_packing_full_model_20261017_review.md']
    reports += [ROOT / f'results/{stem}.md' for stem in STEMS]
    links = []
    for report in reports:
        assert report.is_file(), report
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (report.parent / target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    return dict(status='PASS_FROZEN_PREDECESSOR_RECURSIVE_MANIFESTS_AND_LOCAL_LINKS',
        scope='Archive integrity only. Exact mathematical certificates, independent arithmetic checks, source reviews and bounded numerical searches retain their separate stated scopes. This audit neither proves a geometric statement nor reruns every calculation.',
        checkpoint=CHECKPOINT, previous_archive=PREVIOUS,
        previous_archive_sha256=digest(ROOT / PREVIOUS),
        previous_frozen_sections_preserved=preserved,
        historical_nonbinding_database_metadata=sorted(A.HISTORICAL_DATABASE_METADATA),
        recovered_historical_source_aliases=A.EXTERNAL_SOURCE_ALIASES,
        recovered_provenance_sha256=digest(ROOT / A.RECOVERED_PROVENANCE),
        inherited_bibliographic_entries=bibliography,
        small_certificate_checks=small_checks,
        report_count=len(reports), local_link_count=len(links),
        source_manifest_count=len(manifests), source_reference_count=references,
        nested_record_reference_count=records_count, unique_source_count=len(sources),
        manifests=sorted(manifests, key=lambda r: r['path']),
        source_hashes=dict(sorted(sources.items())),
        top_artifact_hashes={str(p.relative_to(ROOT)): digest(p) for p in top},
        report_hashes={str(p.relative_to(ROOT)): digest(p) for p in reports},
        local_links=links,
        entrypoint_hashes={p: digest(ROOT / p) for p in ('README.md', 'RESEARCH_PLAN.md')},
        entrypoint_scope='Current entrypoints bound separately; previous archive entrypoints remain historical snapshots.',
        auditor_sha256=digest(Path(__file__)))


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--verify', action='store_true')
    args = p.parse_args()
    result = run()
    if args.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], {k: result[k] for k in
        ('report_count', 'local_link_count', 'source_manifest_count',
         'source_reference_count', 'unique_source_count')}, flush=True)
