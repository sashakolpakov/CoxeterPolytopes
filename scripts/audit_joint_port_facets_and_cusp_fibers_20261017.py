#!/usr/bin/env python3
"""Bind joint-port, cusp-fiber and restricted-dimension progress."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/joint_port_facets_and_cusp_fibers_20261017.md'
PREVIOUS = 'results/sharp_finite_edge_flags_and_packing_20261017_audit.json'
OUT = ROOT / 'results/joint_port_facets_and_cusp_fibers_20261017_audit.json'
STEMS = (
    'fresh_sqrt5_restricted_flag_local_cones_20261017',
    'review_sqrt5_restricted_flag_local_cover_20261017',
    'fresh_sqrt5_joint_triangle_square_ports_20261017',
    'fresh_sqrt5_joint_four_port_domain_20261017',
    'fresh_sqrt5_joint_four_port_diagnostic_20261017',
    'fresh_sqrt5_joint_five_port_domain_20261017',
    'fresh_sqrt5_joint_five_port_diagnostic_20261017',
    'review_sqrt5_joint_port_repairs_20261017',
    'fresh_sqrt5_coface_integer_base_scout_20261017',
    'fresh_cusp_finite_port_triangle_capacity_20261017',
    'fresh_simple_cusp_distance_bound_20261017',
    'fresh_simple_cusp_affine_exclusion_20261017',
    'fresh_nonarithmetic_release47_cycle_patch_20261017',
    'fresh_nonarithmetic_release47_replace12_frontier_20261017',
    'fresh_nonarithmetic_release47_nonzero_a2_bridge_20261017',
)



@lru_cache(None)
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
    repair = json.loads((ROOT / 'results/review_sqrt5_joint_port_repairs_20261017.json').read_text())
    assert repair['status'] == 'PASS_INDEPENDENT_FRACTION_JOINT_PORT_REPAIR_UNION'
    stages = repair['stages']
    assert [s['stage'] for s in stages] == [3,4,5]
    assert all(F(s['initial_bound']) >= F(s['bound']) > 12 for s in stages)
    assert F(repair['final_bound']) == F(stages[-1]['bound'])
    simple = json.loads((ROOT / 'results/fresh_simple_cusp_affine_exclusion_20261017.json').read_text())
    facet = json.loads((ROOT / 'results/fresh_simple_cusp_distance_bound_20261017.json').read_text())
    scout = json.loads((ROOT / 'results/fresh_sqrt5_coface_integer_base_scout_20261017.json').read_text())
    assert simple['status'] == 'PASS_EXACT_ALL_SIMPLE_AFFINE_EXCLUSION_DIMENSION61'
    assert simple['simple_cusp_dimension_bound'] == 61
    assert facet['compact_facet_ambient_dimension_bound'] == 127
    assert scout['status'] == 'PASS_EXACT_COFACE_SUMMED_INTEGER_BASE_SCOUT'
    assert scout['counts']['columns'] == 7504 and not scout['violations']
    assert not scout['alpha_changed'] and not scout['mu_changed']
    assert scout['input_bound'] == scout['output_bound'] == repair['final_bound']
    small_checks = dict(scope='Exact scope/count/ratio checks; separate geometric and certificate reviews remain authoritative.',
        final_restricted_joint_bound=repair['final_bound'],
        cumulative_critical_columns=repair['cumulative_verified_critical_columns'],
        simple_cusp_dimension_bound=61, compact_facet_dimension_bound=127,
        coface_summed_integer_base_columns=7504, coface_summed_violations=0)
    reports = [ROOT / CHECKPOINT,
               ROOT / 'results/fresh_sqrt5_joint_five_port_domain_20261017_review.md',
               ROOT / 'results/fresh_simple_cusp_distance_bound_20261017_review.md',
               ROOT / 'results/fresh_simple_cusp_affine_exclusion_20261017_review.md']
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
