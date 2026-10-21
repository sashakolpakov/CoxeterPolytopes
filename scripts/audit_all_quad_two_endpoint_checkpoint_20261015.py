#!/usr/bin/env python3
"""Bind all-quad domains, forest-sensitive weights and the dimension163 candidate."""
import argparse
import json
import re
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/all_quadrilateral_domains_and_two_endpoint_bound_2026_10_15.md'
OUT = ROOT / 'results/all_quadrilateral_domains_and_two_endpoint_bound_2026_10_15_audit.json'
PREVIOUS = 'results/quadratic_five_faces_and_residual_charge_2026_10_15_audit.json'
PREFIXES = (
    'fresh_sqrt5_all_quad_unions_20261015',
    'fresh_sqrt5_all_quad_contacts_20261015',
    'fresh_nikulin_two_budget_frontier_20261015',
    'fresh_silver_full_chamber_20261015',
    'fresh_sqrt5_uniform_face_weights_20261015',
    'fresh_sqrt5_arithmetic_weight_obstructions_20261015',
    'fresh_sqrt5_forest_sensitive_weights_20261015',
    'fresh_sqrt5_sealed_subcore_closure_20261015',
    'fresh_sqrt5_single_chain_completion_20261015',
    'fresh_sqrt5_singular_extension_frontier_20261015',
    'fresh_sqrt5_extended_forest_weights_20261015',
    'fresh_sqrt5_singular_attachment_closure_20261015',
)
EXCLUDED_DIAGNOSTICS = ()

EXTERNAL_SOURCE_ALIASES = {
    '/private/tmp/alexeev_2209_source/all_cases_in_2el_paper.txt':
        'data/references/alexeev_2209_09110v4_all_cases.txt',
}
RECOVERED_PROVENANCE = 'data/references/alexeev_2209_09110v4_provenance.json'


@lru_cache(None)
def digest(path):
    return sha256(path.read_bytes()).hexdigest()

def record_refs(doc):
    if isinstance(doc, dict):
        if isinstance(doc.get('path'), str) and isinstance(doc.get('sha256'), str):
            yield doc['path'], doc['sha256']
        for value in doc.values():
            yield from record_refs(value)
    elif isinstance(doc, list):
        for value in doc:
            yield from record_refs(value)

def run():
    prior = json.loads((ROOT / PREVIOUS).read_text())
    preserved = {}
    for section in ('source_hashes', 'top_artifact_hashes', 'report_hashes'):
        for path, expected in prior[section].items():
            assert digest(ROOT / path) == expected, ('previous frozen archive', path)
        preserved[section] = len(prior[section])
    top = sorted(p for p in (ROOT / 'results').glob('*.json') if p.name.startswith(PREFIXES) and str(p.relative_to(ROOT)) not in EXCLUDED_DIAGNOSTICS)
    pending, seen, sources, manifests = top[:] + [ROOT / RECOVERED_PROVENANCE], set(), {}, []
    references = records_count = 0
    bibliography = []
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
        if 'database_sha256' in doc:
            assert isinstance(doc.get('database'), str), name
            bound.append((doc['database'], doc['database_sha256']))
        for field in ('sources', 'source_hashes'):
            for source, value in doc.get(field, {}).items():
                if isinstance(value, dict):
                    # Some inherited datasets use "sources" for bibliography,
                    # not a path-to-hash manifest. Their containing JSON is bound.
                    assert field == 'sources' and 'title' in value, (name, source)
                    bibliography.append(dict(document=name, citation=source))
                    for key, expected in value.items():
                        if key.endswith('_sha256'):
                            path_key = key[:-7]
                            if isinstance(value.get(path_key), str) and (ROOT / value[path_key]).is_file():
                                bound.append((value[path_key], expected))
                else:
                    bound.append((source, value))
        records = list(record_refs(doc))
        records_count += len(records)
        for source, expected in bound + records:
            assert re.fullmatch('[0-9a-f]{64}', expected), (name, source)
            source = EXTERNAL_SOURCE_ALIASES.get(source, source)
            target = ROOT / source
            assert target.is_file() and digest(target) == expected, (name, source)
            if source in sources:
                assert sources[source] == expected
            sources[source] = expected
            references += 1
            if target.suffix == '.json':
                pending.append(target)
        if bound or records:
            manifests.append(dict(path=name, status=doc.get('status'),
                                  source_references=len(bound), record_references=len(records)))
    reports = [ROOT / CHECKPOINT] + sorted(
        p for p in (ROOT / 'results').glob('*.md') if p.name.startswith(PREFIXES))
    links = []
    for report in reports:
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (report.parent / target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    return dict(
        status='PASS_RECURSIVE_SOURCES_FROZEN_PREDECESSOR_AND_LOCAL_LINKS',
        scope='Archive integrity only. Exact results, independent reviews, bounded numerical diagnostics and unresolved work have separate scopes; this audit alone proves no mathematical claim.',
        checkpoint=CHECKPOINT, previous_archive=PREVIOUS,
        excluded_non_authoritative_producer_checkpoints=list(EXCLUDED_DIAGNOSTICS),
        previous_archive_sha256=digest(ROOT / PREVIOUS),
        previous_frozen_sections_preserved=preserved,
        previous_entrypoints='Historical snapshots; updated entrypoints are bound separately below.',
        report_count=len(reports), local_link_count=len(links), source_manifest_count=len(manifests),
        source_reference_count=references, nested_record_reference_count=records_count,
        unique_source_count=len(sources), manifests=sorted(manifests, key=lambda row: row['path']),
        inherited_bibliographic_entries=bibliography,
        recovered_historical_source_aliases=EXTERNAL_SOURCE_ALIASES,
        recovered_provenance_sha256=digest(ROOT / RECOVERED_PROVENANCE),
        source_hashes=dict(sorted(sources.items())),
        top_artifact_hashes={str(p.relative_to(ROOT)): digest(p) for p in top},
        report_hashes={str(p.relative_to(ROOT)): digest(p) for p in reports}, local_links=links,
        entrypoint_hashes={p: digest(ROOT / p) for p in ('README.md', 'RESEARCH_PLAN.md')},
        auditor_sha256=digest(Path(__file__)),
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--verify', action='store_true')
    args = parser.parse_args()
    result = run()
    if args.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result, indent=2) + '\n')
    print(result['status'], {k: result[k] for k in (
        'report_count', 'local_link_count', 'source_manifest_count',
        'source_reference_count', 'unique_source_count')})
