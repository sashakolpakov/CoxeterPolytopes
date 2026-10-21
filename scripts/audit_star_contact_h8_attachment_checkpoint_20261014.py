#!/usr/bin/env python3
"""Bind star contact families, source couplings and H8 attachment obstructions."""
import argparse
import json
import re
from functools import lru_cache
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/star_contact_families_and_h8_attachment_obstructions_2026_10_14.md'
OUT = ROOT / 'results/star_contact_families_and_h8_attachment_obstructions_2026_10_14_audit.json'
PREVIOUS = 'results/general_endpoint_caps_and_compact_case5_2026_10_14_audit.json'
PREFIXES = (
    'fresh_h9_star_zero_20261014',
    'fresh_h9_h8_facet_cover_20261014',
    'fresh_h9_alternative_mixed_couplings_20261014',
)


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
    top = sorted(p for p in (ROOT / 'results').glob('*.json') if p.name.startswith(PREFIXES))
    pending, seen, sources, manifests = top[:], set(), {}, []
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
        bound = [(p, h) for field in ('sources', 'source_hashes')
                 for p, h in doc.get(field, {}).items()]
        records = list(record_refs(doc))
        records_count += len(records)
        for source, expected in bound + records:
            assert re.fullmatch('[0-9a-f]{64}', expected), (name, source)
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
        previous_archive_sha256=digest(ROOT / PREVIOUS),
        previous_frozen_sections_preserved=preserved,
        previous_entrypoints='Historical snapshots; updated entrypoints are bound separately below.',
        report_count=len(reports), local_link_count=len(links), source_manifest_count=len(manifests),
        source_reference_count=references, nested_record_reference_count=records_count,
        unique_source_count=len(sources), manifests=sorted(manifests, key=lambda row: row['path']),
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
