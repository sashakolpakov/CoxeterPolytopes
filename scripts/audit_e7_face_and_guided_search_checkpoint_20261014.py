#!/usr/bin/env python3
"""Bind new exact E7-face work and guided Gale cases; preserve the old archive."""
import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/e7_face_transfer_and_obstruction_guided_search_2026_10_14.md'
OUT = ROOT/'results/e7_face_transfer_and_obstruction_guided_search_2026_10_14_audit.json'
PREVIOUS = 'results/e8_coupled_caps_and_rank5_gale_2026_10_14_audit.json'
PREFIXES = ('fresh_h13_a11a2_completion_20261014',
            'fresh_h13_a11a2_lattice_gate_20261014',
            'fresh_h20_e7_intrinsic_caps_20261014',
            'fresh_h20_e7_global_escape_gate_20261014',
            'fresh_h9_obstruction_guided_gale_20261014')


@lru_cache(None)
def digest(path): return sha256(path.read_bytes()).hexdigest()


def record_refs(doc):
    if isinstance(doc, dict):
        if isinstance(doc.get('path'), str) and isinstance(doc.get('sha256'), str):
            yield doc['path'], doc['sha256']
        for v in doc.values(): yield from record_refs(v)
    elif isinstance(doc, list):
        for v in doc: yield from record_refs(v)


def run():
    prior = json.loads((ROOT/PREVIOUS).read_text())
    preserved = {}
    for section in ['source_hashes', 'top_artifact_hashes', 'report_hashes']:
        for p, h in prior[section].items():
            assert digest(ROOT/p) == h, ('previous frozen archive', p)
        preserved[section] = len(prior[section])
    top = sorted(p for p in (ROOT/'results').rglob('*.json')
                 if p.name.startswith(PREFIXES) or any(q.startswith(PREFIXES) for q in p.parts[:-1]))
    pending, seen, sources, manifests = list(top), set(), {}, []
    references = records_count = 0
    while pending:
        path = pending.pop()
        name = str(path.relative_to(ROOT))
        if name in seen: continue
        seen.add(name)
        doc = json.loads(path.read_text())
        if not isinstance(doc, dict): continue
        bound = list(doc.get('sources', {}).items())
        records = list(record_refs(doc))
        records_count += len(records)
        for source, expected in bound+records:
            assert re.fullmatch('[0-9a-f]{64}', expected), (name, source)
            target = ROOT/source
            assert target.is_file() and digest(target) == expected, (name, source)
            if source in sources: assert sources[source] == expected
            sources[source] = expected
            references += 1
            if target.suffix == '.json': pending.append(target)
        if bound or records:
            manifests.append(dict(path=name, status=doc.get('status'),
                                  source_references=len(bound), record_references=len(records)))
    reports = [ROOT/CHECKPOINT]+sorted(p for p in (ROOT/'results').glob('*.md') if p.name.startswith(PREFIXES))
    links = []
    for report in reports:
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'): continue
            resolved = (report.parent/target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    return dict(status='PASS_RECURSIVE_SOURCES_CASE_RECORDS_AND_LOCAL_LINKS',
        scope='Archive integrity only. Mathematical replays and independent-review scopes remain those of the individual reports.',
        checkpoint=CHECKPOINT, previous_archive=PREVIOUS,
        previous_archive_sha256=digest(ROOT/PREVIOUS),
        previous_frozen_sections_preserved=preserved,
        previous_entrypoints='Prior README/RESEARCH_PLAN hashes are historical snapshots; current entrypoints are bound below.',
        report_count=len(reports), local_link_count=len(links),
        source_manifest_count=len(manifests), source_reference_count=references,
        nested_record_reference_count=records_count, unique_source_count=len(sources),
        manifests=sorted(manifests, key=lambda v:v['path']),
        source_hashes=dict(sorted(sources.items())),
        top_artifact_hashes={str(p.relative_to(ROOT)):digest(p) for p in top},
        report_hashes={str(p.relative_to(ROOT)):digest(p) for p in reports},
        local_links=links,
        entrypoint_hashes={p:digest(ROOT/p) for p in ('README.md','RESEARCH_PLAN.md')},
        auditor_sha256=digest(Path(__file__)))


if __name__ == '__main__':
    p = argparse.ArgumentParser(); p.add_argument('--verify', action='store_true')
    args = p.parse_args(); result = run()
    if args.verify: assert json.loads(OUT.read_text()) == result
    else: OUT.write_text(json.dumps(result, indent=2)+'\n')
    print(result['status'], {k:result[k] for k in ('report_count','local_link_count',
        'source_manifest_count','source_reference_count','nested_record_reference_count','unique_source_count')})
