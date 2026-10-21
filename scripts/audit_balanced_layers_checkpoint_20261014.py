#!/usr/bin/env python3
"""Check frozen checkpoint links and recursively declared source provenance.

This archive audit does not substitute for mathematical certificate replays.
"""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/balanced_layers_and_compact_carrier_closure_2026_10_14.md'
OUT = ROOT/'results/balanced_layers_and_compact_carrier_closure_2026_10_14_audit.json'
PREFIXES = ('fresh_h20_balanced_boundary_20261014',
            'fresh_h20_balanced_golden_20261014',
            'fresh_h20_balanced_integer',
            'fresh_h20_balanced_half_20261014_review',
            'fresh_h9_bugaenko_closure_20261014')


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run():
    top = sorted(p for p in (ROOT/'results').glob('*.json')
                 if p.name.startswith(PREFIXES))
    todo, seen, sources, manifests = list(top), set(), {}, []
    references = 0
    while todo:
        path = todo.pop()
        name = str(path.relative_to(ROOT))
        if name in seen:
            continue
        seen.add(name)
        doc = json.loads(path.read_text())
        bound = doc.get('sources', {})
        if not bound:
            continue
        for source, expected in bound.items():
            assert re.fullmatch('[0-9a-f]{64}', expected), (name, source)
            target = ROOT/source
            assert target.is_file() and digest(target) == expected, (name, source)
            if source in sources:
                assert sources[source] == expected
            sources[source] = expected
            references += 1
            if target.suffix == '.json':
                todo.append(target)
        manifests.append(dict(path=name, status=doc.get('status'),
                              source_references=len(bound)))
    reports = [ROOT/CHECKPOINT]+sorted(p for p in (ROOT/'results').glob('*.md')
                                     if p.name.startswith(PREFIXES))
    links = []
    for report in reports:
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (report.parent/target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    return dict(status='PASS_RECURSIVE_SOURCE_HASHES_AND_CHECKPOINT_LINKS',
                scope='Provenance and local links only. Mathematical replays and independent reviews are recorded separately; numerical scouts remain numerical.',
                checkpoint=CHECKPOINT,
                report_count=len(reports), local_link_count=len(links),
                source_manifest_count=len(manifests),
                source_reference_count=references, unique_source_count=len(sources),
                manifests=sorted(manifests, key=lambda x: x['path']),
                source_hashes=dict(sorted(sources.items())),
                top_certificate_hashes={str(p.relative_to(ROOT)): digest(p) for p in top},
                report_hashes={str(p.relative_to(ROOT)): digest(p) for p in reports},
                local_links=links,
                entrypoint_hashes={p: digest(ROOT/p) for p in ('README.md', 'RESEARCH_PLAN.md')},
                auditor_sha256=digest(Path(__file__)))


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--verify', action='store_true')
    args = ap.parse_args()
    out = run()
    if args.verify:
        assert json.loads(OUT.read_text()) == out
    else:
        OUT.write_text(json.dumps(out, indent=2)+'\n')
    print(out['status'], {k: out[k] for k in ('report_count', 'local_link_count',
          'source_manifest_count', 'source_reference_count', 'unique_source_count')})
