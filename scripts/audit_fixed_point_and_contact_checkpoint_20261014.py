#!/usr/bin/env python3
"""Check recursive provenance, per-case records and local checkpoint links."""
import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re
ROOT=Path(__file__).resolve().parents[1]
CHECKPOINT='results/fixed_point_symmetry_and_contact_replacements_2026_10_14.md'
OUT=ROOT/'results/fixed_point_symmetry_and_contact_replacements_2026_10_14_audit.json'
PREFIXES=('fresh_h9_fifteen_completion_20261014','fresh_h20_contact_closure_20261014',
          'fresh_h20_structured_extension_20261014','fresh_h9_cubic_corner_20261014',
          'fresh_h12_a12_source_20261014','fresh_h20_e8_cusp_completion_20261014')
@lru_cache(None)
def digest(path):return sha256(path.read_bytes()).hexdigest()
def record_refs(doc):
    if isinstance(doc,dict):
        if isinstance(doc.get('path'),str)and isinstance(doc.get('sha256'),str):
            yield doc['path'],doc['sha256']
        for value in doc.values():yield from record_refs(value)
    elif isinstance(doc,list):
        for value in doc:yield from record_refs(value)
def run():
    top=sorted(p for p in (ROOT/'results').rglob('*.json')
               if p.name.startswith(PREFIXES)or any(part.startswith(PREFIXES)for part in p.parts[:-1]))
    todo=list(top);seen=set();sources={};manifests=[];references=0;record_count=0
    while todo:
        path=todo.pop();name=str(path.relative_to(ROOT))
        if name in seen:continue
        seen.add(name);doc=json.loads(path.read_text())
        if not isinstance(doc,dict):continue
        bound=list(doc.get('sources',{}).items());records=list(record_refs(doc));record_count+=len(records)
        for source,expected in bound+records:
            assert re.fullmatch('[0-9a-f]{64}',expected),(name,source)
            target=ROOT/source
            assert target.is_file()and digest(target)==expected,(name,source)
            if source in sources:assert sources[source]==expected
            sources[source]=expected;references+=1
            if target.suffix=='.json':todo.append(target)
        if bound or records:manifests.append(dict(path=name,status=doc.get('status'),source_references=len(bound),record_references=len(records)))
    reports=[ROOT/CHECKPOINT]+sorted(p for p in (ROOT/'results').glob('*.md')if p.name.startswith(PREFIXES))
    links=[]
    for report in reports:
        for target in re.findall(r'\]\(([^)]+)\)',report.read_text()):
            if '://'in target or target.startswith('#'):continue
            resolved=(report.parent/target.split('#')[0]).resolve()
            assert resolved.is_file()or resolved==OUT,(report,target)
            links.append(dict(report=str(report.relative_to(ROOT)),target=target))
    return dict(status='PASS_RECURSIVE_SOURCES_CASE_RECORDS_AND_LOCAL_LINKS',
        scope='Archive integrity only; mathematical replay and independent-review scopes remain those of the individual reports.',
        checkpoint=CHECKPOINT,report_count=len(reports),local_link_count=len(links),
        source_manifest_count=len(manifests),source_reference_count=references,
        nested_record_reference_count=record_count,unique_source_count=len(sources),
        manifests=sorted(manifests,key=lambda x:x['path']),source_hashes=dict(sorted(sources.items())),
        top_artifact_hashes={str(p.relative_to(ROOT)):digest(p)for p in top},
        report_hashes={str(p.relative_to(ROOT)):digest(p)for p in reports},local_links=links,
        entrypoint_hashes={p:digest(ROOT/p)for p in('README.md','RESEARCH_PLAN.md')},
        auditor_sha256=digest(Path(__file__)))
if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();out=run()
    if args.verify:assert json.loads(OUT.read_text())==out
    else:OUT.write_text(json.dumps(out,indent=2)+'\n')
    print(out['status'],{k:out[k]for k in('report_count','local_link_count','source_manifest_count','source_reference_count','nested_record_reference_count','unique_source_count')})
