#!/usr/bin/env python3
"""Check checkpoint provenance and local links; not a mathematical re-verifier."""
import argparse,hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/noncompact_section_exclusion_and_changed_chamber_search_2026_10_14_audit.json'
CHECKPOINT='results/noncompact_section_exclusion_and_changed_chamber_search_2026_10_14.md'
PREFIXES=('fresh_h20_rank3_global_sections_20261014','fresh_h20_cusp_coverage_20261014',
 'fresh_h20_quartic_orbit_completion_20261014','fresh_h9_nonpolygon_construction_20261014',
 'fresh_h20_quartic_orbit_cap_gate_20261014','fresh_h20_quartic_rational_cap_gate_20261014',
 'fresh_h20_compact_section_coverage_20261014','fresh_cusp_coverage_calibration_20261014')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(ROOT))
def run():
 manifests=[];references=0;sources={};certificates={}
 for p in sorted((ROOT/'results').glob('*.json')):
  if not p.name.startswith(PREFIXES):continue
  x=json.loads(p.read_text());certificates[rel(p)]=sha(p)
  if not x.get('sources'):continue
  for name,digest in x['sources'].items():
   assert isinstance(digest,str) and re.fullmatch('[0-9a-f]{64}',digest),(p,name)
   target=ROOT/name;assert target.is_file() and sha(target)==digest,(p,name)
   if name in sources:assert sources[name]==digest
   sources[name]=digest;references+=1
  manifests.append({'path':rel(p),'status':x.get('status'),'source_references':len(x['sources'])})
 reports=[ROOT/CHECKPOINT]+sorted(p for p in (ROOT/'results').glob('*.md')if p.name.startswith(PREFIXES))
 links=[];report_hashes={}
 for p in reports:
  report_hashes[rel(p)]=sha(p)
  for destination in re.findall(r'\]\(([^)]+)\)',p.read_text()):
   if '://' in destination or destination.startswith('#'):continue
   target=(p.parent/destination.split('#')[0]).resolve()
   assert target.is_file() or target==OUT,(p,destination)
   links.append({'report':rel(p),'destination':destination})
 return {'status':'PASS_CHECKPOINT_SOURCE_HASHES_AND_LOCAL_REPORT_LINKS',
  'scope':'Checks archive provenance and local report links only. It does not rerun or independently prove the mathematical certificates. Numerical scouts remain numerical.',
  'checkpoint':CHECKPOINT,'report_count':len(reports),'local_link_count':len(links),
  'source_manifest_count':len(manifests),'source_reference_count':references,'unique_source_count':len(sources),
  'manifests':manifests,'source_hashes':sources,'certificate_hashes':certificates,
  'report_hashes':report_hashes,'local_links':links,
  'entrypoint_hashes':{name:sha(ROOT/name)for name in['README.md','RESEARCH_PLAN.md']},
  'auditor_sha256':sha(Path(__file__))}
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args();x=run()
 if args.verify:assert json.loads(OUT.read_text())==x
 else:OUT.write_text(json.dumps(x,indent=2)+'\n')
 print(x['status'],{k:x[k]for k in['report_count','local_link_count','source_manifest_count','source_reference_count','unique_source_count']})
