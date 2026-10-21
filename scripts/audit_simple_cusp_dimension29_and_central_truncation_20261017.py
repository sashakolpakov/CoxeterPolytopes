#!/usr/bin/env python3
"""Bind the all-simple dimension29 proof and the central-truncation obstruction."""
import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/simple_cusp_dimension29_and_central_truncation_20261017.md'
PREVIOUS = 'results/orbit_cap_incompatibilities_and_affine_tree_extensions_20261017_audit.json'
OUT = ROOT / 'results/simple_cusp_dimension29_and_central_truncation_20261017_audit.json'
STEMS = (
    'fresh_nikulin_two_budget_frontier_20261015_dimension163',
    'fresh_simple_cusp_affine_tree_completion_20261017',
    'fresh_simple_cusp_affine_tree_terminal_audit_20261017',
    'fresh_all_simple_dimension29_20261017',
    'fresh_nonarithmetic_central_truncation_domain_20261017',
)

@lru_cache(None)
def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def run():
    prior = json.loads((ROOT / PREVIOUS).read_text())
    preserved = {}
    frozen_hashes = {}
    for section in ('source_hashes', 'top_artifact_hashes', 'report_hashes'):
        for name, expected in prior[section].items():
            assert digest(ROOT / name) == expected, ('frozen predecessor', name)
            assert name not in frozen_hashes or frozen_hashes[name] == expected
            frozen_hashes[name] = expected
        preserved[section] = len(prior[section])
    previous_auditor = 'scripts/audit_' + Path(PREVIOUS).name.removesuffix('_audit.json') + '.py'
    assert digest(ROOT / previous_auditor) == prior['auditor_sha256']
    frozen_hashes[previous_auditor] = prior['auditor_sha256']
    frozen_hashes[PREVIOUS] = digest(ROOT / PREVIOUS)
    top = [ROOT / f'results/{stem}.json' for stem in STEMS]
    pending = top[:] + [ROOT / A.RECOVERED_PROVENANCE]
    # Carry every preserved file forward, including predecessors no longer
    # reached from a new package's own dependency graph.
    sources, manifests, seen, bibliography = dict(frozen_hashes), [], set(), []
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
    def read(stem):
        return json.loads((ROOT / f'results/{stem}.json').read_text())
    dimension = read('fresh_all_simple_dimension29_20261017')
    tree = read('fresh_simple_cusp_affine_tree_completion_20261017')
    terminal = read('fresh_simple_cusp_affine_tree_terminal_audit_20261017')
    construction = read('fresh_nonarithmetic_central_truncation_domain_20261017')
    dimension163 = read('fresh_nikulin_two_budget_frontier_20261015_dimension163')
    assert dimension163['combined_candidate_dimension_at_most'] == 163
    assert dimension['dimension_bound'] == 29
    assert dimension['contradiction_from_dimension'] == 30
    assert dimension['minimum_dimension_for_combined_face_rule'] == 18
    assert dimension['arithmetic_hypothesis'] is False
    assert dimension['threshold_polynomials_times8_ascending'] == dict(even=[0,60,4],odd=[36,64,4])
    assert len(dimension['geometric_inputs']) == 7
    assert tree['minimum_dimension'] == 18
    assert terminal['fork_E8_exception']
    counts = construction['counts']
    assert counts['rings'] == 20 and counts['omitted_corner_cases'] == 80
    assert counts['inherited_low_contact_quadruples'] == 35840
    assert counts['one_high_real_branches'] == 1408
    assert counts['one_high_cosine_gaps'] == 41 and counts['two_high_rectangles'] == 6
    assert construction['unchanged_control']['walls'] == 64
    assert construction['unchanged_control']['S6_open_ends'] == 15
    for stem in STEMS:
        assert read(stem)['status'].startswith('PASS'), stem
    small_checks = dict(
        scope='Exact composition/count checks only; written geometric arguments and arithmetic reviews retain their separately stated scope.',
        unrestricted_dimension_bound=163,
        all_simple_finite_volume_dimension_bound=29,
        first_excluded_all_simple_dimension=30,
        nonsimple_ideal_vertices_covered=False,
        arithmetic_hypothesis=False,
        central_truncation_rings=20, central_truncation_one_high_branches=1408,
        central_truncation_two_high_rectangles=6,
        new_example=None)
    reports = [ROOT / CHECKPOINT,
        ROOT / 'results/dimension163_proof_status_20261017.md',
        ROOT / 'results/fresh_simple_finite_quad_vinberg_transfer_20261017.md',
        ROOT / 'results/fresh_simple_finite_quad_vinberg_transfer_20261017_review.md',
        ROOT / 'results/fresh_all_simple_dimension29_20261017_root_review.md',
        ROOT / 'results/fresh_simple_dimension29_literature_scope_20261017.md']
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
    return dict(status='PASS_FROZEN_PREDECESSORS_SIMPLE_DIMENSION29_AND_LOCAL_LINKS',
        scope='Archive integrity and final rational composition checks only. This does not formalize the geometric arguments or rerun every historical calculation.',
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
        manifests=sorted(manifests,key=lambda r:r['path']),
        source_hashes=dict(sorted(sources.items())),
        source_hash_scope='Cumulative preserved predecessor files and current recursive sources; historical entrypoint snapshots excluded.',
        top_artifact_hashes={str(p.relative_to(ROOT)):digest(p) for p in top},
        report_hashes={str(p.relative_to(ROOT)):digest(p) for p in reports},
        local_links=links,
        entrypoint_hashes={p:digest(ROOT/p) for p in ('README.md','RESEARCH_PLAN.md')},
        entrypoint_scope='Current entrypoints only; prior snapshots remain historical.',
        auditor_sha256=digest(Path(__file__)))


if __name__ == '__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--verify',action='store_true');args=ap.parse_args()
    result=run()
    if args.verify:
        assert result == json.loads(OUT.read_text())
    else:
        OUT.write_text(json.dumps(result,indent=2)+'\n')
    print(result['status'],{k:result[k] for k in
        ('report_count','local_link_count','source_manifest_count','source_reference_count','unique_source_count')},flush=True)
