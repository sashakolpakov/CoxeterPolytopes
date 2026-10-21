#!/usr/bin/env python3
"""Bind cap orbits, cross-profile constraints, local covers and affine extensions."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/orbit_cap_incompatibilities_and_affine_tree_extensions_20261017.md'
PREVIOUS = 'results/cross_coface_caps_and_affine_corner_payments_20261017_audit.json'
OUT = ROOT / 'results/orbit_cap_incompatibilities_and_affine_tree_extensions_20261017_audit.json'
STEMS = (
    'fresh_nikulin_two_budget_frontier_20261015_dimension163',
    'fresh_sqrt5_cross_coface_cap_dual_support_20261017',
    'fresh_sqrt5_orbit_integral_cap_capacity_20261017',
    'fresh_sqrt5_cross_orientation_cap_capacity_20261017',
    'fresh_sqrt5_cap_profile_cliques_20261017',
    'fresh_sqrt5_orbit_and_orientation_caps_20261017_review',
    'fresh_sqrt5_orbit_cap_restricted_local_cones_20261017',
    'fresh_sqrt5_tightened_cap_local_cover_update_20261017',
    'fresh_sqrt5_profile_clique_local_cover_update_20261017',
    'fresh_sqrt5_profile_clique_joint_ports_20261017',
    'fresh_simple_cusp_affine_tree_short_ends_20261017',
    'fresh_nonarithmetic_two_central_caps_20261017',
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
    orbit = read('fresh_sqrt5_orbit_integral_cap_capacity_20261017')
    orient = read('fresh_sqrt5_cross_orientation_cap_capacity_20261017')
    clique = read('fresh_sqrt5_cap_profile_cliques_20261017')
    local = read('fresh_sqrt5_orbit_cap_restricted_local_cones_20261017')
    repaired = read('fresh_sqrt5_profile_clique_local_cover_update_20261017')
    joint = read('fresh_sqrt5_profile_clique_joint_ports_20261017')
    review = read('fresh_sqrt5_orbit_and_orientation_caps_20261017_review')
    short = read('fresh_simple_cusp_affine_tree_short_ends_20261017')
    construction = read('fresh_nonarithmetic_two_central_caps_20261017')
    dimension163 = read('fresh_nikulin_two_budget_frontier_20261015_dimension163')
    assert dimension163['combined_candidate_dimension_at_most'] == 163
    assert dimension163['first_excluded_dimension'] == 164
    assert orbit['counts']['columns'] == 815 and orbit['counts']['violations'] == 19
    assert orient['counts']['tightened'] == 215 and orient['counts']['violations'] == 42
    assert clique['counts']['columns'] == 851 and clique['counts']['new_clique_columns'] == 36
    assert clique['counts']['violations'] == 46
    assert len(clique['columns']) == 851
    assert repaired['complete_positive_simplex_forest_coverage']
    assert repaired['counts']['forests'] == 648 and repaired['counts']['block_orbits'] == 4169
    assert repaired['counts']['actual_quad_faces'] == 336
    assert repaired['counts']['raised_simplex_masses'] == 27
    assert repaired['counts']['added_singleton_terms'] == 61
    assert repaired['counts']['new_block_orbits'] == 0
    assert repaired['optimizer_invocations'] == 0
    assert 12 < F(repaired['bound']) <= F(local['bound'])
    assert 12 < F(joint['bound']) <= F(repaired['bound'])
    assert joint['face_multipliers'] == repaired['face_multipliers']
    assert joint['counts']['total_inequalities'] == 2711432
    assert joint['counts']['violated_columns'] == 98
    assert joint['counts']['tested_forests'] == 648
    assert [stage['violated_columns'] for stage in joint['complete_family_stages']] == [21,17,60]
    assert review['final_columns'] == 851 and review['final_violations'] == 46
    assert short['status'] == 'PASS_ALL_RANK_AFFINE_TREE_SHORT_OUTER_END_PAYMENT'
    assert len(short['remaining_central_positions_up_to_left_right_symmetry']) == 13
    for stem in STEMS:
        status = read(stem).get('status','')
        if stem == 'fresh_nonarithmetic_two_central_caps_20261017':
            assert status == 'BOUNDED_NUMERICAL_TWO_CENTER_RING_DISK_SCOUT'
        else:
            assert status.startswith('PASS'), (stem, status)
    small_checks = dict(scope='Exact package counts/ratios only; independent mathematical and arithmetic reviews have separate scope.',
        cap_columns=851, tightened_orientation_capacities=215, distinct_profile_cliques=36,
        original_full_dual_violations=46, covered_local_forests=648,
        updated_local_obstruction_bound=repaired['bound'], bound_is_optimum=False,
        final_combined_obstruction_bound=joint['bound'], joint_inequalities=2711432,
        new_joint_critical_columns=98,
        existing_unrestricted_bound=163,
        new_dimension_bound=None, new_example=None)
    reports = [ROOT / CHECKPOINT,
               ROOT / 'results/dimension163_proof_status_20261017.md',
               ROOT / 'results/fresh_simple_cusp_affine_tree_short_ends_20261017_review.md']
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
        source_hash_scope='Cumulative preserved predecessor files plus all current recursive source references; historical entrypoint snapshots excluded.',
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
