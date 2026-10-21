#!/usr/bin/env python3
"""Bind spectator ground-field capacities, cusp-factor payments and changed-wall tests."""
import argparse
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/spectator_field_gates_and_boundary_replacements_20261017.md'
PREVIOUS = 'results/simple_cusp_dimension29_and_central_truncation_20261017_audit.json'
OUT = ROOT / 'results/spectator_field_gates_and_boundary_replacements_20261017_audit.json'
STEMS = ('fresh_sqrt5_multispectator_cap_capacity_20261017', 'fresh_sqrt5_multispectator_cap_capacity_20261017_review', 'fresh_sqrt5_multispectator_cap_local_cover_update_20261017', 'fresh_sqrt5_multispectator_cap_joint_ports_20261017', 'fresh_sqrt5_triple_a1_cap_capacity_20261017', 'fresh_sqrt5_multispectator_groundfield_cap_capacity_20261017', 'fresh_sqrt5_spectator_groundfield_gate_20261017', 'fresh_sqrt5_triple_groundfield_cap_local_cover_update_20261017', 'fresh_sqrt5_triple_groundfield_cap_joint_ports_20261017', 'fresh_nonarithmetic_release49_hexagon_20261017', 'fresh_nonarithmetic_release49_endpoint_gate_20261017', 'fresh_nonarithmetic_release8_ring_deformation_20261017', 'fresh_cusp_factor_payment_framework_20261017', 'fresh_cusp_link_shape_curvature_20261017', 'fresh_cusp_minimum_factor_rank9_dimension85_20261017', 'review_cusp_minimum_factor_rank9_dimension85_20261017', 'fresh_sqrt5_a3_d5_cap_witnesses_20261017', 'fresh_sqrt5_a3_d5_cap_local_cover_update_20261017', 'fresh_sqrt5_a3_d5_cap_joint_ports_20261017', 'fresh_cusp_spectator_role_payments_20261017', 'fresh_nonarithmetic_release49_62_replacement_20261017', 'fresh_nonarithmetic_release49_62_boundary_20261017')

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
    print('ARCHIVE_PREDECESSORS_PASS', len(frozen_hashes), flush=True)
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
    print('ARCHIVE_RECURSIVE_SOURCES_PASS', len(manifests), len(sources), flush=True)
    helper = 'scripts/audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017.py'
    sources[helper] = digest(ROOT / helper)
    def read(stem):
        return json.loads((ROOT / f'results/{stem}.json').read_text())
    from fractions import Fraction as F
    two = read('fresh_sqrt5_multispectator_cap_capacity_20261017')
    triple = read('fresh_sqrt5_triple_a1_cap_capacity_20261017')
    field = read('fresh_sqrt5_multispectator_groundfield_cap_capacity_20261017')
    independent = read('fresh_sqrt5_spectator_groundfield_gate_20261017')
    local = read('fresh_sqrt5_a3_d5_cap_local_cover_update_20261017')
    joint = read('fresh_sqrt5_a3_d5_cap_joint_ports_20261017')
    release = read('fresh_nonarithmetic_release49_endpoint_gate_20261017')
    deformation = read('fresh_nonarithmetic_release8_ring_deformation_20261017')
    cusp = read('fresh_cusp_factor_payment_framework_20261017')
    rank9 = read('fresh_cusp_minimum_factor_rank9_dimension85_20261017')
    shape = read('fresh_cusp_link_shape_curvature_20261017')
    chamber = read('fresh_nonarithmetic_release49_62_replacement_20261017')
    boundary = read('fresh_nonarithmetic_release49_62_boundary_20261017')
    assert two['counts']['columns'] == 555
    assert triple['counts']['columns'] == 93
    assert field['counts']['tightened_columns'] == 187
    assert independent['counts']['total_pair_tests'] == 79
    assert independent['counts']['newly_forbidden'] == 4
    assert local['counts']['cap_checks'] == 2056
    assert local['counts']['forests'] == 648
    assert local['counts']['block_orbits'] == 4169
    assert joint['face_multipliers'] == local['face_multipliers']
    # Independent Fraction-only composition of the saved cap and packing
    # certificates. This does not reconstruct their geometric provenance.
    from collections import defaultdict
    alpha = {i:F(a) for i,a in local['face_multipliers']}
    mu = {i:F(a) for i,a in joint['vertex_multipliers']}
    weights = defaultdict(F)
    for feature in local['positive_simplex_features']:
        assert sum(m for _,m in feature['block_orbits']) == feature['simplex_dimension']+1
        for j,m in feature['block_orbits']:
            weights[j] += m*alpha[feature['face']]
    literal = {c['index']:c for c in local['literal_block_orbits']}
    packing_terms = 0
    for cover in local['local_covers']:
        v = cover['vertex']; covered = defaultdict(F); mass = F()
        for term in cover['inherited_remapped_packing_mixture']+cover['singleton_packing_increments']:
            a = F(term['mass']); assert a > 0
            used = 0
            for j,mask in term['blocks']:
                c = literal[j]
                assert c['vertex'] == v and mask in c['literal_masks']
                assert not used & mask
                used |= mask; covered[j] += a
            mass += a; packing_terms += 1
        assert mass == F(cover['new_cover_mass']) and mu[v] >= mass
        for j,a in cover['required_marginals']:
            assert F(a) == weights[j] and covered[j] >= weights[j]
    assert packing_terms == 3712
    cap_count = 0
    for tag in local['cap_packages']:
        capdata = read(tag)
        if 'columns' in capdata:
            caprows = capdata['columns']
        else:
            assert tag == 'fresh_sqrt5_a3_d5_cap_witnesses_20261017'
            caprows = capdata['witnesses']; assert len(caprows) == 2
        for c in caprows:
            load = sum((F(m)*alpha.get(i,F()) for i,m in c['positive_face_coefficients']),F())
            assert load <= F(c['capacity'])*alpha.get(c['simplex_face'],F())
            cap_count += 1
    assert cap_count == 2056
    print('ARCHIVE_FRACTION_COMPOSITION_PASS', packing_terms, cap_count, flush=True)
    assert joint['counts']['total_inequalities'] == 2711432
    assert joint['counts']['violated_columns'] == 0
    assert joint['counts']['increased_vertex_masses'] == 0
    for stage in joint['complete_family_stages']:
        assert all(mu[v] >= F(a) for v,a in stage['vertex_required_masses'])
    assert sum(mu.values(),F()) == F(joint['normalizer'])
    assert F(joint['bound']) == F(joint['objective_numerator'])/F(joint['normalizer'])
    assert F(12) < F(joint['bound']) <= F(local['bound'])
    assert F(joint['bound']) < F(read('fresh_sqrt5_multispectator_cap_joint_ports_20261017')['bound'])
    assert release['counts']['finite_corner_recipes'] == 18
    assert release['counts']['obtuse_real_branches'] == 36
    assert release['counts']['nonfinite_contact_recipes'] == 878
    assert [r['boundary']['open_ends'] for r in deformation['variants']] == [18,15]
    assert [r['lost_old_finite_count'] for r in deformation['variants']] == [119,4]
    assert not deformation['new_exact_chamber'] and not deformation['new_finite_volume_example']
    assert cusp['two_factor_dimension_bound'] == 117
    assert rank9['dimension_bound'] == 85 and rank9['inherited_unrestricted_bound'] == 163
    assert rank9['enumerated_nonempty_partitions'] == 50378399
    assert rank9['excluded_range_partitions'] == 50330195
    assert rank9['maximum_cusp_factor_count'] is None
    assert shape['status'] == 'PASS_EXACT_CUSP_LINK_SHAPE_CURVATURE'
    assert chamber['total_wall_count'] == 64-len(chamber['deleted_old_walls'])+1 == 63
    assert len(chamber['native_gates']) == 62 and len(chamber['facet_witnesses']) == 63
    assert len(chamber['retained_old_finite_supports']) == 197
    assert len(chamber['new_finite_support_trials']) == 9
    assert len(chamber['new_finite_vertices']) == 2
    assert chamber['complete_finite_vertex_count'] == 199
    assert F(chamber['strict_escape']['norm'][0]) > 0
    assert all(F(a[1]) < 0 for a in chamber['strict_escape']['products'].values())
    assert F(chamber['inherited64_escape_new_product'][0]) > 0
    assert F(chamber['old60_escape_after_replacement']['products']['66'][0]) > 0
    assert boundary['before64']['open_ends'] == 15
    assert boundary['after63']['open_ends'] == 10
    assert len(boundary['after63']['rays']) == 32 and len(boundary['after63']['edges']) == 48
    assert not boundary['compact_side_faces_lost']
    assert not boundary['after63']['singular_triples'] and not boundary['after63']['nonsimple_rays'] and not boundary['after63']['ideal_rays']
    assert cusp['hypothesis'] == 'Every irreducible affine cusp factor has rank at least9.'
    for stem in STEMS:
        status = read(stem)['status']
        assert status.startswith('PASS') or (stem.startswith('fresh_nonarithmetic_') and status.startswith('BOUNDED_NUMERICAL')), (stem,status)
    small_checks = dict(
        scope='Exact counts and composition of saved certificates; no formal verification of the written geometric arguments or numerical construction diagnostics.',
        unrestricted_dimension_bound=163,
        all_simple_finite_volume_dimension_bound=29,
        compact_arithmetic_dimension_bound=21,
        minimum_rank9_dimension_bound=85,
        exact_excluded_range_partitions=50330195,
        arithmetic_unique_selected_cap_columns=1501,
        arithmetic_cap_checks_including_preserved_weaker_columns=2056,
        independent_fraction_packing_terms=packing_terms,
        arithmetic_fixed_field_tightened_columns=187,
        arithmetic_joint_bound=joint['bound'],
        release49_finite_endpoint_recipes=18, release49_excluded_real_branches=36,
        continuous_deformation_open_ends=[18,15],
        construction_control_walls=63,construction_control_finite_vertices=199,
        construction_control_S6_open_ends=10,construction_control_infinite_volume=True,
        new_example=None)
    reports = [ROOT / CHECKPOINT,
        ROOT / 'results/fresh_cusp_factor_payment_framework_20261017_root_review.md',
        ROOT / 'results/fresh_cusp_link_shape_curvature_20261017_root_review.md',
        ROOT / 'results/fresh_nonarithmetic_release8_ring_deformation_20261017_review.md',
        ROOT / 'results/fresh_nonarithmetic_release49_62_replacement_20261017_root_review.md']
    reports += [ROOT / f'results/{stem}.md' for stem in STEMS if not stem.startswith('review_cusp_')]
    reports.append(ROOT / 'results/fresh_cusp_minimum_factor_rank9_dimension85_20261017_review.md')
    links = []
    for report in reports:
        assert report.is_file(), report
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (report.parent / target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    return dict(status='PASS_FROZEN_PREDECESSORS_SPECTATOR_FIELDS_BOUNDARY_TESTS_AND_LOCAL_LINKS',
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
