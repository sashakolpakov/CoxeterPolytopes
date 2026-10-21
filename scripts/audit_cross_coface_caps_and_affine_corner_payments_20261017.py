#!/usr/bin/env python3
"""Bind cross-coface capacities, affine corner proofs and construction exclusions."""
import argparse
from fractions import Fraction as F
from functools import lru_cache
from hashlib import sha256
import json
from pathlib import Path
import re

import audit_finite_edge_coefficient10_refined_arithmetic_checkpoint_20261017 as A

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/cross_coface_caps_and_affine_corner_payments_20261017.md'
PREVIOUS = 'results/joint_port_facets_and_cusp_fibers_20261017_audit.json'
OUT = ROOT / 'results/cross_coface_caps_and_affine_corner_payments_20261017_audit.json'
STEMS = (
    'fresh_sqrt5_labelled_coface_base_packing_20261017',
    'fresh_simplex_base_packing_cone_20261017',
    'review_simplex_base_packing_cone_20261017',
    'fresh_sqrt5_complete_coface_base_cones_20261017',
    'fresh_sqrt5_singleton_coface_covers_20261017',
    'fresh_sqrt5_cross_coface_cap_capacity_20261017',
    'review_sqrt5_cross_coface_full_base_witness_20261017',
    'fresh_sqrt5_cross_coface_cap_capacity_20261017_review',
    'fresh_sqrt5_cross_coface_cap_repair_20261017',
    'fresh_sqrt5_cross_coface_cap_full_model_20261017',
    'fresh_simple_cusp_vinberg_weight_frontier_20261017',
    'fresh_simple_cusp_affine_a_quadrilateral_20261017',
    'review_simple_cusp_affine_a_quadrilateral_20261017',
    'fresh_simple_cusp_affine_tree_corridor_20261017',
    'fresh_nonarithmetic_release47_central_patch_20261017',
    'fresh_nonarithmetic_central_ring_contact_gate_20261017',
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
    cap = read('fresh_sqrt5_cross_coface_cap_capacity_20261017')
    repair = read('fresh_sqrt5_cross_coface_cap_repair_20261017')
    model = read('fresh_sqrt5_cross_coface_cap_full_model_20261017')
    cones = read('fresh_sqrt5_complete_coface_base_cones_20261017')
    singles = read('fresh_sqrt5_singleton_coface_covers_20261017')
    acycle = read('fresh_simple_cusp_affine_a_quadrilateral_20261017')
    corridor = read('fresh_simple_cusp_affine_tree_corridor_20261017')
    rings = read('fresh_nonarithmetic_central_ring_contact_gate_20261017')
    assert cap['status'] == 'PASS_EXACT_CROSS_COFACE_A1_CAP_CAPACITIES'
    assert cap['counts']['columns'] == 241 and cap['counts']['violations'] == 53
    assert repair['status'] == 'PASS_EXACT_CROSS_COFACE_CAP_DUAL_REPAIR'
    assert F(repair['uniform_singleton_bound']) <= F(repair['bound']) < 12
    assert F(repair['initial_bound']) > 12
    assert repair['counts']['inherited_faces'] == 82806
    assert repair['counts']['inherited_signed_columns'] == 94481
    assert repair['counts']['new_cap_columns'] == 241
    assert model['solver_invocations'] == 1
    assert model['prepared']['counts']['faces'] == 82806
    assert model['prepared']['counts']['vertex_forests'] == 6558
    assert model['prepared']['counts']['variables'] == 428768
    assert model['prepared']['counts']['new_signed_variables'] == 241
    assert model['prepared']['counts']['new_finite_positive_cuts'] == 217
    certs = {c['kind']: c for c in model['certificates']}
    lower = F(certs['EXACT_PARTIAL_MODEL_DUAL']['bound'])
    upper = F(certs['EXACT_PARTIAL_MODEL_PRIMAL']['maximum'])
    assert lower <= upper
    assert cones['status'] == 'PASS_EXACT_COMPLETE_COFACE_BASE_PACKING_CONES'
    assert cones['counts']['complete_inequalities'] == 510656 and not cones['violations']
    assert singles['counts']['slices'] == 34520 and singles['counts']['positive_slices'] == 870
    assert acycle['status'] == 'PASS_ALL_RANK_AFFINE_A_QUADRILATERAL_PAYMENT'
    assert acycle['counts']['rectangle_corners'] == 256
    assert corridor['status'] == 'PASS_ALL_RANK_AFFINE_TREE_CENTRAL_CORRIDOR_PAYMENT'
    assert rings['counts'] == dict(rings=20, low_contact_quadruples=35840, high_real_branches=64)
    small_checks = dict(scope='Exact package-scope/count/ratio checks; independent proof and arithmetic reviews remain separate.',
        cross_coface_columns=241, previous_witness_violations=53,
        repaired_fixed_witness_bound=repair['bound'], repaired_bound_is_optimum=False,
        enlarged_finite_model_lower_bound=str(lower), enlarged_finite_model_upper_bound=str(upper),
        enlarged_model_solver_invocations=1,
        singleton_coface_slices=34520, affine_a_rectangle_checks=256,
        long_tree_corridor_closed=True, full_quadrilateral_payment_open=True,
        exact_central_ring_exclusions=20, new_dimension_bound=None, new_example=None)
    reports = [ROOT / CHECKPOINT,
               ROOT / 'results/fresh_sqrt5_cross_coface_cap_repair_20261017_review.md',
               ROOT / 'results/fresh_sqrt5_cross_coface_cap_full_model_20261017_review.md',
               ROOT / 'results/fresh_sqrt5_complete_coface_base_cones_20261017_review.md',
               ROOT / 'results/fresh_simple_cusp_affine_tree_corridor_20261017_review.md',
               ROOT / 'results/fresh_nonarithmetic_central_ring_contact_gate_20261017_review.md']
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
