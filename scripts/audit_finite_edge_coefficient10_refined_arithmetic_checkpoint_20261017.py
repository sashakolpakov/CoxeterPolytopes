#!/usr/bin/env python3
"""Preserve the archive and bind coefficient10, refined arithmetic and boundary packages."""
import argparse
import json
import re
from functools import lru_cache
from fractions import Fraction
from hashlib import sha256
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT = 'results/finite_edge_coefficient10_and_refined_arithmetic_20261017.md'
OUT = ROOT / 'results/finite_edge_coefficient10_and_refined_arithmetic_20261017_audit.json'
PREVIOUS = 'results/finite_edge_triangles_forced_bases_and_boundary_2026_10_17_audit.json'
STEMS = (
    'fresh_finite_edge_seven_simplex_20261017',
    'fresh_seven_simplex_affine_and_d8_20261017',
    'fresh_finite_edge_cusp_tradeoff10_20261017',
    'fresh_all_ideal_seven_simplex_realization_20261017',
    'fresh_all_ideal_seven_simplex_literature_20261017',
    'fresh_sqrt5_actual_provenance_dual_20261017',
    'fresh_sqrt5_vertex_refined_base_capacity_20261017',
    'fresh_sqrt5_integer_simplex_port_packing_20261017',
    'fresh_sqrt5_vertex_refined_integer_union_20261017',
    'fresh_sqrt5_labelled_packing_full_model_20261017',
    'fresh_sqrt5_labelled_packing_dual_support_20261017',
    'fresh_nonarithmetic_b_release_costs_20261017',
    'fresh_nonarithmetic_b_release_incidence_20261017',
    'fresh_nonarithmetic_b_release52_patch_20261017',
    'fresh_nonarithmetic_b_release52_continuous_20261017',
    'fresh_nonarithmetic_b_release52_left_cluster_20261017',
    'fresh_nonarithmetic_b_release30_replacement_20261017',
)


EXCLUDED_DIAGNOSTICS = ()
HISTORICAL_DATABASE_METADATA = {
    'data/enumerations/fresh_sqrt5_nonsingular_attachment_search_20261015_batch610.json',
    'data/enumerations/fresh_sqrt5_nonsingular_attachment_search_20261015_batch1097.json',
}
# These immutable snapshots explicitly retain historical active-DB metadata;
# their frozen reports say the DB is not a dependency. Snapshot bytes and all
# source/map/domain certificates remain bound. No current DB is exempted.

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
    top = sorted(ROOT / f'results/{stem}.json' for stem in STEMS)
    assert all(p.is_file() for p in top)
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
        if 'database_sha256' in doc and name not in HISTORICAL_DATABASE_METADATA:
            assert isinstance(doc.get('database'), str), name
            bound.append((doc['database'], doc['database_sha256']))
        # The full augmented LP stores its immutable preparation manifest
        # under prepared.sources. This is a named manifest, not solver data.
        manifest_docs = [doc]
        if name == 'results/fresh_sqrt5_labelled_packing_full_model_20261017.json':
            assert isinstance(doc.get('prepared', {}).get('sources'), dict)
            manifest_docs.append(doc['prepared'])
        for field, mapping in ((field, md.get(field, {}))
                               for md in manifest_docs for field in ('sources', 'source_hashes')):
            for source, value in mapping.items():
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
    reports = [ROOT / CHECKPOINT,
        ROOT / 'results/fresh_vertex_weighted_prism_base_capacity_20261017.md',
        ROOT / 'results/fresh_sqrt5_integer_simplex_port_packing_20261017_review.md',
        ROOT / 'results/fresh_nonarithmetic_b_release_scouts_20261017.md'] + sorted(
        ROOT / f'results/{stem}.md' for stem in STEMS if (ROOT / f'results/{stem}.md').is_file())
    links = []
    for report in reports:
        for target in re.findall(r'\]\(([^)]+)\)', report.read_text()):
            if '://' in target or target.startswith('#'):
                continue
            resolved = (report.parent / target.split('#')[0]).resolve()
            assert resolved.is_file() or resolved == OUT, (report, target)
            links.append(dict(report=str(report.relative_to(ROOT)), target=target))
    small_certificate_checks = {}
    lp_path = ROOT / 'results/fresh_sqrt5_labelled_packing_full_model_20261017.json'
    if lp_path in top:
        lp = json.loads(lp_path.read_text())
        assert lp['status'] == 'PASS_EXACT_PARTIAL_MODEL_CERTIFICATES'
        assert lp['solver_invocations'] == 1
        assert lp['prepared']['dimension'] == 13 and lp['prepared']['target'] == 12
        primal, dual = lp['certificates']
        assert primal['kind'] == 'EXACT_PARTIAL_MODEL_PRIMAL'
        assert dual['kind'] == 'EXACT_PARTIAL_MODEL_DUAL'
        alpha, mu = dual['face_multipliers'], dual['vertex_multipliers']
        assert len(alpha) == 5638 and len(mu) == 758
        for rows, limit in ((alpha, lp['prepared']['counts']['faces']),
                            (mu, lp['prepared']['counts']['vertex_forests'])):
            assert len({i for i, _ in rows}) == len(rows)
            assert all(0 <= i < limit and Fraction(v) >= 0 for i, v in rows)
        assert Fraction(dual['normalizer']) > 0
        assert Fraction(dual['bound']) == Fraction(dual['objective_numerator']) / Fraction(dual['normalizer'])
        assert dual['excludes_target12'] is True and Fraction(dual['bound']) > 12
        assert primal['meets_target12'] is False and Fraction(primal['maximum']) > 12
        assert Fraction(primal['maximum']) == max(Fraction(v) for _, v in primal['vertex_budgets'])
        assert Fraction(primal['repair_delta']) >= 0 and primal['integer_denominator'] > 0
        assert all(n >= 0 for _, n in primal['rounded_variable_numerators'])
        assert Fraction(dual['bound']) <= Fraction(primal['maximum'])
        small_certificate_checks = dict(
            scope='Independent small sign, indexing, ratio and target checks; full matrix arithmetic remains the separately replayed producer certificate.',
            nonnegative_face_multipliers=len(alpha), nonnegative_vertex_multipliers=len(mu),
            exact_dual_bound=dual['bound'], exact_primal_maximum=primal['maximum'],
            target12_excluded_within_partial_model=True)
        support_path = ROOT / 'results/fresh_sqrt5_labelled_packing_dual_support_20261017.json'
        assert support_path in top
        support = json.loads(support_path.read_text())
        assert support['status'] == 'PASS_EXACT_ENLARGED_MODEL_DUAL_SUPPORT'
        restriction = support['exported_only_repaired_dual']
        assert not support['exported_repair_added_simplex_rows']
        row_flags = {r['face']: r['has_exported_mixed_feature'] for r in support['positive_face_support']}
        assert len(restriction['face_multipliers']) == 4037
        assert len(restriction['vertex_multipliers']) == 687
        assert all(row_flags[i] and Fraction(v) >= 0 for i, v in restriction['face_multipliers'])
        assert all(Fraction(v) >= 0 for _, v in restriction['vertex_multipliers'])
        assert Fraction(restriction['normalizer']) > 0
        assert Fraction(restriction['bound']) == Fraction(restriction['objective_numerator']) / Fraction(restriction['normalizer']) > 12
        small_certificate_checks['exported_only'] = dict(
            faces=4037, forests=687, exact_bound=restriction['bound'],
            no_outside_repair_rows=True,
            scope='Saved-support membership and small exact sign/ratio checks, not a duplicate full load recount.')
    return dict(
        status='PASS_RECURSIVE_SOURCES_FROZEN_PREDECESSOR_AND_LOCAL_LINKS',
        scope='Archive integrity only. Exact results, independent reviews, bounded numerical diagnostics and unresolved work have separate scopes; this audit alone proves no mathematical claim.',
        checkpoint=CHECKPOINT, previous_archive=PREVIOUS,
        excluded_non_authoritative_producer_checkpoints=list(EXCLUDED_DIAGNOSTICS),
        historical_nonbinding_database_metadata=sorted(HISTORICAL_DATABASE_METADATA),
        previous_archive_sha256=digest(ROOT / PREVIOUS),
        previous_frozen_sections_preserved=preserved,
        previous_entrypoints='Historical snapshots; updated entrypoints are bound separately below.',
        report_count=len(reports), local_link_count=len(links), source_manifest_count=len(manifests),
        source_reference_count=references, nested_record_reference_count=records_count,
        unique_source_count=len(sources), manifests=sorted(manifests, key=lambda row: row['path']),
        inherited_bibliographic_entries=bibliography,
        recovered_historical_source_aliases=EXTERNAL_SOURCE_ALIASES,
        recovered_provenance_sha256=digest(ROOT / RECOVERED_PROVENANCE),
        small_certificate_checks=small_certificate_checks,
        source_hashes=dict(sorted(sources.items())),
        top_artifact_hashes={str(p.relative_to(ROOT)): digest(p) for p in top},
        report_hashes={str(p.relative_to(ROOT)): digest(p) for p in reports}, local_links=links,
        entrypoint_hashes={p: digest(ROOT / p) for p in ('README.md', 'RESEARCH_PLAN.md')},
        entrypoint_scope='Final entrypoints bound separately from the preserved historical archive.',
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
