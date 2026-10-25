#!/usr/bin/env python3
"""SOURCE ONLY: exact newer whole-provider links, with universal names separated.

Eight collections retain their original source identifiers. Explicit current
high-row IDs generate full finite completion keys. Two named universal models
and three old ideal17 updates keep their own exact receiver interfaces rather
than receiving fabricated high-row IDs.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from hashlib import sha256
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
HELPER = 'scripts/index_finite57_funded_selected_whole_providers_20261025.py'
HELPER_SHA = '0a245deec21b283f0d47af0323cf9557fd6ba4465f6e25a2b7c47897f2a46e86'
MAX_OUTPUT = 500 * 1024
VARIABLE = 'fresh_rank13_variable_context_retention_20261018'
A14 = 'fresh_a14_endpoint_a1_context_half_transfer_20261018'
COUPLED = 'fresh_d17_leaf_a16_context_coupled_transfer_20261018'


def need(ok, message):
    if not ok:
        raise ValueError(message)


def raw(relative):
    path = (ROOT / relative).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded newer-provider source path')
    with path.open('rb') as stream:
        data = stream.read(2 * 1024**2 + 1)
    need(0 < len(data) <= 2 * 1024**2, 'bounded newer-provider source read')
    return data


def prepare():
    data = raw(HELPER)
    need(sha256(data).hexdigest() == HELPER_SHA, 'newer whole helper identity')
    helper = types.ModuleType('finite57_newer_whole_helper')
    helper.__file__ = str(ROOT / HELPER)
    exec(compile(data, str(ROOT / HELPER), 'exec'), helper.__dict__)
    r, c = helper.prepare()
    return helper, r, c


def source_rows(c):
    for collection, (stem, path, single) in enumerate(c.NEWER):
        found = c.field(stem, path)
        for ordinal, proof in enumerate([found] if single else found):
            identifier = proof.get('index', proof.get('case', proof.get('name', 'A14')))
            yield collection, stem, path, ordinal, identifier, proof


def current_receipt(c, stem, identifier, proof, location):
    gram, white, gp, wp = c.old_context(proof, path=location)
    need(gram is not None and white is not None, 'newer exact physical context')
    receipt = c.receipt(stem, identifier, gram, white, proof['corners'], dict(gram=gp, white=wp))
    need(receipt['direct_current_receiver_condition'] and not receipt['legacy_pointwise_bounds_used'],
         'newer provider uses DIRECT current receiver only')
    return gram, white, receipt


def run(shard, helper, r, c):
    high, cores, partition = r.domain(c)
    reference = c.field(c.B13, ('current_interfaces', 'newer_corners'))
    wanted = [(p['source'], c.canonical(p['index']), tuple(p['marks'])) for p in reference]
    domain = []
    all_rows = list(source_rows(c))
    for collection, stem, path, ordinal, identifier, proof in all_rows:
        domain.extend((stem, c.canonical(identifier), tuple(q['marks'])) for q in proof['corners'])
    need(len(domain) == len(set(domain)) == len(wanted) and set(domain) == set(wanted),
         'ENTIRE eight-collection source/index/ordered-mark domain')
    original17 = c.field(c.OLD17, ('order17', 'replacement_cases'))
    old_by_case = {p['case']: p for p in original17}
    need(len(old_by_case) == len(original17) == 9, 'unique nine original ideal17 receivers')
    results = []
    for collection, stem, path, ordinal, identifier, proof in all_rows:
        if collection % 4 != shard:
            continue
        single = c.NEWER[collection][2]
        location = '/'.join(path) + ('' if single else '/' + str(ordinal))
        gram, white, receipt = current_receipt(c, stem, identifier, proof, location)
        common = dict(source=stem, source_path=location, source_identifier=identifier,
                      actual_source_selector=list(path), collection_is_single_object=single,
                      direct_current_receipt=receipt)
        if stem == VARIABLE and identifier in ('order16', 'cycle14'):
            results.append(dict(**common, domain='UNIVERSAL_NAMED_CURRENT_OUTPUT',
                receiver_key=dict(source=stem, name=identifier),
                high_row_id_assigned=False,
                required_geometry='Accepted complete named context/full-core seal and the earlier universal theorem using this exact model; no fabricated selected-row correspondence.'))
            continue
        if stem == COUPLED and path == ('local_proof', 'old_ideal17_seals'):
            need(identifier in (662, 675, 828), 'exact three coupled old ideal17 source cases')
            old = old_by_case[identifier]
            n = len(gram)
            old_white = sorted(set(range(n)) - set(old['black']))
            old_gram = cores['(25, 25, ((9, 9),))']['gram']
            need(c.E.matrix(gram) == c.E.matrix(old_gram) and white == old_white,
                 'coupled update exact old ideal17 physical receiver')
            need({tuple(q['marks']) for q in proof['corners']} ==
                 {tuple(q['marks']) for q in old['corners']}, 'coupled old receiver complete ordered-corner domain')
            results.append(dict(**common, domain='ORIGINAL_IDEAL17_CURRENT_UPDATE',
                original_receiver_key=dict(source=c.OLD17, case=identifier),
                high_row_id_assigned=False,
                required_geometry='Exact old-d17 positive-kernel/radical seal and its old/new correspondence provider; the source case number is not an authoritative new high index.'))
            continue
        if stem == VARIABLE:
            mapping = {'D15_1147': 1147, 'D15_1269': 1269}
            need(identifier in mapping, 'exact two selected D15 named seals')
            index = mapping[identifier]
        elif stem == A14:
            need(proof['affected_contexts'] == [1639] and identifier == 'A14',
                 'exact sole A14 selected context; original symbolic ID preserved')
            index = 1639
        else:
            need(type(proof.get('index')) is int and identifier == proof['index'],
                 'explicit current high index; no automatic case/name conversion')
            index = proof['index']
        need(0 <= index < len(high), 'newer authoritative high index range')
        case = high[index]
        core = cores[case['core_key']]
        provider = helper.record(c, r, stem, '/'.join(path), ordinal, proof, case, core)
        need(provider is not None and provider['full_completion_keys'],
             'selected newer model has exact current-funded whole finite completion keys')
        need(all(q['key']['target_order'] <= case['component_order_bound']
                 and len(case['marked_component']) <= q['key']['target_order']
                 for q in provider['full_completion_keys']), 'whole newer provider respects actual target eligibility')
        results.append(dict(**common, domain='CURRENT_SELECTED_WHOLE_FINITE',
                            authoritative_high_index=index, whole_provider=provider,
                            model_locator_convention='Helper source_row0 for a singleton is its explicit one-element wrapper, not a literal /0 JSON key.' if single else 'Actual source list ordinal.'))
    return dict(shard=shard, shards=4, entire_source_identifier_corner_domain=domain,
                exact_newer_provider_rows=results,
                scope='Eight newer collections only. Universal names and old ideal17 IDs remain distinct from authoritative high IDs; final G6 union and their source theorem applicability remain separate.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', choices=tuple(f'newer-{i}' for i in range(4)), required=True)
    parser.add_argument('--parent-lock', action='store_true', required=True)
    args = parser.parse_args()
    relative = str(Path(__file__).resolve().relative_to(ROOT))
    before = raw(relative)
    helper, r, c = prepare()
    out = run(int(args.batch[-1]), helper, r, c)
    for path in c.READ:
        need(c.digest(path) == c.PINS[path], 'newer provider consumed input/source changed')
    need(raw(relative) == before and sha256(raw(HELPER)).hexdigest() == HELPER_SHA, 'newer/helper source changed')
    need(sha256(helper.raw(helper.ROUTES)).hexdigest() == helper.ROUTES_SHA, 'native route definitions changed')
    need(sha256(helper.raw(helper.COLLECTOR)).hexdigest() == helper.COLLECTOR_SHA, 'direct collector changed')
    need(sha256(c.raw(c.MANIFEST)).hexdigest() == c.MANIFEST_SHA, 'base110 manifest changed')
    payload = json.dumps(dict(status='NEWER_CURRENT_WHOLE_PROVIDER_INDEX_ONLY', batch=args.batch,
        source_sha256=sha256(before).hexdigest(), current_terminal_bridge_complete=False,
        consumed_inputs={p: c.PINS[p] for p in sorted(c.READ)}, fuel_used=c.FUEL[0], result=out), sort_keys=True)
    need(len(payload.encode()) <= MAX_OUTPUT, 'newer provider output cap; no truncated domains')
    print(payload)


if __name__ == '__main__':
    main()
