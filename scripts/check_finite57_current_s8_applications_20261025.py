#!/usr/bin/env python3
"""Bind exact current deep/tree applications to the accepted conditional DAG.

This adapter does not prove its ordinary geometric or computed assumptions.
The frozen G4 controller, literal library and exact engine are its only loaded
repository modules. Legacy mathematical workers are read, never imported.
"""
if not __debug__:
    raise SystemExit('Optimization is forbidden')
import os
import sys
if sys.flags.optimize or 'PYTHONOPTIMIZE' in os.environ:
    raise SystemExit('Optimization environment is forbidden')

import argparse
from copy import deepcopy
import fcntl
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import types

ROOT = Path(__file__).resolve().parents[1]
G4 = 'scripts/check_finite57_g4_lower_cycle_terminal_bindings_20261025.py'
G4_SHA = 'f17014594587a68ef766ae245db372550fe89adc02bedf9e44c477a8d79eb242'
REPORT = 'results/review_finite57_current_S8_conditional_bridge_fresh_20261025.md'
REPORT_SHA = '2c97525c6005627fb78e627dbca31895271112fd02556f72a6de7034f9d23f48'
ACCEPTANCE = 'results/finite57_current_S8_conditional_bridge_root_acceptance_20261025.json'
ACCEPTANCE_SHA = '7cfc187fe9a25e488e0952d063ab7771702a2732baeafaea6077c8ef628cff4f'
LIB = 'scripts/finite57_literal_weights_audit_20261024.py'
LIB_SHA = 'f6eba7af77af50b821e9e99c589388c1352f9ba5a99052aaeee64393f5ada1e0'
TEX_SHA = 'c300db69c58211ffdf69c22970b289325224ad284d9e63a3c28d88821fca82a0'
WORKING_SHA = '0a01a79659c1c522f587315252e7dd1019fcca3281913cce089c4d4e154b7900'
INPUTS = 'results/finite57_current_s8_application_source_inputs_20261025.json'
INPUTS_SHA = 'eff4692873f4bd854e4dea9b6ecb1ede64cce253b77783d120e457af59397fb1'
PROFILE = 'fresh_distance8_ideal_corner_replacement_20261017'
EXTRA_PINS = {
    G4: G4_SHA,
    REPORT: REPORT_SHA,
    ACCEPTANCE: ACCEPTANCE_SHA,
    'manuscripts/revisions/finite_volume57_consolidated.tex': WORKING_SHA,
    'scripts/' + PROFILE + '.py': '858a000e6d3ed04fe19604c04b977c385cfe79dcc3001c4815bfd82bd73ab813',
    'results/' + PROFILE + '.md': '83b6c9eb437a0eb39ac3caf1fbd482fba313cae120824b311006100273d1417a',
    'results/' + PROFILE + '.json': 'dc51b741b002d6596277f102fc79997f7a4be3c479369eec973e1e5d680e4057',
}
PROFILES = ('isolated', 'outer_end3', 'outer_end4', 'outer_next3',
            'inner_end3', 'inner_next3', 'both_ends3', 'inner_end4')
DEEP = {
    6: 'fresh_depth6_finite6_zero_transfer_20261018',
    7: 'fresh_depth6_finite7_adaptive_transfer_20261018',
    8: 'fresh_depth6_finite8_zero_transfer_20261018',
}
KEY_FIELDS = {'high_index', 'core_key', 'base_nodes', 'selected_marks',
              'target_kind', 'target_order', 'affine', 'induced_map', 'marked_image'}
REQUEST_FIELDS = {'provider_id', 'ambient_dimension', 'factor_kind', 'factor_rank',
                  'canonical_marks', 'marked_distance', 'feature_depths',
                  'current_intrinsic_weight', 'source_family', 'status'}
MAX_OUTPUT = 500 * 1024
MAX_APPLICATIONS = 10000
RULE = dict(literal_source=LIB, literal_source_sha256=LIB_SHA,
            canonical_source_sha256=TEX_SHA, working_source_sha256=WORKING_SHA,
            weight_section_sha256='a895621180ab57f3d52d814e652283be469f9800c8af814fcfeeed8ef43b56e8')
PROOF = dict(path=REPORT, sha256=REPORT_SHA, acceptance=ACCEPTANCE,
             acceptance_sha256=ACCEPTANCE_SHA,
             interpretation='Source-bound accepted conditional mathematical argument; hashes identify bytes, not truth.')
DAG = ('S5', 'F6', 'S6', 'F7', 'S7', 'F8', 'S8')
BACK_EDGES = ('protected74 application provider', 'adaptive9 high-core provider',
              'large finiteD leaf9 half provider', 'affineD rank11/order12 opposite-leaf provider')


def need(ok, message):
    if not ok:
        raise ValueError(message)


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(',', ':'), ensure_ascii=True)


def raw(name, cap=4 * 1024**2):
    path = (ROOT / name).resolve()
    need(path.is_relative_to(ROOT) and path.is_file(), 'bounded repository source path')
    with path.open('rb') as stream:
        before = os.fstat(stream.fileno())
        data = stream.read(cap + 1)
        after = os.fstat(stream.fileno())
    identity = lambda st: (st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns, st.st_ctime_ns)
    need(identity(before) == identity(after) == identity(path.stat()), 'source changed during read')
    need(0 < len(data) <= cap, 'bounded source bytes')
    return data


def integers(values, bound, length=None):
    need(type(values) is list and (length is None or len(values) == length), 'integer-list schema')
    need(all(type(i) is int and 0 <= i < bound for i in values), 'integer-list domain')
    need(len(values) == len(set(values)), 'integer-list injection')
    return values


def selected_key(key):
    """Same nine-field finite domain as frozen union435e, preserving order."""
    need(type(key) is dict and set(key) == KEY_FIELDS, 'exact selected full-key fields')
    need(type(key['high_index']) is int and 0 <= key['high_index'] < 2614, 'selected high-index domain')
    need(type(key['core_key']) is str and 0 < len(key['core_key']) <= 256, 'selected physical-core key')
    base = integers(key['base_nodes'], 64)
    marks = integers(key['selected_marks'], 64, 2)
    need(10 <= len(base) <= 17 and set(marks) <= set(base), 'selected base/mark domain')
    n = key['target_order']
    need(type(n) is int and len(base) <= n <= 17 and type(key['affine']) is bool,
         'selected complete target order/affine type')
    need(key['target_kind'] in (('A', 'B', 'C', 'D') if key['affine'] else ('A', 'B', 'D')),
         'selected classical target family')
    mapping = integers(key['induced_map'], n, len(base))
    need(key['marked_image'] == sorted(mapping[base.index(m)] for m in marks),
         'selected full-map/mark association')
    return deepcopy(key)


def target_hypotheses(W, kind, order, affine, marks):
    """Recompute the entire typed hypothesis, without geometric inference."""
    need(type(affine) is bool and type(order) is int and 9 <= order <= 58,
         'ambient58 whole-target order/affine domain')
    need(kind in (('A', 'B', 'C', 'D') if affine else ('A', 'B', 'D')), 'classical target family')
    integers(marks, order, 2)
    need(marks == sorted(marks), 'canonical sorted marked pair')
    graph = W.classical_graph(kind, order, affine=affine)
    ds = [W._distances(graph, mark) for mark in marks]
    features = [i for i, row in enumerate(graph)
                if (len(row) >= 3 if affine else len(row) != 2)
                or any(label != 3 for label in row.values())]
    depths = [min(d[f] for f in features) for d in ds] if features else [None, None]
    value = W.affine_weight(kind, order, *marks) if affine else W.finite_same_weight((kind, order), *marks)
    # The incomplete-context D13 negative guard is conservative and explicit.
    if not affine and (kind, order) == ('D', 13) and tuple(marks) in ((0, 11), (0, 12)):
        value = min(value, F(1, 2))
    return dict(target=dict(affine=affine, kind=kind, order=order,
                            rank=order - int(affine), canonical_marks=list(marks)),
                marked_distance=ds[0][marks[1]], feature_depths=depths,
                feature_vertices=features, current_intrinsic_weight=str(value),
                positive_unmarked_or_context_credits_used=False)


def assumptions(stage):
    """Ordinary premises, not axioms inferred from matching file identities."""
    return [
        dict(id='actual_star', kind='geometric', status='ASSUMED',
             proposition='One actual quadrilateral W+a+b+c+d, opposite pairs a/c,b/d; norm2 roots, PD W and each actual edge, PSD corners, indefinite opposite pairs, negative index<=1; every original wall retained.'),
        dict(id='whole_factor', kind='geometric', status='ASSUMED',
             proposition='Whenever a complete canonical target is keyed, it is the WHOLE varying marked factor in that same actual face with the stated complete physical map; affine kernel/heights and arbitrary orthogonal finite spectators satisfy factor-fiber. A physical-path-only key asserts only its original actual roots and white-interior path, with no whole-target inference.'),
        dict(id='complete_classification', kind='external_mathematical', status='ASSUMED',
             proposition='Complete finite/affine and minimal-indefinite classifications and all-native support reductions, without arithmeticity or a finite57 dimension theorem.'),
        dict(id='complete_computed_premises', kind='source_bound_computed', status='ASSUMED',
             source=PROOF, stage=stage,
             proposition='The COMPLETE profile/inverse/Schur-polynomial/rectangle/spectrum/minimal-core propositions and all-white whole-factor seals used in accepted report Sections2–3; no sampled rank, hash, count or PASS string substitutes for these propositions.'),
    ]


def application(theorem, key, stage, predecessors):
    return dict(status='BOUND_CONDITIONAL_WHOLE_FACE_APPLICATION', theorem_id=theorem,
                application_key=key, rule_identity=RULE, proof_reference=PROOF,
                assumptions=assumptions(stage),
                required_predecessor_theorems=list(predecessors),
                conclusion=dict(quantity='total_current_charge_of_same_actual_face', lower='1'),
                intrinsic_corner_unit=False, theorem_execution=False)


def apply_deep(request, W, selected_full_key=None):
    """Bind one *complete* G4 request, retaining any nine-field selected key."""
    need(type(request) is dict and set(request) == REQUEST_FIELDS, 'complete original deep request fields')
    distance = request['marked_distance']
    need(type(distance) is int and distance in DEEP, 'exact F6/F7/F8 distance')
    need(request['status'] == 'UNBOUND_G6_CURRENT_OUTPUT_BRIDGE' and
         type(request['ambient_dimension']) is int and request['ambient_dimension'] == 58,
         'original unbound ambient58 request')
    theorem = 'finite57.current.deep_finite' + str(distance) + '.v1'
    need(request['provider_id'] == theorem and request['source_family'] == DEEP[distance],
         'provider/distance/source association')
    need(type(request['feature_depths']) is list and len(request['feature_depths']) == 2 and
         all(type(d) is int and d >= 6 for d in request['feature_depths']), 'exact two integer feature depths')
    current = target_hypotheses(W, request['factor_kind'], request['factor_rank'], False,
                                request['canonical_marks'])
    need(current['marked_distance'] == distance and
         current['feature_depths'] == request['feature_depths'] and
         all(type(d) is int and d >= 6 for d in current['feature_depths']) and
         request['current_intrinsic_weight'] == current['current_intrinsic_weight'] == '0',
         'exact distance, BOTH depths>=6 and current intrinsic zero')
    full = None if selected_full_key is None else selected_key(selected_full_key)
    if full is not None:
        need(not full['affine'] and full['target_kind'] == request['factor_kind'] and
             full['target_order'] == request['factor_rank'] and full['marked_image'] == request['canonical_marks'],
             'same complete selected target and deep request')
    key = dict(selected_full_key=full, original_deep_request=deepcopy(request), **current)
    predecessor = 'finite57.current.S' + str(distance - 1) + '.whole_face.v1'
    return application(theorem, key, 'F' + str(distance), (predecessor,))


def apply_s8(full_key, W):
    """Exact selected-map application; independent S8 execution remains needed."""
    full = selected_key(full_key)
    current = target_hypotheses(W, full['target_kind'], full['target_order'], full['affine'], full['marked_image'])
    need(current['marked_distance'] <= 8, 'S8 applies only to the actual marked white-interior path<=8')
    out = application('finite57.current.S8.whole_face.v1', dict(selected_full_key=full, **current),
                      'S8', ('finite57.current.deep_finite8.v1', 'finite57.current.S7.whole_face.v1'))
    out['required_actual_runtime_provider'] = 'finite57.current.S8.whole_face.v1'
    return out


def independent_s8_application(full_key, W):
    """Associate a provider's explicit short-path escape with its SAME face.

    The selected pair itself may be far apart. This function neither proves
    that an escape occurs nor silently substitutes that pair for its path.
    """
    full = selected_key(full_key)
    current = target_hypotheses(W, full['target_kind'], full['target_order'], full['affine'], full['marked_image'])
    out = application('finite57.current.S8.whole_face.v1', dict(selected_full_key=full, **current),
                      'S8-provider-escape', ('finite57.current.deep_finite8.v1', 'finite57.current.S7.whole_face.v1'))
    out['antecedent'] = dict(kind='accepted_provider_proof_alternative', status='ASSUMED_IF_THIS_BRANCH_OCCURS',
        proposition='The exact matched provider proof produces an actual adjacent white-interior path<=8 in this SAME original actual face. The selected marked pair need not be that path.')
    out['provider_composition_premise'] = dict(status='REQUIRES_EXACT_MATCHED_ACCEPTED_PROVIDER_IMPLICATION',
        same_full_key=deepcopy(full),
        logical_form='The accepted source proof exhaustively gives exclusion OR whole-current-face payment OR an actual shortpath<=8 in this SAME face; equivalently its exclusion/payment proof may assume NO such shortpath.',
        provenance='Caller must associate the complete provider receipt/key and source-bound accepted proof. A requiresS8 flag alone proves no disjunction and no actual escape.')
    out['required_actual_runtime_provider'] = 'finite57.current.S8.whole_face.v1'
    return out


def apply_physical_s8(high_row, core, engine, fuel):
    """Bind a high-prefix path in the supplied source-bound physical core."""
    need(type(high_row) is dict and type(core) is dict, 'actual high row/core dictionaries')
    index = high_row['index']; name = high_row['core_key']
    need(type(index) is int and 0 <= index < 2614 and type(name) is str and 0 < len(name) <= 256,
         'complete actual high-index/core-key association')
    nodes = high_row['marked_component']; marks = high_row['preferred_corner']
    gram = engine.matrix(core['gram']); n = len(gram)
    integers(nodes, n); integers(marks, n, 2)
    omitted = integers(high_row['omitted_black'], n, 2)
    need(not set(omitted) & set(nodes) and not set(omitted) & set(marks), 'actual selected/omitted black association')
    need(set(marks) <= set(nodes), 'actual preferred marks in the entire displayed marked component')
    corner_nodes = [i for i in range(n) if i not in omitted]
    parts = engine.components(gram, corner_nodes, fuel)
    need(any(set(part) == set(nodes) and set(marks) <= set(part) for part in parts),
         'same complete displayed marked component at the actual selected corner')
    whites = set(range(n)) - set(marks) - set(omitted)
    labels = engine.labels(gram, nodes)
    start, end = (nodes.index(mark) for mark in marks)
    parents = {start: None}; todo = [start]
    for v in todo:
        engine.charge(fuel)
        for w, label in enumerate(labels[v]):
            if label and w not in parents:
                parents[w] = v; todo.append(w)
    need(end in parents, 'connected actual preferred marks')
    local = [end]
    while local[-1] != start:
        local.append(parents[local[-1]])
    local_path = list(reversed(local)); path = [nodes[i] for i in local_path]
    need(1 <= len(path) - 1 <= 8 and set(path[1:-1]) <= whites,
         'actual physical white-interior marked path<=8')
    key = dict(high_index=index, core_key=name, ordered_base_nodes=list(nodes),
               selected_marks=list(marks), source_method=high_row['method'],
               omitted_black=list(omitted), original_physical_white_nodes=sorted(whites),
               physical_core_sha256=sha256(canonical(core).encode()).hexdigest(),
               physical_white_path=dict(vertices=path, length=len(path) - 1,
                                        native_labels=[labels[a][b] for a, b in zip(local_path, local_path[1:])]),
               actual_whole_target='UNSPECIFIED: the displayed component is not asserted whole')
    out = application('finite57.current.S8.whole_face.v1', key, 'S8-physical-prefix',
                      ('finite57.current.deep_finite8.v1', 'finite57.current.S7.whole_face.v1'))
    out['assumptions'].append(dict(id='exact_high_core_source_association', kind='source_bound_geometric', status='ASSUMED',
        proposition='Caller binds this complete high row to its exact frozen core_key and original actual-star occurrence. Its ordered marked_component/marks are actual roots, and all nonmarked path vertices are white in the selected adjacent corner. The core hash identifies bytes only.'))
    out['required_actual_runtime_provider'] = 'finite57.current.S8.whole_face.v1'
    return out


def tree_application(row):
    """A symbolic physical path in the actual spine; no target is invented."""
    need(type(row) is dict and set(row) == {'P', 'c', 'd', 'kind', 'length', 'finite_reason', 'ideal_reason'},
         'complete immediate tree source row')
    P = row['P']; c = row['c']; d = row['d']
    need(type(P) is int and P in (6, 7, 8) and c in PROFILES and d in PROFILES and
         'inner_next3' in (c, d) and row['kind'] == 'immediate_payment' and
         type(row['length']) is int and row['length'] == P - 1 and
         row['finite_reason'] == 'distance<=8' and
         row['ideal_reason'] == ('distance<=5' if P == 6 else 'persistent feature within1 of new mark'),
         'exact immediate source branch and current prior-stage path')
    # Choose one source-guaranteed side, even when both are inner_next3.
    side = 'c' if c == 'inner_next3' else 'd'
    reverse = side == 'd'
    white_positions = list(range(2, P)) if not reverse else list(range(P - 2, 0, -1))
    path = [side] + ['spine_' + str(i) for i in white_positions] + [('a' if reverse else 'b')]
    need(len(path) == len(set(path)) == P and len(white_positions) == P - 2,
         'complete physical white-interior output path')
    branch = dict(P=P, c=c, d=d)
    key = dict(source_branch=branch, complete_source_row=deepcopy(row), chosen_side=side,
               original_path=dict(marks=['a', 'b'], vertices=['a'] + ['spine_' + str(i) for i in range(1, P)] + ['b'],
                                  white_positions=list(range(1, P)), length=P),
               adjacent_corner_marks=[side, 'a' if reverse else 'b'],
               physical_white_path=dict(vertices=path, white_positions=white_positions, length=P - 1,
                                        attachment_position=P - 2 if reverse else 2),
               actual_whole_target='UNSPECIFIED: do not infer a complete factor from this path',
               persistent_feature=dict(vertex=path[1], distance_from_new_mark=1,
                                       original_retained_spine_neighbors=[P - 3, P - 1] if reverse else [1, 3]))
    prior = 'finite57.current.S' + str(P - 1) + '.whole_face.v1'
    predecessors = ('finite57.current.cycle10.whole_face.v1', 'finite57.current.cycle11.whole_face.v1') if P == 6 else (
        'finite57.current.S' + str(P - 2) + '.whole_face.v1', 'finite57.current.deep_finite' + str(P - 1) + '.v1')
    out = application(prior, key, 'tree-immediate-P' + str(P), predecessors)
    out['assumptions'].append(dict(id='actual_inner_next3_path', kind='geometric_source_implication', status='ASSUMED',
        proposition='Accepted whole-spine profile theorem: spine0=a,...,spineP=b are actual original roots with each inner vertex white; chosen opposite root meets the second white node simply. The two spine neighbors of this attachment remain in the indicated adjacent actual corner, so its degree>=3 persists. No induced path is asserted to be a whole factor.'))
    return out


def prepare():
    lock = ROOT / 'formal/.proof-batches/exclusive.lock'
    need(lock.is_file(), 'shared proof lock exists')
    with lock.open('r+') as stream:
        try:
            fcntl.flock(stream.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(stream.fileno(), fcntl.LOCK_UN)
            raise ValueError('shared parent proof lock must already be held')
    encoded = raw(G4); need(sha256(encoded).hexdigest() == G4_SHA, 'frozen G4f170 source')
    c = types.ModuleType('finite57_current_s8_frozen_g4'); c.__file__ = str(ROOT / G4)
    exec(compile(encoded, c.__file__, 'exec'), c.__dict__)
    manifest = c.raw(c.MANIFEST)
    need(sha256(manifest).hexdigest() == c.MANIFEST_SHA, 'unchanged old64 register')
    rows = c.decode(manifest)['inputs']
    need(type(rows) is list and len(rows) == 64, 'exact old64 input register')
    for r in rows:
        need(r['path'] not in c.PINS, 'duplicate frozen source input')
        c.PINS[r['path']] = r['sha256']; c.SIZES[r['path']] = r['bytes']
    c.PINS[c.MANIFEST] = c.MANIFEST_SHA; c.SIZES[c.MANIFEST] = len(manifest); c.READ[c.MANIFEST] = c.MANIFEST_SHA
    input_bytes = raw(INPUTS); need(sha256(input_bytes).hexdigest() == INPUTS_SHA, 'exact successor application register')
    inputs = c.decode(input_bytes)
    need(inputs['schema'] == 'finite57-current-s8-application-source-inputs-v1' and
         inputs['base_register'] == dict(path=c.MANIFEST, bytes=len(manifest), sha256=c.MANIFEST_SHA, entries=64),
         'unchanged old64 base register association')
    extras = inputs['additional_inputs']
    need(type(extras) is list and len(extras) == len(EXTRA_PINS) and
         {r['path']: r['sha256'] for r in extras} == EXTRA_PINS,
         'exact additional source input domain; no duplicates or replacement')
    c.PINS[INPUTS] = INPUTS_SHA; c.SIZES[INPUTS] = len(input_bytes); c.READ[INPUTS] = INPUTS_SHA
    for name, digest in EXTRA_PINS.items():
        data = raw(name); need(sha256(data).hexdigest() == digest, 'accepted application source: ' + name)
        need(len(data) == next(r['bytes'] for r in extras if r['path'] == name), 'additional source byte length')
        need(name not in c.PINS or c.PINS[name] == digest, 'extra/old64 source collision')
        c.PINS[name] = digest; c.SIZES[name] = len(data); c.READ[name] = digest
    acceptance = c.decode(c.captured(ACCEPTANCE))
    need(acceptance['report_sha256'] == REPORT_SHA and
         acceptance['accepted_conditional_dag'] == ['direct literal short4 and independent small ideal priorities', *DAG] and
         tuple(acceptance['absent_back_edges']) == BACK_EDGES and
         acceptance['executed_S8_provider'] is False and acceptance['current_terminal_bridge_complete'] is False,
         'precise accepted conditional DAG and runtime boundary')
    canonical_tex = c.captured('manuscripts/finite_volume57.tex')
    working_tex = c.captured('manuscripts/revisions/finite_volume57_consolidated.tex')
    canonical_weights = b''.join(canonical_tex.splitlines(keepends=True)[1250:1560])
    working_weights = b''.join(working_tex.splitlines(keepends=True)[1319:1629])
    need(canonical_weights == working_weights and sha256(canonical_weights).hexdigest() == RULE['weight_section_sha256'],
         'accepted actual canonical/working weight sections byte-identical')
    c.E = c.load('finite57_current_s8_engine', c.ENGINE, c.captured(c.ENGINE))
    c.W = c.load('finite57_current_s8_literal', c.LIB, c.captured(c.LIB))
    need(c.W.TEX_SHA256 == TEX_SHA and c.PINS[LIB] == LIB_SHA, 'literal current rule identity')
    return c, manifest


def run(batch, c):
    if batch.startswith('deep-'):
        shard = int(batch.rsplit('-', 1)[1]); need(0 <= shard < 16, 'deep rank-mod16 shard')
        for stem in DEEP.values():
            for folder, ext in (('scripts', 'py'), ('results', 'md'), ('results', 'json')):
                c.captured(f'{folder}/{stem}.{ext}')
        requests = []; ranks = [rank for rank in range(9, 59) if rank % 16 == shard]
        # This is exactly G4's distance<=8 deep-zero subdomain, repartitioned
        # by rank modulo16 for complete bounded output. No physical search.
        for rank, kind in product(ranks, ('A', 'B', 'D')):
            graph = c.W.classical_graph(kind, rank)
            distances = [c.W._distances(graph, i) for i in range(rank)]
            features = c.W._features(graph, finite=True)
            for a in range(rank):
                for b in range(a + 1, rank):
                    c.E.charge(c.FUEL)
                    distance = distances[a][b]
                    if distance not in DEEP or c.W.finite_same_weight((kind, rank), a, b) != 0:
                        continue
                    depths = [min(distances[m][f] for f in features) for m in (a, b)]
                    need(min(depths) >= 6, 'all current deep zeros have BOTH required depths')
                    requests.append(dict(provider_id='finite57.current.deep_finite' + str(distance) + '.v1',
                        ambient_dimension=58, factor_kind=kind, factor_rank=rank, canonical_marks=[a, b],
                        marked_distance=distance, feature_depths=depths, current_intrinsic_weight='0',
                        source_family=DEEP[distance], status='UNBOUND_G6_CURRENT_OUTPUT_BRIDGE'))
        need(0 < len(requests) <= MAX_APPLICATIONS, 'complete nonempty deep application domain')
        applications = [compact(apply_deep(p, c.W)) for p in requests]
        keys = [canonical(p['application_key']['original_deep_request']) for p in applications]
        need(len(keys) == len(set(keys)) == len(requests), 'each complete outgoing key bound exactly once')
        return dict(deep_applications=applications,
                    rank_partition=dict(modulus=16, residue=shard, ranks=ranks),
                    domain='Exactly EVERY original G4 deep-zero request at these ranks; all16 shards partition the old4 domain. LargeD half requests are outside this application domain.')
    need(batch == 'tree-immediate', 'named application batch')
    doc = c.packet('cycle18-tree')
    rows = c.domain(doc['tree']['rows'], lambda r: (r['P'], r['c'], r['d']),
                    list(product((6, 7, 8, 9), PROFILES, PROFILES)))
    expected = {(P, a, b) for P, a, b in product((6, 7, 8), PROFILES, PROFILES) if 'inner_next3' in (a, b)}
    immediate = [r for r in rows if r['kind'] == 'immediate_payment']
    c.domain(immediate, lambda r: (r['P'], r['c'], r['d']), expected)
    need(all((r['kind'] == 'immediate_payment') == ((r['P'], r['c'], r['d']) in expected) for r in rows),
         'complete source partition/branch association')
    need(all(r['kind'] == ('sealed_finite_A7' if (r['P'], r['c'], r['d']) == (6, 'inner_end3', 'inner_end3')
                          else 'positive_gap') for r in rows if r['kind'] != 'immediate_payment'),
         'exact remaining positive-gap and sole whole-A7 source branches')
    return dict(tree_applications=[compact(tree_application(r)) for r in immediate],
                source_partition_keys=[dict(P=r['P'], c=r['c'], d=r['d'], kind=r['kind']) for r in rows],
                positive_gap_inequalities_recomputed=False,
                domain='All45 exact source immediate branches, each with a chosen actual side and complete symbolic physical white path; no whole target is invented.')


def compact(app):
    """Keep every key; share only the identical conditional proof contract."""
    key = app['application_key']
    if 'original_deep_request' in key:
        # The original request already retains distance, BOTH depths, literal
        # zero, family/rank and marks. These were separately recomputed above;
        # sharing their agreed values avoids printing each value twice.
        key = dict(original_deep_request=key['original_deep_request'],
                   selected_full_key=key['selected_full_key'], verified_canonical_target=key['target'])
    out = dict(contract_id=app['theorem_id'], application_key=key)
    extra = [p for p in app['assumptions'] if p['id'] == 'actual_inner_next3_path']
    if extra:
        out['additional_assumptions'] = extra
    return out


def contracts():
    out = {}
    for k in (5, 6, 7):
        theorem = 'finite57.current.S' + str(k) + '.whole_face.v1'
        predecessors = ('finite57.current.cycle10.whole_face.v1', 'finite57.current.cycle11.whole_face.v1') if k == 5 else (
            'finite57.current.S' + str(k - 1) + '.whole_face.v1', 'finite57.current.deep_finite' + str(k) + '.v1')
        out[theorem] = application(theorem, {}, 'S' + str(k), predecessors)
    for k in DEEP:
        theorem = 'finite57.current.deep_finite' + str(k) + '.v1'
        out[theorem] = application(theorem, {}, 'F' + str(k), ('finite57.current.S' + str(k - 1) + '.whole_face.v1',))
        out[theorem]['recomputed_hypotheses'] = dict(
            distance='application_key.original_deep_request.marked_distance',
            both_feature_depths='application_key.original_deep_request.feature_depths',
            intrinsic_zero='application_key.original_deep_request.current_intrinsic_weight',
            target='application_key.verified_canonical_target',
            exact_source_association='application_key.original_deep_request.source_family')
    for value in out.values():
        del value['application_key']
    return out


def association_controls(c):
    """Small adversarial associations, not a second corridor enumeration."""
    def request(k):
        return dict(provider_id='finite57.current.deep_finite' + str(k) + '.v1', ambient_dimension=58,
                    factor_kind='A', factor_rank=k + 13, canonical_marks=[6, k + 6], marked_distance=k,
                    feature_depths=[6, 6], current_intrinsic_weight='0', source_family=DEEP[k],
                    status='UNBOUND_G6_CURRENT_OUTPUT_BRIDGE')
    positives = [apply_deep(request(k), c.W)['application_key'] for k in DEEP]
    bad = []
    q = request(6); q['source_family'] = DEEP[7]; bad.append(('wrong-source-for-distance', q))
    q = request(6); q['canonical_marks'] = [5, 11]; bad.append(('same-distance-near-feature-not-zero', q))
    q = request(6); q['feature_depths'] = [6, 7]; bad.append(('only-one-depth-corrupted', q))
    q = request(6); q['canonical_marks'] = [6, 13]; bad.append(('wrong-distance-and-target-mark', q))
    rejected = []
    for name, q in bad:
        try:
            apply_deep(q, c.W)
        except ValueError:
            rejected.append(name)
        else:
            raise ValueError('association control accepted: ' + name)
    row = dict(P=8, c='inner_next3', d='isolated', kind='immediate_payment', length=7,
               finite_reason='distance<=8', ideal_reason='persistent feature within1 of new mark')
    tree = tree_application(row)
    q = dict(row, c='isolated')
    try:
        tree_application(q)
    except ValueError:
        rejected.append('immediate-side-with-no-inner-next3')
    else:
        raise ValueError('tree source-side mutation accepted')
    full = dict(high_index=0, core_key='CONTROL_ONLY_NO_REALIZATION', base_nodes=list(range(10)),
                selected_marks=[0, 9], target_kind='A', target_order=10, affine=False,
                induced_map=list(range(10)), marked_image=[0, 9])
    escape = independent_s8_application(full, c.W)
    need(escape['application_key']['marked_distance'] == 9 and
         escape['antecedent']['status'] == 'ASSUMED_IF_THIS_BRANCH_OCCURS',
         'a selected far pair is kept separate from its assumed provider escape')
    q = deepcopy(full); q['marked_image'] = [0, 8]
    try:
        independent_s8_application(q, c.W)
    except ValueError:
        rejected.append('wrong-full-map-marked-image')
    else:
        raise ValueError('false full-map/image association accepted')
    return dict(association_control_positives=positives + [tree['application_key'], escape['application_key']],
                rejected_associations=rejected,
                scope_of_controls='Three exact deep source roles, one physical tree branch, and one explicitly hypothetical far-pair provider escape; six rejected false associations. Not a geometry test or domain replay.')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--batch', choices=tuple('deep-' + format(i, '02d') for i in range(16)) + ('tree-immediate', 'association-controls'), required=True)
    p.add_argument('--parent-lock', action='store_true', required=True)
    args = p.parse_args(); own = str(Path(__file__).resolve().relative_to(ROOT)); source = raw(own)
    c, manifest = prepare(); result = association_controls(c) if args.batch == 'association-controls' else run(args.batch, c)
    for name, digest in c.READ.items():
        need(sha256(raw(name)).hexdigest() == digest, 'consumed input changed')
    need(raw(G4) and sha256(raw(G4)).hexdigest() == G4_SHA and raw(c.MANIFEST) == manifest and raw(own) == source,
         'controller/register source bytes changed')
    status = 'PASS_CURRENT_S8_APPLICATION_ASSOCIATION_CONTROLS' if args.batch == 'association-controls' else 'PASS_CURRENT_S8_CONDITIONAL_APPLICATION_BINDINGS'
    out = dict(status=status, batch=args.batch,
               source_sha256=sha256(source).hexdigest(), consumed_inputs=c.READ,
               rule_identity=RULE, proof_reference=PROOF, application_contracts=contracts(), accepted_conditional_dag=list(DAG),
               absent_back_edges=list(BACK_EDGES), **result,
               conditional_applications_checked=True, assumptions_discharged=False,
               independent_S8_verified=False, executed_S8_provider=False,
               current_terminal_bridge_complete=False, global_finite57_verified=False,
               scope='Exact current typed/source application bindings at explicit ordinary hypotheses; neither intrinsic unit restoration nor a runtime S8 conjunction.')
    payload = canonical(out)
    need(len(payload.encode()) + 1 <= MAX_OUTPUT, 'complete application output cap; refuse, never truncate')
    print(payload)


if __name__ == '__main__':
    main()
