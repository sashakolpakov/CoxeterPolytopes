#!/usr/bin/env python3
"""Replay the two remaining eager finite57 producers by bounded exact inputs.

--verify --batch b13:support (inside the unchanged bounded wrapper)
--ledger results/finite57_loader_replay_new.json --plan
--ledger results/finite57_loader_replay_new.json --max-jobs 1

No generation mode. Complete producer equivalence requires every listed batch.
Only the coordinator writes its explicitly named, new-prefix runtime ledger.
"""
import argparse
import ast
import fcntl
import gc
from hashlib import sha256
import importlib
import json
import mmap
import os
from pathlib import Path
import signal
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]
NAMES = dict(b13='fresh_b13_endpoint_retention_20261018',
             frontier='fresh_d17_a16_coupled_ambient58_frontier_20261018')
PROJECTOR = 'check_d13_transfer_batches_20261020'
RUNNER = 'scripts/run_bounded_research_check_20261018.py'
DRIVER = 'scripts/verify_manuscript.py'
FIXED = {
    'scripts/'+NAMES['b13']+'.py': '5d0a782872fbb8270fea1706d95bfbf072d46c8f8a99cdae9a62330ac9cb9d53',
    'results/'+NAMES['b13']+'.json': '2bb463c567b393d7c95f211f8199c013bc0c34716a4dc98c98251fad7db2b2da',
    'scripts/'+NAMES['frontier']+'.py': '7b3f966491b8780e68b64ad6a92a805708c94a1915c14a6085a78c25e20f4d1c',
    'results/'+NAMES['frontier']+'.json': '0fb359b3b07d6fbabb31adbf611d769f41ebe855075076301975e086055b1223',
    'scripts/'+PROJECTOR+'.py': 'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    RUNNER: '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
    DRIVER: 'f6a0b39f981c41f3c7d75b184ee3d274e309616bcfa72b693a3b018bfc7394a1',
    'scripts/build_illustrated_ra_manuscript_20261018.py': '27f14501b52a38c92121266aceb58eb2e1d36bfc21c2615ab0b4e9be4f9dee15',
    'scripts/run_lean_batches_signal_safe_20261018.py': 'c0ddc98cd98e55f669252e4f138345665c2ce2b0e4adf12d8d4ce96aed65af4d',
    'scripts/run_lean_batches_20261018.py': '97822ecf4d9aef2bf04aa429b582c4cc36cf6d1b6f0364579251278f66fae72c',
}
CLOSURE_COUNT = 132
CLOSURE_SHA = '04aa42fe551f300dc6385e518766742a25c6a41c07310ace6a3b6be46902420e'
STEM_COUNT = 158
STEM_SHA = 'e393a2ff27ac8004c9dfbfc07c8f5a6522ef359218e2f1121199c3f9e8d59efb'
B13_PARTS = ('metadata', 'support', 'repairs', 'literal',
             'later-0', 'later-1', 'later-2', 'later-3',
             'newer', 'variants', 'd12', 'universal')
FRONT_PARTS = ('metadata', 'cover', 'increasing', 'complement')
BATCHES = ('parser',) + tuple(t+':'+b for t, parts in
    (('b13', B13_PARTS), ('frontier', FRONT_PARTS))
    for b in ('pins-0', 'pins-1', 'pins-2', 'pins-3')+parts)
B13_TOP = set('status sources baseline old_B13_cost new_B13_cost changed_pair weighted_pairs complete_selected_support whole_component_repairs current_interfaces universal_retained_roles ideal_weights_changed added_coupons triangle_minimum quadrilateral_minimum no_new_forest_count new_dimension_bound scope'.split())
B13_INTERFACE = set('literals literal_changes later_corners later_corner_count newer_corners D14_variants D12_positive_domains D12_B13_outputs'.split())
FRONT_TOP = set('status sources local_theorem scale representative_cover new_context_family finite exact_actual_finite_maximum positive_ideal_profiles actual_cusp_surplus finite_curvature strict_global_constant necessary_ambient58_inequality minimum_finite_vertices completed_dimension_bound new_dimension_bound scope'.split())
FINITE_KEYS = set('classes actual_fixed_fibers fresh_actual_residual_tables explicit_options Bellman_inequalities attaining_predecessors exact_maximum classes_above_curvature exact_obstructing_forests all_other_38_classes_at_most_curvature'.split())
PROJ = None


def digest(relative):
    p = (ROOT/relative).resolve()
    assert p.is_relative_to(ROOT) and p.is_file(), relative
    h = sha256()
    with p.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def source_pins():
    pending, pins = list(NAMES.values()), {}
    while pending:
        name = pending.pop()
        rel = 'scripts/'+name+'.py'
        if rel in pins:
            continue
        pins[rel] = digest(rel)
        tree = ast.parse((ROOT/rel).read_text(), filename=rel)
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                names = [a.name for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.module:
                assert not node.level
                names = [node.module]
            else:
                names = []
            for imported in names:
                if (ROOT/'scripts'/(imported+'.py')).is_file():
                    pending.append(imported)
    assert len(pins) == CLOSURE_COUNT
    assert sha256(json.dumps(pins, sort_keys=True, separators=(',', ':')).encode()).hexdigest() == CLOSURE_SHA
    return pins


def document(stem, *paths):
    return PROJ.document(stem, *paths)


def merge_pins(pins, extra):
    for p, h in extra.items():
        assert isinstance(p, str) and isinstance(h, str) and len(h) == 64
        assert p not in pins or pins[p] == h, ('inconsistent inherited pin', p)
        pins[p] = h


def ancestor_pins(projector, stem, outer):
    # Both frozen original run() source domains are exactly own script plus
    # script/JSON/report triples for STEMS. Derive that exact set without a
    # mathematical-module import; a frozen metadata list hash also binds it.
    stems = sorted(Path(p).stem for p in outer
                   if Path(p).parent == Path('results') and Path(p).suffix == '.json')
    assert len(stems) == STEM_COUNT and len(set(stems)) == STEM_COUNT
    assert sha256(json.dumps(stems, separators=(',', ':')).encode()).hexdigest() == STEM_SHA
    paths = {'scripts/'+stem+'.py'} | {f'{folder}/{s}.{ext}' for s in stems
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))}
    assert paths == set(outer), 'original complete source-domain schema'
    pins = dict(outer)
    for s in stems:
        doc = projector.document(s, ('status',), ('sources',))
        assert doc['status'].startswith('PASS'), s
        assert isinstance(doc['sources'], dict)
        merge_pins(pins, doc['sources'])
    return stems, pins


def equal(native, expected, encode):
    # The original enc -> JSON tuple normalization, one bounded object at a time.
    assert json.loads(json.dumps(encode(native))) == expected, 'saved output differs'


class ArrayView:
    """Read-only exact-path array adapter using the already audited lexer.

    No entire array is decoded. Visited objects reject duplicate keys;
    skipped values retain the original JsonProjection frozen-hash boundary.
    Each decoded array item has the SAME24MiB cap as JsonProjection.
    """
    def __init__(self, path):
        assert 0 < path.stat().st_size <= PROJ.MAX_FILE
        self.stream = path.open('rb')
        self.data = mmap.mmap(self.stream.fileno(), 0, access=mmap.ACCESS_READ)
        root_end = PROJ.JsonProjection.end_value(self.data, self.ws(0))
        assert self.ws(root_end) == len(self.data), 'trailing JSON'

    def close(self):
        self.data.close()
        self.stream.close()

    def ws(self, i):
        return PROJ.SPACE.match(self.data, i).end()

    def entries(self, start):
        d, seen = self.data, set()
        i = self.ws(start)
        assert d[i] == 123, 'expected object'
        i = self.ws(i+1)
        if d[i] == 125:
            return
        while True:
            token = PROJ.STRING.match(d, i)
            assert token and token.end()-i <= PROJ.MAX_KEY
            key = json.loads(d[i:token.end()])
            assert key not in seen, 'duplicate key'
            seen.add(key)
            i = self.ws(token.end())
            assert d[i] == 58
            lo = self.ws(i+1)
            hi = PROJ.JsonProjection.end_value(d, lo)
            yield key, lo, hi
            i = self.ws(hi)
            if d[i] == 125:
                return
            assert d[i] == 44
            i = self.ws(i+1)

    def span(self, path):
        lo, hi = self.ws(0), len(self.data)
        for key in path:
            fields = list(self.entries(lo))  # Keys/offsets only, never values.
            hits = [(a, b) for k, a, b in fields if k == key]
            assert len(hits) == 1, ('missing exact path', path)
            lo, hi = hits[0]
        return lo, hi

    def keys(self, path=()):
        return {k for k, _, _ in self.entries(self.span(path)[0])}

    def array(self, path):
        lo, hi = self.span(path)
        assert self.data[lo] == 91
        ranges, i = [], self.ws(lo+1)
        if self.data[i] == 93:
            assert i+1 == hi
            return LazyArray(self, ranges)
        while True:
            end = PROJ.JsonProjection.end_value(self.data, i)
            ranges.append((i, end))
            i = self.ws(end)
            if self.data[i] == 93:
                assert i+1 == hi
                return LazyArray(self, ranges)
            assert self.data[i] == 44
            i = self.ws(i+1)


class LazyArray:
    def __init__(self, view, ranges):
        self.view, self.ranges, self.cached = view, ranges, (-1, None)

    def __len__(self):
        return len(self.ranges)

    def __getitem__(self, index):
        assert isinstance(index, int) and 0 <= index < len(self.ranges)
        if self.cached[0] != index:
            lo, hi = self.ranges[index]
            assert hi-lo <= PROJ.MAX_FIELD, 'array item exceeds selected-value cap'
            self.cached = index, json.loads(self.view.data[lo:hi])
        return self.cached[1]


def parser_fixture():
    # Independent stdlib full decodes are confined to tiny owned fixtures.
    examples = [dict(a=dict(rows=[])), dict(a=dict(rows=[0, None, [1, 2], {'x': '[\\\"]'}])),
                dict(noise=[{'rows': 'wrong path'}], a=dict(rows=[{'k': 'é'}]))]
    with tempfile.TemporaryDirectory(prefix='finite57_array_fixture_') as folder:
        path = Path(folder)/'fixture.json'
        for obj in examples:
            path.write_text(json.dumps(obj))
            view = ArrayView(path)
            try:
                rows = view.array(('a', 'rows'))
                assert [rows[i] for i in range(len(rows))] == obj['a']['rows']
                assert view.keys() == set(obj)
            finally:
                view.close()
        for raw in ('{"a":0,"a":{"rows":[]}}', '{"a":{"rows":[1,]}}', '{"a":{"rows":[]}}x'):
            path.write_text(raw)
            view = None
            try:
                view = ArrayView(path)
                rows = view.array(('a', 'rows'))
                for i in range(len(rows)):
                    rows[i]
            except (AssertionError, ValueError, IndexError):
                pass
            else:
                raise AssertionError('invalid array fixture accepted')
            finally:
                if view is not None:
                    view.close()
    return dict(valid=3, rejected=3, mathematical_replay=False)


def b13_replay(n, saved, part):
    expected = saved['current_interfaces']
    if part == 'support':
        equal(n.support(document(n.A.OVER, ('isolation', 'cases'), ('isolation', 'cores'))), saved['complete_selected_support'], n.enc)
    elif part == 'repairs':
        docs = {n.O.TAG: document(n.O.TAG, ('conditional_long_orbit', 'cases')),
                n.COUPON: document(n.COUPON, ('additional_rules',))}
        equal(n.repairs(docs, document(n.A.OVER, ('isolation', 'cases'), ('isolation', 'cores'))), saved['whole_component_repairs'], n.enc)
    elif part == 'literal':
        docs = PROJ.literal_documents(n)
        old = document(n.S.TAG, ('current_interfaces', 'literals'))['current_interfaces']['literals']
        rows = []
        for prev, (st, case, g, w, cs) in zip(old, n.T.literal_inputs(docs), strict=True):
            assert (prev['source'], prev['case']) == (st, n.enc(case))
            corners = []
            for before, (marks, *_) in zip(prev['corners'], cs, strict=True):
                assert before['marks'] == marks
                r = n.corner(g, w, marks)
                assert not r['debit']
                corners.append(dict(**r, old_lower=before['lower'], new_lower=before['lower']))
            rows.append(dict(source=st, case=case, corners=corners, total=prev['total']))
        assert len(rows) == 216
        equal(rows, expected['literals'], n.enc)
    elif part.startswith('later-'):
        stems = PROJ.chunks(tuple(n.B.FINITE))[int(part[-1])]
        rows = []
        for st in stems:
            doc = document(st, *((key,) for key in n.B.FINITE[st]))
            for key in n.B.FINITE[st]:
                for path, p, g, w, _ in n.B.corners(doc[key], key):
                    r = n.corner(n.A.matrix(g), w, p['marks'])
                    assert not r['debit']
                    rows.append(dict(source=st, path=path, **r))
            del doc
            gc.collect()
        equal(rows, [r for r in expected['later_corners'] if r['source'] in stems], n.enc)
    elif part == 'newer':
        s = n.S
        collections = ((s.R.TAG, ('local_proof', 'whole_proofs'), False),
            (s.H.TAG, ('local_proof', 'whole_proofs'), False),
            (s.S.TAG, ('local_proof', 'whole_proofs'), False),
            (s.N.TAG, ('local_proof', 'proofs'), False),
            (s.R.R.TAG, ('local_seals',), False),
            (s.H.Z.U.TAG, ('local_proof',), True),
            (s.TAG, ('local_proof', 'old_ideal17_seals'), False),
            (s.TAG, ('local_proof', 'A16_seals'), False))
        rows = []
        for st, path, single in collections:
            doc = document(st, path)
            ps = doc
            for key in path:
                ps = ps[key]
            for p in [ps] if single else ps:
                physical = p.get('old_complete_seal', p)
                for c in p['corners']:
                    r = n.corner(n.A.matrix(physical['gram']), physical['white'], c['marks'])
                    assert not r['debit']
                    rows.append(dict(source=st, index=p.get('index', p.get('case', p.get('name', 'A14'))), **r))
            del doc, ps
            gc.collect()
        equal(rows, expected['newer_corners'], n.enc)
    elif part == 'variants':
        rows = []
        for p in document(n.S.L.C.R.S.TAG, ('selected_proofs',))['selected_proofs']:
            for v in p.get('variants', []):
                rr = [n.corner(n.A.matrix(v['gram']), v['white'], c['marks']) for c in v['corners']]
                assert not any(r['debit'] for r in rr)
                rows.append(dict(index=p['index'], extra_contact=v['extra_contact'], corners=rr))
        assert len(rows) == 24
        equal(rows, expected['D14_variants'], n.enc)
    elif part == 'd12':
        rows, domains = [], []
        for st, count in ((n.B.PART, 280), (n.B.REST, 52)):
            doc = document(st, ())
            ps = [(path, p) for path, p in n.W.walk(doc) if n.F(p['payment']) > 0]
            assert len(ps) == count
            for path, p in ps:
                if any(b['type'] == 'B13' for b in p['roles']):
                    assert list(n.PAIR) not in p['images']
                    rows.append(dict(source=st, path=path, marks=p['marks'], images=p['images'], payment=p['payment']))
            domains.append(dict(source=st, positive_chosen_outputs=count))
            del doc, ps
        equal(rows, expected['D12_B13_outputs'], n.enc)
        equal(domains, expected['D12_positive_domains'], n.enc)
    elif part == 'universal':
        docs = {n.FORK: document(n.FORK, ('cases',)),
                n.BC.TAG: document(n.BC.TAG, ('whole_finite_corners',)),
                n.B10.TAG: document(n.B10.TAG, ('finite_corners',))}
        equal(n.universal_roles(docs), saved['universal_retained_roles'], n.enc)
    else:
        raise AssertionError('unknown B13 part')


def metadata(n, saved, target, snapshot):
    if target == 'b13':
        paths = {f'scripts/{n.TAG}.py', f'scripts/{n.V.TAG}.py', f'results/{n.V.TAG}.md'}
        paths |= {f'{folder}/{s}.{ext}' for s in n.STEMS for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))}
        assert paths == set(saved['sources'])
        weights = [dict(pair=p, weight='0' if p == n.PAIR else '1') for p in n.combinations(range(13), 2)]
        assert n.F(n.S.actual_base(('B', 13)), 64) == 78 and n.F(n.actual_base(('B', 13)), 64) == 77
        assert sum(n.F(r['weight']) for r in weights) == 77
        out = dict(status='PASS_B13_ENDPOINT_EXISTING_COUPON_RETENTION',
            sources={p: snapshot[p] for p in sorted(paths)}, baseline=n.S.TAG,
            old_B13_cost='78', new_B13_cost='77', changed_pair=n.PAIR, weighted_pairs=weights,
            ideal_weights_changed=False, added_coupons=[], triangle_minimum='2', quadrilateral_minimum='1',
            no_new_forest_count=True, new_dimension_bound=None,
            scope='Whole B13 endpoint unit removed using two full old affine/radical seals and an existing distinct A1/D12 shortleaf unit; complete current chosen-output retention. No forest DP or new bound.')
        assert set(saved['current_interfaces']) == B13_INTERFACE
        assert saved['current_interfaces']['literal_changes'] == 0
        assert saved['current_interfaces']['later_corner_count'] == len(saved['current_interfaces']['later_corners']) == 564
        assert all(r['source'] in n.B.FINITE for r in saved['current_interfaces']['later_corners'])
    else:
        paths = {f'scripts/{n.TAG}.py'} | {f'{folder}/{s}.{ext}' for s in n.STEMS for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))}
        assert paths == set(saved['sources'])
        assert document(n.P.TAG, ('ideal_weights_changed',))['ideal_weights_changed'] is False
        d = document(n.S.TAG, ('positive_ideal_profiles',), ('actual_cusp_surplus',))
        assert d == dict(positive_ideal_profiles=2552, actual_cusp_surplus='36')
        out = dict(status='PASS_D17_A16_COUPLED_AMBIENT58_FRONTIER',
            sources={p: snapshot[p] for p in sorted(paths)}, local_theorem=n.P.TAG, scale=64,
            exact_actual_finite_maximum='378', positive_ideal_profiles=2552, actual_cusp_surplus='36',
            finite_curvature='377', strict_global_constant='4205', necessary_ambient58_inequality='v_f >= 4205 + 36*v_i',
            minimum_finite_vertices=4205, completed_dimension_bound=58, new_dimension_bound=None,
            scope='Complete actual frontier under the five sealed exchanges and A13 guard. Fresh actual residual certificates, two complete increasing-family bounds, no new global representative DP or ideal convolution. Eight B13/D13 plus three rank15 forests remain.')
    equal(out, {k: saved[k] for k in out}, n.enc)


def frontier_replay(n, view, part):
    if part == 'cover':
        docs = {n.C.TAG: document(n.C.TAG, ('finite', 'expanded_alternatives')),
                n.N.TAG: document(n.N.TAG, ('finite', 'classes'), ('increasing_family_upper', 'rows'))}
        expected = document(n.TAG, ('representative_cover',))['representative_cover']
        equal(n.representative_cover(docs), expected, n.enc)
    elif part == 'increasing':
        docs = {n.C.TAG: document(n.C.TAG, ('finite', 'tables')),
                n.N.TAG: document(n.N.TAG, ('finite', 'classes'))}
        equal(n.increasing_family(docs), document(n.TAG, ('new_context_family',))['new_context_family'], n.enc)
    elif part == 'complement':
        assert view.keys(('finite',)) == FINITE_KEYS
        prior = view.array(('finite', 'classes'))
        assert len(prior) == 39
        docs = {n.N.TAG: document(n.N.TAG, ('finite', 'classes'))}
        # Unchanged complete native complement, including every table/option.
        # Only the expected classes are lazy; its integer mathematics is untouched.
        actual = n.complement(docs, {'classes': prior})
        assert set(actual) == FINITE_KEYS and len(actual['classes']) == 39
        expected = document(n.TAG, *(('finite', k) for k in sorted(FINITE_KEYS-{'classes'})))['finite']
        equal({k: v for k, v in actual.items() if k != 'classes'}, expected, n.enc)
        for i, row in enumerate(actual['classes']):
            equal(row, prior[i], n.enc)
    else:
        raise AssertionError('unknown frontier part')


def run_batch(batch):
    global PROJ
    assert __debug__, 'assertions must remain enabled'
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 1 and os.getppid() == owner
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('r+') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise AssertionError('bounded --parent-lock required')
    signal.alarm(55)
    own = digest(Path(__file__).relative_to(ROOT))
    for p, h in FIXED.items():
        assert digest(p) == h, ('fixed input changed', p)
    PROJ = importlib.import_module(PROJECTOR)  # stdlib-only definitions.
    assert Path(PROJ.__file__).resolve() == ROOT/'scripts'/f'{PROJECTOR}.py'
    if batch == 'parser':
        details = parser_fixture()
    else:
        target, part = batch.split(':')
        native_pins = source_pins()
        cert = ROOT/'results'/f'{NAMES[target]}.json'
        view = ArrayView(cert)
        try:
            assert view.keys() == (B13_TOP if target == 'b13' else FRONT_TOP)
            excluded = ({'complete_selected_support', 'whole_component_repairs', 'universal_retained_roles'}
                        if target == 'b13' else {'representative_cover', 'new_context_family', 'finite'})
            saved = document(NAMES[target], *((k,) for k in sorted((B13_TOP if target == 'b13' else FRONT_TOP)-excluded)))
            pins = dict(FIXED)
            stems, inherited = ancestor_pins(PROJ, NAMES[target], saved['sources'])
            merge_pins(pins, native_pins)
            merge_pins(pins, inherited)
            assert {p: digest(p) for p in pins} == pins, 'saved input/source pins'
            n = importlib.import_module(NAMES[target])  # No original run() call.
            assert Path(n.__file__).resolve() == ROOT/'scripts'/f'{NAMES[target]}.py'
            assert len(n.STEMS) == len(set(n.STEMS))
            assert set(n.STEMS) == set(stems)
            if part.startswith('pins-'):
                for st in PROJ.chunks(n.STEMS)[int(part[-1])]:
                    doc = document(st, ('status',), ('sources',))
                    assert doc['status'].startswith('PASS')
                    assert isinstance(doc['sources'], dict)
                    for p, h in doc['sources'].items():
                        assert digest(p) == h, (st, p)
                details = dict(nested_pin_group=int(part[-1]))
            elif part == 'metadata':
                metadata(n, saved, target, pins)
                details = dict(metadata=True)
            elif target == 'b13':
                expected = document(n.TAG, *((k,) for k in ('complete_selected_support', 'whole_component_repairs', 'universal_retained_roles')))
                saved.update(expected)
                b13_replay(n, saved, part)
                details = dict(interface=part)
            else:
                frontier_replay(n, view, part)
                details = dict(finite_phase=part)
            assert {p: digest(p) for p in pins} == pins, 'input/source changed'
            assert source_pins() == native_pins
        finally:
            view.close()
    assert os.getppid() == owner and digest(Path(__file__).relative_to(ROOT)) == own
    for p, h in FIXED.items():
        assert digest(p) == h
    signal.alarm(0)
    print('PASS_FINITE57_LOADER_BATCH '+json.dumps(dict(batch=batch, worker_sha256=own,
          complete_producer_replay_requires_all=list(BATCHES), details=details), sort_keys=True), flush=True)


def coordinator(args):
    # Definition-only existing driver; it launches only the unchanged guarded runner.
    for p, h in FIXED.items():
        assert digest(p) == h
    current_pins = dict(FIXED)
    merge_pins(current_pins, source_pins())  # Static hashes, no mathematical imports.
    projector = importlib.import_module(PROJECTOR)
    for stem in NAMES.values():
        outer = projector.document(stem, ('sources',))['sources']
        _, inherited = ancestor_pins(projector, stem, outer)
        merge_pins(current_pins, inherited)
    driver = importlib.import_module('verify_manuscript')
    own = digest(Path(__file__).relative_to(ROOT))
    path = (ROOT/args.ledger).resolve()
    assert path.parent == ROOT/'results' and path.name.startswith('finite57_loader_replay_') and path.suffix == '.json'
    lock = path.with_suffix('.lock')
    with lock.open('a') as handle:
        fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        if path.exists():
            assert path.stat().st_size < 2*1024**2
            data = json.loads(path.read_text())
            assert data['worker_sha256'] == own and data['fixed_pins'] == FIXED
        else:
            data = dict(worker_sha256=own, fixed_pins=FIXED, jobs=[])

        def missing():
            # Cached pin-unit successes never substitute for current mandatory
            # nested source/status inputs: all union files are checked anew.
            assert {p: digest(p) for p in current_pins} == current_pins, 'current complete inherited pin union'
            good = set()
            for r in data['jobs']:
                if r['status'] != 'PASS_DECLARED_SCOPE' or r['exit_code'] != 0 or r['source_sha256'] != own:
                    continue
                raw = r.get('output', '')
                assert sha256(raw.encode()).hexdigest() == r['stdout_sha256']
                lines = [x for x in raw.splitlines() if x.startswith('PASS_FINITE57_LOADER_BATCH ')]
                assert len(lines) == 1
                answer = json.loads(lines[0].split(' ', 1)[1])
                assert answer['worker_sha256'] == own and answer['batch'] in BATCHES
                assert r['id'] == answer['batch']
                assert r['execution_runner_sha256'] == FIXED[RUNNER]
                assert answer['complete_producer_replay_requires_all'] == list(BATCHES)
                guard = json.loads(raw.splitlines()[-1])
                assert guard['status'] == 'PASS_BOUNDED_RESEARCH_CHECK'
                assert guard['checker_sha256'] == own and guard['runner_sha256'] == FIXED[RUNNER]
                assert guard['command'] == [str(Path(__file__).relative_to(ROOT)), '--verify', '--batch', answer['batch']]
                assert guard['wall_limit_seconds'] == 59 and guard['rss_limit_bytes'] == 1024**3
                assert guard['threads'] == 1 and guard['lock_owner'] == 'wrapper'
                assert 0 <= guard['elapsed_seconds'] <= 59  # The runner rounds milliseconds.
                assert 0 < guard['peak_process_tree_rss_bytes'] < 1024**3
                assert guard['initial_resources']['free_disk_gib'] >= 40
                assert guard['initial_resources']['available_memory_gib'] >= 6
                good.add(answer['batch'])
            return [b for b in BATCHES if b not in good]

        def report():
            todo = missing()
            print(json.dumps(dict(status='COMPLETE_RECORDED_COVERAGE' if not todo else 'INCOMPLETE',
                  plan_only=args.plan, missing=todo, ledger=args.ledger), sort_keys=True), flush=True)

        if args.plan:
            report()
            return
        assert 1 <= args.max_jobs <= len(BATCHES)
        for _ in range(args.max_jobs):
            todo = missing()
            if not todo:
                break
            assert digest(Path(__file__).relative_to(ROOT)) == own
            for p, h in FIXED.items():
                assert digest(p) == h
            batch = todo[0]
            job = dict(id=batch, kind='script', source=str(Path(__file__).relative_to(ROOT)),
                command=[sys.executable, '-B', str(ROOT/RUNNER), '--parent-lock', '--',
                         str(Path(__file__).relative_to(ROOT)), '--verify', '--batch', batch],
                coverage='Only '+batch+'; complete equivalence requires all named batches')
            result = driver.run_job(job)
            data['jobs'].append(result)  # Failures/deferred attempts remain in the ledger.
            data['status'] = 'COMPLETE_RECORDED_COVERAGE' if not missing() else 'INCOMPLETE'
            raw = (json.dumps(data, indent=2)+'\n').encode()
            assert len(raw) < 2*1024**2
            with tempfile.NamedTemporaryFile(dir=path.parent, prefix='.finite57_loader_', delete=False) as stream:
                tmp = Path(stream.name)
                stream.write(raw)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(tmp, path)
            if result['status'] != 'PASS_DECLARED_SCOPE':
                break
        report()


def main():
    assert __debug__
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true')
    parser.add_argument('--batch', choices=BATCHES)
    parser.add_argument('--ledger')
    parser.add_argument('--plan', action='store_true')
    parser.add_argument('--max-jobs', type=int, default=0)
    args = parser.parse_args()
    if args.batch:
        assert args.verify and not args.ledger and not args.plan and not args.max_jobs
        run_batch(args.batch)
    else:
        assert args.ledger and not args.verify and (args.plan or args.max_jobs)
        coordinator(args)


if __name__ == '__main__':
    main()
