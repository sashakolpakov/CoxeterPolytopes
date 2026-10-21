#!/usr/bin/env python3
"""Bounded, replay-only partition of the frozen ambient58 D13 checker.

Run ONE named batch through run_bounded_research_check_20261018.py
--parent-lock.  No certificate is written.  The conjunction of all BATCHES
replays every original output field and every original input status/pin.
The original source and certificate remain immutable.

JsonProjection is intentionally importable without importing the native
mathematical modules.  It projects exact object paths, never marker matches,
and closes its mmap before returning the bounded decoded projection.
"""
import argparse
import fcntl
import gc
from hashlib import sha256
import importlib
import json
import mmap
import os
from pathlib import Path
import re
import signal

ROOT = Path(__file__).resolve().parents[1]
OLD = 'fresh_d13_three_rank15_context_transfer_20261018'
SOURCE = f'scripts/{OLD}.py'
CERT = f'results/{OLD}.json'
WRAPPER = 'scripts/run_bounded_research_check_20261018.py'
FIXED = {
    SOURCE: 'a95ea679854ba907ead8640964a258c3baeb3e998052a3e2e3ff20f230cce486',
    CERT: '2c7d615cd38123bef20cffad68de8ab30baa91826c147a6bb084892590e8ca86',
    WRAPPER: '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}
MAX_CERT = 2 * 1024**2
MAX_FILE = 40 * 1024**2
MAX_FIELD = 24 * 1024**2
MAX_KEY = 4096
BATCHES = (
    'pins-0', 'pins-1', 'pins-2', 'pins-3', 'metadata',
    'support', 'rank15', 'one-sided', 'quiet', 'literal',
    'later-0', 'later-1', 'later-2', 'later-3',
    'newer', 'variants', 'd12', 'universal', 'finite',
)
TOP = {
    'status', 'sources', 'ambient_dimension', 'changed_pairs',
    'selected_support', 'selected_proofs', 'rank15_contact_classification',
    'one_sided_whole_models', 'quiet_repairs', 'current_interfaces',
    'universal_output_inventory', 'finite', 'scope',
}
INTERFACE = {
    'literal_corner_rows', 'literal_configurations', 'later_corners',
    'newer_corners', 'D14_variants', 'D12_positive_domain',
    'D12_D13_outputs', 'no_changed_same_component_role_in_these_domains',
    'universal_unstored_exception',
}
SPACE = re.compile(rb'\s*')
STRING = re.compile(rb'"(?:[^"\\]|\\.)*"')
STRUCTURAL = re.compile(rb'"(?:[^"\\]|\\.)*"|[{}\[\]]')
SCALAR = re.compile(rb'[^,\]}\s]+')


def digest(relative):
    """Streaming digest: never read a scientific input into one byte string."""
    path = (ROOT / relative).resolve()
    assert path.is_relative_to(ROOT) and path.is_file(), relative
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def bounded_read(relative, limit):
    with (ROOT / relative).open('rb') as stream:
        data = stream.read(limit + 1)
    assert len(data) <= limit, ('bounded read', relative, limit)
    return data


class JsonProjection:
    """Exact object-path projection with bounded allocation.

    Unselected composites are scanned structurally on the read-only mmap;
    their arrays/scalars are not decoded or copied.  Strings, including
    escaped quotes and brackets in strings, are one lexical token.  The
    visited objects reject duplicate keys.  The frozen file digest, checked
    before and after each batch, binds skipped JSON to the original input.
    Thus this is an input projection, not a new validator of every skipped
    JSON scalar's grammar.  All selected values use strict json.loads.
    """
    def __init__(self, relative):
        self.relative = relative
        self.path = ROOT / relative

    @staticmethod
    def end_value(data, start):
        assert start < len(data), 'missing JSON value'
        first = data[start]
        if first == 34:
            token = STRING.match(data, start)
            assert token is not None, 'unterminated JSON string'
            return token.end()
        if first in (123, 91):
            stack = [125 if first == 123 else 93]
            for token in STRUCTURAL.finditer(data, start + 1):
                at = token.start()
                char = data[at]
                if char == 34:
                    continue
                if char in (123, 91):
                    stack.append(125 if char == 123 else 93)
                else:
                    assert stack and stack.pop() == char, 'JSON nesting'
                    if not stack:
                        return token.end()
            raise AssertionError('unterminated JSON composite')
        token = SCALAR.match(data, start)
        assert token is not None, 'missing JSON scalar'
        return token.end()

    def load(self, paths, optional=()):
        paths = tuple(tuple(p) for p in paths)
        assert paths and len(paths) == len(set(paths)), 'duplicate selectors'
        assert all(all(isinstance(k, str) for k in p) for p in paths)
        optional = set(optional)
        assert all((key,) in paths for key in optional)
        trie = {}
        for path in paths:
            here = trie
            for key in path:
                assert None not in here, 'overlapping selectors'
                here = here.setdefault(key, {})
            assert not here, 'overlapping selectors'
            here[None] = True
        assert 0 < self.path.stat().st_size <= MAX_FILE, self.relative
        with self.path.open('rb') as stream:
            with mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
                def ws(i):
                    return SPACE.match(data, i).end()

                def take(i, wanted, prefix=()):
                    i = ws(i)
                    if None in wanted:
                        end = self.end_value(data, i)
                        assert end - i <= MAX_FIELD, ('selected field cap', self.relative)
                        return json.loads(data[i:end]), end
                    assert data[i] == 123, ('expected selected object', self.relative)
                    i = ws(i + 1)
                    out, seen = {}, set()
                    if data[i] == 125:
                        assert set(wanted) <= (optional if not prefix else set()), ('missing selectors', self.relative)
                        return out, i + 1
                    while True:
                        token = STRING.match(data, i)
                        assert token is not None and token.end() - i <= MAX_KEY
                        key = json.loads(data[i:token.end()])
                        assert key not in seen, ('duplicate key', self.relative, key)
                        seen.add(key)
                        i = ws(token.end())
                        assert data[i] == 58, 'expected colon'
                        i = ws(i + 1)
                        if key in wanted:
                            out[key], i = take(i, wanted[key], prefix + (key,))
                        else:
                            i = self.end_value(data, i)
                        i = ws(i)
                        if data[i] == 125:
                            i += 1
                            break
                        assert data[i] == 44, 'expected comma'
                        i = ws(i + 1)
                    assert set(wanted) - set(out) <= (optional if not prefix else set()), ('missing selectors', self.relative)
                    return out, i

                result, end = take(0, trie)
                assert ws(end) == len(data), 'trailing JSON data'
                return result


def document(stem, *paths):
    return JsonProjection(f'results/{stem}.json').load(paths)


def equal(native, expected, encode):
    # Exactly the frozen run's enc -> JSON normalization (tuples to lists).
    actual = json.loads(json.dumps(encode(native)))
    assert actual == expected, 'saved output projection differs'


def chunks(sequence, count=4):
    width = (len(sequence) + count - 1) // count
    groups = tuple(tuple(sequence[i*width:(i+1)*width]) for i in range(count))
    assert tuple(x for group in groups for x in group) == tuple(sequence)
    assert sum(map(len, groups)) == len(set(sequence))
    return groups


def pin_batch(old, saved, number):
    """Original run's all-STEMS status and inherited sources-map loop."""
    stems = chunks(old.STEMS)[number]
    cache = dict(saved['sources'])
    for stem in stems:
        doc = JsonProjection(f'results/{stem}.json').load(
            (('status',), ('sources',)), optional=('sources',))
        assert doc['status'].startswith('PASS'), stem
        assert isinstance(doc.get('sources', {}), dict)
        for relative, expected in doc.get('sources', {}).items():
            if relative not in cache:
                cache[relative] = digest(relative)
            assert cache[relative] == expected, (stem, relative)
        del doc
    # Any inherited pins outside the outer certificate are checked again.
    for relative, expected in cache.items():
        assert digest(relative) == expected, ('input changed', relative)


def literal_documents(old):
    """The complete selectors used by frozen T.literal_inputs, unchanged."""
    t = old.T
    plan = {}

    def add(stem, *paths):
        plan.setdefault(stem, []).extend(paths)

    add(t.A.OVER, ('isolation', 'cores', '(25, 25, ((9, 9),))'))
    add(t.A.OLD16, ('exceptional_D15',))
    add(t.A.OLD17, ('order17',))
    add(t.A.D13, ('extreme_equality', 'rows'))
    add(t.A.D14, ('equality', 'rows'))
    add(t.T.PARTIAL, ('cases',))
    add(t.D.TAG, ('cases',))
    add(t.T.NEW_D, ('two_A1_A3_spectators', 'rows'))
    return {stem: document(stem, *paths) for stem, paths in plan.items()}


def replay(old, saved, name, snapshot):
    s, b, t = old.S, old.B, old.T
    expect = saved['current_interfaces']
    if name.startswith('pins-'):
        pin_batch(old, saved, int(name[-1]))
    elif name == 'metadata':
        output = dict(
            status='PASS_AMBIENT58_D13_THREE_RANK15_CONTEXT_LOCAL_TRANSFER_AND_FINITE_CLOSURE',
            sources={p: snapshot[p] for p in sorted(saved['sources'])},
            ambient_dimension=58, changed_pairs=sorted(old.PAIRS),
            scope='Dimension-specific local rule and exact finite-count composition; no new global dimension theorem is claimed by this packet alone.')
        equal(output, {k: saved[k] for k in output}, old.enc)
        constants = dict(
            literal_configurations=216, D12_positive_domain=332,
            no_changed_same_component_role_in_these_domains=True,
            universal_unstored_exception='Affine-D12 quiet and arbitrary one-sided outputs, repaired separately')
        equal(constants, {k: expect[k] for k in constants}, old.enc)
        sources = set(b.FINITE)
        assert len(expect['later_corners']) == 564
        assert all(row['source'] in sources for row in expect['later_corners'])
    elif name == 'support':
        over = document(old.A.OVER, ('isolation', 'cases'), ('isolation', 'cores'))
        support, hits = old.support(over)
        equal(support, saved['selected_support'], old.enc)
        selected = [old.selected_proof(c, core, maps) for c, core, maps in hits]
        equal(selected, saved['selected_proofs'], old.enc)
    elif name == 'rank15':
        equal(old.rank15_spectators(), saved['rank15_contact_classification'], old.enc)
    elif name == 'one-sided':
        equal(old.one_sided_models(), saved['one_sided_whole_models'], old.enc)
    elif name == 'quiet':
        docs = {
            old.U.TAG: document(old.U.TAG, ('nonextreme_zero_cases',)),
            old.Z.O.TAG: document(old.Z.O.TAG, ('target',)),
        }
        equal(old.quiet_models(docs), saved['quiet_repairs'], old.enc)
    elif name == 'literal':
        docs = literal_documents(old)
        rows = []
        for stem, case, gram, white, corners in t.literal_inputs(docs):
            for marks, *_ in corners:
                images = old.long_images(gram, white, marks)
                assert not set(images) & old.PAIRS
                rows.append(dict(source=stem, case=case, marks=marks, images=images))
        equal(rows, expect['literal_corner_rows'], old.enc)
    elif name.startswith('later-'):
        stems = chunks(tuple(b.FINITE))[int(name[-1])]
        rows = []
        for stem in stems:
            keys = b.FINITE[stem]
            doc = document(stem, *((key,) for key in keys))
            for key in keys:
                for path, corner, gram, white, _ in b.corners(doc[key], key):
                    images = old.long_images(old.A.matrix(gram), white, corner['marks'])
                    assert not set(images) & old.PAIRS
                    rows.append(dict(source=stem, path=path, marks=corner['marks'], images=images))
            del doc
            gc.collect()
        equal(rows, [r for r in expect['later_corners'] if r['source'] in stems], old.enc)
    elif name == 'newer':
        # Same eight collections, same order as frozen current_interfaces.
        collections = (
            (s.R.TAG, ('local_proof', 'whole_proofs'), False),
            (s.H.TAG, ('local_proof', 'whole_proofs'), False),
            (s.S.TAG, ('local_proof', 'whole_proofs'), False),
            (s.N.TAG, ('local_proof', 'proofs'), False),
            (s.R.R.TAG, ('local_seals',), False),
            (s.H.Z.U.TAG, ('local_proof',), True),
            (s.TAG, ('local_proof', 'old_ideal17_seals'), False),
            (s.TAG, ('local_proof', 'A16_seals'), False),
        )
        rows = []
        for stem, path, single in collections:
            doc = document(stem, path)
            proofs = doc
            for key in path:
                proofs = proofs[key]
            if single:
                proofs = [proofs]
            for proof in proofs:
                physical = proof.get('old_complete_seal', proof)
                for corner in proof['corners']:
                    images = old.long_images(old.A.matrix(physical['gram']), physical['white'], corner['marks'])
                    assert not set(images) & old.PAIRS
                    rows.append(dict(source=stem, marks=corner['marks'], images=images))
            del doc, proofs
            gc.collect()
        equal(rows, expect['newer_corners'], old.enc)
    elif name == 'variants':
        rows = []
        doc = document(s.L.C.R.S.TAG, ('selected_proofs',))
        for proof in doc['selected_proofs']:
            for variant in proof.get('variants', []):
                images = [old.long_images(old.A.matrix(variant['gram']), variant['white'], corner['marks'])
                          for corner in variant['corners']]
                assert not any(set(im) & old.PAIRS for im in images)
                rows.append(dict(index=proof['index'], images=images))
        assert len(rows) == 24
        equal(rows, expect['D14_variants'], old.enc)
    elif name == 'd12':
        rows = []
        for stem, count in ((b.PART, 280), (b.REST, 52)):
            # These two complete documents are each below1.3MiB; walk's
            # traversal/path semantics are exactly those of the original.
            doc = document(stem, ())
            positive = [(path, p) for path, p in s.W.walk(doc) if old.F(p['payment']) > 0]
            assert len(positive) == count
            for path, proof in positive:
                if any(role['type'] == 'D13' for role in proof['roles']):
                    assert not set(map(tuple, proof['images'])) & old.PAIRS
                    rows.append(dict(source=stem, path=path, images=proof['images'], payment=proof['payment']))
            del doc, positive
        equal(rows, expect['D12_D13_outputs'], old.enc)
    elif name == 'universal':
        stem = 'fresh_d12_low_protected_adaptive_interfaces_20261018'
        docs = {stem: document(stem, ('low',), ('protected',), ('adaptive',))}
        equal(old.universal_roles(docs), saved['universal_output_inventory'], old.enc)
    elif name == 'finite':
        docs = {old.C.TAG: document(old.C.TAG, ('exact_actual_finite_maximum',),
                                  ('finite', 'exact_obstructing_forests'))}
        equal(old.finite(docs), saved['finite'], old.enc)
    else:
        raise AssertionError('unreachable batch')


def parent_lock():
    owner = int(os.environ.get('COXETER_CHECK_OWNER', '0'))
    assert owner > 1 and os.getppid() == owner, 'bounded parent wrapper required'
    with (ROOT / 'formal/.proof-batches/exclusive.lock').open('a') as probe:
        try:
            fcntl.flock(probe, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            pass
        else:
            fcntl.flock(probe, fcntl.LOCK_UN)
            raise AssertionError('--parent-lock required')
    return owner


def main():
    if not __debug__:
        raise RuntimeError('assertions must remain enabled')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', required=True)
    parser.add_argument('--batch', choices=BATCHES, required=True)
    args = parser.parse_args()
    owner = parent_lock()
    signal.alarm(55)  # The common wrapper retains its59s/1GiB limits.
    self_path = Path(__file__).relative_to(ROOT).as_posix()
    self_hash = digest(self_path)
    for path, expected in FIXED.items():
        assert digest(path) == expected, ('fixed pin', path)
    saved_bytes = bounded_read(CERT, MAX_CERT)
    assert sha256(saved_bytes).hexdigest() == FIXED[CERT]
    saved = json.loads(saved_bytes)
    assert set(saved) == TOP and set(saved['current_interfaces']) == INTERFACE
    assert len(saved['sources']) == 484
    before = {p: digest(p) for p in saved['sources']}
    assert before == saved['sources'], 'outer input pins'
    # Definition imports only: the frozen module's run() is never called.
    old = importlib.import_module(OLD)
    assert len(old.STEMS) == len(set(old.STEMS)) == 161
    paths = {SOURCE} | {
        f'{folder}/{stem}.{ext}' for stem in old.STEMS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))
    }
    assert paths == set(saved['sources']), 'original complete source domain'
    assert len(old.B.FINITE) == 22
    replay(old, saved, args.batch, before)
    assert os.getppid() == owner, 'parent exited'
    assert {p: digest(p) for p in before} == before, 'input pins changed'
    assert bounded_read(CERT, MAX_CERT) == saved_bytes, 'saved certificate changed'
    assert digest(self_path) == self_hash, 'worker changed'
    for path, expected in FIXED.items():
        assert digest(path) == expected, ('fixed pin changed', path)
    print('PASS_D13_FROZEN_TRANSFER_REPLAY_BATCH', args.batch,
          'outer_pins', len(before), 'saved_certificate', FIXED[CERT],
          'complete_replay_requires_all', len(BATCHES), flush=True)


if __name__ == '__main__':
    main()
