#!/usr/bin/env python3
"""Bounded exact replay of the immutable distance-nine/dimension108 packet.

No output certificate is written. Run under the bounded parent-lock wrapper.
Every interval in every phase is required; the contract explains why their
conjunction is the original complete run, including output order and counts.
No mathematical producer is run merely by importing this module.
"""
import argparse
from collections import Counter
import fcntl
import gc
from hashlib import sha256
import importlib
from functools import lru_cache
from itertools import combinations, product
import json
import mmap
from math import comb
import os
from pathlib import Path
import signal

from check_d13_transfer_batches_20261020 import JsonProjection

ROOT = Path(__file__).resolve().parents[1]
OLD = 'fresh_selective_distance9_dimension108_20261017'
SOURCE = f'scripts/{OLD}.py'
CERT = f'results/{OLD}.json'
FIXED = {
    SOURCE: 'f5cd7d0181cf9d5f96b10b6d4d6e7654caa3369b18ffa7473a16033f3124c81a',
    CERT: 'e7c6d1a1e56d8788e7d7a7d02bd5dd71837e3a3aef4ae768c1efc0a187bcacd9',
    'scripts/check_d13_transfer_batches_20261020.py':
        'be8a42d1d3687204af4d559c6bd8a28779dbfac68aeccdd95d41a7172645d555',
    'scripts/run_bounded_research_check_20261018.py':
        '4d8a2a74fa9fddcac796ca60f59505e633eb1f9912d4418a4e563debf4594132',
}
LIMITS = {
    'metadata': 1, 'domains': 1, 'overlap-domain': 384,
    'disjoint': 384, 'inertia': 85, 'overlaps': 384,
    'isolation': 949, 'profiles': 8,
}
MAX_ROW = 1024**2
MAX_FILE = 32 * 1024**2
TOP = {
    'status', 'sources', 'dimension_bound', 'inherited_dimension_bound',
    'unrestricted_finite_volume', 'arithmetic_hypothesis', 'cusp_restriction',
    'excluded_dimensions', 'weight_rule', 'finite_budget', 'ideal_budget',
    'census', 'disjoint_obstructions', 'low_core_lemma', 'isolation',
    'reduced_profiles', 'dimension_minima', 'boundary_at_n108',
    'compact_facet_ambient_bound', 'compact_facet_parity_slack_times8', 'scope',
}
CENSUS_KEYS = {
    'marked_contexts', 'ports', 'allowed_ports', 'overlap_generated',
    'consistent_overlaps', 'overlap_counts', 'overlap_records',
}
JOIN_KEYS = {
    'candidate_joins', 'white_rejections', 'remaining_joins',
    'remaining_context_pair_multiplicities', 'rows',
    'underlying_inertia_types', 'multiple_edge_obstructions', 'bond_forcing',
}
ISO_KEYS = {
    'cores', 'cases', 'component_order_cap', 'extra_cross_shortcuts',
    'extra_cross_pd_rejections', 'residual_tail_lemma',
}


def digest(relative):
    path = (ROOT / relative).resolve()
    assert path.is_relative_to(ROOT) and path.is_file(), relative
    h = sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def document(relative, *paths, optional=()):
    return JsonProjection(relative).load(paths, optional=optional)


@lru_cache(maxsize=128)
def saved(*path):
    out = document(CERT, path)
    for key in path:
        out = out[key]
    return out


def equal(actual, expected, label):
    # The frozen run normalizes tuples and integer object keys through JSON.
    assert json.loads(json.dumps(actual)) == expected, ('output differs', label)


class Arrays:
    """Read one array row at a time; no 21MiB certificate-wide json.loads.

    The object-path locator rejects duplicate keys on visited objects. Skipped
    values are bound by the before/after immutable certificate digest. Each
    selected array item is decoded with json.loads and a one-MiB item cap.
    """
    @staticmethod
    def ws(data, i):
        return _projection.SPACE.match(data, i).end()

    @classmethod
    def object(cls, data, start):
        i = cls.ws(data, start)
        assert data[i] == 123, 'expected object'
        i = cls.ws(data, i + 1)
        entries = {}
        if data[i] == 125:
            return entries
        while True:
            token = _projection.STRING.match(data, i)
            assert token and token.end() - i <= 4096
            key = json.loads(data[i:token.end()])
            assert key not in entries, 'duplicate object key'
            i = cls.ws(data, token.end())
            assert data[i] == 58
            i = cls.ws(data, i + 1)
            end = JsonProjection.end_value(data, i)
            entries[key] = (i, end)
            i = cls.ws(data, end)
            if data[i] == 125:
                return entries
            assert data[i] == 44
            i = cls.ws(data, i + 1)

    @classmethod
    def locate(cls, data, path):
        start = cls.ws(data, 0)
        end = JsonProjection.end_value(data, start)
        assert cls.ws(data, end) == len(data), 'trailing data after root'
        for key in path:
            entries = cls.object(data, start)
            assert key in entries, ('missing array path', path)
            start, end = entries[key]
        return start, end

    @classmethod
    def keys(cls, path=()):
        with (ROOT / CERT).open('rb') as stream:
            with mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
                start, _ = cls.locate(data, path)
                return set(cls.object(data, start))

    @classmethod
    def walk(cls, path, visit):
        assert 0 < (ROOT / CERT).stat().st_size <= MAX_FILE
        with (ROOT / CERT).open('rb') as stream:
            with mmap.mmap(stream.fileno(), 0, access=mmap.ACCESS_READ) as data:
                start, end = cls.locate(data, path)
                assert data[start] == 91
                i, count = cls.ws(data, start + 1), 0
                if data[i] == 93:
                    assert i + 1 == end
                    return 0
                while True:
                    stop = JsonProjection.end_value(data, i)
                    assert stop - i <= MAX_ROW, ('array item cap', path)
                    visit(count, json.loads(data[i:stop]))
                    count += 1
                    i = cls.ws(data, stop)
                    if data[i] == 93:
                        assert i + 1 == end
                        return count
                    assert data[i] == 44
                    i = cls.ws(data, i + 1)

    @classmethod
    def select(cls, path, predicate=lambda index, row: True):
        out = []
        total = cls.walk(path, lambda i, row: out.append(row) if predicate(i, row) else None)
        return out, total


# Import only the lightweight reusable scanner; its guarded main does not run.
_projection = importlib.import_module('check_d13_transfer_batches_20261020')


def context_key(context):
    return context['diagram'], tuple(context['black'])


def row_context(row, side='left'):
    return row[side][0], tuple(row[side][1])


def setup(old):
    # contexts_and_ports reads only this one field of its input dictionary.
    data = {old.CENSUS: document(f'results/{old.CENSUS}.json', ('diagrams',))}
    return old.contexts_and_ports(data)


def join_domain(old, grams, ports):
    """Original ordered port Cartesian product and white Schur classification.

    No obstruction is taken from saved rows. Only matrix-witness production is
    deferred to its own batches. Whole context groups retain all their ports.
    """
    allowed = [p for p in ports if p['allowed']]
    inverse = {}

    def diagonal(port):
        di, black = port['diagram'], tuple(port['black'])
        key = di, black
        if key not in inverse:
            ids = [i for i in range(grams[di].rows) if i not in black]
            a, scales = old.rationalize(grams[di].extract(ids, ids))
            old.D.pd(a)
            iv = a.inv()
            inverse[key] = {i: scales[k] * iv[k, k] for k, i in enumerate(ids)}
        return inverse[key][port['node']]

    rows, good, types = [], {}, set()
    for left, right in product(allowed, repeat=2):
        if left['distance'] + right['distance'] < 9:
            continue
        x, y = diagonal(left), diagonal(right)
        rows.append((left, right, x, y))
        if x * y < 1:
            di, dj = left['diagram'], right['diagram']
            key = di, tuple(left['black']), dj, tuple(right['black'])
            good.setdefault(key, []).append((left['node'], right['node']))
            types.add(tuple(sorted(((di, left['node']), (dj, right['node'])))))
    assert len(rows) == 4442 and sum(1 for _, _, x, y in rows if x*y >= 1) == 1103
    assert Counter(map(len, good.values())) == {1: 2007, 2: 348, 3: 180, 4: 24}
    assert len(types) == 85
    return rows, good, sorted(types)


def overlap_domain(old, grams, contexts, lo, hi, stats):
    """Exact frozen overlap loops, before inverse/PD-obstruction decisions.

    The only localization is the LEFT context range. Right contexts and every
    candidate matching size remain unchanged; matches spanning ports stay whole.
    """
    high = [r for r in contexts if max(r['distances']) >= 3]
    assert all(sum(d >= 3 for d in r['distances']) <= 5 for r in high)
    selected = {context_key(r) for r in contexts[lo:hi]}
    for left, right in product(high, repeat=2):
        if context_key(left) not in selected:
            continue
        g, h = grams[left['diagram']], grams[right['diagram']]
        candidates = [(i, j) for i in range(g.rows) for j in range(h.rows)
                      if left['distances'][i] + right['distances'][j] >= 10]
        bound = min(sum(d >= 3 for d in left['distances']),
                    sum(d >= 3 for d in right['distances']))
        for size in range(1, bound + 1):
            for match in combinations(candidates, size):
                stats['generated'] += 1
                ls, rs = [i for i, _ in match], [j for _, j in match]
                if len(set(ls)) < size or len(set(rs)) < size:
                    continue
                if any(g[i, k] != h[j, l] for i, j in match for k, l in match):
                    continue
                stats['consistent'] += 1
                yield left, right, match


def raw_overlap_count(contexts, lo=0, hi=384):
    """Count the unchanged combinations domain, without Gram pruning."""
    high = [r for r in contexts if max(r['distances']) >= 3]
    owners = {context_key(c) for c in contexts[lo:hi]}
    total = 0
    for left, right in product(high, repeat=2):
        if context_key(left) not in owners:
            continue
        count = sum(a+b >= 10 for a in left['distances'] for b in right['distances'])
        bound = min(sum(d >= 3 for d in left['distances']),
                    sum(d >= 3 for d in right['distances']))
        total += sum(comb(count, size) for size in range(1, min(count, bound)+1))
    return total


def survivor_rows():
    rows, total = Arrays.select(('census', 'overlap_records'),
                                lambda i, row: 'surviving_type' in row)
    assert total == 20718 and len(rows) == 949
    return rows


def input_pins(stems, checked):
    for stem in stems:
        doc = document(f'results/{stem}.json', ('status',), ('sources',), optional=('sources',))
        assert doc['status'].startswith('PASS'), stem
        for path, expected in doc.get('sources', {}).items():
            if path not in checked:
                assert digest(path) == expected, ('inherited source pin', stem, path)
            else:
                assert checked[path] == expected, ('inconsistent inherited pin', stem, path)
            checked[path] = expected


def metadata(old, contexts):
    expected = dict(
        status='PASS_SELECTIVE_DISTANCE9_CONTROLLED_EXTENSIONS_DIMENSION108',
        dimension_bound=108, inherited_dimension_bound=116,
        unrestricted_finite_volume=True, arithmetic_hypothesis=False,
        cusp_restriction=None, excluded_dimensions=list(range(109, 117)),
        weight_rule='1 if marked distance<=9 OR both marks lie in one component of order<=17;0 otherwise',
        finite_budget='9n', ideal_budget='sum min(9,r_i/2)', low_core_lemma=old.EXT,
        compact_facet_ambient_bound=79,
        scope='Exact finite classification, controlled relevant-component extensions and scalar inequalities with explicit inherited geometry and complete small-order lemma. No construction, external refereeing, formalization or priority claim.')
    equal(expected, {key: saved(key) for key in expected}, 'metadata constants')
    assert Arrays.keys() == TOP
    assert Arrays.keys(('census',)) == CENSUS_KEYS
    assert Arrays.keys(('disjoint_obstructions',)) == JOIN_KEYS
    assert Arrays.keys(('isolation',)) == ISO_KEYS
    census, relevant, shorts, extra, cap, corekeys = Counter(), 0, 0, 0, 0, set()
    survivors = survivor_rows()

    def case(index, row):
        nonlocal relevant, shorts, extra, cap
        source = survivors[index]
        assert (row['type'], row['left'], row['right'], row['match']) == (
            source['surviving_type'], source['left'], source['right'], source['match'])
        di, _ = source['left']; dj, _ = source['right']
        key = str((di, dj, tuple(map(tuple, source['match']))))
        assert row['core_key'] == key
        corekeys.add(key)
        relevant += row['method'] == 'relevant_component'
        shorts += row['extra_cross_bond_checks']['shortcut_rejections']
        extra += len(row['extra_cross_bond_checks']['edge_pd_z_rejections'])
        cap = max(cap, row['component_order_bound'])
    assert Arrays.walk(('isolation', 'cases'), case) == 949
    assert len(corekeys) == 54 and Arrays.keys(('isolation', 'cores')) == corekeys
    assert (relevant, shorts, extra, cap) == (69, 60831, 888, 17)
    equal(dict(component_order_cap=17, extra_cross_shortcuts=60831,
               extra_cross_pd_rejections=888,
               residual_tail_lemma='Every connected PD norm2 Coxeter Gram of rankm has inverse diagonal>=m/(m+1). Across all contacted components this gives totalrank<=floor(rho/(1-rho)) forrho<1. When an outside core-contacting root coexists, the exact two-residual Schur tests exclude or seal every tail leading from the preferred marked component.'),
          {key: saved('isolation', key) for key in ISO_KEYS - {'cores', 'cases'}},
          'isolation global scalars')
    context_order = {context_key(c): i for i, c in enumerate(contexts)}
    last_owner = -1
    def tally(index, row):
        nonlocal last_owner
        owner = row_context(row)
        assert owner in context_order and row_context(row, 'right') in context_order
        assert context_order[owner] >= last_owner, 'global left-context order'
        last_owner = context_order[owner]
        census['inverse_rejection' if 'inverse_rejection' in row else
               'edge_pd_rejection' if 'edge_pd_rejection' in row else row['surviving_type']] += 1
    assert Arrays.walk(('census', 'overlap_records'), tally) == 20718
    equal(dict(census), saved('census', 'overlap_counts'), 'all overlap counters')
    assert census['inverse_rejection'] == 18721 and census['edge_pd_rejection'] == 1048
    assert sum(v for k, v in census.items() if k not in ('inverse_rejection', 'edge_pd_rejection')) == 949
    assert raw_overlap_count(contexts) == saved('census', 'overlap_generated') == 312316
    assert saved('census', 'consistent_overlaps') == 20718
    margins, profile_keys = [], []
    def profile(index, row):
        margins.append(old.F(row['margin']))
        profile_keys.append((row['link_dimension'], row['ranks']))
    assert Arrays.walk(('reduced_profiles',), profile) == 756
    complete_profile_keys = [(d, ranks) for d in range(108, 116) for ranks in old.profiles(d)]
    assert len(complete_profile_keys) == 756
    equal(complete_profile_keys, json.loads(json.dumps(profile_keys)), 'complete profile domain/order')
    assert min(margins) == old.F(45873, 2476099)
    boundary = old.budget((18,)*5 + (17,))
    assert old.F(boundary['margin']) == old.F(-4659, 260642)
    equal(boundary, saved('boundary_at_n108'), 'boundary')
    t = old.sp.Symbol('t'); facet = []
    for start in (79, 80):
        d = start + 2*t
        beta = (d-1)*(d-5)/8 if start % 2 else d*(d-6)/8
        polynomial = old.sp.Poly(old.sp.expand(8*(beta-9*(d+1))), t)
        facet.append([int(polynomial.nth(k)) for k in range(3)])
    assert facet == [[12, 160, 4], [88, 164, 4]]
    equal(facet, saved('compact_facet_parity_slack_times8'), 'facet polynomial')


def domains(old, grams, contexts, ports):
    equal(ports, saved('census', 'ports'), 'complete port census')
    assert saved('census', 'marked_contexts') == len(contexts) == 384
    assert saved('census', 'allowed_ports') == 378
    joins, good, types = join_domain(old, grams, ports)
    expected, total = Arrays.select(('disjoint_obstructions', 'rows'))
    assert total == 4442
    actual_keys = [([a['diagram'], a['black'], a['node']],
                    [b['diagram'], b['black'], b['node']]) for a, b, _, _ in joins]
    equal(actual_keys, [[row['left'], row['right']] for row in expected], 'ordered join domain')
    fields = dict(candidate_joins=4442, white_rejections=1103, remaining_joins=3339,
                  remaining_context_pair_multiplicities=dict(Counter(map(len, good.values()))),
                  bond_forcing='Every label>=4 is excluded by the actual-edge normalized inverse test or the supporting-white Schur product. All nonempty subsets of up tofour remaining simple cross edges are exhausted.')
    equal(fields, {key: saved('disjoint_obstructions', key) for key in fields}, 'join counters')
    mult = []
    for (di, u, dj, v), edges in sorted(good.items()):
        for size in range(2, len(edges)+1):
            for chosen in combinations(edges, size):
                mult.append(([di, list(u)], [dj, list(v)], list(map(list, chosen))))
    expected, total = Arrays.select(('disjoint_obstructions', 'multiple_edge_obstructions'))
    assert total == len(mult) == 1332
    equal(mult, [[row['left'], row['right'], row['edges']] for row in expected], 'all multi-edge subset keys/order')
    expected, total = Arrays.select(('disjoint_obstructions', 'underlying_inertia_types'))
    assert total == len(types) == 85
    equal(types, [row['ports'] for row in expected], 'sorted recomputed inertia types')


def overlap_keys(old, grams, contexts, lo, hi):
    stats = Counter(); keys = []
    for left, right, match in overlap_domain(old, grams, contexts, lo, hi, stats):
        keys.append(([left['diagram'], left['black']], [right['diagram'], right['black']], list(map(list, match))))
    owners = {context_key(c) for c in contexts[lo:hi]}
    expected, total = Arrays.select(('census', 'overlap_records'), lambda i, row: row_context(row) in owners)
    assert stats['generated'] == raw_overlap_count(contexts, lo, hi)
    assert stats['consistent'] == len(expected) and total == 20718
    equal(keys, [[row['left'], row['right'], row['match']] for row in expected], 'complete ordered overlap domain')


def disjoint(old, grams, contexts, ports, lo, hi):
    joins, good, _ = join_domain(old, grams, ports)
    chosen_contexts = {context_key(c) for c in contexts[lo:hi]}
    rows = []
    for left, right, x, y in joins:
        if context_key(left) not in chosen_contexts:
            continue
        di, dj = left['diagram'], right['diagram']
        u, v = left['black'], right['black']; i, j = left['node'], right['node']
        core, mp = old.D.union(grams[di], grams[dj], ())
        core[i, mp[j]] = core[mp[j], i] = -1
        white = [k for k in range(core.rows) if k not in u + [mp[k] for k in v]]
        row = dict(left=[di, u, i], right=[dj, v, j], white_inverse_product=str(x*y))
        if x*y >= 1:
            a, scales = old.rationalize(core.extract(white, white))
            obstruction = old.z_obstruction(a); assert obstruction
            row['white_rejection'] = dict(nodes=white, gauge_squared=scales, **obstruction)
        else:
            maxima = [max(map(old.sp.Rational, p['inverse_diagonals'])) for p in (left, right)]
            assert max(maxima) >= 1 or 2*x*y >= 1
            row['larger_bond_obstruction'] = dict(endpoint_maxima=list(map(str, maxima)), twice_white_inverse_product=str(2*x*y))
            row['inertia_type'] = [list(p) for p in sorted(((di, i), (dj, j)))]
        rows.append(row)
    expected, _ = Arrays.select(('disjoint_obstructions', 'rows'),
                                lambda i, row: row_context(row) in chosen_contexts)
    equal(rows, expected, 'context-owned single joins')
    multiple = []
    for (di, u, dj, v), edges in sorted(good.items()):
        if (di, u) not in chosen_contexts:
            continue
        base, mp = old.D.union(grams[di], grams[dj], ())
        black = list(u) + [mp[i] for i in v]
        for size in range(2, len(edges)+1):
            for chosen in combinations(edges, size):
                g = base.copy()
                for x, y in chosen:
                    g[x, mp[y]] = g[mp[y], x] = -1
                cert = old.edge_obstruction(g, black); assert cert
                multiple.append(dict(left=[di, list(u)], right=[dj, list(v)],
                                     edges=list(map(list, chosen)), edge_pd_rejection=cert))
    expected, _ = Arrays.select(('disjoint_obstructions', 'multiple_edge_obstructions'),
                                lambda i, row: row_context(row) in chosen_contexts)
    equal(multiple, expected, 'context-owned all multi-edge subsets')


def inertia(old, grams, ports, lo, hi):
    _, _, types = join_domain(old, grams, ports)
    rows = []
    for typ in types[lo:hi]:
        (a, b), (c, d) = typ
        gram, mp = old.D.union(grams[a], grams[c], ())
        gram[b, mp[d]] = gram[mp[d], b] = -1
        h, scales = old.rationalize(gram)
        cert = old.exact_inertia(h); assert cert['signature'][1] == 2
        rows.append(dict(ports=[list(p) for p in typ], gram=old.matrix(h), gauge_squared=scales, **cert))
    expected, total = Arrays.select(('disjoint_obstructions', 'underlying_inertia_types'), lambda i, row: lo <= i < hi)
    assert total == 85
    equal(rows, expected, 'exact inertia types')


def overlaps(old, grams, contexts, values, lo, hi):
    rows = []; stats = Counter()
    for left, right, match in overlap_domain(old, grams, contexts, lo, hi, stats):
        g, h = grams[left['diagram']], grams[right['diagram']]
        ls, rs = [i for i, _ in match], [j for _, j in match]
        row = dict(left=[left['diagram'], left['black']], right=[right['diagram'], right['black']], match=list(map(list, match)))
        bad = None
        for i, j in match:
            for side, x, y, rr, other, shared, other_black in (
                    ('left', i, j, left, h, rs, right['black']),
                    ('right', j, i, right, g, ls, left['black'])):
                ext = [k for k in range(other.rows) if k != y and other[y, k] and k not in shared]
                if ext and max(values(rr['diagram'], tuple(rr['black']), x)) >= 2:
                    assert all(k not in other_black for k in ext)
                    bad = dict(side=side, node=x, external_neighbor=ext[0],
                               inverse_diagonals=list(map(str, values(rr['diagram'], tuple(rr['black']), x))))
                    break
            if bad:
                break
        if bad:
            row['inverse_rejection'] = bad
        else:
            core, mp = old.D.union(g, h, match)
            black = left['black'] + [mp[i] for i in right['black']]
            obstruction = old.edge_obstruction(core, black)
            if obstruction:
                row['edge_pd_rejection'] = obstruction
            else:
                row['surviving_type'] = old.overlap_type(left['diagram'], right['diagram'], match)
        rows.append(row)
    owners = {context_key(c) for c in contexts[lo:hi]}
    expected, _ = Arrays.select(('census', 'overlap_records'), lambda i, row: row_context(row) in owners)
    equal(rows, expected, 'left-context overlap classifications')
    return dict(stats)


def isolation(old, grams, lo, hi):
    survivors = survivor_rows()
    expected, total = Arrays.select(('isolation', 'cases'), lambda i, row: lo <= i < hi)
    assert total == 949
    cases, cache = [], {}
    for row in survivors[lo:hi]:
        di, u = row['left']; dj, v0 = row['right']
        match = tuple(map(tuple, row['match']))
        native, mp = old.D.union(grams[di], grams[dj], match)
        v = [mp[i] for i in v0]
        assert min(old.P.distances(native, u)[j] for j in v) > 9
        cross = old.extra_cross(native, mp, match, u, v, di, dj, grams)
        key = str((di, dj, match))
        if key not in cache:
            g, scales = old.rationalize(native)
            ks = g.nullspace(); iv = None if ks else g.inv()
            cert = old.exact_inertia(g); assert cert['signature'][1] == 1
            core = dict(left_diagram=di, right_diagram=dj, match=list(map(list, match)),
                        gram=old.matrix(g), gauge_squared=scales, **cert)
            if ks:
                assert len(ks) == 1
                core['radical'] = list(map(str, ks[0]))
            else:
                core['inverse'] = old.matrix(iv)
            equal(core, saved('isolation', 'cores', key), 'whole isolation core')
            cache[key] = (g, scales, ks, iv)
        g, scales, ks, iv = cache[key]
        du, dv = old.P.distances(native, u), old.P.distances(native, v)
        corners, choices = [], []
        for p, q in product(u, v):
            corner = old.corner(g, scales, u, v, p, q)
            if corner:
                corners.append(corner)
                choice = old.ordinary_corner(g, scales, iv, ks, corner, du, dv)
                if choice:
                    choices.append(choice)
        if not choices:
            assert not ks
            for corner in corners:
                choice = old.relevant_corner(g, scales, iv, corner, du, dv)
                if choice:
                    choices.append(choice)
        assert choices, 'no controlled preferred corner'
        if choices[0]['method'] == 'relevant_component':
            choice = min(choices, key=lambda c: (c['component_order_bound'], len(c['extension_cases']), c['preferred_corner']))
        else:
            choice = min(choices, key=lambda c: (c['component_order_bound'], len(c['possible_contacts'])))
        cases.append(dict(type=row['surviving_type'], left=[di, u], right=[dj, v0], right_copy_map=mp,
                          match=list(map(list, match)), core_key=key, left_black=u, right_black=v,
                          extra_cross_bond_checks=cross, **choice))
    equal(cases, expected, 'all isolation case fields in survivor order')


def profiles(old, lo, hi):
    rows, minima = [], []
    for d in range(108 + lo, 108 + hi):
        ranks = list(old.profiles(d))
        assert len(ranks) == len(set(tuple(sorted(rs)) for rs in ranks))
        group = [old.budget(rs) for rs in ranks]
        assert all(old.F(row['margin']) > 0 for row in group)
        n = d + 1
        beta = old.F(n*(n-6), 8) if n % 2 == 0 else old.F((n-1)*(n-5), 8)
        assert beta > 9*n
        minima.append(dict(ambient_dimension=n, profiles=len(group),
                           minimum_margin=str(min(old.F(row['margin']) for row in group)), finite_slack=str(beta-9*n)))
        rows.extend(group)
    expected, total = Arrays.select(('reduced_profiles',), lambda i, row: 108+lo <= row['link_dimension'] < 108+hi)
    assert total == 756
    equal(rows, expected, 'all reduced profiles')
    expected, total = Arrays.select(('dimension_minima',), lambda i, row: lo <= i < hi)
    assert total == 8
    equal(minima, expected, 'dimension minima')


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
    parser.add_argument('--phase', choices=tuple(LIMITS), required=True)
    parser.add_argument('--start', type=int, default=0)
    parser.add_argument('--stop', type=int)
    args = parser.parse_args()
    lo, hi = args.start, args.stop if args.stop is not None else args.start + 1
    assert 0 <= lo < hi <= LIMITS[args.phase], 'phase interval outside complete domain'
    owner = parent_lock()
    signal.alarm(45)  # Wrapper remains the authoritative59s/1GiB resource cap.
    self_path = Path(__file__).relative_to(ROOT).as_posix()
    own_hash = digest(self_path)
    for path, expected in FIXED.items():
        assert digest(path) == expected, ('fixed source pin', path)
    snapshot = saved('sources')
    assert isinstance(snapshot, dict) and len(snapshot) == 23
    before = {path: digest(path) for path in snapshot}
    assert before == snapshot, 'outer source pins'
    stems = (
        'fresh_high_rank_quasilanner_census_20261017',
        'fresh_simple_cusp_distance_bound_20261017',
        'fresh_selective_distance10_dimension116_20261017',
        'fresh_selective_distance10_dimension116_20261017_review',
        'fresh_distance9_low_core_20261017',
        'fresh_cusp_factor_payment_framework_20261017',
        'fresh_cusp_link_shape_curvature_20261017',
    )
    checked, detail = dict(before), {}
    # Repeat ALL inherited status/pin checks before each mathematical batch.
    # Projection prevents decoding their scientific tables; digest() streams.
    input_pins(stems, checked)
    old = importlib.import_module(OLD)  # Definition-only: never old.run().
    assert old.INPUTS == stems and len(set(stems)) == 7
    paths = {SOURCE, str(Path(old.D.__file__).relative_to(ROOT))} | {
        f'{folder}/{stem}.{ext}' for stem in old.INPUTS
        for folder, ext in (('scripts', 'py'), ('results', 'json'), ('results', 'md'))}
    assert paths == set(snapshot), 'original full source domain'
    if args.phase == 'profiles':
        profiles(old, lo, hi)
    else:
        _, grams, contexts, ports, values = setup(old)
        assert len({context_key(c) for c in contexts}) == len(contexts) == 384
        if args.phase == 'metadata':
            metadata(old, contexts)
        elif args.phase == 'domains':
            domains(old, grams, contexts, ports)
        elif args.phase == 'overlap-domain':
            overlap_keys(old, grams, contexts, lo, hi)
        elif args.phase == 'disjoint':
            disjoint(old, grams, contexts, ports, lo, hi)
        elif args.phase == 'inertia':
            inertia(old, grams, ports, lo, hi)
        elif args.phase == 'overlaps':
            detail = overlaps(old, grams, contexts, values, lo, hi)
        elif args.phase == 'isolation':
            isolation(old, grams, lo, hi)
        else:
            raise AssertionError('unreachable phase')
    gc.collect()
    assert os.getppid() == owner, 'parent exited'
    assert {path: digest(path) for path in checked} == checked, 'input pins changed'
    assert digest(self_path) == own_hash, 'replay worker changed'
    for path, expected in FIXED.items():
        assert digest(path) == expected, ('fixed pin changed', path)
    print(json.dumps(dict(status='PASS_DISTANCE9_REPLAY_PARTITION', phase=args.phase,
                          start=lo, stop=hi, detail=detail, complete_phase_domains=LIMITS,
                          source=own_hash, certificate=FIXED[CERT]), sort_keys=True), flush=True)


if __name__ == '__main__':
    main()
