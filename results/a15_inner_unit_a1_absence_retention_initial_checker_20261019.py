#!/usr/bin/env python3
"""Solver-free A15(1,13) output-retention audit over saved physical interfaces.

This checks exact finite interfaces. Universal contact-domain completeness,
short-path routing and arbitrary-support consequences are written premises,
identified in the accompanying report. No producer or scientific engine is
imported; no historical census or forest optimization is executed.
"""
import argparse
import ast
from collections import Counter, deque
import fcntl
from fractions import Fraction as F
from functools import lru_cache
from hashlib import file_digest
import json
from pathlib import Path
import resource
import signal
import sys

import run_lean_batches_20261018 as guard

ROOT = Path(__file__).resolve().parents[1]
TAG = 'a15_inner_unit_a1_absence_retention_20261018'
OVER = 'fresh_distance8_high_overlap_20261017'
SUPPORT = 'fresh_finite_rank15_long_pair_support_scout_20261018'
OLD16 = 'fresh_order16_finite_cross_payments_20261017'
OLD17 = 'fresh_ideal7_finite_cross_payments_20261017'
D13 = 'fresh_affine_d13_long_pair_transfer_20261017'
D14 = 'fresh_affine_d14_long_pair_transfer_20261017'
PARTIAL = 'fresh_d17_partial_weight_transfer_20261017'
D16PART = 'fresh_d16_partial_weight_transfer_20261017'
NEWD = 'fresh_affine_d11_leaf_pair_transfer_20261017'
A14 = 'fresh_finite_a14_near_endpoint_zero_transfer_20261018'
CHOSEN = ('fresh_d12_radical_saturated_half_replacements_20261018',
          'fresh_d12_remaining_selected_half_replacements_20261018')
NEWER = (
 ('fresh_rank15_sparse_payment_exchange_20261018', ('selected_proofs',)),
 ('fresh_rank16_a1_trade_composed_ambient58_frontier_20261018', ('local_proof', 'whole_proofs')),
 ('fresh_d16_a1_long_leaf_d17_zero_ambient58_frontier_20261018', ('local_proof', 'whole_proofs')),
 ('fresh_a16_a1_seven_eighth_second_coupon_20261018', ('local_proof', 'whole_proofs')),
 ('fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018', ('local_proof', 'proofs')),
 ('fresh_rank13_variable_context_retention_20261018', ('local_seals',)),
 ('fresh_a14_endpoint_a1_context_half_transfer_20261018', ('local_proof',)),
 ('fresh_d17_leaf_a16_context_coupled_transfer_20261018', ('local_proof', 'old_ideal17_seals')),
 ('fresh_d17_leaf_a16_context_coupled_transfer_20261018', ('local_proof', 'A16_seals')),
 ('fresh_a16_double_large_d_physical_scout_20261018', ('local', 'selected_and_later_seals')),
 ('fresh_d17_node3_large_d_guard_transfer_20261018', ('local_proof', 'proofs')),
 ('fresh_d16_node3_rank16_large_d_guard_scout_20261018', ('geometry',)),
 ('fresh_a17_inner_unit_rank16_rank17_guard_scout_20261018', ('geometry',)),
 ('fresh_order16_cd_sparse_sixteenth_transfer_20261018', ('rules',)),
 ('fresh_order10_d_opposite_fork_fractional_transfer_20261018', ('mixed_whole_corners',)),
 ('fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018', ('local',)),
 ('fresh_order11_fractional_cusp_closure_local_20261018', ('B_corner_templates',)),
 ('fresh_order11_fractional_cusp_closure_local_20261018', ('D_next', 'templates')),
 ('fresh_order11_fractional_cusp_closure_local_20261018', ('B_zero_repair',)),
)
# These later modifications only extend already sealed external witnesses or
# change named typed coupons. Their saved interfaces are nevertheless read.
EXTENSIONS = (
 'fresh_d16_rank17_witness_extension_scout_20261018',
 'fresh_a16_d16_spectator_extension_scout_20261018',
 'fresh_a16_rank17_spectator_extension_scout_20261018',
 'fresh_order10_d_quarter_ideal57_scout_20261018',
 'fresh_order16_d_quarter_epsilon_ideal57_scout_20261018',
)
LARGE_SECTIONS = {
 'fresh_d17_a9_a2_context_zero_ambient58_frontier_20261018': ('local_proof',),
 'fresh_order16_d_quarter_epsilon_ideal57_scout_20261018': ('local',),
 'fresh_order17_d13_shortleaf_seven_thirtyseconds_scout_20261018': ('local',),
}


def path(stem):
    return ROOT / 'results' / (stem + '.json')


def digest(p):
    with p.open('rb') as stream:
        return file_digest(stream, 'sha256').hexdigest()


def memory():
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform != 'darwin':
        rss *= 1024
    assert rss < 1024**3, 'one-GiB checker memory cap'
    return rss


def load(stem):
    p = path(stem)
    if p.stat().st_size > 4 * 1024**2:
        assert stem in LARGE_SECTIONS, (stem, 'declare selected sections; no blind large load')
        result = {key: saved_value(stem, key)
                  for key in ('status', *LARGE_SECTIONS[stem])}
    else:
        with p.open() as stream:
            result = json.load(stream)
    memory()
    assert str(result.get('status', '')).startswith('PASS'), stem
    return result


def finite_inventory():
    source = ROOT / 'scripts/fresh_d12_half_face_composition_20261018.py'
    tree = ast.parse(source.read_text())
    values = [ast.literal_eval(n.value) for n in tree.body
              if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name)
                 and t.id == 'FINITE' for t in n.targets)]
    assert len(values) == 1 and len(values[0]) == 22
    return values[0]


def input_paths():
    stems = {OVER, SUPPORT, OLD16, OLD17, D13, D14, PARTIAL, D16PART, NEWD,
             A14, 'fresh_rank15_sparse_payment_exchange_20261018',
             *CHOSEN, *EXTENSIONS, *finite_inventory(), *(s for s, _ in NEWER)}
    return sorted(path(s) for s in stems)


def streamed_value(stem, marker, back):
    # Marker selectors are fixed, source-bound field names/core identities.
    # No whole-file scientific JSON load is used above the4MiB threshold.
    with path(stem).open('rb') as stream:
        buf = b''
        while True:
            chunk = stream.read(65536)
            if not chunk:
                raise AssertionError(('missing saved marker', stem, marker))
            buf += chunk
            at = buf.find(marker)
            if at < 0:
                buf = buf[-len(marker):]
                continue
            buf = buf[at + len(marker) - back:]
            while True:
                try:
                    return json.JSONDecoder().raw_decode(buf.decode().lstrip())[0]
                except json.JSONDecodeError:
                    chunk = stream.read(65536)
                    assert chunk, 'truncated selected overlap object'
                    buf += chunk
                    assert len(buf) < 4 * 1024**2, 'selected object cap'


def saved_value(stem, key):
    # All selected immutable inputs use an unspaced key/colon separator.
    marker = (json.dumps(key)+':').encode()
    return streamed_value(stem, marker, 0)


def native(index):
    marker = ('{"index":' + str(index) + ',').encode()
    case = streamed_value(OVER, marker, len(marker))
    key = (json.dumps(case['core_key'], separators=(',', ':')) + ':{').encode()
    return case, streamed_value(OVER, key, 1)


@lru_cache(maxsize=128)
def matrix(raw):
    return tuple(tuple(F(x) for x in row) for row in raw)


def mat(raw):
    return matrix(tuple(tuple(row) for row in raw))


def parts(g, ids):
    todo = set(ids)
    result = []
    while todo:
        start = min(todo)
        seen = {start}
        stack = [start]
        while stack:
            i = stack.pop()
            for j in todo - seen:
                if g[i][j]:
                    seen.add(j)
                    stack.append(j)
        todo -= seen
        result.append(sorted(seen))
    return result


def pd(g, ids):
    a = [[g[i][j] for j in ids] for i in ids]
    for k in range(len(a)):
        pivot = a[k][k]
        if pivot <= 0:
            return False
        for i in range(k + 1, len(a)):
            for j in range(i, len(a)):
                a[i][j] -= a[i][k] * a[k][j] / pivot
                a[j][i] = a[i][j]
    return True


def distances(g, starts):
    d = {i: 0 for i in starts}
    queue = deque(starts)
    while queue:
        i = queue.popleft()
        for j in range(len(g)):
            if i != j and g[i][j] and j not in d:
                d[j] = d[i] + 1
                queue.append(j)
    return d


def a15_role(raw, white, marks, ids=None):
    ids = sorted(set(white or ()) | set(marks)) if ids is None else sorted(ids)
    if len(ids) < 15:
        return None
    g = mat(raw)
    blocks = parts(g, ids)
    ns = next((b for b in blocks if set(marks) <= set(b)), None)
    if ns is None or len(ns) != 15:
        return None
    adj = {i: [j for j in ns if i != j and g[i][j]] for i in ns}
    ends = [i for i in ns if len(adj[i]) == 1]
    if len(ends) != 2 or any(len(adj[i]) not in (1, 2) for i in ns):
        return None
    for i in ns:
        if g[i][i] <= 0:
            return None
        for j in adj[i]:
            if g[i][j] >= 0 or 4 * g[i][j]**2 != g[i][i] * g[j][j]:
                return None
    order = [min(ends)]
    while len(order) < 15:
        nxt = [j for j in adj[order[-1]] if j not in order]
        assert len(nxt) == 1
        order.append(nxt[0])
    marked = sorted(order.index(m) for m in marks)
    images = sorted({tuple(marked), tuple(sorted(14-i for i in marked))})
    if images != [(1, 13)]:
        return None
    if not pd(g, ids):
        return dict(nonfinite=True, images=images)
    return dict(nonfinite=False, images=images, nodes=ns,
                displayed_A1=[b[0] for b in blocks if len(b) == 1],
                corner_nodes=ids)


def physical_rows(x, prefix='', g=None, white=None, index=None,
                  require_components=False):
    if isinstance(x, dict):
        inherited = x.get('old_complete_seal', x.get('old_seal', {}))
        g = x.get('completed_gram', x.get('gram', inherited.get('gram', g)))
        white = x.get('completed_white', x.get('white', inherited.get('white', white)))
        if white is None and isinstance(x.get('edge_checks'), dict):
            white = x['edge_checks'].get('white')
        index = x.get('index', x.get('case', index))
        cs = x.get('corners')
        if white is None and g is not None and isinstance(cs, list) and cs:
            black = {m for c in cs if isinstance(c, dict) for m in c.get('marks', [])}
            if len(black) == 4:
                white = sorted(set(range(len(g))) - black)
        blocks = x.get('components', x.get('whole_components'))
        if ('marks' in x and g is not None
                and (white is not None or 'corner_nodes' in x)
                and (blocks is not None or not require_components)):
            yield prefix, x, g, white, index
        if not require_components and g is not None and white is not None:
            for key in ('first_corner', 'second_corner'):
                marks = x.get(key)
                if isinstance(marks, list) and len(marks) == 2 and all(isinstance(i,int) for i in marks):
                    yield prefix+'/'+key, dict(marks=marks), g, white, index
        for key, value in x.items():
            yield from physical_rows(value, prefix+'/'+key, g, white, index,
                                     require_components)
    elif isinstance(x, list):
        for i, value in enumerate(x):
            yield from physical_rows(value, prefix+'/'+str(i), g, white, index,
                                     require_components)


def d13_gram(row):
    g = [[2 * int(i == j) for j in range(19)] for i in range(19)]
    left, right = row['leaf_attachments']
    edges = [(i, i+1) for i in range(9)] + [(0,10),(0,11),(9,12),(9,13)]
    edges += [(14,left),(15,right),(16,row['A3_attached_to']),(16,17),(17,18)]
    for i, j in edges:
        g[i][j] = g[j][i] = -1
    return g


def literal_rows(common_gram):
    for stem in (OLD16, OLD17):
        doc = load(stem)
        saved = doc['exceptional_D15']['cases'] if stem == OLD16 else doc['order17']['replacement_cases']
        g = doc['exceptional_D15']['gram'] if stem == OLD16 else common_gram
        for row in saved:
            white = sorted(set(range(len(g))) - set(row['black']))
            yield stem, row['case'], g, white, row['corners']
    for stem, key in ((D13, 'extreme_equality'), (D14, 'equality')):
        for i, row in enumerate(load(stem)[key]['rows']):
            yield stem, i, d13_gram(row) if stem == D13 else row['gram'], row['white'], row['corners']
    for row in load(PARTIAL)['cases']:
        yield PARTIAL, row['index'], row['gram'], row['white'], row['corners']
    for row in load(D16PART)['cases']:
        for i, completed in enumerate(row['completions']):
            if 'corners' in completed:
                yield D16PART, [row['index'],i], completed['gram'], completed['white'], completed['corners']
    for i, row in enumerate(load(NEWD)['two_A1_A3_spectators']['rows']):
        yield NEWD, i, row['gram'], row['white'], row['corners']


def core_proofs():
    c, core = native(1142)
    g = mat(core['gram'])
    z = list(map(F, core['radical'][0]))
    assert len(g) == 19 and core['gauge'] == [1]*19
    assert all(sum(g[i][j]*z[j] for j in range(19)) == 0 for i in range(19))
    u = [max(F(), -x) for x in z]
    u[9] += F(1,2)
    timelike = sum(u[i]*g[i][j]*u[j] for i in range(19) for j in range(19))
    assert timelike < 0
    allowed = c['possible_contacts']
    assert allowed == [0,1,2,3,10,11,12,13] and all(z[i] for i in allowed)
    left = c['left'][1]
    right = [c['right_copy_map'][i] for i in c['right'][1]]
    white = sorted(set(range(19)) - set(left+right))
    assert pd(g, white)
    assert all(pd(g, sorted(white+[b])) for b in left+right)
    dl, dr = distances(g, left), distances(g, right)
    bounds = [min(dl[i]+dr[j], dl[j]+dr[i])+2 for i in allowed for j in allowed if z[i]*z[j] < 0]
    assert bounds and max(bounds) <= 4
    selected = a15_role(core['gram'], white, c['preferred_corner'])
    assert selected['displayed_A1'] == [3,13]
    # All15 target labels and every reversal map are checked without graph search.
    support = next(r for r in load(SUPPORT)['records'] if r['kind'] == 'A')
    assert support['complete_eligible_patterns'] == 536 and support['complete_maps'] == 482
    hits = [r for r in support['selected_patterns'] if [1,13] in r['marked_images']]
    assert len(hits) == 1 and hits[0]['index'] == 1142 and hits[0]['displayed_order'] == 15
    assert len(hits[0]['maps']) == 2
    for mp in hits[0]['maps']:
        assert mp['marked_images'] == [1,13]
        ns = c['marked_component']
        assert len(mp['images']) == len(ns) == len(set(mp['images'])) == 15
        assert all((g[a][b] == -1) == (abs(mp['images'][i]-mp['images'][j]) == 1)
                   for i,a in enumerate(ns) for j,b in enumerate(ns) if i != j)
    proofs = [dict(index=1142, source='selected_A15', A1=[3,13],
                   radical_relation=True, timelike_square=str(timelike),
                   complete_allowed_contacts=allowed, short_bound=max(bounds),
                   arbitrary_support='written radical-sign and first-contact proof')]
    consumers = {}
    for index, marks, a1 in ((1157,[0,12],3),(1617,[2,10],13)):
        case, data = native(index)
        gg = mat(data['gram']); J = mat(data['inverse']); n = len(gg)
        assert n == 18 and data['gauge'] == [1]*n
        assert all(sum(gg[i][k]*J[k][j] for k in range(n)) == int(i == j)
                   for i in range(n) for j in range(n))
        ids = case['marked_component']; remain = sorted(set(case['possible_contacts'])-set(ids))
        assert remain == [2,3,12,13] and all(J[i][i] >= 4 for i in remain)
        ends = case['end_sets']
        assert set(ends[0]).isdisjoint(ends[1]) and set(ends[0]+ends[1]) == set(case['possible_contacts'])
        assert all(J[i][j] > 0 for end in ends for i in end for j in end)
        l = case['left'][1]; r = [case['right_copy_map'][i] for i in case['right'][1]]
        dl, dr = distances(gg,l), distances(gg,r)
        bound = max(dl[i] for i in ends[0]) + max(dr[i] for i in ends[1]) + 2
        assert bound <= 8
        w = sorted(set(range(n))-set(l+r))
        out = a15_role(data['gram'], w, marks)
        assert out and out['displayed_A1'] == [a1]
        for corner in case['all_corners']:
            assert pd(gg, sorted(set(w)|set(corner['preferred_corner'])))
        consumers[index] = dict(gram=data['gram'], white=w, marks=marks, A1=a1)
        proofs.append(dict(index=index, source='whole_A14_projection_consumer', A1=[a1],
                     remaining_contacts=remain, inverse_diagonals=[str(J[i][i]) for i in remain],
                     short_bound=bound, arbitrary_support='written positive-end/projection and first-contact proof'))
    return core['gram'], proofs, consumers


def run():
    common, proofs, consumers = core_proofs()
    counts = Counter(); hits = []; gaps = []; nonfinite = 0
    def one(stem, where, raw, white, row, index=None):
        nonlocal nonfinite
        counts[stem] += 1
        result = a15_role(raw, white, row['marks'], row.get('corner_nodes'))
        if not result:
            return
        if result['nonfinite']:
            nonfinite += 1
            return
        hit = dict(source=stem, path=where, index=index, marks=row['marks'], **result)
        if (index == 1142 and result['displayed_A1'] == [3,13]
                and row['marks'] == [0,10]
                and sorted(white) == sorted(set(range(19))-{0,2,10,12})
                and mat(raw) == mat(common)):
            hit['whole_A1_proof'] = 'selected1142_full_radical'
        elif isinstance(index, int) and index in consumers and (mat(raw) == mat(consumers[index]['gram'])
                and sorted(white) == consumers[index]['white']
                and row['marks'] == consumers[index]['marks']
                and result['displayed_A1'] == [consumers[index]['A1']]):
            # Only the original selected-A14 route has the original whole-A14
            # premise. A matching subcore in another route is NOT enough.
            allowed = {A14, 'fresh_finite_rank14_existing_payment_transfer_20261018'}
            if stem in allowed:
                hit['whole_A1_proof'] = 'original_whole_A14_full_core_seal'
            else:
                gaps.append(hit)
        else:
            gaps.append(hit)
        hits.append(hit)
    legacy = 0
    for stem, keys in finite_inventory().items():
        doc = load(stem)
        for key in keys:
            for where,row,g,w,index in physical_rows(doc[key],key,require_components=True):
                one(stem,where,g,w,row,index);legacy += 1
        del doc
    assert legacy == 564, ('legacy interface coverage mismatch',legacy)
    literal = 0
    for stem,index,g,w,cs in literal_rows(common):
        for i,row in enumerate(cs):
            one(stem,'literal/'+str(index)+'/'+str(i),g,w,row,index)
        literal += 1
    assert literal == 216
    # Read one file at a time. Some sections of one file are repeated by design.
    newer_counts = []
    for stem, keys in NEWER:
        doc = load(stem); section = doc
        for key in keys:
            section = section[key]
        before = counts[stem]
        for where,row,g,w,index in physical_rows(section,'/'.join(keys)):
            one(stem,where,g,w,row,index)
        assert counts[stem] > before, ('uncovered physical section',stem,keys)
        newer_counts.append(dict(source=stem,section=list(keys),corners=counts[stem]-before))
        del doc, section
    chosen = []
    def role_rows(x):
        if isinstance(x,dict):
            if 'payment' in x and isinstance(x.get('roles'),list):
                yield x
            for v in x.values():
                yield from role_rows(v)
        elif isinstance(x,list):
            for v in x:
                yield from role_rows(v)
    for stem, expected in zip(CHOSEN,(280,52)):
        rs = [r for r in role_rows(load(stem)) if F(r['payment']) > 0]
        assert len(rs) == expected
        assert all(not any(b.get('type') == 'A15' for b in r['roles']) for r in rs)
        chosen.append(dict(source=stem,positive_outputs=expected,A15_outputs=0))
    extension_records = []
    for stem in EXTENSIONS:
        doc = load(stem)
        # These records may contain additional copied physical corners.
        before = counts[stem]
        for where,row,g,w,index in physical_rows(doc):
            one(stem,where,g,w,row,index)
        extension_records.append(dict(source=stem,corners=counts[stem]-before,
            rule=doc.get('rule',doc.get('local',{}).get('new_node11_coupon'))))
        del doc
    if gaps:
        raise AssertionError(json.dumps(dict(status='UNRESOLVED_ACTUAL_A1_SEAL',gaps=gaps[:12])))
    assert hits and any(h['source'] == A14 for h in hits), 'expected A14 consumers missing'
    files = input_paths()
    sources = {str(p.relative_to(ROOT)):digest(p) for p in files}
    for p in (Path(__file__), ROOT/'results/fresh_a15_inner_unit_a1_absence_research_20261018.md',
              ROOT/'scripts/fresh_d12_half_face_composition_20261018.py'):
        sources[str(p.relative_to(ROOT))] = digest(p)
    return dict(status='PASS_A15_INNER_UNIT_A1_ABSENCE_SAVED_OUTPUT_RETENTION',
        sources=sources,input_files=len(files),input_bytes=sum(p.stat().st_size for p in files),
        selected_patterns=536,selected_maps=482,selected_target_contexts=[1142],
        physical_whole_A1_seals=proofs,later_packets=22,later_corners=legacy,
        literal_configurations=literal,corner_counts=dict(counts),target_hits=hits,
        newer_sections=newer_counts,chosen_D12=chosen,extensions=extension_records,
        excluded_nonfinite_A15_displays=nonfinite,new_positive_finite_families=0,
        exact_ordinary_delta='-a15*[a1=0]',new_dimension_bound=None,
        universal_geometry_formalized=False,
        scope='Exact saved output-role audit plus finite core algebra. Universal selected/nonselected coverage, complete contact domains and actual whole-factor consequences are the written proofs, not conclusions from hashes or displayed subcores.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='Create new certificate only after every check passes.')
    parser.add_argument('--verify',action='store_true',help='Reconstruct and compare the saved certificate.')
    args = parser.parse_args()
    assert not (args.write and args.verify)
    target = ROOT/'results'/(TAG+'.json')
    if args.write:
        assert not target.exists(), 'refuse to overwrite evidence'
    guard.resources(); memory(); signal.alarm(59)
    with (ROOT/'formal/.proof-batches/exclusive.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        data = run(); guard.resources(start=False); peak = memory()
        if args.write:
            target.write_text(json.dumps(data,indent=2)+'\n')
        if args.verify:
            assert data == json.loads(target.read_text())
    signal.alarm(0)
    print(json.dumps(dict(status=data['status'],input_files=data['input_files'],
        input_bytes=data['input_bytes'],later_corners=data['later_corners'],
        literal_configurations=data['literal_configurations'],target_hits=len(data['target_hits']),
        peak_rss_bytes=peak,new_dimension_bound=None)),flush=True)


if __name__ == '__main__':
    main()
